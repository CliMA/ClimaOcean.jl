# NEMOTKEVerticalDiffusivity closure type, dispatch interface, and closure-field allocation.

"""
    NEMOTKEVerticalDiffusivity{TD, P, FT}

NEMO 3.6 TKE vertical-mixing closure (`zdftke` + `zdfevd`).
References: Blanke & Delecluse 1993; Gaspar et al. 1990; Madec et al. 2017.
"""
struct NEMOTKEVerticalDiffusivity{TD, P, FT} <: AbstractScalarDiffusivity{TD, VerticalFormulation, 1}
    parameters          :: P
    maximum_viscosity   :: FT
    maximum_diffusivity :: FT
end

"""
    NEMOTKEVerticalDiffusivity(time_discretization = VerticallyImplicitTimeDiscretization(),
                               FT = Float64;
                               parameters = nothing,
                               maximum_viscosity = 1.0,
                               maximum_diffusivity = 1.0,
                               kwargs...)

Construct a `NEMOTKEVerticalDiffusivity` with NEMO 3.6 ORCA1 OMIP-2 defaults.
Pass extra keyword args (e.g., `Cᴸ = 0.2`) to override individual `NEMOTKEParameters`.
If `parameters` is supplied, individual parameter kwargs are not allowed.
"""
function NEMOTKEVerticalDiffusivity(time_discretization = VerticallyImplicitTimeDiscretization(),
                                    FT::DataType        = Float64;
                                    parameters          = nothing,
                                    maximum_viscosity   = 1.0,
                                    maximum_diffusivity = 1.0,
                                    kwargs...)

    if !isnothing(parameters) && !isempty(kwargs)
        throw(ArgumentError("Pass either `parameters` or individual parameter kwargs, not both."))
    end

    p = isnothing(parameters) ? NEMOTKEParameters(FT; kwargs...) : parameters
    TD = typeof(time_discretization)
    P  = typeof(p)
    return NEMOTKEVerticalDiffusivity{TD, P, FT}(p,
                                                 convert(FT, maximum_viscosity),
                                                 convert(FT, maximum_diffusivity))
end

NEMOTKEVerticalDiffusivity(FT::DataType; kw...) =
    NEMOTKEVerticalDiffusivity(VerticallyImplicitTimeDiscretization(), FT; kw...)

Adapt.adapt_structure(to, clo::NEMOTKEVerticalDiffusivity{TD, P, FT}) where {TD, P, FT} =
    NEMOTKEVerticalDiffusivity{TD, P, FT}(adapt(to, clo.parameters),
                                          clo.maximum_viscosity,
                                          clo.maximum_diffusivity)

#####
##### Type aliases and dispatch
#####

const NEMOTKEVD       = NEMOTKEVerticalDiffusivity
const NEMOTKEVDArray  = AbstractArray{<:NEMOTKEVD}
const FlavorOfNEMOTKE = Union{NEMOTKEVD, NEMOTKEVDArray}

@inline viscosity_location(::FlavorOfNEMOTKE)   = (Center(), Center(), Face())
@inline diffusivity_location(::FlavorOfNEMOTKE) = (Center(), Center(), Face())

@inline viscosity(::FlavorOfNEMOTKE, K)       = K.κu
@inline diffusivity(::FlavorOfNEMOTKE, K, id) = K.κc

with_tracers(tracers, closure::FlavorOfNEMOTKE) = closure

@inline time_discretization(::NEMOTKEVerticalDiffusivity{TD}) where TD = TD()

#####
##### Closure-field allocation
#####
#
# Storage layout:
#   e           — TKE (Center, Center, Center). Persistent prognostic field.
#   eⁿ          — snapshot of e at start of each outer RK3 step.
#   ℓ           — mixing length (Center, Center, Center).
#   κu          — momentum diffusivity (Center, Center, Face).
#   κc          — tracer diffusivity   (Center, Center, Face).
#   N²          — Brunt-Väisälä² at faces, cached for re-use across phases.
#   u★²         — column scalar (Center, Center, Nothing).
#   τx, τy      — column wind-stress components (Center, Center, Nothing).
#   ℵ           — sea-ice fraction (Center, Center, Nothing).

function build_closure_fields(grid, clock, tracer_names, bcs, closure::FlavorOfNEMOTKE)
    e   = Field((Center(), Center(), Center()), grid)
    eⁿ  = Field((Center(), Center(), Center()), grid)
    γ   = Field((Center(), Center(), Center()), grid)
    ℓ   = Field((Center(), Center(), Center()), grid)
    κu  = Field((Center(), Center(), Face()),   grid)
    κc  = Field((Center(), Center(), Face()),   grid)
    N²  = Field((Center(), Center(), Face()),   grid)

    u★²    = Field{Center, Center, Nothing}(grid)
    τx     = Field{Center, Center, Nothing}(grid)
    τy     = Field{Center, Center, Nothing}(grid)
    ℵ      = Field{Center, Center, Nothing}(grid)

    # Initialize e to the floor so the first sqrt is well-defined. eⁿ and γ
    # are overwritten by step_closure_prognostics! before any read.
    parent(e) .= closure.parameters.minimum_TKE
    parent(ℓ) .= closure.parameters.minimum_mixing_length

    return (; κu, κc, e, eⁿ, ℓ, γ, N², u★², τx, τy, ℵ)
end

#####
##### Show
#####

function Base.summary(closure::NEMOTKEVerticalDiffusivity)
    TD = nameof(typeof(time_discretization(closure)))
    return string("NEMOTKEVerticalDiffusivity{", TD, "}")
end

function Base.show(io::IO, closure::NEMOTKEVerticalDiffusivity)
    p = closure.parameters
    print(io, summary(closure), '\n',
              "├── Cᴷ:    ", prettysummary(p.Cᴷ),    " (rn_ediff)\n",
              "├── Cᴰ:    ", prettysummary(p.Cᴰ),    " (rn_ediss)\n",
              "├── Cᵇ:    ", prettysummary(p.Cᵇ),    " (rn_ebb)\n",
              "├── Cᴸ:    ", prettysummary(p.Cᴸ),    " (rn_lc)\n",
              "├── Cᶠ:    ", prettysummary(p.Cᶠ),    " (rn_efr)\n",
              "├── Cˢ:    ", prettysummary(p.Cˢ),    " (Stokes-drift coeff)\n",
              "├── κᶜⁿᵛ:  ", prettysummary(p.κᶜⁿᵛ),  " (rn_avevd)\n",
              "├── νᵇ:    ", prettysummary(p.νᵇ),    " (avmb)\n",
              "├── κᵇ:    ", prettysummary(p.κᵇ),    " (avtb)\n",
              "├── minimum_TKE:                  ", prettysummary(p.minimum_TKE),           " (rn_emin)\n",
              "├── minimum_surface_TKE:          ", prettysummary(p.minimum_surface_TKE),   " (rn_emin0)\n",
              "├── minimum_mixing_length:        ", prettysummary(p.minimum_mixing_length), " (rn_mxl0)\n",
              "├── mixing_length_formulation:        ", p.mixing_length_formulation,        " (nn_mxl)\n",
              "├── wave_penetration_formulation:     ", p.wave_penetration_formulation,     " (nn_etau)\n",
              "├── surface_length_scale_formulation: ", p.surface_length_scale_formulation, " (nn_htau)\n",
              "├── apply_langmuir_circulation:        ", p.apply_langmuir_circulation, '\n',
              "├── apply_wave_penetration:            ", p.apply_wave_penetration, '\n',
              "├── apply_enhanced_vertical_diffusion: ", p.apply_enhanced_vertical_diffusion, '\n',
              "├── apply_evd_to_momentum:             ", p.apply_evd_to_momentum, '\n',
              "├── apply_prandtl_richardson:          ", p.apply_prandtl_richardson, '\n',
              "├── maximum_viscosity:   ", prettysummary(closure.maximum_viscosity), '\n',
              "└── maximum_diffusivity: ", prettysummary(closure.maximum_diffusivity))
end
