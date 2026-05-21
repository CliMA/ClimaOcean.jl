# KPPVerticalDiffusivity closure type, dispatch interface, and closure-field allocation.

# Wrapper that opts the per-tracer surface-BC NamedTuple out of recursive
# `fill_halo_regions!` traversal in `update_state!`.
struct KPPTracerBCs{T}
    bcs :: T
end

@inline Base.getindex(k::KPPTracerBCs, i) = k.bcs[i]

Adapt.adapt_structure(to, k::KPPTracerBCs) = KPPTracerBCs(adapt(to, k.bcs))

fill_halo_regions!(::KPPTracerBCs, args...; kwargs...) = nothing

"""
    KPPVerticalDiffusivity{TD, P, FT}

K-Profile Parameterization vertical mixing closure (Large, McWilliams, & Doney 1994).
"""
struct KPPVerticalDiffusivity{TD, P, FT} <: AbstractScalarDiffusivity{TD, VerticalFormulation, 1}
    parameters          :: P
    maximum_viscosity   :: FT
    maximum_diffusivity :: FT
end

"""
    KPPVerticalDiffusivity(time_discretization = VerticallyImplicitTimeDiscretization(),
                           FT = Float64;
                           parameters = KPPParameters(FT),
                           maximum_viscosity = 1.0,
                           maximum_diffusivity = 1.0)

Construct a `KPPVerticalDiffusivity` with MITgcm defaults. The radiation profile
(for SW penetration) is auto-detected at compute time via `get_radiative_forcing(model)`.
"""
function KPPVerticalDiffusivity(time_discretization = VerticallyImplicitTimeDiscretization(),
                                FT::DataType        = Float64;
                                parameters          = KPPParameters(FT),
                                maximum_viscosity   = 1.0,
                                maximum_diffusivity = 1.0)

    TD = typeof(time_discretization)
    P  = typeof(parameters)
    return KPPVerticalDiffusivity{TD, P, FT}(parameters,
                                             convert(FT, maximum_viscosity),
                                             convert(FT, maximum_diffusivity))
end

KPPVerticalDiffusivity(FT::DataType; kw...) = KPPVerticalDiffusivity(VerticallyImplicitTimeDiscretization(), FT; kw...)

Adapt.adapt_structure(to, clo::KPPVerticalDiffusivity{TD, P, FT}) where {TD, P, FT} =
    KPPVerticalDiffusivity{TD, P, FT}(adapt(to, clo.parameters),
                                      clo.maximum_viscosity,
                                      clo.maximum_diffusivity)

#####
##### Type aliases and dispatch
#####

const KPPVD       = KPPVerticalDiffusivity
const KPPVDArray  = AbstractArray{<:KPPVD}
const FlavorOfKPP = Union{KPPVD, KPPVDArray}

@inline viscosity_location(::FlavorOfKPP)   = (Center(), Center(), Face())
@inline diffusivity_location(::FlavorOfKPP) = (Center(), Center(), Face())

@inline viscosity(::FlavorOfKPP, K)       = K.κu
@inline diffusivity(::FlavorOfKPP, K, id) = K.κc

with_tracers(tracers, closure::FlavorOfKPP) = closure

@inline time_discretization(::KPPVerticalDiffusivity{TD}) where TD = TD()

#####
##### Closure-field allocation
#####

function build_closure_fields(grid, clock, tracer_names, bcs, closure::FlavorOfKPP)
    # Output diffusivities + nonlocal-transport coefficient.
    κu = Field((Center(), Center(), Face()), grid)
    κc = Field((Center(), Center(), Face()), grid)
    γ  = Field((Center(), Center(), Face()), grid)

    # Cached interior K (ν, κ) at every face — filled once per step, reused by
    # the column sweep and by the per-interface kernel.
    νᵢ = Field((Center(), Center(), Face()), grid)
    κᵢ = Field((Center(), Center(), Face()), grid)

    # Column-level scalars (one value per (i, j) column).
    hbl  = Field{Center, Center, Nothing}(grid)
    u★   = Field{Center, Center, Nothing}(grid)
    Bo   = Field{Center, Center, Nothing}(grid)
    α    = Field{Center, Center, Nothing}(grid)
    G1u  = Field{Center, Center, Nothing}(grid)
    dG1u = Field{Center, Center, Nothing}(grid)
    G1s  = Field{Center, Center, Nothing}(grid)
    dG1s = Field{Center, Center, Nothing}(grid)

    top_tracer_bcs = KPPTracerBCs(NamedTuple(name => bcs[name].top for name in tracer_names))

    return (; κu, κc, γ, νᵢ, κᵢ, hbl, u★, Bo, α, G1u, dG1u, G1s, dG1s, top_tracer_bcs)
end

#####
##### Show
#####

function Base.summary(closure::KPPVerticalDiffusivity)
    TD = nameof(typeof(time_discretization(closure)))
    return string("KPPVerticalDiffusivity{", TD, "}")
end

function Base.show(io::IO, closure::KPPVerticalDiffusivity)
    p = closure.parameters
    print(io, summary(closure), '\n',
              "├── Riᶜ:  ", prettysummary(p.Riᶜ),  '\n',
              "├── Ri∞:  ", prettysummary(p.Ri∞),  '\n',
              "├── κᵥ:   ", prettysummary(p.κᵥ),   '\n',
              "├── ν₀ˢʰ: ", prettysummary(p.ν₀ˢʰ), '\n',
              "├── κ₀ˢʰ: ", prettysummary(p.κ₀ˢʰ), '\n',
              "├── C★:   ", prettysummary(p.C★),   '\n',
              "├── maximum_viscosity:   ", prettysummary(closure.maximum_viscosity),   '\n',
              "└── maximum_diffusivity: ", prettysummary(closure.maximum_diffusivity))
end
