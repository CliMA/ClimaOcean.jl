using Printf
using KernelAbstractions: @index, @kernel
using Oceananigans.Operators: Δzᶜᶜᶜ, ℑxᶠᵃᵃ, ℑyᵃᶠᵃ, ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using Dates: month
using Oceananigans.Grids: λnode, φnode, znode, λnodes, φnodes, Center
using Oceananigans.Architectures: on_architecture, architecture
using Oceananigans.DistributedComputations: @root, Distributed
using Oceananigans.BoundaryConditions: DiscreteBoundaryFunction, getbc, fill_halo_regions!
using Oceananigans.Fields: Field, CenterField, interior
using Oceananigans.Advection: CWENOZ
using Oceananigans.ImmersedBoundaries: bottom_height_field, mask_immersed_field!
using Oceananigans.Utils: launch!
using Adapt: Adapt
using ClimaSeaIce
using ClimaSeaIce.Rheologies: ElastoViscoPlasticRheology
using ClimaSeaIce.SeaIceThermodynamics: LinearLiquidus
using NumericalEarth: DataWrangling
using NumericalEarth.Bathymetry: remove_minor_basins!, atlantic_ocean_basin, pacific_ocean_basin
using NumericalEarth.Oceans: MultipleFluxes, FreshwaterExchange, extract_freshwater_flux, freshwater_exchange
using NumericalEarth.EarthSystemModels.InterfaceComputations: computed_fluxes,
                                                              ConservativeIceFreshwater,
                                                              ScaledIceFreshwater,
                                                              VirtualSaltFluxIceFreshwater,
                                                              ZeroHeatContentMeltwater,
                                                              InterfaceTemperatureMeltwater
using SeawaterPolynomials.TEOS10: Sᴬ_from_Sᴾ, Θ_from_T
using Oceananigans.TurbulenceClosures: IsopycnalSkewSymmetricDiffusivity,
                                       TriadIsopycnalSkewSymmetricDiffusivity,
                                       ConvectiveAdjustmentVerticalDiffusivity,
                                       AdvectiveFormulation, DiffusiveFormulation
using Oceananigans.Utils: NormalDivision
using NumericalEarth.EarthSystemModels.InterfaceComputations: COARELogarithmicSimilarityProfile,
                                                              WindDependentWaveFormulation,
                                                              MomentumRoughnessLength,
                                                              TemperatureDependentAirViscosity,
                                                              ScalarRoughnessLength,
                                                              atmosphere_sea_ice_stability_functions,
                                                              MomentumBasedFrictionVelocity,
                                                              LargeYeagerTransferCoefficients,
                                                              FixedIterations,
                                                              large_yeager_stability_functions,
                                                              RelativeVelocity,
                                                              WindVelocity,
                                                              ConvectiveGustiness

#####
##### Flux configurations
#####

"""
    corrected_atmosphere_ocean_fluxes(FT = Float64)

COARE 3.6-consistent atmosphere-ocean flux formulation with:
- Wind-dependent Charnock parameter (Edson et al. 2013, eq. 13)
- COARE logarithmic similarity profile (no ψ(ℓ/L) term)
- Minimum gustiness = 0.5 m/s (CICE / NCAR CORE-II convention)
- Temperature-dependent air viscosity
"""
function corrected_atmosphere_ocean_fluxes(FT = Float64;
                                           subgrid_velocities = ConvectiveGustiness{FT}(minimum_gustiness = FT(0.5)))
    air_kinematic_viscosity = TemperatureDependentAirViscosity(FT)
    return SimilarityTheoryFluxes(FT;
                                  similarity_form              = COARELogarithmicSimilarityProfile(),
                                  subgrid_velocities           = subgrid_velocities,
                                  momentum_roughness_length    = MomentumRoughnessLength(FT;
                                  wave_formulation             = WindDependentWaveFormulation(FT),
                                  air_kinematic_viscosity      = TemperatureDependentAirViscosity(FT)),
                                  temperature_roughness_length = ScalarRoughnessLength(FT; air_kinematic_viscosity),
                                  water_vapor_roughness_length = ScalarRoughnessLength(FT; air_kinematic_viscosity))
end

"""
    corrected_atmosphere_sea_ice_fluxes(FT = Float64; momentum_roughness_length, scalar_roughness_length)

Atmosphere-sea ice flux formulation with:
- SHEBA/Paulson+Grachev stability functions (existing default, correct)
- Fixed momentum roughness z0 = 5e-4 m (CICE/SHEBA standard; Andreas et al. 2010)
- Fixed scalar roughness z0t = z0q = 5e-5 m (Andreas 1987: z0t ≈ z0/10 at R*≈7)
- COARE logarithmic similarity profile
- Minimum gustiness = 0.2 m/s

Both roughnesses are geometric constants rather than wind-dependent, because sea ice carries no gravity
waves; the roughness is set by ridges, floe edges and sastrugi. The SHEBA value describes multiyear pack,
and smoother first-year ice sits nearer 1e-4 m, so `momentum_roughness_length` is exposed to let the drift
speed be varied over its observed range.
"""
corrected_atmosphere_sea_ice_fluxes(FT = Float64;
                                    momentum_roughness_length = 5e-4,
                                    scalar_roughness_length = 5e-5) =
    SimilarityTheoryFluxes(FT;
                           stability_functions          = atmosphere_sea_ice_stability_functions(FT),
                           similarity_form              = COARELogarithmicSimilarityProfile(),
                           subgrid_velocities           = ConvectiveGustiness{FT}(minimum_gustiness = FT(0.2)),
                           momentum_roughness_length    = FT(momentum_roughness_length),
                           temperature_roughness_length = FT(scalar_roughness_length),
                           water_vapor_roughness_length = FT(scalar_roughness_length))

"""
    corrected_ice_ocean_heat_flux()

Three-equation ice-ocean heat flux with momentum-based friction velocity
computed from actual ice-ocean stress (McPhee 1992, 2008; SHEBA median u*≈0.01 m/s).

`heat_transfer_coefficient` defaults to McPhee's Stanton number 0.0057, which is the value
consistent with a *computed* friction velocity. The 0.0095 of Shi et al. (2021) is calibrated
against a fixed u* = 0.002 m/s; pairing it with u* = sqrt(|tau|/rho), which is ~0.006 m/s under the
Arctic pack, inflates the exchange velocity alpha_h u* by a factor of three.
The salt transfer coefficient follows at the standard ratio R = alpha_h / alpha_s = 35.
"""
corrected_ice_ocean_heat_flux(; heat_transfer_coefficient = 0.0057) =
    ThreeEquationHeatFlux(; heat_transfer_coefficient,
                            salt_transfer_coefficient = heat_transfer_coefficient / 35,
                            friction_velocity = MomentumBasedFrictionVelocity())

"""
    ncar_atmosphere_ocean_fluxes(FT = Float64)

OMIP-2 standard atmosphere-ocean flux formulation using the Large & Yeager
(2004, 2009) bulk algorithm. Iterates directly on transfer coefficients (Cd, Ch, Ce),
NOT on roughness lengths. Uses 5 fixed iterations with Paulson stability functions.
"""
ncar_atmosphere_ocean_fluxes(FT = Float64) =
    CoefficientBasedFluxes(FT;
                           transfer_coefficients = LargeYeagerTransferCoefficients(FT),
                           solver_stop_criteria = FixedIterations(5))

"""
    ncar_atmosphere_sea_ice_fluxes(FT = Float64)

NCAR/CORE atmosphere-sea ice flux formulation with full Monin-Obukhov
similarity theory and stability corrections:
- Paulson (1970) + linear stable (-5ζ) stability functions (same as NCAR ocean)
- Fixed z0 = z0t = z0q = 5e-4 m (CICE default; SHEBA standard)
- Wind speed floor at 0.5 m/s
- COARE logarithmic similarity profile (no ψ(ℓ/L) term)

Over ice the roughness lengths are fixed geometric constants (not wind-dependent),
so the standard MOST roughness-length iteration is consistent here (unlike the
ocean case where the NCAR polynomial Cd requires its own solver).
"""
ncar_atmosphere_sea_ice_fluxes(FT = Float64) =
    SimilarityTheoryFluxes(FT;
                           stability_functions          = large_yeager_stability_functions(FT),
                           similarity_form              = COARELogarithmicSimilarityProfile(),
                           subgrid_velocities           = ConvectiveGustiness{FT}(gustiness_parameter = FT(0),
                                                                                  minimum_gustiness   = FT(0.5)),
                           momentum_roughness_length    = FT(5e-4),
                           temperature_roughness_length = FT(5e-4),
                           water_vapor_roughness_length = FT(5e-4))

"""
    build_coupled_model(ocean, sea_ice, atmosphere, radiation, land, flux_configuration;
                        velocity_formulation = :relative)

Build the `OceanSeaIceModel` with the specified flux configuration.
Options for `flux_configuration`: `:default`, `:corrected`, `:ncar`.
Options for `velocity_formulation`:  `:relative`, `:wind`
"""
function build_coupled_model(ocean, sea_ice, atmosphere, radiation, land, flux_configuration;
                             velocity_formulation::Symbol = :relative,
                             sea_ice_ocean_heat_transfer_coefficient = 0.0057,
                             sea_ice_momentum_roughness_length = 5e-4,
                             ice_freshwater_delivery = ConservativeIceFreshwater(),
                             ice_meltwater_enthalpy = ZeroHeatContentMeltwater())
    FT = eltype(ocean.model.grid)
    if flux_configuration == :default
        interfaces = ComponentInterfaces(atmosphere, ocean, sea_ice; radiation, land,
                                         ice_freshwater_delivery, ice_meltwater_enthalpy)
        return OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation, land, interfaces)
    end

    velocity_difference_obj = velocity_formulation == :relative ? RelativeVelocity() :
                              velocity_formulation == :wind     ? WindVelocity()     :
                              error("Unknown velocity_formulation: $velocity_formulation. Options: :relative, :wind")

    if flux_configuration == :corrected
        interfaces = ComponentInterfaces(atmosphere, ocean, sea_ice;
                                         radiation,
                                         land,
                                         atmosphere_ocean_fluxes   = corrected_atmosphere_ocean_fluxes(FT),
                                         atmosphere_sea_ice_fluxes = corrected_atmosphere_sea_ice_fluxes(FT; momentum_roughness_length = sea_ice_momentum_roughness_length),
                                         sea_ice_ocean_heat_flux   = corrected_ice_ocean_heat_flux(; heat_transfer_coefficient = sea_ice_ocean_heat_transfer_coefficient),
                                         ice_freshwater_delivery,
                                         ice_meltwater_enthalpy,
                                         atmosphere_ocean_velocity_difference   = velocity_difference_obj,
                                         atmosphere_sea_ice_velocity_difference = velocity_difference_obj)
    elseif flux_configuration == :ncar
        interfaces = ComponentInterfaces(atmosphere, ocean, sea_ice;
                                         radiation,
                                         land,
                                         atmosphere_ocean_fluxes   = ncar_atmosphere_ocean_fluxes(FT),
                                         atmosphere_sea_ice_fluxes = ncar_atmosphere_sea_ice_fluxes(FT),
                                         sea_ice_ocean_heat_flux   = corrected_ice_ocean_heat_flux(; heat_transfer_coefficient = sea_ice_ocean_heat_transfer_coefficient),
                                         ice_freshwater_delivery,
                                         ice_meltwater_enthalpy,
                                         atmosphere_ocean_velocity_difference   = velocity_difference_obj,
                                         atmosphere_sea_ice_velocity_difference = velocity_difference_obj)
    else
        error("Unknown flux_configuration: $flux_configuration. Options: :default, :corrected, :ncar")
    end

    return OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation, land, interfaces)
end

#####
##### Conservative salinity restoring
#####
#
# The surface-salinity restoring must inject zero net salt globally (standard OMIP practice,
# e.g. NorESM/BLOM). `ConservativeSurfaceFluxRestoring` wraps a bare `SurfaceFluxRestoring` and
# stores a pre-corrected, zero-wet-mean flux (`corrected_flux = raw - ⟨raw⟩`, immersed cells
# excluded from the mean) that the salinity top BC reads directly via `getbc`. Because the stored
# field is always its own zero-mean field, the applied restoring integrates to zero over the wet
# ocean exactly, independently of when `update_restoring_flux!` refreshes it from the model state.

# Materializes the restoring flux into a field the boundary condition reads, which is what lets the
# open-water weight and the zero-mean correction be applied at all. The two are independent: either can
# be active alone.
struct ConservativeSurfaceFluxRestoring{R, F, W, M, N} <: Function
    flux            :: R   # wrapped raw surface-flux restoring (getbc-compatible)
    corrected_flux  :: F   # 2D field storing the applied flux; read by the boundary condition
    open_water      :: W   # 2D weight the restoring acts through: 1 everywhere, or 1 - ℵ under ice
    mean_flux       :: M   # host-side scratch reductions
    mean_open_water :: N
    normalize       :: Bool
end

function ConservativeSurfaceFluxRestoring(flux, grid; normalize = true)
    corrected_flux = Field{Center, Center, Nothing}(grid)
    open_water = Field{Center, Center, Nothing}(grid)
    return ConservativeSurfaceFluxRestoring(flux, corrected_flux, open_water,
                                            Field(Average(corrected_flux, dims=(1, 2))),
                                            Field(Average(open_water, dims=(1, 2))),
                                            normalize)
end

# The boundary condition reads only the pre-corrected field; the weight and the reductions over it are
# host-side scratch and are dropped on the device.
Adapt.adapt_structure(to, sf::ConservativeSurfaceFluxRestoring) =
    ConservativeSurfaceFluxRestoring(Adapt.adapt(to, sf.flux),
                                     Adapt.adapt(to, sf.corrected_flux),
                                     Adapt.adapt(to, sf.open_water),
                                     nothing, nothing, sf.normalize)

@inline Oceananigans.BoundaryConditions.getbc(sf::ConservativeSurfaceFluxRestoring, i, j, grid, clock, fields) =
    @inbounds sf.corrected_flux[i, j, 1]

@inline open_water_fraction(::Nothing, i, j, grid) = one(grid)
@inline open_water_fraction(ℵ, i, j, grid) = @inbounds one(grid) - ℵ[i, j, 1]

@kernel function _materialize_surface_flux!(buffer, weight, flux, grid, clock, fields, ice_concentration)
    i, j = @index(Global, NTuple)
    w = open_water_fraction(ice_concentration, i, j, grid)
    @inbounds weight[i, j, 1] = w
    @inbounds buffer[i, j, 1] = w * getbc(flux, i, j, grid, clock, fields)
end

# Refresh the stored restoring flux from the current model state: materialize the wrapped raw flux,
# weight it by the open-water fraction, and remove a multiple of that same weight so the applied flux
# both integrates to zero over the wet ocean and vanishes wherever the weight does. Passing
# `ice_concentration = nothing` leaves the weight at one and recovers `raw - ⟨raw⟩` exactly.
function update_restoring_flux!(sf::ConservativeSurfaceFluxRestoring, model, ice_concentration = nothing)
    grid   = model.grid
    arch   = architecture(grid)
    fields = merge(model.velocities, model.tracers)

    launch!(arch, grid, :xy, _materialize_surface_flux!, sf.corrected_flux, sf.open_water,
            sf.flux, grid, model.clock, fields, ice_concentration)

    if sf.normalize
        compute!(sf.mean_flux)
        compute!(sf.mean_open_water)
        interior(sf.corrected_flux) .-= interior(sf.mean_flux) ./ interior(sf.mean_open_water) .* interior(sf.open_water)
    end

    return nothing
end

# Pull the `ConservativeSurfaceFluxRestoring` off the salinity top BC, or `nothing` when the
# restoring is not conservative-corrected (a bare `SurfaceFluxRestoring` or a plain field).
conservative_restoring(bc::DiscreteBoundaryFunction)         = conservative_restoring(bc.func)
conservative_restoring(mf::MultipleFluxes)                   = conservative_restoring(mf.additional_fluxes)
conservative_restoring(fe::FreshwaterExchange)               = conservative_restoring(fe.additional)
conservative_restoring(sf::ConservativeSurfaceFluxRestoring) = sf
conservative_restoring(other)                                = nothing

struct RefreshSalinityRestoring{R}
    restoring          :: R
    mask_under_sea_ice :: Bool
end

@inline restoring_ice_concentration(::Nothing) = nothing
@inline restoring_ice_concentration(sea_ice) = sea_ice.model.ice_concentration

function (r::RefreshSalinityRestoring)(sim)
    ℵ = r.mask_under_sea_ice ? restoring_ice_concentration(sim.model.sea_ice) : nothing
    update_restoring_flux!(r.restoring, sim.model.ocean.model, ℵ)
    return nothing
end

#####
##### Global freshwater-flux normalization
#####
#
# A forced configuration cannot close its surface water budget: precipitation and runoff are prescribed
# while evaporation is computed from the model state, so the global mean net freshwater flux is non-zero
# and the ocean drifts in volume. Most OMIP-2 models remove that mean each step (Tsujino et al. 2020);
# ICON-O normalizes the global mean sea level daily for the same reason. Only the atmospheric part is
# corrected: the sea-ice exchange moves water between reservoirs rather than adding it, so normalizing it
# would break the ocean + sea-ice water budget. Its mean is recovered as ⟨Jʷᵃᵒ⟩ = ⟨Jʷ⟩ - ⟨Jʷⁱᵒ⟩, since the
# assembled flux carries both. Subtracting in place also fixes the salinity dilution, because
# `carried_tracer_flux` reads the same field.

# The salinity content flux is a `ZeroField` (incoming freshwater carries no salt), so salt follows the
# correction through `Sᴺ · carrying_flux` on its own. The heat content flux `Jᴴ = Tᵒᶜ Jʷᵃᵒ` does not: since
# the salinity/temperature top BCs evaluate `cᴺ · carrying_flux - content_flux`, correcting only the
# carrying flux would leave a spurious surface heat flux `-ρ cₚ Tᵒᶜ ⟨Jʷᵃᵒ⟩`. Removing `Tᵒᶜ ⟨Jʷᵃᵒ⟩` from the
# content flux as well leaves the carried heat flux unchanged: the water removed leaves at the local
# surface temperature, so closing the water budget stays heat-neutral.

# `:none` disables the correction, `:timestep` removes the instantaneous global mean, and `:annual`
# removes a running mean relaxed over a year. `Bool` is accepted for `:none` / `:timestep`.
freshwater_normalization(mode::Symbol) = mode
freshwater_normalization(mode::Bool) = ifelse(mode, :timestep, :none)

function freshwater_averaging_timescale(mode)
    normalization = freshwater_normalization(mode)
    normalization ∈ (:none, :timestep, :annual) ||
        throw(ArgumentError("normalize_freshwater must be :none, :timestep or :annual, got $normalization"))
    return normalization === :annual ? 365.25days : nothing
end

"""
    struct NormalizeTotalWater{S, O, I, N, A, FT}

Removes the drift in total ocean + sea-ice + snow water by adjusting the free surface, rescaling
tracer concentrations by the same factor so their content is unchanged. The reductions that measure
the total are materialized once here and recomputed in place, so stepping allocates nothing. The
reference, total and correction stay 1×1×1 arrays that are only ever combined by broadcasting:
reading one as a scalar would index a GPU array from the host.
"""
struct NormalizeTotalWater{S, O, I, N, A, FT}
    surface_height  :: S     # workspace holding `η` on the ocean grid (see `compute_total_water!`)
    ocean_water     :: O     # ∫η dA
    sea_ice_water   :: I     # ∫hⁱ ℵ dA, or `nothing` without sea ice
    snow_water      :: N     # ∫hˢ ℵ dA, or `nothing` without snow
    reference_water :: A     # total water to hold to, in freshwater-equivalent volume
    total_water     :: A     # workspace for the current total
    correction      :: A     # workspace for δ, the sea-level adjustment
    density_ratios  :: NTuple{2, FT}   # ρⁱ/ρᵒ and ρˢ/ρᵒ
    wet_area        :: FT    # ∫dA over wet columns, the area the correction is spread over
    relaxation      :: FT    # fraction of the excess removed each step
end

function NormalizeTotalWater(coupled_model, Δt, averaging_timescale)
    grid = coupled_model.ocean.model.grid
    FT = eltype(grid)

    unit_area = Field{Center, Center, Nothing}(grid)
    set!(unit_area, 1)
    wet_area = Array(interior(compute!(Field(Integral(unit_area)))))[1]

    surface_height = Field{Center, Center, Nothing}(grid)
    ocean_water = Field(Integral(surface_height))

    sea_ice = coupled_model.sea_ice
    ρᵒ = coupled_model.interfaces.ocean_properties.reference_density

    sea_ice_water, snow_water, density_ratios = if isnothing(sea_ice)
        nothing, nothing, (zero(FT), zero(FT))
    else
        si = sea_ice.model
        hˢ = si.snow_thickness
        (Field(Integral(si.ice_thickness * si.ice_concentration)),
         isnothing(hˢ) ? nothing : Field(Integral(hˢ * si.ice_concentration)),
         (convert(FT, uniform_value(si.sea_ice_density) / ρᵒ),
          convert(FT, uniform_value(si.snow_density) / ρᵒ)))
    end

    reference_water, total_water, correction = (similar(interior(ocean_water)) for _ in 1:3)
    relaxation = isnothing(averaging_timescale) ? one(FT) : convert(FT, Δt / averaging_timescale)

    normalization = NormalizeTotalWater(surface_height, ocean_water, sea_ice_water, snow_water,
                                        reference_water, total_water, correction, density_ratios,
                                        convert(FT, wet_area), relaxation)

    compute_total_water!(normalization, coupled_model)
    reference_water .= total_water

    return normalization
end

# Remove the volume drift by adjusting the free surface. Avoids problems with spurious convection 
# appearing where it should not when using a surface restoring.
@kernel function _correct_sea_level!(η, T, S, grid, δ, k_top)
    i, j = @index(Global, NTuple)
    h = Oceananigans.Grids.static_column_depthᶜᶜᵃ(i, j, grid)
    @inbounds ηⁿ = η[i, j, k_top]
    H = h + ηⁿ
    @inbounds d = δ[1, 1, 1]
    correctable = (h > 0) & (H - d > 0)
    f = ifelse(correctable, H / (H - d), one(grid))

    for k in 1:size(grid, 3)
        @inbounds T[i, j, k] *= f
        @inbounds S[i, j, k] *= f
    end

    @inbounds η[i, j, k_top] = ηⁿ - ifelse(correctable, d, zero(grid))
end

function correct_sea_level!(ocean_model, δ)
    grid = ocean_model.grid
    Nx, Ny, Nz = size(grid)
    η = ocean_model.free_surface.displacement

    launch!(architecture(grid), grid, Oceananigans.Utils.KernelParameters(1:Nx, 1:Ny), _correct_sea_level!,
            η, ocean_model.tracers.T, ocean_model.tracers.S, grid, δ, Nz + 1)

    # Refresh `σⁿ` from the corrected `η`.
    launch!(architecture(grid), grid, Oceananigans.Models.surface_kernel_parameters(grid),
            Oceananigans.Models.HydrostaticFreeSurfaceModels._update_zstar_scaling!, η, grid)

    return nothing
end

@inline uniform_value(ρ) = ρ
@inline uniform_value(ρ::Oceananigans.Fields.ConstantField) = ρ.constant

"""
    compute_total_water!(n::NormalizeTotalWater, coupled_model)

Accumulate into `n.total_water` the water held by the ocean, sea ice and snow in freshwater-equivalent
volume, up to the fixed reference volume `V₀` that cancels when differences are taken. `∫η dA` is a
reduction over wet columns rather than over every cell, which keeps the summation roundoff orders of
magnitude below the drift being removed. Recomputes `n`'s reductions in place and allocates nothing.
"""
function compute_total_water!(n::NormalizeTotalWater, coupled_model)
    # `Integral` over the free surface's own grid over-counts: with a `SplitExplicitFreeSurface` the
    # displacement lives on a grid with extended halos, and reducing there returns an area 1.5× the
    # wet area. Copying onto the ocean grid, where the same reduction gives ∫dA correctly, costs one
    # 2-D copy per step against three reductions of the same size.
    interior(n.surface_height) .= interior(coupled_model.ocean.model.free_surface.displacement)
    n.total_water .= interior(compute!(n.ocean_water))

    ρⁱ, ρˢ = n.density_ratios
    isnothing(n.sea_ice_water) || (n.total_water .+= ρⁱ .* interior(compute!(n.sea_ice_water)))
    isnothing(n.snow_water)    || (n.total_water .+= ρˢ .* interior(compute!(n.snow_water)))

    return n.total_water
end

# Checkpoint the reference state
Oceananigans.Simulations.callback_state(n::NormalizeTotalWater) =
    (; reference_water = Array(n.reference_water)[1])

function Oceananigans.Simulations.restore_callback_state!(n::NormalizeTotalWater, state)
    n.reference_water .= state.reference_water
    return n
end

function (n::NormalizeTotalWater)(sim)
    compute_total_water!(n, sim.model)
    n.correction .= n.relaxation .* (n.total_water .- n.reference_water) ./ n.wet_area
    correct_sea_level!(sim.model.ocean.model, n.correction)
    return nothing
end

#####
##### Main simulation builder
#####

# Every stage of `omip_simulation` below is collective, and at eddying resolution each one takes
# minutes, so a rank that stalls in a mismatched collective is invisible: the job simply stops
# producing output until the scheduler kills it. Tag each completed stage with the rank and the
# elapsed time, so a hang is localized both to a stage and to the ranks that never reached it.
function log_setup_stage(arch, stage, t₀)
    rank = arch isa Distributed ? arch.local_rank : 0
    @info @sprintf("omip_simulation setup [rank %d] %-28s %8.1f s", rank, stage, time() - t₀)
    flush(stderr)
    return nothing
end

"""
    omip_simulation(config::Symbol = :halfdegree; kwargs...)

Create a fully coupled ocean--sea-ice--atmosphere OMIP simulation.

The single positional argument selects the grid configuration:

- `:halfdegree`    -- 720x360   `TripolarGrid`
- `:quarterdegree` -- NEMO eORCA025 (1/4ᵒ) mesh
- `:twelfthdegree` -- NEMO eORCA12 (1/12ᵒ) mesh
- `:orca`          -- NEMO eORCA1 (~1ᵒ) mesh
- `:test`          -- NEMO eORCA1 (~1ᵒ) mesh, locally-runnable preset for reproducing the
                      quarter-degree spurious high-latitude ice + surface salinity drift.
                      Overrides `Nz = 15`, `Δz_top = 1.5` m, `Δt = 45minutes`, and a `10days`
                      biharmonic-viscosity timescale. GM/Redi is disabled (`κ_skew = κ_symmetric =
                      nothing`) to keep short test runs cheap; momentum advection follows `:orca`.

Returns a `Simulation` wrapping an `OceanSeaIceModel`. The simulation
already has a progress callback attached, and (when `diagnostics=true`)
the OMIP-protocol output writers from [`add_omip_diagnostics!`](@ref).

To restart from a previous run, simply call

    run!(sim; pickup = true)

which uses Oceananigans' built-in `Checkpointer` machinery — no extra
plumbing is needed because `NumericalEarth.EarthSystemModels` provides
`prognostic_state` / `restore_prognostic_state!` for the coupled model.

# Keyword arguments

- `arch`: architecture (`CPU()` or `GPU()`). Default: `CPU()`.
- `Nz::Int`: number of vertical levels. Per-config default: `15` for `:test`, `100` otherwise.
- `depth`: maximum ocean depth in metres. Default: `5500`.
- `Δz_top`: target surface-cell thickness in metres (sets the exponential vertical scale). Per-config
  default: `1.5` for `:quarterdegree`/`:twelfthdegree`/`:test`, `nothing` (scale derived from
  `depth`/`Nz`) otherwise.
- `κ_skew`, `κ_symmetric`: GM/Redi diffusivities. Per-config defaults: `nothing` (no isopycnal
  diffusivity) for the eddy-resolving `:quarterdegree`/`:twelfthdegree` and for `:test`, `800` for
  `:orca`, `250` for `:halfdegree`. Either may instead be `:nemo`, which selects NEMO's Treguier
  et al. (1997) coefficient (`nn_aei_ijk_t = nn_aht_ijk_t = 21`, the setting CMCC's ORCA1 uses):
  the internal Rossby radius squared times the baroclinic growth rate, recomputed every step and
  held depth-uniform. GM tapers to zero equatorward of 20°, Redi rises to its reference value there
  and carries a floor of one fifth of it. See [`NEMOEddyCoefficients`](@ref).
- `isopycnal_formulation`: which closure carries GM/Redi. `:standard` (default) is
  `IsopycnalSkewSymmetricDiffusivity`; `:triad` is `TriadIsopycnalSkewSymmetricDiffusivity`, which
  evaluates the tensor on the Griffies triad stencil. The triad closure carries no
  `skew_flux_formulation`, so `:triad` requires the default `:diffusive`.
- `skew_flux_formulation`: how the GM skew transport is applied. `:diffusive` (default) adds it to
  the tracer flux; `:advective` builds the eddy-induced velocity and advects with it, which also
  makes the bolus transport available as a model field. Those two are equivalent continuously, not
  discretely. `:boundary_value` is a different transport, not a different discretization of the same
  one: the eddy transport in each column solves the boundary-value problem of Ferrari et al. (2010),
  `(c² ∂²/∂z² − N²) Υ = −N² Υᴳᴹ` with `Υ = 0` at the surface and the bottom, which low-passes the
  baroclinic modes, satisfies the boundary conditions without tapering, and interpolates through
  weakly stratified layers with no floor on `N²`. It is applied advectively and requires a
  depth-independent `κ_skew` (a number or `:nemo`). See [`BoundaryValueTransport`](@ref).
  Ignored when `κ_skew` is `nothing`.
- `boundary_value_mode_number`, `boundary_value_minimum_speed`: `M` and `c_min` setting the speed
  `c = max(c_min, (M π)⁻¹ ∫ N dz)` that weights the second-order operator, used only when
  `skew_flux_formulation = :boundary_value`. Defaults: `2` and `0.1` m s⁻¹. Larger `M` filters less
  and gives a larger transport; `M = 1` is the first baroclinic mode, whose amplitude is about half
  the truncated GM transport.
- `biharmonic_timescale`: horizontal biharmonic-viscosity timescale. Per-config default: `nothing`
  (no biharmonic viscosity) for `:quarterdegree`/`:twelfthdegree`, `10days` for `:test`, `50days`
  otherwise.
- `forcing_dir`: directory for JRA55 forcing data. Default: `"forcing_data"`.
- `restoring_dir`: directory for restoring/IC climatology. Default: `"climatology"`.
- `piston_velocity`: surface salinity restoring piston velocity in m/day. Default: `1/6`.
  Restoring is applied uniformly over the ocean surface, including under sea ice unless
  `restoring_under_sea_ice = false`.
- `restoring_under_sea_ice`: whether the surface-salinity restoring acts under sea ice. Default:
  `true`, the OMIP-2 convention. Set `false` to weight it by the open-water fraction `1 - ℵ`, since
  WOA is poorly constrained beneath ice and the restoring there works against the ice--ocean salt
  flux. When `normalize_salinity` is also on, the zero-global-mean correction is spread over the
  open-water weight alone, so the applied flux still injects no net salt while vanishing under ice.
  The two are independent — either may be used without the other.
- `normalize_freshwater`: removal of the drift in total ocean + sea-ice + snow water, held to its
  initial value by lowering or raising the free surface. Tracer concentrations are rescaled by the
  same factor within each column, so heat and salt content are unchanged and the stratification is
  scaled rather than shifted — a surface-flux correction instead lands on the top cell alone, which
  crosses convection thresholds in marginally-stable columns. Because the controller measures the
  total rather than predicting it from a flux, it is self-correcting, ignores ocean--sea-ice exchange
  (which leaves the total untouched), and still catches snowfall intercepted by the ice and
  sublimation off it. Options:
    * `:none` (or `false`, the default): no correction.
    * `:timestep` (or `true`): remove the whole excess each step. Pins the total exactly, but also
      removes the physical seasonal cycle of land water storage, which is comparable in size to the
      secular imbalance.
    * `:annual`: remove `Δt / 365.25days` of the excess each step, which removes the drift while
      leaving the seasonal cycle in place.
- `river_spread_radius`: radius in degrees over which each river mouth's discharge is divided equally
  among the wet cells around its landing cell. A geographic footprint keeps the freshwater flux per
  unit area put as the grid refines; a fixed cell count concentrates it by the square of the
  resolution ratio and drives coastal salinity to zero. Per-config default: `1.2` for the refined
  grids, `nothing` for `:orca`/`:test`, which pins them to the historical cell-count footprint so
  their existing integrations stay reproducible.
- `river_spread_cells`: number of cells in that footprint, nearest first — a cap when
  `river_spread_radius` is set, and the footprint itself when it is `nothing`. Per-config default:
  `nothing` (uncapped) for the refined grids, `8` for `:orca`/`:test`.
- `atlantic_runoff_diversion`: fraction of the river and iceberg discharge landing in the Atlantic
  that is delivered to the Pacific instead, at the latitude it was diverted from. Conserves the
  global freshwater input. Default: 0.
- `initial_condition_blend_depth`: depth in metres above which the initial T and S come from WOA
  **Monthly** for the month of `start_date` rather than WOA Annual, tapering linearly to the annual
  field at the monthly climatology's own reach (~1500 m). WOA Annual averages a seasonal cycle whose
  winter half has kilometre-deep mixed layers, so a January start initialized from it begins with a
  seasonal thermocline the season does not have. `nothing` (default) uses WOA Annual throughout and
  reproduces the previous model exactly.
- `northern_sea_ice_initial_date`, `southern_sea_ice_initial_date`: the ECCO4Monthly dates the initial
  sea-ice thickness and concentration are taken from, north and south of the equator. The two
  hemispheres reach their minimum six months apart, so the single date both default to,
  `DateTime(1993, 1, 1)`, starts the Arctic at its seasonal maximum while the Antarctic is at its
  minimum; `DateTime(1993, 9, 1)` in the north gives the summer pack in both. Sea ice never reaches
  the equator, so the two fields are stitched there with no taper. Defaults: `DateTime(1993, 1, 1)`
  for both, which reproduces the previous model exactly.
- `sea_ice_ocean_drag_reference_depth`: depth in metres over which the ocean velocity is averaged to
  give the reference of the ice-ocean drag. McPhee's `Cᵢₒ = 5.5e-3` is defined against the under-ice
  boundary layer, tens of metres thick, while the topmost cell is 1.5 m and is dragged along by the ice
  itself, so referencing the drag there under-brakes the pack. `nothing` uses the topmost cell and
  reproduces the previous model exactly. Default: `nothing`.
- `ice_salinity`: bulk salinity of the sea ice in psu, carried as a `ConstantField`. It sets the salt
  returned per unit of melt, `Jˢ = Eᵢ Sˢⁱ / ρᵒᶜ`, so the freshwater a melting cell delivers goes as
  `(Sᴺ - Sˢⁱ)/Sᴺ`: 0.885 at 4 psu against 0.828 at 6, for `Sᴺ = 34.9`. Multi-year Arctic ice is 2-4 psu
  and first-year ice 5-8. Default: `4`.
- `ice_freshwater_fraction`: the fraction of the sea ice-ocean mass exchange delivered to the ocean,
  volume and salt alike. The withheld water leaves the ocean + ice + snow total and
  `normalize_freshwater` returns it globally through the free surface, so the global budget closes
  while the local delivery is scaled. Default: `1`, the full exchange.
- `ice_melt_mixing`, `ice_melt_mixing_κ`, `ice_melt_mixing_depth`, `ice_melt_mixing_threshold`: extra
  vertical tracer diffusivity over the top `ice_melt_mixing_depth` metres wherever the sea ice is
  melting into the ocean faster than `ice_melt_mixing_threshold` (m s⁻¹). The ice-ocean exchange is
  delivered to the surface cell, so a melt event leaves a lid one cell thick that the vertical closure
  must erode; in reality it is stirred through the ice-ocean boundary layer under keels of 1–3 m. Same
  device as `river_mixing`, and unlike it the footprint follows the melt from step to step. ⚠ `κ` is
  *not* the river value — 0.1 m² s⁻¹ would mix ≈131 m in a day, which is convective adjustment, while
  5e-4 gives ≈9 m day⁻¹. Defaults: `false`, `5e-4` m² s⁻¹, `10` m, `1e-9` m s⁻¹.
- `under_ice_viscosity`, `under_ice_viscosity_depth`: floor on the vertical viscosity over the top
  `under_ice_viscosity_depth` metres, scaled by the ice concentration, the closure otherwise untouched.
  Under ice CATKE's momentum diffusivity falls to molecular values below ≈18 m, so the ice stress reaches
  an Ekman layer of 3–16 m against 11–196 m in open water; this bounds how deep it reaches. McPhee's u★ℓ
  scaling with u★ ≈ 1 cm s⁻¹ and ℓ ≈ 1 m gives 1e-2 m² s⁻¹. Defaults: `nothing` (off), `20` m.
- `ice_arch_region`, `ice_arch_stress`, `ice_arch_months`: a seasonal ice arch over the longitude–latitude
  box `ice_arch_region = (λ₁, λ₂, φ₁, φ₂)`, as a basal stress of `ice_arch_stress` N m⁻² arresting the ice
  from the first to the last month of `ice_arch_months` inclusive, wrapping the year. `nothing` switches it
  off. Two boxes are in use. Nares Strait `(-78, -58, 77.5, 82.5)` over `(12, 7)` is the observed arch season
  (Kwok 2005, 2010), which holds the annual export to ≈ 130 km³ against the 450 the model exports without it.
  Baffin Bay and Davis Strait `(-80, -48, 65, 80)` over `(1, 12)` arrests the Davis delivery outright: the
  ψ response is linear in the removed Davis export and blind to Fram (C21-3, C21-5), and no run has yet
  varied Davis while holding the global ice state fixed. ⚠ Davis Strait is far too deep for a basal stress
  to be physical there — read the second box as a mechanism probe, like `with_ice_dynamics = false`, not as
  a tuning. Defaults: `nothing`, `100`, `(12, 7)`.
- `ice_meltwater_at_interface_temperature`: deliver ice meltwater at the interface temperature `Tᵦ`
  rather than at Conservative Temperature 0. `ClimaSeaIce` produces basal meltwater at `Tᵦ`, so the
  default hands the ocean roughly 0.23 W m⁻² per m yr⁻¹ of basal melt that it never paid for.
  ⚠ The correction is applied to the whole ice mass flux, over-correcting the top-melt fraction, which
  is produced near 0; read it as an upper bound. Default: `false`.
- `ice_virtual_salt_flux`: deliver that exchange as a salt flux at fixed ocean volume — the classical
  virtual salt flux `Jˢ = Jʷ (Sᴺ − Sˢⁱ)` — instead of as a real volume flux, which isolates the volume
  pathway from the freshwater amount. Exact only for `Sᴺ` uniform over the column, and it does not
  conserve total salt. Overrides `ice_freshwater_fraction`. Default: `false`.
- `river_mixing`, `river_mixing_κ`, `river_mixing_depth`: extra vertical tracer diffusivity applied over
  the whole spread footprint (cf. NEMO `rn_avt_rnf` over `rn_hrnf`). Defaults: `true`, `0.1` m² s⁻¹,
  `10` m.
- `start_date`, `end_date`: bracket for forcing/restoring metadata. Defaults: 1958-01-01 .. 2018-01-01.
- `Δt`: simulation time step. Per-config default: `5minutes` for `:twelfthdegree`, `20minutes` for
  `:quarterdegree`, `30minutes` otherwise.
- `stop_time`: stop time for the wrapping `Simulation`. Default: `Inf`.
- `thickness_categories`: number of equal-area sub-grid ice thickness categories used for the
  effective conductivity of [Fichefet and Morales Maqueda (1997)](@cite fichefet1997sensitivity).
  Conducting through the cell-mean thickness underestimates growth because the conductive flux goes
  as ``1/h``; placing `N` categories at ``(2i-1) h / N`` multiplies the conductivity by
  ``\\sum_i 1/(2i-1)`` — 1.53 for `N = 3`, 1.79 for `N = 5`. Default: `1`, which conducts through
  the mean and applies no enhancement.
- `implicit_bottom_drag::Bool`: if `true` (default), the bottom and immersed quadratic drag are affine
  fluxes `J = λ φᵦ` with `λ = -μ |u|` in the vertical solver's diagonal. `false` applies both
  explicitly, the treatment this replaced, so the pair can be A/B'd. The two differ most in shallow
  fast cells such as narrow straits, where the explicit stability number `μ |u| Δt / Δz` is order ten.
- `bottom_drag_background_velocity`: unresolved velocity `uᵦ` added in quadrature to the resolved speed
  in the quadratic bottom drag, `τ = μ u √(uᵦ² + |u|²)`, standing for tides and other motions the grid
  does not carry. It only bites where the resolved flow is weak: at `|u| = 1 cm/s` and `uᵦ = 10 cm/s`
  the stress is ten times the resolved one, at `50 cm/s` within a percent of it. Default `0`; GFDL's
  OM4 uses `0.1`.
- `barotropic_substeps`: number of split-explicit substeps the free surface takes per `Δt`. The
  barotropic gravity wave must stay inside a substep, so a refined grid or a longer `Δt` needs more;
  too few blows the free surface up on the first step. Per-config default: `200` for
  `:quarterdegree`/`:twelfthdegree`, `100` otherwise. A warning names the count the grid needs.
- `flux_configuration`: surface flux formulation. Options:
   * `:default` — current defaults (Edson/COARE with constant Charnock 0.02)
   * `:corrected` — COARE 3.6 with wind-dependent Charnock, fixed ice roughness, momentum-based u*
   * `:ncar` — OMIP-2 standard Large & Yeager (2004) bulk formulae
- `sea_ice_momentum_roughness_length`: aerodynamic roughness z₀ of the ice surface, m, used by
  `:corrected`. 5e-4 is the SHEBA multiyear-pack value; smooth first-year ice is nearer 1e-4, which cuts
  the neutral drag coefficient by about a quarter and the free-drift speed by about a seventh.
- `vertical_closure::Symbol`: ocean vertical-mixing closure. Options:
   * `:catke` — CATKE TKE-based scheme (default).
   * `:simple` — `ConvectiveAdjustmentVerticalDiffusivity(convective_κz=1)` plus a
     depth-step background `VerticalScalarDiffusivity` (κ=10⁻², ν=10⁻² in upper
     100 m; κ=10⁻⁵, ν=10⁻⁴ below). For diagnostic A/B tests vs CATKE.
   * `:nori` — NORi base Richardson-number closure
     (xkykai/NORiOceanParameterization.jl, vendored as
     `nori_base_closure.jl`). Calibrated defaults; no `Cᵇ` parameter.
   * `:rbvd` — Oceananigans' built-in `RiBasedVerticalDiffusivity`
     (Richardson-number-based, with κ-clip and time-averaged smoothing).
     A battle-tested comparison point for `:nori`; no `Cᵇ` parameter.
   * `:kpp` — KPP boundary-layer scheme (Large 1994 / MITgcm), vendored
     in `KPP/`. Includes nonlocal tracer flux + SW-aware Bf. No `Cᵇ`.
   * `:nemo_tke` — NEMO 3.6 TKE scheme (Blanke & Delecluse 1993; Gaspar et al.
     1990; Madec et al. 2017), vendored in `NEMOTKE/`. OMIP-2 ORCAOne preset:
     prognostic e, gradient-limited length scale, Langmuir + Mellor-Blumberg
     wave penetration + EVD on static instability. No `Cᵇ`.
- `background_vertical_diffusivity`: interior background tracer diffusivity κ added underneath the
  primary vertical closure (`:catke` and `:rbvd` only — the other closures set their own interior
  background internally and reject this keyword). Options:
   * `:henyey` (default) — the latitude-dependent internal-wave scaling of Henyey et al. (1986),
     κ = max(2×10⁻⁶, 10⁻⁵ |sin φ|), i.e. 2×10⁻⁶ m² s⁻¹ at the equator rising to 10⁻⁵ m² s⁻¹ at the poles.
   * `:bryan_lewis` — the Bryan & Lewis (1979) depth profile,
     κ = 0.8×10⁻⁴ + (1.05×10⁻⁴/π) atan[4.5×10⁻³ (|z| − 2500)] m² s⁻¹, i.e. 3×10⁻⁵ in the upper ocean
     rising to 1.3×10⁻⁴ in the abyss. Buys the deep upwelling without diffusing the thermocline the
     way a uniform 10⁻⁴ does.
   * `:abyssal_henyey` — `:henyey` in the thermocline plus an arctangent enhancement reaching
     5×10⁻⁵ m² s⁻¹ by 5000 m, so κ ≈ 2×10⁻⁶ at the equatorial surface and ≈ 5×10⁻⁵ in the abyss.
     The upper-ocean value sets the global heat uptake and the abyssal value the deep ventilation;
     this holds the former at Henyey and raises only the latter.
   * a number — a uniform κ in m² s⁻¹ (e.g. `3e-5`, `1e-4`). The background diffusivity sets the
     diapycnal upwelling that closes the AMOC's lower limb, so raising it strengthens the
     overturning (Bryan 1987 gives AMOC ∝ κ^(2/3) in the diffusive limit); it also deepens the
     thermocline, so watch the tropical SST and mixed-layer depth alongside the AMOC.
- `background_vertical_viscosity`: the matching background momentum viscosity ν in m² s⁻¹, subject to
  the same per-closure restriction. `nothing` (default) uses 10⁻⁴ m² s⁻¹ for both `:catke` and
  `:rbvd`. Set it alongside the diffusivity to hold the background Prandtl number fixed while
  scanning κ.
- `implicit_vertical_advection::Bool`: if `true` (default), tracer and momentum vertical advection use
  `AdaptiveVerticallyImplicitDiscretization(cfl=0.5)` (switches the vertical advective flux to implicit
  where the vertical Courant number is large — e.g. in thin near-surface cells). If `false`, fully
  explicit `WENO`/`WENOVectorInvariant`. Use `false` to isolate adaptive-implicit advection effects.
- `boundary_scheme::Symbol`: the reconstruction the WENO buffer chain terminates in, used in the one cell
  whose stencil no longer fits — the domain buffer and, on an `ImmersedBoundaryGrid`, any cell adjacent to
  an inactive node. Options:
   * `:default` — `Centered(order=2)` for tracers and `UpwindBiased(order=1)` for momentum.
   * `:upwind` — first-order upwind, monotone, in exactly those cells while the interior keeps the full order.
   * `:cwenoz` — the third-order central-WENO reconstruction of Semplice, Travaglia and Puppo (2022),
     whose stencil extends only inwards.
- `tracer_boundary_scheme::Symbol`, `momentum_boundary_scheme::Symbol`: the same choice made separately for the
  tracer and the momentum reconstructions, both defaulting to `boundary_scheme`. Setting one of them alone
  isolates which of the two the boundary treatment acts through.
- `temperature_reference_variation`, `salinity_reference_variation`, `momentum_reference_variation`:
  the cell-to-cell variation of the reconstructed field below which CWENOZ reads the stencil as noise, setting
  `ϵ = reference_variation²`. It carries the units of the field and no grid spacing, so one value serves every
  direction and every cell thickness. The constant candidate takes over above a variation of about
  `8.6 × reference_variation` between adjacent averages. All three default to `0`, which reads `ϵ` off the
  stencil: a nonzero value acts only on the horizontal, where it is worth some tens of m²/s against a skew
  diffusivity of order 1e3, and is inert on the vertical, where the boundary reconstruction competes against a
  background diffusivity of order 1e-5. The horizontal momentum terms reconstruct a vorticity, a divergence flux
  and a squared velocity, so they take no variation; `momentum_reference_variation` is a speed in m/s and applies
  to the vertical reconstruction alone.
- `velocity_formulation::Symbol`: Δu used by the bulk formula. Options:
   * `:relative` — `Δu = u_atm − u_ocean` (OMIP-2 α=1, default).
   * `:wind` — `Δu = u_atm` (ignores ocean current). For isolating bulk-formula
     response from current feedback (e.g. when an over-strong ACC self-reinforces).
- `diagnostics::Bool`: whether to attach OMIP diagnostics. Default: `true`.
- `surface_averaging_interval`, `field_averaging_interval`: averaging windows.
- `checkpoint_interval`: interval between checkpoint writes.
- `output_dir`, `filename_prefix`, `file_splitting_interval`: output configuration.
"""
function omip_simulation(config::Symbol = :halfdegree;
                         arch = CPU(),
                         Nz = ConfigDefault(),
                         depth = 5500,
                         Δz_top = ConfigDefault(),
                         κ_skew = ConfigDefault(),
                         κ_symmetric = ConfigDefault(),
                         skew_flux_formulation = :diffusive,
                         isopycnal_formulation = :standard,
                         Cᵇ = 0.28,
                         Cᵘⁿᵇ = 0.0,
                         Cᶠ = 1.0,
                         Cᶠ⁰ = 1e9,
                         Cᶠᵟ = 0.75,
                         Cᵉc = 0.112,
                         biharmonic_timescale = ConfigDefault(),
                         biharmonic_viscosity = nothing,
                         forcing_dir = joinpath(get(ENV, "DATA", ""), "forcing_data"),
                         staging_dir = nothing,
                         backend_size = 50,
                         restoring_dir = "climatology",
                         piston_velocity = 1 / 6, # m / day
                         initial_condition_blend_depth = nothing,
                         start_date = DateTime(1958, 1, 1),
                         end_date = DateTime(2018, 1, 1),
                         Δt = ConfigDefault(),
                         stop_time = Inf,
                         flux_configuration = :default,
                         vertical_closure = :catke,
                         boundary_value_mode_number = 2,
                         boundary_value_minimum_speed = 0.1,
                         background_vertical_diffusivity = :henyey,
                         background_vertical_viscosity = nothing,
                         implicit_vertical_advection = true,
                         tracer_advection_order = 7,
                         boundary_scheme = :default,
                         tracer_boundary_scheme = boundary_scheme,
                         momentum_boundary_scheme = boundary_scheme,
                         temperature_reference_variation = 0,
                         salinity_reference_variation = 0,
                         momentum_reference_variation = 0,
                         implicit_bottom_drag = true,
                         bottom_drag_background_velocity = 0,
                         velocity_formulation = :relative,
                         Cᵂu★ = nothing,
                         with_snow = false,
                         with_ice_dynamics = true,
                         with_ocean_surface_tilt = false,
                         with_landfast_basal_stress = true,
                         sea_ice_ocean_heat_transfer_coefficient = 0.0057,
                         sea_ice_momentum_roughness_length = 5e-4,
                         sea_ice_lateral_boundary_condition = :no_slip,
                         sea_ice_ocean_drag_coefficient = 5.5e-3,
                         sea_ice_ocean_drag_reference_depth = 6,
                         ice_compressive_strength = 27500,
                         ice_salinity = 4,
                         northern_sea_ice_initial_date = DateTime(1993, 1, 1),
                         southern_sea_ice_initial_date = DateTime(1993, 1, 1),
                         partial_cell_bathymetry = false,
                         mixed_layer_tapering = false,
                         normalize_salinity = true,
                         restoring_under_sea_ice = true,
                         normalize_freshwater = false,
                         river_mixing = true,
                         river_mixing_κ = 0.1,
                         river_mixing_depth = 10,
                         river_spread_radius = ConfigDefault(),
                         river_spread_cells = ConfigDefault(),
                         atlantic_runoff_diversion = 0,
                         ice_freshwater_fraction = 1,
                         ice_virtual_salt_flux = false,
                         ice_meltwater_at_interface_temperature = true,
                         sea_ice_liquidus = :teos10,
                         ice_melt_mixing = false,
                         ice_melt_mixing_κ = 5e-4,
                         ice_melt_mixing_depth = 10,
                         ice_melt_mixing_threshold = 1e-9,
                         under_ice_viscosity = nothing,
                         under_ice_viscosity_depth = 20,
                         ice_arch_region = nothing,
                         ice_arch_stress = 100,
                         ice_arch_months = (12, 7),
                         barotropic_substeps = ConfigDefault(),
                         chlorophyll = :seawifs,
                         thickness_categories = 1,
                         snow_thickness_categories = thickness_categories,
                         itd_shape = nothing,
                         bbl_diffusivity = nothing,
                         bbl_transport_coefficient = nothing,
                         overflow_restoring_timescale = nothing,
                         labrador_restoring_timescale = nothing,
                         diagnostics = true,
                         field_mean_interval = 5days,
                         surface_averaging_interval = 5days,
                         field_averaging_interval = 15days,
                         checkpoint_interval = 360days,
                         output_dir = ".",
                         filename_prefix = string(config),
                         file_splitting_interval = 360days)

    cfg = Val(config)

    # Resolve resolution-sensitive parameters to their per-configuration defaults unless the
    # user passed an explicit value (see `config_*` below).
    Nz                   = resolve_config_default(Nz,                   config_Nz(cfg))
    Δz_top               = resolve_config_default(Δz_top,               config_Δz_top(cfg))
    κ_skew               = resolve_config_default(κ_skew,               config_κ_skew(cfg))
    κ_symmetric          = resolve_config_default(κ_symmetric,          config_κ_symmetric(cfg))
    biharmonic_timescale = resolve_config_default(biharmonic_timescale, config_biharmonic_timescale(cfg))
    river_spread_radius  = resolve_config_default(river_spread_radius,  config_river_spread_radius(cfg))
    river_spread_cells   = resolve_config_default(river_spread_cells,   config_river_spread_cells(cfg))
    barotropic_substeps  = resolve_config_default(barotropic_substeps,  config_barotropic_substeps(cfg))
    Δt                   = resolve_config_default(Δt,                   config_Δt(cfg))

    check_depth_independent_skew_coefficient(κ_skew, skew_flux_formulation)
    check_isopycnal_formulation(isopycnal_formulation, skew_flux_formulation)

    setup_t₀ = time()
    log_setup_stage(arch, "start", setup_t₀)

    grid = build_grid(cfg, arch, Nz, depth; Δz_top, partial_cell_bathymetry)
    log_setup_stage(arch, "grid", setup_t₀)

    # When staging_dir is provided, JRA55 data is read from fast scratch
    # with symlink fallback to the slow source directory.
    if !isnothing(staging_dir)
        setup_staging_directory(forcing_dir, staging_dir)
        atmosphere_dir = staging_dir
    else
        atmosphere_dir = forcing_dir
    end

    # Build the land before the ocean so its river routing can seed enhanced vertical mixing at
    # river mouths — an extra closure that keeps concentrated runoff from freshening a cell to zero.
    # Closing the shallow Ob/Yenisei gulfs relocates their mouths ~2–3° onto the deeper Kara Sea shelf,
    # so the routing search must reach that far. A fixed geographic reach keeps it resolution-independent.
    Nx, Ny, _ = size(grid)
    maximum_search_radius = max(5, ceil(Int, 3 / ((360 / Nx + 180 / Ny) / 2)))
    # Basin masks are only needed when discharge is being moved between them, and the flood fill
    # behind them is serial, so build them on request.
    flux_diversion = if atlantic_runoff_diversion > 0
        basin_mask(basin) = Bool.(dropdims(Array(interior(basin.mask)); dims = 3))
        (; fraction = atlantic_runoff_diversion,
           from = basin_mask(atlantic_ocean_basin(grid)),
           to   = basin_mask(pacific_ocean_basin(grid)))
    else
        nothing
    end

    land = JRA55PrescribedLand(grid; dir = atmosphere_dir, dataset = MultiYearJRA55(),
                               start_date, end_date, time_indices_in_memory = backend_size, prefetch = true,
                               maximum_search_radius,
                               spread_radius = river_spread_radius,
                               n_spread_cells = river_spread_cells,
                               flux_diversion)
    log_setup_stage(arch, "land", setup_t₀)

    # Built here because the ocean closure is constructed before the coupled model that owns the real
    # ice-ocean flux; `RefreshIceMeltDiffusivity` fills this field once the simulation exists.
    ice_melt_diffusivity = ice_melt_mixing ? CenterField(grid) : nothing
    ice_melt_mask = ice_melt_mixing ?
        ice_melt_mixing_mask(grid; κ = ice_melt_mixing_κ, mixing_depth = ice_melt_mixing_depth) :
        nothing
    ice_melt_κ_closure = ice_melt_mixing ? ice_melt_vertical_diffusivity(ice_melt_diffusivity) : nothing

    under_ice_ν_field = isnothing(under_ice_viscosity) ? nothing : CenterField(grid)
    under_ice_ν_mask = isnothing(under_ice_viscosity) ? nothing :
        ice_melt_mixing_mask(grid; κ = under_ice_viscosity, mixing_depth = under_ice_viscosity_depth)
    under_ice_ν_closure = isnothing(under_ice_viscosity) ? nothing : under_ice_vertical_viscosity(under_ice_ν_field)

    river_κ = river_mixing ?
        river_mouth_vertical_diffusivity(grid, land.river_routing; κ = river_mixing_κ, mixing_depth = river_mixing_depth) :
        nothing

    nemo_eddy_coefficients = uses_nemo_eddy_coefficients(κ_skew, κ_symmetric) ?
        NEMOEddyCoefficients(grid) : nothing

    cesm_eddy_coefficients = uses_cesm_eddy_coefficients(κ_skew, κ_symmetric) ?
        CESMEddyCoefficients(grid) : nothing

    hybrid_eddy_coefficients = uses_hybrid_eddy_coefficients(κ_skew, κ_symmetric) ?
        HybridEddyCoefficients(grid) : nothing

    eddy_slope_limiter = mixed_layer_tapering ? MixedLayerTapering(grid) : nothing

    diffusive_forcing, bottom_boundary_layer = bottom_boundary_layer_forcing(grid, bbl_diffusivity)
    advective_forcing, advective_bottom_boundary_layer =
        advective_bottom_boundary_layer_forcing(grid, bbl_transport_coefficient)

    restoring_forcing = overflow_restoring_forcing(grid, overflow_restoring_timescale)
    labrador_forcing  = labrador_restoring_forcing(grid, labrador_restoring_timescale; restoring_dir)

    ocean_forcing = merge_tracer_forcings(
                        merge_tracer_forcings(merge_tracer_forcings(diffusive_forcing, advective_forcing),
                                              restoring_forcing),
                        labrador_forcing)

    ocean = build_ocean(cfg, grid;
                        forcing = ocean_forcing,
                        κ_skew, κ_symmetric, Cᵇ, Cᵘⁿᵇ, Cᶠ, Cᶠ⁰, Cᶠᵟ, Cᵉc,
                        barotropic_substeps, Δt,
                        nemo_eddy_coefficients,
                        cesm_eddy_coefficients,
                        hybrid_eddy_coefficients,
                        eddy_slope_limiter,
                        boundary_value_mode_number,
                        boundary_value_minimum_speed,
                        biharmonic_timescale,
                        biharmonic_viscosity,
                        vertical_closure,
                        background_vertical_diffusivity,
                        background_vertical_viscosity,
                        implicit_vertical_advection,
                        tracer_advection_order,
                        boundary_scheme,
                        tracer_boundary_scheme,
                        momentum_boundary_scheme,
                        temperature_reference_variation,
                        salinity_reference_variation,
                        momentum_reference_variation,
                        implicit_bottom_drag,
                        bottom_drag_background_velocity,
                        skew_flux_formulation,
                        isopycnal_formulation,
                        restoring_under_sea_ice,
                        Cᵂu★,
                        restoring_dir, piston_velocity, chlorophyll,
                        initial_condition_blend_depth,
                        normalize_salinity,
                        additional_tracer_closure = filter(!isnothing, (river_κ, ice_melt_κ_closure, under_ice_ν_closure)),
                        start_date, end_date)
    log_setup_stage(arch, "ocean", setup_t₀)

    snow_thermodynamics = with_snow ?
        NumericalEarth.SeaIces.default_snow_thermodynamics(grid; thickness_categories = snow_thickness_categories) : nothing
    sea_ice = build_sea_ice(cfg, grid, ocean; restoring_dir, snow_thermodynamics, with_ice_dynamics,
                            with_ocean_surface_tilt, sea_ice_liquidus,
                            with_landfast_basal_stress, sea_ice_lateral_boundary_condition,
                            sea_ice_ocean_drag_coefficient, sea_ice_ocean_drag_reference_depth,
                            ice_compressive_strength, ice_salinity,
                            northern_sea_ice_initial_date, southern_sea_ice_initial_date,
                            thickness_categories, itd_shape,
                            ice_arch_region, ice_arch_stress, ice_arch_months)
    log_setup_stage(arch, "sea ice", setup_t₀)

    atmosphere, radiation = omip_forcing(arch, sea_ice;
                                         forcing_dir = atmosphere_dir,
                                         start_date,
                                         end_date,
                                         backend_size)
    log_setup_stage(arch, "atmosphere", setup_t₀)

    ice_freshwater_delivery = if ice_virtual_salt_flux
        VirtualSaltFluxIceFreshwater()
    elseif ice_freshwater_fraction == 1
        ConservativeIceFreshwater()
    else
        ScaledIceFreshwater(convert(eltype(grid), ice_freshwater_fraction))
    end

    ice_meltwater_enthalpy = ice_meltwater_at_interface_temperature ?
        InterfaceTemperatureMeltwater() : ZeroHeatContentMeltwater()

    coupled = build_coupled_model(ocean, sea_ice, atmosphere, radiation, land, flux_configuration;
                                  velocity_formulation, sea_ice_ocean_heat_transfer_coefficient,
                                  sea_ice_momentum_roughness_length,
                                  ice_freshwater_delivery, ice_meltwater_enthalpy)
    log_setup_stage(arch, "coupled model", setup_t₀)

    simulation = Simulation(coupled; Δt, stop_time)
    log_setup_stage(arch, "simulation", setup_t₀)

    # Only rank 0 creates dirs; others barrier inside @root and proceed once
    # the dirs exist. mkpath is idempotent so a race-free retry would also
    # work, but @root keeps the pattern symmetric with the staging code.
    @root for dir in [forcing_dir, restoring_dir, output_dir]
        if !isdir(dir)
            mkdir(dir)
        end
    end

    # Stage JRA55 data from slow disk to fast scratch
    if !isnothing(staging_dir)
        staging_callback = JRA55DataStagingCallback(; source_dir = forcing_dir,
                                                      staging_dir,
                                                      start_date)
        # Run monthly (≈1440 iterations at Δt=30min) — well ahead of year boundaries.
        # The callback only copies files at year transitions; otherwise it returns immediately.
        add_callback!(simulation, staging_callback, IterationInterval(1440))
    end

    # Keep the surface-salinity restoring globally salt-conserving: refresh its stored
    # zero-mean flux from the ocean state each step. Primed once here so the first step
    # already sees a valid corrected flux.
    salt_restoring = conservative_restoring(ocean.model.tracers.S.boundary_conditions.top.condition)
    if !isnothing(salt_restoring)
        refresh_restoring = RefreshSalinityRestoring(salt_restoring, !restoring_under_sea_ice)
        refresh_restoring(simulation)
        add_callback!(simulation, refresh_restoring, IterationInterval(1))
    end

    # The ice-melt diffusivity is rebuilt from the ice-ocean freshwater flux each step; prime it so the
    # first step already sees the real field rather than zeros.
    if !isnothing(ice_melt_diffusivity)
        refresh_ice_melt = RefreshIceMeltDiffusivity(ice_melt_diffusivity, ice_melt_mask,
                                                     convert(eltype(grid), ice_melt_mixing_threshold))
        refresh_ice_melt(simulation)
        add_callback!(simulation, refresh_ice_melt, IterationInterval(1))
    end

    if !isnothing(under_ice_ν_field)
        refresh_under_ice_ν = RefreshUnderIceViscosity(under_ice_ν_field, under_ice_ν_mask)
        refresh_under_ice_ν(simulation)
        add_callback!(simulation, refresh_under_ice_ν, IterationInterval(1))
    end

    if !isnothing(ice_arch_region)
        refresh_arch = RefreshArchStress(start_date)
        refresh_arch(simulation)
        add_callback!(simulation, refresh_arch, IterationInterval(1))
    end

    # NEMO recomputes its Treguier coefficient every step from the current stratification. Primed here
    # so the first step sees a valid field rather than zeros.
    if !isnothing(nemo_eddy_coefficients)
        compute_nemo_eddy_coefficients!(nemo_eddy_coefficients, ocean.model)
        add_callback!(simulation, RefreshNEMOEddyCoefficients(nemo_eddy_coefficients), IterationInterval(1))
    end

    # The face transports must exist before the first tendency evaluation, so seed them here as well
    # as refreshing them each step.
    if !isnothing(bottom_boundary_layer)
        update_bottom_boundary_layer!(simulation, bottom_boundary_layer)
        add_callback!(simulation, BottomBoundaryLayerUpdate(bottom_boundary_layer), IterationInterval(1))
    end

    if !isnothing(advective_bottom_boundary_layer)
        update_advective_bottom_boundary_layer!(simulation, advective_bottom_boundary_layer)
        add_callback!(simulation, AdvectiveBottomBoundaryLayerUpdate(advective_bottom_boundary_layer),
                      IterationInterval(1))
    end

    # Same for CESM's stratification-dependent coefficient.
    if !isnothing(cesm_eddy_coefficients)
        compute_cesm_eddy_coefficients!(cesm_eddy_coefficients, ocean.model)
        add_callback!(simulation, RefreshCESMEddyCoefficients(cesm_eddy_coefficients), IterationInterval(1))
    end

    # Same for the Treguier × Danabasoglu-Marshall hybrid.
    if !isnothing(hybrid_eddy_coefficients)
        compute_hybrid_eddy_coefficients!(hybrid_eddy_coefficients, ocean.model)
        add_callback!(simulation, RefreshHybridEddyCoefficients(hybrid_eddy_coefficients), IterationInterval(1))
    end

    # The mixed-layer taper reads a depth field refreshed from the model state each step.
    if !isnothing(eddy_slope_limiter)
        compute_tapering_mixed_layer_depth!(eddy_slope_limiter, ocean.model)
        add_callback!(simulation, RefreshMixedLayerTapering(eddy_slope_limiter), IterationInterval(1))
    end

    # Hold the global ocean volume fixed by removing the global mean of the atmospheric freshwater
    # flux. Registered after the flux assembly at the end of `update_state!`, so the correction is in
    # place for the ocean step that consumes it.
    if freshwater_normalization(normalize_freshwater) !== :none
        averaging_timescale = freshwater_averaging_timescale(normalize_freshwater)
        add_callback!(simulation, NormalizeTotalWater(coupled, Δt, averaging_timescale), IterationInterval(1))
    end


    wall_time = Ref(time_ns())
    add_callback!(simulation, omip_progress_callback(wall_time), IterationInterval(1))

    if diagnostics
        add_omip_diagnostics!(simulation;
                              surface_averaging_interval,
                              field_averaging_interval,
                              field_mean_interval,
                              checkpoint_interval,
                              output_dir,
                              filename_prefix,
                              file_splitting_interval)

        # Dispatches to the active method only for cfg == Val(:twelfthdegree);
        # other configurations get the no-op fallback.
        add_ke_spectrum_diagnostic!(simulation, cfg;
                                     output_dir,
                                     filename_prefix,
                                     flush_interval = field_averaging_interval)
    end

    return simulation
end

#####
##### WOA → TEOS-10 conversion utilities
#####
##### WOA's `t_an` is sea_water_temperature (in-situ, °C) and `s_an` is
##### sea_water_practical_salinity (PSS-78). Oceananigans' default
##### `TEOS10EquationOfState` expects Conservative Temperature (Θ) and
##### Absolute Salinity (S_A). The functions below convert WOA fields to the
##### TEOS-10 conventions in place, using SeawaterPolynomials (CPU only —
##### the SAAR atlas read is host-resident and the loop body is scalar).
#####

# Approximate hydrostatic pressure in dbar from depth z [m] (cell-center, negative for ocean).
@inline approx_pressure_dbar(z) = max(zero(z), -z)

"""
    woa_to_teos10!(T_field, S_field)

Convert WOA in-situ temperature `t [°C]` and Practical Salinity `S_P` to TEOS-10 Conservative Temperature `Θ`
and Absolute Salinity `S_A`, in place. Both fields must live on the same grid. The conversion runs on the host;
data is copied to/from the device automatically.
"""
function woa_to_teos10!(T_field, S_field)
    grid = T_field.grid
    cpu_arch = Oceananigans.DistributedComputations.cpu_architecture(architecture(grid))
    cpu_grid = on_architecture(cpu_arch, grid)
    Nx, Ny, Nz = size(grid)
    T_h = Array(interior(T_field))
    S_h = Array(interior(S_field))
    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        t  = T_h[i, j, k]
        SP = S_h[i, j, k]
        (isnan(t) || isnan(SP)) && continue
        λ = λnode(i, j, k, cpu_grid, Center(), Center(), Center())
        φ = φnode(i, j, k, cpu_grid, Center(), Center(), Center())
        z = znode(i, j, k, cpu_grid, Center(), Center(), Center())
        p = approx_pressure_dbar(z)
        SA = Sᴬ_from_Sᴾ(SP, p, λ, φ)
        Θ  = Θ_from_T(SA, t, p)
        T_h[i, j, k] = Θ
        S_h[i, j, k] = SA
    end
    copyto!(interior(T_field), T_h)
    copyto!(interior(S_field), S_h)
    return T_field, S_field
end

"""
    blend_monthly_initial_condition!(T_field, S_field, date, restoring_dir, blend_depth)

Blend WOA Monthly for the month of `date` into `T_field` and `S_field` above `blend_depth`, leaving the
WOA Annual values they already carry below it.

WOA Annual is a mean over a seasonal cycle whose winter half has kilometre-deep mixed layers, so a
1 January start initialized from it begins with a seasonal thermocline the season does not have. The
monthly climatology has the right surface structure but reaches only ~1500 m and is noisier where the
observations are sparse, so it is used where it is informative and the annual field is kept at depth.

The weight is 1 above `blend_depth` and tapers linearly to 0 at the monthly dataset's own vertical
reach, taken from its native grid rather than assumed; cells outside that reach keep the annual value
exactly. The monthly field is converted to Θ and Sᴬ cell by cell before it is blended, so the result
mixes TEOS-10 quantities and never mixed conventions.
"""
function blend_monthly_initial_condition!(T_field, S_field, date, restoring_dir, blend_depth)
    grid = T_field.grid
    cpu_arch = Oceananigans.DistributedComputations.cpu_architecture(architecture(grid))
    cpu_grid = on_architecture(cpu_arch, grid)
    Nx, Ny, Nz = size(grid)

    Tmeta = Metadatum(:temperature; dir = restoring_dir, dataset = WOAMonthly(), date)
    Smeta = Metadatum(:salinity;    dir = restoring_dir, dataset = WOAMonthly(), date)
    Tnative = Field(Tmeta, architecture(grid))
    Snative = Field(Smeta, architecture(grid))

    # `set!(field, metadata)` refuses a target grid deeper than the dataset, and the monthly
    # climatology reaches only ~1500 m; interpolating directly fills deeper cells with the nearest
    # interior value, which the taper below discards.
    Tm = CenterField(grid)
    Sm = CenterField(grid)
    DataWrangling.interpolate_physical!(Tm, Tnative, Tmeta)
    DataWrangling.interpolate_physical!(Sm, Snative, Smeta)

    deepest = Tnative.grid.Lz
    deepest > blend_depth ||
        throw(ArgumentError("initial_condition_blend_depth = $blend_depth m is below the monthly \
                             climatology's reach of $(round(deepest)) m"))

    Ta = Array(interior(T_field)); Sa = Array(interior(S_field))
    Tb = Array(interior(Tm));      Sb = Array(interior(Sm))

    # The monthly field is still in-situ T and Practical Salinity, and carries sentinel fills that are
    # not NaN outside its coverage — `woa_to_teos10!`'s `isnan` guard lets those through into
    # `Θ_from_T`, which throws out of `sqrt`. Convert here instead, behind a physical-range test.
    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        t = Tb[i, j, k]; SP = Sb[i, j, k]
        (isfinite(t) && isfinite(SP) && 1 < SP < 45 && -5 < t < 45) || continue
        isfinite(Ta[i, j, k]) && isfinite(Sa[i, j, k]) || continue
        d = -znode(i, j, k, cpu_grid, Center(), Center(), Center())
        w = clamp((deepest - d) / (deepest - blend_depth), 0, 1)
        w > 0 || continue
        λ = λnode(i, j, k, cpu_grid, Center(), Center(), Center())
        φ = φnode(i, j, k, cpu_grid, Center(), Center(), Center())
        p = approx_pressure_dbar(-d)
        Sᴬᵐ = Sᴬ_from_Sᴾ(SP, p, λ, φ)
        Θᵐ  = Θ_from_T(Sᴬᵐ, t, p)
        Ta[i, j, k] = w * Θᵐ  + (1 - w) * Ta[i, j, k]
        Sa[i, j, k] = w * Sᴬᵐ + (1 - w) * Sa[i, j, k]
    end

    copyto!(interior(T_field), Ta)
    copyto!(interior(S_field), Sa)

    @info "blended WOA Monthly $(month(date)) into the initial condition above $(blend_depth) m, " *
          "tapering to the monthly climatology's reach of $(round(Int, deepest)) m"

    return T_field, S_field
end

"""
    woa_salinity_fts_to_teos10!(fts)

Convert each time slice of a WOA Practical Salinity `FieldTimeSeries` to TEOS-10
Absolute Salinity, in place. Requires that all time indices be in memory
(use `time_indices_in_memory = length(metadata)`).
"""
function woa_salinity_fts_to_teos10!(fts)
    grid = fts.grid
    cpu_arch = Oceananigans.DistributedComputations.cpu_architecture(architecture(grid))
    cpu_grid = on_architecture(cpu_arch, grid)
    Nx, Ny, Nz = size(grid)
    Nt = length(fts.times)
    for t_idx in 1:Nt
        S_int = interior(fts[t_idx])
        S_h   = Array(S_int)
        for k in 1:Nz, j in 1:Ny, i in 1:Nx
            SP = S_h[i, j, k]
            isnan(SP) && continue
            λ = λnode(i, j, k, cpu_grid, Center(), Center(), Center())
            φ = φnode(i, j, k, cpu_grid, Center(), Center(), Center())
            z = znode(i, j, k, cpu_grid, Center(), Center(), Center())
            p = approx_pressure_dbar(z)
            S_h[i, j, k] = Sᴬ_from_Sᴾ(SP, p, λ, φ)
        end
        copyto!(S_int, S_h)
    end
    return fts
end

#####
##### Shared closure utilities
#####

@inline νhb(i, j, k, grid, ℓx, ℓy, ℓz, clock, fields, λ) = Oceananigans.Operators.Az(i, j, k, grid, ℓx, ℓy, ℓz)^2 / λ

# Background tracer diffusivity following Henyey et al. (1986).
@inline henyey_diffusivity(x, y, z, t) = max(2e-6, 1e-5 * abs(sind(y)))

# Bryan & Lewis (1979) depth-dependent background diffusivity, in the form GFDL models carry:
# κ = 0.8×10⁻⁴ + (1.05×10⁻⁴/π) atan[4.5×10⁻³ (|z| − 2500)], i.e. 3×10⁻⁵ m² s⁻¹ in the upper ocean
# rising across a ~2500 m transition to 1.3×10⁻⁴ m² s⁻¹ in the abyss. It buys the deep diapycnal
# upwelling that a uniform κ only reaches by also diffusing the thermocline.
@inline bryan_lewis_diffusivity(x, y, z, t) = 0.8e-4 + (1.05e-4 / π) * atan(4.5e-3 * (-z - 2500))

# Resolve the `background_vertical_diffusivity` option into something `VerticalScalarDiffusivity`
# accepts: `:henyey` keeps the latitude-dependent internal-wave scaling above (2×10⁻⁶ at the
# equator rising to 10⁻⁵ at the poles), `:bryan_lewis` the depth profile above, and a number sets a
# uniform interior κ instead. The background is the diapycnal diffusivity that feeds the upwelling
# closing the AMOC's lower limb, so it is the knob for the Bryan (1987) κ^(2/3) sensitivity test.
# Henyey in the thermocline, enhanced in the abyss. The two roles of the background diffusivity
# separate cleanly by depth: the *upper-ocean* value sets the global heat uptake — a uniform 3×10⁻⁵
# and Bryan & Lewis, which share that value but differ fourfold in the abyss, drift identically —
# while the *abyssal* value sets the diapycnal upwelling that keeps bottom water ventilated. Holding
# the thermocline at Henyey and raising only the abyss buys the second without paying for the first.
#
# The shape is Bryan & Lewis's arctangent, shifted to vanish at the surface and normalized so the
# enhancement reaches `ABYSSAL_ENHANCEMENT` at 5000 m.
const ABYSSAL_TRANSITION_DEPTH = 2500     # m, the Bryan & Lewis (1979) inflection
const ABYSSAL_TRANSITION_SCALE = 4.5e-3   # m⁻¹
const ABYSSAL_ENHANCEMENT = 5e-5          # m² s⁻¹ added to Henyey by 5000 m

@inline function abyssal_enhancement(z)
    shape = atan(ABYSSAL_TRANSITION_SCALE * (-z - ABYSSAL_TRANSITION_DEPTH)) +
            atan(ABYSSAL_TRANSITION_SCALE *       ABYSSAL_TRANSITION_DEPTH)
    return ABYSSAL_ENHANCEMENT * shape / 2atan(ABYSSAL_TRANSITION_SCALE * ABYSSAL_TRANSITION_DEPTH)
end

@inline abyssal_henyey_diffusivity(x, y, z, t) = henyey_diffusivity(x, y, z, t) + abyssal_enhancement(z)

resolve_background_diffusivity(κ::Number) = κ
resolve_background_diffusivity(κ::Symbol) =
    κ === :henyey         ? henyey_diffusivity :
    κ === :bryan_lewis    ? bryan_lewis_diffusivity :
    κ === :abyssal_henyey ? abyssal_henyey_diffusivity :
    throw(ArgumentError("background_vertical_diffusivity must be :henyey, :bryan_lewis, :abyssal_henyey or a number, got :$κ"))

# Default background momentum viscosity, shared by the closures that carry an explicit background.
# `nothing` keeps it, a number overrides it.
const default_background_viscosity = 1e-4

resolve_background_viscosity(ν) = isnothing(ν) ? default_background_viscosity : ν

# The Richardson-number and TKE closures carry their own interior background (NEMO's avtb/avmb,
# KPP's κⁱʷ/νⁱʷ, NORi's calibrated floor), so an external one would double-count rather than
# replace it.
function check_no_background_mixing(κ, ν, vertical_closure)
    if κ !== :henyey || !isnothing(ν)
        throw(ArgumentError("background_vertical_diffusivity / background_vertical_viscosity are not \
                             supported for the :$vertical_closure closure, which sets its own interior \
                             background internally"))
    end
    return nothing
end

# Step-function background diffusivity for the :simple closure.
# Strong mixing in the upper 100 m, weak interior diffusivity below.
@inline ν_step_simple(x, y, z, t) = ifelse(z >= -100, 1e-2, 1e-4)
@inline κ_step_simple(x, y, z, t) =
      z >= -10  ? 5e-2 :       # mimic BL mixing
      z >= -100 ? 1e-2 :
                  1e-5

# GM discretization. The first two forms are equivalent continuously but not discretely: `:diffusive`
# adds a skew flux to the tracer equation, `:advective` builds an eddy-induced velocity and advects
# with it. Only the latter puts the bolus transport into a velocity field the rest of the model
# (and `bolus_meridional_volume_flux_operation`) can read directly. `:boundary_value` is a different
# transport altogether — see `BoundaryValueTransport` — and is advective for the same reason.
gm_skew_flux_formulation(formulation::Symbol) =
    formulation === :diffusive ? DiffusiveFormulation() :
    formulation === :advective ? AdvectiveFormulation() :
    throw(ArgumentError("skew_flux_formulation must be :diffusive or :advective, got :$formulation"))

# `κ_skew` must be depth-independent under the boundary-value problem: it parameterizes stirring by
# the barotropic eddy velocity, and the vertical structure of the transport comes from the column
# problem instead (Ferrari et al. 2010, Section 4.1). The CESM and hybrid coefficients carry their
# own vertical shape, which would apply a vertical structure twice.
# The triad closure carries no `skew_flux_formulation`: its skew flux is always the diffusive form on
# the triad stencil, so pairing it with an advective or boundary-value transport would silently drop one
# of the two choices.
function check_isopycnal_formulation(isopycnal_formulation, skew_flux_formulation)
    isopycnal_formulation in (:standard, :triad) ||
        throw(ArgumentError("isopycnal_formulation must be :standard or :triad, got :$isopycnal_formulation"))
    if isopycnal_formulation === :triad && skew_flux_formulation !== :diffusive
        throw(ArgumentError("isopycnal_formulation = :triad supports only skew_flux_formulation = :diffusive, \
                             got :$skew_flux_formulation"))
    end
    return nothing
end

function check_depth_independent_skew_coefficient(κ_skew, skew_flux_formulation)
    if skew_flux_formulation === :boundary_value && (κ_skew === :cesm || κ_skew === :hybrid)
        throw(ArgumentError("skew_flux_formulation = :boundary_value requires a depth-independent \
                             κ_skew (a number or :nemo); the :$κ_skew coefficient has vertical \
                             structure, which the boundary-value problem supplies itself"))
    end
    return nothing
end

# Build a vertical-mixing closure tuple. The eddy and horizontal
# components are common to every option; the primary vertical closure
# and any background κ/ν are selected by `vertical_closure`.
function omip_closure(vertical_closure::Symbol;
                      κ_skew, κ_symmetric, 
                      Cᵇ = 0.28, 
                      Cᵘⁿᵇ = 0.0,
                      Cᶠ = 1.0,
                      Cᶠ⁰ = 1e9,
                      Cᶠᵟ = 0.75,
                      Cᵉc = 0.112,
                      biharmonic_timescale,
                      biharmonic_viscosity = nothing,
                      skew_flux_formulation = :diffusive,
                      isopycnal_formulation = :standard,
                      eddy_slope_limiter = nothing,
                      boundary_value_mode_number = 2,
                      boundary_value_minimum_speed = 0.1,
                      background_vertical_diffusivity = :henyey,
                      background_vertical_viscosity = nothing,
                      Cᵂu★ = nothing)

    background_κ = resolve_background_diffusivity(background_vertical_diffusivity)
    background_ν = resolve_background_viscosity(background_vertical_viscosity)

    primary, background = if vertical_closure == :catke
        # CATKEMixingLength is @kwdef over a single FT, so a mixed Int/Float keyword set has no method
        mixing_length = CATKEMixingLength(; Cᵇ = Float64(Cᵇ), Cᵘⁿᵇ = Float64(Cᵘⁿᵇ), Cᵉc = Float64(Cᵉc),
                                            Cᶠ = Float64(Cᶠ), Cᶠ⁰ = Float64(Cᶠ⁰), Cᶠᵟ = Float64(Cᶠᵟ))
        tke_eq = isnothing(Cᵂu★) ? CATKEEquation() : CATKEEquation(; Cᵂu★)
        catke = CATKEVerticalDiffusivity(VerticallyImplicitTimeDiscretization();
                                         mixing_length,
                                         maximum_viscosity=3,
                                         maximum_tracer_diffusivity=3,
                                         maximum_tke_diffusivity=3,
                                         negative_tke_damping_time_scale=10, # (seconds)
                                         turbulent_kinetic_energy_equation = tke_eq)
        catke, VerticalScalarDiffusivity(κ=background_κ, ν=background_ν)
    elseif vertical_closure == :simple
        check_no_background_mixing(background_vertical_diffusivity, background_vertical_viscosity, vertical_closure)
        convective = ConvectiveAdjustmentVerticalDiffusivity(VerticallyImplicitTimeDiscretization();
                                                             convective_κz = 1.0,
                                                             convective_νz = 1.0)
        background = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(); κ=κ_step_simple, ν=ν_step_simple)
        convective, background
    elseif vertical_closure == :nori
        check_no_background_mixing(background_vertical_diffusivity, background_vertical_viscosity, vertical_closure)
        NORiBaseVerticalDiffusivity(), nothing
    elseif vertical_closure == :rbvd
        convective = RiBasedVerticalDiffusivity(; horizontal_Ri_filter = Oceananigans.TurbulenceClosures.FivePointHorizontalFilter())
        background = VerticalScalarDiffusivity(κ=background_κ, ν=background_ν)
        convective, background
    elseif vertical_closure == :kpp
        check_no_background_mixing(background_vertical_diffusivity, background_vertical_viscosity, vertical_closure)
        KPPVerticalDiffusivity(), nothing
    elseif vertical_closure == :nemo_tke
        check_no_background_mixing(background_vertical_diffusivity, background_vertical_viscosity, vertical_closure)
        NEMOTKEVerticalDiffusivity(), nothing
    else
        error("Unknown vertical_closure: $vertical_closure. Options: :catke, :simple, :nori, :rbvd, :kpp, :nemo_tke")
    end

    # The boundary-value problem replaces the skew transport only, so Redi mixing rides along in a
    # companion closure carrying `κ_symmetric` alone.
    eddy = if isnothing(κ_skew) | isnothing(κ_symmetric)
        ()
    elseif skew_flux_formulation === :boundary_value
        limiter = isnothing(eddy_slope_limiter) ? FluxTapering(1e-2) : eddy_slope_limiter
        transport = BoundaryValueTransport(; κ_skew, slope_limiter = limiter,
                                           mode_number = boundary_value_mode_number,
                                           minimum_speed = boundary_value_minimum_speed)
        redi = IsopycnalSkewSymmetricDiffusivity(; κ_skew = nothing, κ_symmetric,
                                                 slope_limiter = limiter)
        (transport, redi)
    elseif isopycnal_formulation === :triad
        limiter = isnothing(eddy_slope_limiter) ? FluxTapering(1e-2) : eddy_slope_limiter
        # The triad constructor defaults to an explicit discretization, unlike every other closure here.
        # κ_symmetric S² over Δz_top = 1.5 m is stable only below Δt ≈ 15 s, so the vertical component
        # must be implicit. `TriadSlopeTapering` limits each triad on its own slope, which the
        # per-cell factor the closure applies otherwise does not bound.
        (TriadIsopycnalSkewSymmetricDiffusivity(VerticallyImplicitTimeDiscretization();
                                                κ_skew, κ_symmetric,
                                                slope_limiter = TriadSlopeTapering(limiter)),)
    else
        limiter = isnothing(eddy_slope_limiter) ? FluxTapering(1e-2) : eddy_slope_limiter
        (IsopycnalSkewSymmetricDiffusivity(; κ_skew, κ_symmetric, slope_limiter = limiter,
                                           skew_flux_formulation = gm_skew_flux_formulation(skew_flux_formulation)),)
    end

    horizontal_viscosity = if !isnothing(biharmonic_viscosity)
        HorizontalScalarBiharmonicDiffusivity(ν=biharmonic_viscosity)
    elseif !isnothing(biharmonic_timescale)
        HorizontalScalarBiharmonicDiffusivity(ν=νhb,
                                              discrete_form=true,
                                              parameters=biharmonic_timescale)
    else
        nothing
    end

    return filter(!isnothing, (primary, eddy..., horizontal_viscosity, background))
end

# Enhanced vertical mixing at river mouths (cf. NEMO `rn_avt_rnf` over `rn_hrnf`): an extra tracer
# diffusivity `κ` over the top `mixing_depth` metres at the routed river-mouth cells, mixing the
# fresh plume downward so a coastal surface cell cannot be freshened to zero. Added to the closure.
@inline river_mouth_κ(i, j, k, grid, clock, fields, mask) = @inbounds mask[i, j, k]

# The same idea as `river_mouth_κ`, applied where the sea ice is *melting* into the ocean. The ice-ocean
# freshwater flux is delivered to the top cell, so a melt event freshens one 1.5 m cell and CATKE has to
# erode the resulting lid; in reality the meltwater is stirred through the ice-ocean boundary layer under
# keels of 1-3 m. The melt condition is dynamic, but the kernel stays a bare lookup and the melt test is
# applied on the host by `RefreshIceMeltDiffusivity`, which rewrites the field every step. A compound
# parameter here (the flux, a mask and a threshold) widens the tendency kernel's argument types far
# enough that inference gives up on unrelated arguments and the GPU compile fails.
@inline ice_melt_κ(i, j, k, grid, clock, fields, κ) = @inbounds κ[i, j, k]

"""
    ice_melt_mixing_mask(grid; κ = 5e-4, mixing_depth = 10)

The depth taper of the ice-melt diffusivity: `κ` over the top `mixing_depth` metres and zero below.
`RefreshIceMeltDiffusivity` copies it into the live diffusivity wherever the ice is melting.

⚠ `κ` is **not** the river value: 0.1 m² s⁻¹ mixes √(2κt) ≈ 131 m in a day, which is convective
adjustment rather than boundary-layer stirring. The default 5e-4 gives ≈9 m day⁻¹, reaching the
default `mixing_depth` in about ten days — the right order for an under-ice boundary layer.
"""
function ice_melt_mixing_mask(grid; κ = 5e-4, mixing_depth = 10)
    zc = Array(znodes(grid, Center()))
    Nz = size(grid, 3)
    mask_data = zeros(eltype(grid), size(grid)...)
    for k in 1:Nz
        zc[k] > -mixing_depth && (mask_data[:, :, k] .= κ)
    end

    mask = CenterField(grid)
    set!(mask, mask_data)

    return mask
end

"""
    ice_melt_vertical_diffusivity(diffusivity)

Extra vertical tracer diffusivity read from the live `diffusivity` field, which
`RefreshIceMeltDiffusivity` rewrites each step from the ice-ocean freshwater flux.
"""
ice_melt_vertical_diffusivity(diffusivity) =
    VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization();
                              κ = ice_melt_κ, discrete_form = true,
                              loc = (Center, Center, Center), parameters = diffusivity)

"""
The ocean closure is built before the coupled model exists, so the ice-melt diffusivity is given its own
field and this callback rewrites it each step: the depth taper wherever the ice-ocean freshwater flux
exceeds `threshold`, zero elsewhere. The test is against a small positive `threshold` rather than zero so
a column hovering near no net exchange does not flicker the diffusivity on and off between steps;
`threshold` is in m s⁻¹, and 1e-9 is about 0.03 m yr⁻¹ of meltwater.
"""
struct RefreshIceMeltDiffusivity{K, M, FT}
    diffusivity :: K
    mask :: M
    threshold :: FT
end

function (r::RefreshIceMeltDiffusivity)(sim)
    io = sim.model.interfaces.sea_ice_ocean_interface
    isnothing(io) && return nothing
    Jʷ = parent(io.fluxes.freshwater)
    parent(r.diffusivity) .= ifelse.(Jʷ .> r.threshold, parent(r.mask), 0)
    return nothing
end

@inline under_ice_ν(i, j, k, grid, clock, fields, ν) = @inbounds ν[i, j, k]

"""
    under_ice_vertical_viscosity(viscosity)

Extra vertical viscosity read from the live `viscosity` field, which `RefreshUnderIceViscosity` rewrites
each step as the depth taper scaled by the ice concentration.
"""
under_ice_vertical_viscosity(viscosity) =
    VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization();
                              ν = under_ice_ν, discrete_form = true,
                              loc = (Center, Center, Center), parameters = viscosity)

struct RefreshUnderIceViscosity{V, M}
    viscosity :: V
    mask :: M
end

function (r::RefreshUnderIceViscosity)(sim)
    sea_ice = sim.model.sea_ice
    isnothing(sea_ice) && return nothing
    ℵ = parent(sea_ice.model.ice_concentration)
    parent(r.viscosity) .= parent(r.mask) .* ℵ
    return nothing
end

"""
A basal stress that adds a seasonal ice arch to `landfast`: `stress` (N m⁻²) over the cells inside the
longitude–latitude box `bounds = (λ₁, λ₂, φ₁, φ₂)` while the month lies in `months` (first to last inclusive,
wrapping the year), zero otherwise. `region` and `arch` are built on the ice velocity grid by
`materialize_basal_stress`; `RefreshArchStress` rewrites `arch` each step. `minimum_speed` regularizes τ/u
exactly as in `LandfastBasalStress`.
"""
struct SeasonalArchStress{L, R, A, FT, M, B}
    landfast :: L
    region :: R
    arch :: A
    stress :: FT
    minimum_speed :: FT
    months :: M
    bounds :: B
end

seasonal_arch_stress(landfast, FT; bounds, stress = 100, months = (12, 7)) =
    SeasonalArchStress(landfast, nothing, nothing, convert(FT, stress), convert(FT, 5e-5), months,
                       convert.(FT, Tuple(bounds)))

in_arch_bounds(λ, φ, (λ₁, λ₂, φ₁, φ₂)) = λ₁ <= (λ > 180 ? λ - 360 : λ) <= λ₂ && φ₁ <= φ <= φ₂

function ClimaSeaIce.SeaIceDynamics.materialize_basal_stress(b::SeasonalArchStress, grid)
    isnothing(b.region) || return b
    landfast = ClimaSeaIce.SeaIceDynamics.materialize_basal_stress(b.landfast, grid)
    λ = Array(λnodes(grid, Center(), Center(), Center()))
    φ = Array(φnodes(grid, Center(), Center(), Center()))
    Nx, Ny = size(grid, 1), size(grid, 2)
    region_data = [in_arch_bounds(λ[i, j], φ[i, j], b.bounds) ? one(eltype(grid)) : zero(eltype(grid))
                   for i in 1:Nx, j in 1:Ny]
    region = Field{Center, Center, Nothing}(grid)
    set!(region, region_data)
    arch = Field{Center, Center, Nothing}(grid)
    @root @info "Ice arch $(b.bounds): $(Int(sum(region_data))) cells arrested at $(b.stress) N m⁻², months $(b.months)"
    return SeasonalArchStress(landfast, region, arch, b.stress, b.minimum_speed, b.months, b.bounds)
end

@inline landfast_magnitude(i, j, k, grid, ::Nothing, fields) = zero(grid)
@inline landfast_magnitude(i, j, k, grid, b, fields) = ClimaSeaIce.SeaIceDynamics.basal_stress_magnitude(i, j, k, grid, b, fields)

@inline ClimaSeaIce.SeaIceDynamics.basal_stress_magnitude(i, j, k, grid, b::SeasonalArchStress, fields) =
    landfast_magnitude(i, j, k, grid, b.landfast, fields) + @inbounds b.arch[i, j, 1]

# `basal_τ{x,y}_coefficient` dispatch on the concrete stress type, so a new stress needs both, not just
# `basal_stress_magnitude`. Bodies mirror the `LandfastBasalStress` ones.
@inline function ClimaSeaIce.SeaIceDynamics.basal_τx_coefficient(i, j, k, grid, b::SeasonalArchStress, fields)
    kᵇ = ℑxᶠᵃᵃ(i, j, 1, grid, ClimaSeaIce.SeaIceDynamics.basal_stress_magnitude, b, fields)
    u  = @inbounds fields.u[i, j, k]
    v  = ℑxyᶠᶜᵃ(i, j, k, grid, fields.v)
    return kᵇ / (sqrt(u^2 + v^2) + b.minimum_speed)
end

@inline function ClimaSeaIce.SeaIceDynamics.basal_τy_coefficient(i, j, k, grid, b::SeasonalArchStress, fields)
    kᵇ = ℑyᵃᶠᵃ(i, j, 1, grid, ClimaSeaIce.SeaIceDynamics.basal_stress_magnitude, b, fields)
    u  = ℑxyᶜᶠᵃ(i, j, k, grid, fields.u)
    v  = @inbounds fields.v[i, j, k]
    return kᵇ / (sqrt(u^2 + v^2) + b.minimum_speed)
end

Adapt.adapt_structure(to, b::SeasonalArchStress) =
    SeasonalArchStress(Adapt.adapt(to, b.landfast), Adapt.adapt(to, b.region), Adapt.adapt(to, b.arch),
                       b.stress, b.minimum_speed, b.months, b.bounds)

Base.summary(b::SeasonalArchStress) =
    "SeasonalArchStress($(b.stress) N m⁻², months $(b.months), bounds $(b.bounds)) over $(summary(b.landfast))"

struct RefreshArchStress{D}
    start_date :: D
end

function (r::RefreshArchStress)(sim)
    b = sim.model.sea_ice.model.dynamics.basal_stress
    m = month(r.start_date + Second(round(Int, sim.model.clock.time)))
    m₁, m₂ = b.months
    in_season = m₁ <= m₂ ? (m₁ <= m <= m₂) : (m >= m₁ || m <= m₂)
    parent(b.arch) .= (in_season ? b.stress : zero(b.stress)) .* parent(b.region)
    return nothing
end

function river_mouth_vertical_diffusivity(grid, river_routing; κ = 0.1, mixing_depth = 10)
    zc = Array(znodes(grid, Center()))
    Nz = size(grid, 3)
    mask_data = zeros(eltype(grid), size(grid)...)

    for routing in values(river_routing)
        ti = Array(routing.target_i)
        tj = Array(routing.target_j)
        for n in eachindex(ti), k in 1:Nz
            zc[k] > -mixing_depth && (mask_data[ti[n], tj[n], k] = κ)
        end
    end

    mask = CenterField(grid)
    set!(mask, mask_data)

    return VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization();
                                     κ = river_mouth_κ, discrete_form = true,
                                     loc = (Center, Center, Center), parameters = mask)
end

#####
##### Salinity restoring (shared by both configurations)
#####

# Surface-only restoring, applied uniformly in space (no ice mask).
# Wrapped in a `ConservativeSurfaceFluxRestoring` so it rides on the ocean's top-flux BC
# via the `additional_surface_fluxes` kwarg of `ocean_simulation` while injecting zero net
# salt globally (OMIP zero-global-mean convention). The stored corrected flux is refreshed
# each step by `update_restoring_flux!` (registered as a callback in `omip_simulation`).
# WOA Practical Salinity is converted to TEOS-10 Absolute Salinity at setup so
# the restoring target matches the ocean prognostic-S convention.
function salinity_surface_restoring(grid, dataset;
                                    restoring_dir,
                                    piston_velocity,
                                    conservative = true,
                                    mask_under_sea_ice = false)

    Nz = size(grid, 3)
    Δz_surface = CUDA.@allowscalar Δzᶜᶜᶜ(1, 1, Nz, grid)

    rate = piston_velocity / (Δz_surface * days)

    Smetadata = Metadata(:salinity; dir = restoring_dir, dataset)

    restoring = DatasetRestoring(Smetadata, Oceananigans.Architectures.architecture(grid);
                                 rate,
                                 time_indices_in_memory = length(Smetadata))

    woa_salinity_fts_to_teos10!(restoring.field_time_series)

    surface_restoring = SurfaceFluxRestoring(restoring)

    # Either option needs the flux materialized into a field; without both, the bare inline restoring
    # is cheaper and behaves identically.
    materialize = conservative | mask_under_sea_ice
    return materialize ? ConservativeSurfaceFluxRestoring(surface_restoring, grid; normalize = conservative) :
                         surface_restoring
end

#####
##### Grid builder
#####

function find_exponential_scale(Nz, depth, Δzᵀ; tolerance = 1e-7, maxiter = 200)
    Δzᵁ = depth / Nz
    Δzᵀ < Δzᵁ || throw(ArgumentError("Δzᵀ = $Δzᵀ must be < depth/Nz = $Δzᵁ"))
    Δzᵀ > 0   || throw(ArgumentError("Δzᵀ = $Δzᵀ must be positive"))

    Δz_at_scale(h) = depth * expm1(Δzᵁ / h) / expm1(depth / h)

    h⁻ = Δzᵁ / 1000
    h⁺ = 1000 * depth

    for _ in 1:maxiter
        h = (h⁻ + h⁺) / 2
        Δz = Δz_at_scale(h)
        abs(Δz - Δzᵀ) <= tolerance * Δzᵀ && return h
        Δz < Δzᵀ ? (h⁻ = h) : (h⁺ = h)
    end
    error("Could not converge to scale matching Δz_top = $Δzᵀ within relative tolerance $tolerance")
end

exponential_scale(Nz, depth, ::Nothing) = 1300
exponential_scale(Nz, depth, Δz_top)    = find_exponential_scale(Nz, depth, Δz_top)

# Partial bottom cells resolve sill depths and slopes continuously instead of in full-cell
# steps. Documented benefits: mean-circulation and boundary-current realism (Gulf Stream
# separation, NAC path — Barnier et al. 2006) and reduced staircase entrainment of downslope
# overflows (Winton et al. 1998). They mitigate but do not cure the too-shallow NADW that
# every configuration shares (zero-crossing ~2900 m vs ~4300 m in RAPID); the documented full
# fix in z-coordinate models is a dedicated overflow parameterization (Legg et al. 2009;
# Danabasoglu et al. 2010).
bottom_immersed_boundary(bottom_height, partial_cell_bathymetry) =
    partial_cell_bathymetry ? PartialCellBottom(bottom_height) : GridFittedBottom(bottom_height)

function build_grid(config, arch, Nz, depth; Δz_top = nothing, partial_cell_bathymetry = false)

    Nx = config == Val(:halfdegree) ? 720 : throw("Configuration $(config) does not exist")

    Ny = Nx ÷ 2

    scale = exponential_scale(Nz, depth, Δz_top)
    z_faces = ExponentialDiscretization(Nz, -depth, 0; scale, mutable=true)

    base_grid = TripolarGrid(arch;
                             size = (Nx, Ny, Nz),
                             z = z_faces,
                             halo = (8, 8, 8))

    bottom_height = regrid_bathymetry(base_grid;
                                    minimum_depth = 20,
                                    major_basins = 1,
                                    interpolation_passes = 25)

    return ImmersedBoundaryGrid(base_grid, bottom_immersed_boundary(bottom_height, partial_cell_bathymetry); active_cells_map = true)
end

build_grid(::Val{:orca}, arch, Nz, depth; Δz_top = nothing, partial_cell_bathymetry = false)          = build_grid(ORCAOne(),     arch, Nz, depth; Δz_top, partial_cell_bathymetry)
build_grid(::Val{:quarterdegree}, arch, Nz, depth; Δz_top = nothing, partial_cell_bathymetry = false) = build_grid(ORCAQuarter(), arch, Nz, depth; Δz_top, partial_cell_bathymetry)
build_grid(::Val{:twelfthdegree}, arch, Nz, depth; Δz_top = nothing, partial_cell_bathymetry = false) = build_grid(ORCATwelfth(), arch, Nz, depth; Δz_top, partial_cell_bathymetry)

# The Gulf of Ob and the Yenisei Gulf are ~5 m deep for hundreds of kilometres, so their full river
# discharge lands in a single 1.5 m top cell with no water column to mix into and the salinity collapses.
# Closing the sub-`minimum_depth` cells of each gulf turns them to land; the river routing then relocates
# the discharge onto the deeper Kara Sea shelf just north, where the top ~10 m spans five cells.
# Boxes are (λ_min, λ_max, φ_min, φ_max) in degrees.
const kara_river_closures = ((68.0, 77.0, 66.0, 72.6),   # Gulf of Ob
                             (77.0, 85.0, 70.0, 73.8))    # Yenisei Gulf

@kernel function _close_shallow_regions!(bottom_height, grid, regions, minimum_depth)
    i, j = @index(Global, NTuple)
    λ = λnode(i, j, 1, grid, Center(), Center(), Center())
    φ = φnode(i, j, 1, grid, Center(), Center(), Center())
    @inbounds z = bottom_height[i, j, 1]
    shallow = (z < 0) & (z > -minimum_depth)
    closed = false
    for (λ₀, λ₁, φ₀, φ₁) in regions
        closed = closed | (shallow & (λ ≥ λ₀) & (λ ≤ λ₁) & (φ ≥ φ₀) & (φ ≤ φ₁))
    end
    @inbounds bottom_height[i, j, 1] = ifelse(closed, oftype(z, 100), z)
end

function close_shallow_river_regions(grid; regions = kara_river_closures, minimum_depth = 10, partial_cell_bathymetry = false)
    arch      = architecture(grid)
    underlying = grid.underlying_grid
    bottom    = bottom_height_field(grid)
    launch!(arch, underlying, :xy, _close_shallow_regions!, bottom, underlying, regions,
            convert(eltype(grid), minimum_depth))
    fill_halo_regions!(bottom)
    remove_minor_basins!(bottom, 1)
    return ImmersedBoundaryGrid(underlying, bottom_immersed_boundary(bottom, partial_cell_bathymetry); active_cells_map = true)
end

function build_grid(dataset::ORCADataset, arch, Nz, depth; Δz_top = nothing, partial_cell_bathymetry = false)

    scale = exponential_scale(Nz, depth, Δz_top)
    z_faces = ExponentialDiscretization(Nz, -depth, 0; scale, mutable=true)

    grid = ORCAGrid(arch;
                    dataset,
                    Nz,
                    z = z_faces,
                    halo = (8, 8, 8),
                    with_bathymetry = true,
                    partial_cell_bathymetry,
                    major_basins = 1,
                    active_cells_map = true)

    return grid # close_shallow_river_regions(grid)
end

# Locally-runnable testing configuration: the NEMO eORCA1 (~1ᵒ) mesh, used to reproduce the
# quarter-degree spurious high-latitude ice + surface salinity drift at a fraction of the cost.
build_grid(::Val{:test}, arch, Nz, depth; Δz_top = nothing, partial_cell_bathymetry = false) = build_grid(Val(:orca), arch, Nz, depth; Δz_top, partial_cell_bathymetry)

#####
##### ORCA builder
#####

using Oceananigans.TimeSteppers: AdaptiveVerticallyImplicitDiscretization, ExplicitTimeDiscretization
using Oceananigans.Utils: NormalDivision

# `nothing` keeps the `WENOVectorInvariant` per-term defaults: vorticity_order = 9, everything else 5.
config_momentum_advection_order(::Val{:orca})          = 5
config_momentum_advection_order(::Val{:test})          = 5
config_momentum_advection_order(::Val{:halfdegree})    = 5
config_momentum_advection_order(::Val{:quarterdegree}) = nothing
config_momentum_advection_order(::Val{:twelfthdegree}) = nothing

#####
##### Boundary reconstructions
#####

# The reconstruction a WENO buffer chain terminates in, used in the one cell whose stencil no longer fits:
# the domain buffer and, on an `ImmersedBoundaryGrid`, any cell adjacent to an `inactive_node`
# (`Advection/immersed_advective_fluxes.jl`). The interior keeps the full order.
#
#   :default  `Centered(order=2)` for tracers and `UpwindBiased(order=1)` for momentum, the Oceananigans defaults
#   :upwind   first-order upwind, monotone, in exactly those cells and nowhere else
#   :cwenoz   the third-order central-WENO reconstruction of Semplice, Travaglia and Puppo (2022), whose stencil
#             extends only inwards, blending an inward parabola, a linear polynomial and a constant with Z-weights
#
# `reference_variation` sets the oscillation scale ϵ = reference_variation² below which CWENOZ reads the stencil as
# noise and keeps third order. It carries the units of the reconstructed field and no grid spacing, so one value serves
# every direction; zero reads ϵ off the stencil, which is a pure shape measure and limits at any amplitude.
tracer_boundary_reconstruction(::Val{:default}, reference_variation) = nothing
tracer_boundary_reconstruction(::Val{:upwind},  reference_variation) = UpwindBiased(order=1)
tracer_boundary_reconstruction(::Val{:cwenoz},  reference_variation) = CWENOZ(; reference_variation)

momentum_boundary_reconstruction(::Val{:default}, reference_variation) = UpwindBiased(order=1)
momentum_boundary_reconstruction(::Val{:upwind},  reference_variation) = UpwindBiased(order=1)
momentum_boundary_reconstruction(::Val{:cwenoz},  reference_variation) = CWENOZ(; reference_variation)

function boundary_scheme_value(boundary_scheme)
    boundary_scheme ∈ (:default, :upwind, :cwenoz) ||
        throw(ArgumentError("boundary_scheme must be :default, :upwind or :cwenoz, got $boundary_scheme"))

    return Val(boundary_scheme)
end

"""
    split_tracer_advection(order, time_discretization, boundary_scheme, reference_variation)

Tracer advection of order `order` whose reconstructions terminate in a boundary scheme carrying
`reference_variation`. One value serves both directions: the variation is in units of the tracer and carries no
grid spacing. The split is only so that the vertical direction takes `time_discretization`, which is where the
adaptive-implicit treatment applies.
"""
function split_tracer_advection(order, time_discretization, boundary_scheme, reference_variation)

    tracer_boundary_scheme = tracer_boundary_reconstruction(boundary_scheme, reference_variation)

    horizontal = WENO(; order, boundary_scheme = tracer_boundary_scheme)
    vertical   = WENO(; order, time_discretization, boundary_scheme = tracer_boundary_scheme)

    return FluxFormAdvection(horizontal, horizontal, vertical)
end

"""
    split_momentum_advection(order, time_discretization, boundary_scheme, reference_variation)

Vector-invariant momentum advection whose four reconstructions terminate in a boundary scheme chosen per direction,
reproducing `WENOVectorInvariant` in every other respect. The vorticity, divergence and kinetic-energy-gradient terms
are the horizontal ones, and they reconstruct a vorticity, a divergence flux and a squared velocity: three different
units, so no single variation carries them and they always read the oscillation scale off the stencil. The vertical
term reconstructs velocity, so `reference_variation` is a speed in m/s and applies there alone.
"""
function split_momentum_advection(order, time_discretization, boundary_scheme, reference_variation)

    vorticity_order, remaining_order = isnothing(order) ? (9, 5) : (order, order)

    horizontal_boundary_scheme = momentum_boundary_reconstruction(boundary_scheme, 0)
    vertical_boundary_scheme   = momentum_boundary_reconstruction(boundary_scheme, reference_variation)

    vorticity_scheme               = WENO(order=vorticity_order, boundary_scheme=horizontal_boundary_scheme)
    divergence_scheme              = WENO(order=remaining_order, boundary_scheme=horizontal_boundary_scheme)
    kinetic_energy_gradient_scheme = WENO(order=remaining_order, boundary_scheme=horizontal_boundary_scheme)
    vertical_advection_scheme      = WENO(; order = remaining_order, time_discretization,
                                            boundary_scheme = vertical_boundary_scheme)

    return VectorInvariant(; vorticity_scheme,
                             vertical_advection_scheme,
                             divergence_scheme,
                             kinetic_energy_gradient_scheme)
end

struct ConfigDefault end

@inline resolve_config_default(value, default)           = value
@inline resolve_config_default(::ConfigDefault, default) = default

config_Nz(::Val)        = 100
config_Nz(::Val{:test}) = 15

config_κ_skew(::Val{:orca})          = 800
config_κ_skew(::Val{:halfdegree})    = 250
config_κ_skew(::Val{:quarterdegree}) = nothing
config_κ_skew(::Val{:twelfthdegree}) = nothing
config_κ_skew(::Val{:test})          = nothing

config_κ_symmetric(::Val{:orca})          = 800
config_κ_symmetric(::Val{:halfdegree})    = 250
config_κ_symmetric(::Val{:quarterdegree}) = nothing
config_κ_symmetric(::Val{:twelfthdegree}) = nothing
config_κ_symmetric(::Val{:test})          = nothing

# River-mouth footprint. `:orca` and its `:test` preset keep the historical eight-cell footprint so
# their long integrations stay reproducible; the refined grids use a geographic radius, because a fixed
# cell count shrinks the footprint as the square of the resolution ratio and concentrates the
# freshwater flux until coastal salinity hits zero.
config_river_spread_radius(::Val)                 = 1.2
config_river_spread_radius(::Val{:orca})          = nothing
config_river_spread_radius(::Val{:test})          = nothing

config_river_spread_cells(::Val)                 = nothing
config_river_spread_cells(::Val{:orca})          = 8
config_river_spread_cells(::Val{:test})          = 8

# Split-explicit barotropic substeps. The free surface is stable only while the barotropic gravity
# wave stays inside a substep, √(gH) Δτ ≲ 0.7 Δx with Δτ ∝ Δt / substeps, so the count a grid needs
# grows as it is refined: eORCA025 carries a CFL of 0.92 at Δt = 20 minutes with 100, and 1.38 — 30 000
# cells past unity, and a first-step blow-up — at Δt = 30 minutes. `validate_barotropic_substeps`
# reports the count the configuration actually needs.
#
# ⚠ The substep count and Δt must move TOGETHER. `orca` runs at Δt = 5400 s, where the generic 100
# gives a gravity-wave Courant number of 1.16 and the free surface NaNs on the first step; 300 is what
# the whole Option A family runs at (M★ = 216, halo 218, inside Ny = 297; the ceiling is 409).
config_barotropic_substeps(::Val)                 = 100
config_barotropic_substeps(::Val{:orca})          = 300
config_barotropic_substeps(::Val{:quarterdegree}) = 200
config_barotropic_substeps(::Val{:twelfthdegree}) = 200

config_biharmonic_timescale(::Val)                 = 50days
config_biharmonic_timescale(::Val{:quarterdegree}) = nothing
config_biharmonic_timescale(::Val{:twelfthdegree}) = nothing
config_biharmonic_timescale(::Val{:test})          = 10days

config_Δt(::Val)                 = 30minutes
config_Δt(::Val{:orca})          = 5400        # 90 minutes; needs 300 barotropic substeps, see above
config_Δt(::Val{:quarterdegree}) = 20minutes
config_Δt(::Val{:twelfthdegree}) = 5minutes
config_Δt(::Val{:test})          = 45minutes

config_Δz_top(::Val)                 = nothing
config_Δz_top(::Val{:quarterdegree}) = 1.5
config_Δz_top(::Val{:twelfthdegree}) = 1.5
config_Δz_top(::Val{:test})          = 1.5

# Buoyancy gradients are only needed by the GM/Redi closures; the eddy-resolving
# configurations run without isopycnal diffusivities, so skip materializing them.
config_materialize_buoyancy_gradients(::Val)                 = true
config_materialize_buoyancy_gradients(::Val{:quarterdegree}) = false
config_materialize_buoyancy_gradients(::Val{:twelfthdegree}) = false
config_materialize_buoyancy_gradients(::Val{:test})          = false


"""
    omip_radiative_forcing(grid, chlorophyll, restoring_dir)

Return the `TwoColorRadiation` the ocean is forced with. `chlorophyll` is `:seawifs` for the SeaWiFS
monthly climatology, cycled every year, or anything `TwoColorRadiation` accepts directly — a number for
globally uniform optics, a surface `Field`, or a `FieldTimeSeries`.

`:none` returns `nothing`, which routes the whole shortwave into the surface heat flux instead of the
interior: `shortwave_radiative_forcing` falls back to returning the flux to the boundary condition
rather than stashing it for a penetrating scheme. The vertical closure's surface buoyancy flux is
built from the tracer boundary conditions, so this is also the only setting under which it sees the
shortwave at all.
"""
function omip_radiative_forcing(grid, chlorophyll, restoring_dir)
    chlorophyll === :none && return nothing

    if chlorophyll === :seawifs
        dates = (DateTime(2000, 1, 1), DateTime(2000, 12, 1))
        metadata = Metadata(:chlorophyll; dataset = SeaWiFSMonthly(), dates, dir = restoring_dir)
        chlorophyll = FieldTimeSeries(metadata, grid)
    end

    return TwoColorRadiation(grid; chlorophyll)
end

@inline function barotropic_courant_number(i, j, k, grid, Δτ, g)
    H  = static_column_depthᶜᶜᵃ(i, j, grid)
    Δx = Δxᶜᶜᶜ(i, j, k, grid)
    Δy = Δyᶜᶜᶜ(i, j, k, grid)
    return sqrt(g * max(0, H)) * Δτ * sqrt(1 / Δx^2 + 1 / Δy^2)
end

# The free surface is stable only while the barotropic gravity wave stays inside a substep. A fixed
# substep count silently violates that when `Δt` is raised or the grid is refined, and the free surface
# then blows up on the very first step with nothing to show for it but a NaN. eORCA025 at Δt = 30 minutes
# with 100 substeps carries 1.38 in the Ross Sea and 30 000 cells past unity. The Courant number is taken
# cell by cell, because the grid's finest spacing and its greatest depth are nowhere near each other and
# pairing them condemns configurations that run.
function barotropic_free_surface(grid, substeps, Δt; cfl = 0.7)
    free_surface = SplitExplicitFreeSurface(grid; substeps)
    Δτ = free_surface.substepping.fractional_step_size * Δt
    g  = free_surface.gravitational_acceleration

    courant = KernelFunctionOperation{Center, Center, Center}(barotropic_courant_number, grid, Δτ, g)
    maximum_courant_number = maximum(courant)

    if maximum_courant_number > cfl
        needed = ceil(Int, substeps * maximum_courant_number / cfl)
        @warn string(substeps, " barotropic substeps carry a gravity-wave Courant number of ",
                     round(maximum_courant_number, digits=2), " at Δt = ", prettytime(Δt),
                     ", above ", cfl, "; this grid needs ", needed,
                     ". Raise `barotropic_substeps` or lower `Δt`.")
    end

    return free_surface
end

function build_ocean(config, grid;
                     κ_skew, κ_symmetric, Cᵇ = 0.28, Cᵘⁿᵇ = 0.0, Cᵉc = 0.112,
                     Cᶠ = 1.0, Cᶠ⁰ = 1e9, Cᶠᵟ = 0.75,
                     barotropic_substeps = 100,
                     Δt,
                     restoring_dir, piston_velocity,
                     initial_condition_blend_depth = nothing,
                     chlorophyll = :seawifs,
                     biharmonic_timescale,
                     biharmonic_viscosity = nothing,
                     vertical_closure = :catke,
                     implicit_vertical_advection = true,
                     tracer_advection_order = 7,
                     boundary_scheme = :default,
                     tracer_boundary_scheme = boundary_scheme,
                     momentum_boundary_scheme = boundary_scheme,
                     temperature_reference_variation = 0,
                     salinity_reference_variation = 0,
                     momentum_reference_variation = 0,
                     implicit_bottom_drag = true,
                     bottom_drag_background_velocity = 0,
                     skew_flux_formulation = :diffusive,
                     isopycnal_formulation = :standard,
                     nemo_eddy_coefficients = nothing,
                     cesm_eddy_coefficients = nothing,
                     hybrid_eddy_coefficients = nothing,
                     eddy_slope_limiter = nothing,
                     boundary_value_mode_number = 2,
                     boundary_value_minimum_speed = 0.1,
                     background_vertical_diffusivity = :henyey,
                     background_vertical_viscosity = nothing,
                     restoring_under_sea_ice = true,
                     Cᵂu★ = nothing,
                     normalize_salinity = true,
                     additional_tracer_closure = nothing,
                     forcing = NamedTuple(),
                     start_date, end_date)

    κ_skew      = resolve_nemo_coefficient(κ_skew,      nemo_eddy_coefficients, :skew_coefficient)
    κ_symmetric = resolve_nemo_coefficient(κ_symmetric, nemo_eddy_coefficients, :symmetric_coefficient)
    κ_skew      = resolve_cesm_coefficient(κ_skew,      cesm_eddy_coefficients, :skew_coefficient)
    κ_symmetric = resolve_cesm_coefficient(κ_symmetric, cesm_eddy_coefficients, :symmetric_coefficient)
    κ_skew      = resolve_hybrid_coefficient(κ_skew,      hybrid_eddy_coefficients, :skew_coefficient)
    κ_symmetric = resolve_hybrid_coefficient(κ_symmetric, hybrid_eddy_coefficients, :symmetric_coefficient)

    if !isnothing(κ_skew) && !isnothing(κ_symmetric)
        κ_skew, κ_symmetric = fold_safe_constant_coefficients(grid, κ_skew, κ_symmetric)
    end

    additional_surface_fluxes = if piston_velocity == 0
        NamedTuple()
    else
        salt_restoring = salinity_surface_restoring(grid, WOAMonthly(); restoring_dir, piston_velocity,
                                                    conservative = normalize_salinity,
                                                    mask_under_sea_ice = !restoring_under_sea_ice)
        (; S = salt_restoring)
    end

    closure = omip_closure(vertical_closure;
                           κ_skew, κ_symmetric, Cᵇ, Cᵘⁿᵇ, Cᶠ, Cᶠ⁰, Cᶠᵟ, Cᵉc,
                           biharmonic_timescale, biharmonic_viscosity,
                           skew_flux_formulation,
                           isopycnal_formulation,
                           eddy_slope_limiter,
                           boundary_value_mode_number,
                           boundary_value_minimum_speed,
                           background_vertical_diffusivity,
                           background_vertical_viscosity,
                           Cᵂu★)
    extra_closures = additional_tracer_closure isa Tuple ? additional_tracer_closure :
                     isnothing(additional_tracer_closure) ? () : (additional_tracer_closure,)
    closure = (closure..., extra_closures...)
    coriolis = HydrostaticSphericalCoriolis(scheme = Oceananigans.Coriolis.EnstrophyConserving())

    time_discretization = implicit_vertical_advection ?
        AdaptiveVerticallyImplicitDiscretization(cfl=0.5) : ExplicitTimeDiscretization()

    tracer_boundary_scheme   = boundary_scheme_value(tracer_boundary_scheme)
    momentum_boundary_scheme = boundary_scheme_value(momentum_boundary_scheme)

    momentum_advection = split_momentum_advection(config_momentum_advection_order(config),
                                                  time_discretization, momentum_boundary_scheme,
                                                  momentum_reference_variation)

    # Turbulent kinetic energy keeps the `ocean_simulation` default: its variation scale is ~1e-3 m²/s².
    tracer_advection = (T = split_tracer_advection(tracer_advection_order, time_discretization,
                                                   tracer_boundary_scheme, temperature_reference_variation),
                        S = split_tracer_advection(tracer_advection_order, time_discretization,
                                                   tracer_boundary_scheme, salinity_reference_variation))

    ocean = ocean_simulation(grid;
                             Δt = 1minutes,
                             radiative_forcing = omip_radiative_forcing(grid, chlorophyll, restoring_dir),
                             momentum_advection,
                             tracer_advection,
                             coriolis,
                             implicit_bottom_drag,
                             bottom_drag_background_velocity,
                             timestepper = :SplitRungeKutta3,
                             materialize_buoyancy_gradients = config_materialize_buoyancy_gradients(config),
                             free_surface = barotropic_free_surface(grid, barotropic_substeps, Δt),
                             additional_surface_fluxes,
                             forcing,
                             closure)

    # Load WOA Annual T (in-situ, °C) and S (Practical) onto the model grid,
    # convert to TEOS-10 Conservative T and Absolute Salinity in place, then
    # initialize the prognostic ocean state from the converted fields.
    T_init = CenterField(grid)
    S_init = CenterField(grid)
    set!(T_init, Metadatum(:temperature; dir=restoring_dir, dataset=WOAAnnual()))
    set!(S_init, Metadatum(:salinity;    dir=restoring_dir, dataset=WOAAnnual()))
    woa_to_teos10!(T_init, S_init)

    isnothing(initial_condition_blend_depth) ||
        blend_monthly_initial_condition!(T_init, S_init, start_date, restoring_dir,
                                         initial_condition_blend_depth)

    set!(ocean.model, T=T_init, S=S_init)

    return ocean
end

#####
##### Sea Ice builder
#####

# `:teos10` is the relation fitted to the TEOS-10 freezing point in Conservative Temperature, which is
# what the ocean carries; `:linear` restores ClimaSeaIce's own (0, 0.054) default, up to 0.032 K warmer.
resolve_liquidus(::Val{:teos10}, FT) = NumericalEarth.SeaIces.conservative_temperature_liquidus(FT)
resolve_liquidus(::Val{:linear}, FT) = LinearLiquidus(FT)
resolve_liquidus(name::Symbol, FT) = resolve_liquidus(Val(name), FT)

function build_sea_ice(config, grid, ocean; restoring_dir, snow_thermodynamics = nothing,
                       with_ice_dynamics = true,
                       with_ocean_surface_tilt = false,
                       sea_ice_liquidus = :teos10,
                       with_landfast_basal_stress = true,
                       sea_ice_lateral_boundary_condition = :no_slip,
                       sea_ice_ocean_drag_coefficient = 5.5e-3,
                       sea_ice_ocean_drag_reference_depth = 6,
                       ice_compressive_strength = 27500,
                       ice_salinity = 4,
                       northern_sea_ice_initial_date = DateTime(1993, 1, 1),
                       southern_sea_ice_initial_date = DateTime(1993, 1, 1),
                       thickness_categories = 1,
                       itd_shape = nothing,
                       ice_arch_region = nothing,
                       ice_arch_stress = 100,
                       ice_arch_months = (12, 7))

    basal_stress = with_landfast_basal_stress ? LandfastBasalStress(eltype(grid)) : nothing
    isnothing(ice_arch_region) ||
        (basal_stress = seasonal_arch_stress(basal_stress, eltype(grid); bounds = ice_arch_region,
                                             stress = ice_arch_stress, months = ice_arch_months))

    dynamics = if with_ice_dynamics
        rheology = ElastoViscoPlasticRheology(eltype(grid); ice_compressive_strength)
        NumericalEarth.SeaIces.sea_ice_dynamics(grid, ocean; basal_stress, sea_ice_ocean_drag_coefficient,
                                                sea_ice_ocean_drag_reference_depth, rheology,
                                                with_ocean_surface_tilt)
    else
        nothing
    end

    sea_ice = sea_ice_simulation(grid, ocean;
                                 advection = ClimaSeaIce.IncrementalRemapping(),
                                 lateral_boundary_condition = sea_ice_lateral_boundary_condition,
                                 dynamics,
                                 liquidus = resolve_liquidus(sea_ice_liquidus, eltype(grid)),
                                 ice_salinity,
                                 thickness_categories, itd_shape,
                                 snow_thermodynamics)

    set_sea_ice_initial_condition!(sea_ice.model, restoring_dir,
                                   northern_sea_ice_initial_date, southern_sea_ice_initial_date)

    return sea_ice
end

"""
    set_sea_ice_initial_condition!(model, restoring_dir, northern_date, southern_date)

Initialize `h` and `ℵ` from ECCO4Monthly, taking the northern hemisphere from `northern_date` and the
southern hemisphere from `southern_date`.

The two hemispheres reach their ice minimum six months apart, so one date for the whole globe starts
one of them at its seasonal maximum: an ORCA run beginning on 1 January carries an Arctic pack at its
March extent, which the first spring melts into the Labrador. Sea ice never reaches the equator, so
the two fields are stitched at ``φ = 0`` with no taper.
"""
function set_sea_ice_initial_condition!(model, restoring_dir, northern_date, southern_date)
    thickness(date)     = Metadatum(:sea_ice_thickness;     dir=restoring_dir, dataset=ECCO4Monthly(), date)
    concentration(date) = Metadatum(:sea_ice_concentration; dir=restoring_dir, dataset=ECCO4Monthly(), date)

    set!(model, h = thickness(southern_date), ℵ = concentration(southern_date))

    northern_date == southern_date && return model

    h = model.ice_thickness
    ℵ = model.ice_concentration
    hˢ = Array(interior(h)); ℵˢ = Array(interior(ℵ))

    set!(model, h = thickness(northern_date), ℵ = concentration(northern_date))

    hⁿ = Array(interior(h)); ℵⁿ = Array(interior(ℵ))

    grid = h.grid
    cpu_arch = Oceananigans.DistributedComputations.cpu_architecture(architecture(grid))
    cpu_grid = on_architecture(cpu_arch, grid)
    Nx, Ny = size(grid)[1:2]

    for j in 1:Ny, i in 1:Nx
        φnode(i, j, 1, cpu_grid, Center(), Center(), Center()) < 0 || continue
        hⁿ[i, j, 1] = hˢ[i, j, 1]
        ℵⁿ[i, j, 1] = ℵˢ[i, j, 1]
    end

    copyto!(interior(h), hⁿ)
    copyto!(interior(ℵ), ℵⁿ)

    @info "sea ice initialized from ECCO4Monthly $(northern_date) north of the equator and " *
          "$(southern_date) south of it"

    return model
end

#####
##### Progress callback
#####

# Per-step resource probe, on unless `OMIP_PROBE=0`. It separates the candidate causes of the
# multi-second stalls that appear part-way through a run while the timestep cost itself stays flat:
# a Julia heap at its ceiling shows up as `gc`/`pauses`, CUDA.jl's reclaim path forces a `full` sweep on
# every allocation it cannot serve, a background thread parked inside a blocking NetCDF read shows up as
# `safepoint`, and a filling device pool shows up as `gpu_free` falling with `gpu_pool` pinned at `gpu_total`.
const resource_probe_enabled = get(ENV, "OMIP_PROBE", "1") == "1"

# `total_time_to_safepoint` postdates some supported Julia versions; resolved once so the per-step read
# stays a plain field access.
const gc_reports_safepoint_time = hasfield(Base.GC_Num, :total_time_to_safepoint)

@inline safepoint_time_ns(gc) =
    gc_reports_safepoint_time ? Int64(getfield(gc, :total_time_to_safepoint)) : zero(Int64)

device_memory_report(arch) = nothing

function device_memory_report(::GPU)
    CUDA.functional() || return nothing
    try
        reserved = CUDA.cached_memory()
        return (; free = CUDA.free_memory(),
                  total = CUDA.total_memory(),
                  reserved = reserved === missing ? 0 : Int(reserved))
    catch
        return nothing
    end
end

function omip_progress_callback(wall_time)
    initial_gc = Base.gc_num()
    previous_gc_time        = Ref(Int64(Base.gc_time_ns()))
    previous_collections    = Ref(Int64(initial_gc.pause))
    previous_full_sweeps    = Ref(Int64(initial_gc.full_sweep))
    previous_safepoint_time = Ref(safepoint_time_ns(initial_gc))
    previous_allocated      = Ref(Int64(Base.gc_bytes()))

    function progress(sim)
        sea_ice = sim.model.sea_ice
        ocean   = sim.model.ocean

        hmax = maximum(sea_ice.model.ice_thickness)
        ℵmax = maximum(sea_ice.model.ice_concentration)
        Tmax = maximum(ocean.model.tracers.T)
        Tmin = minimum(ocean.model.tracers.T)
        Smax = maximum(ocean.model.tracers.S)
        Smin = minimum(ocean.model.tracers.S)
        umax = maximum(ocean.model.velocities.u)
        vmax = maximum(ocean.model.velocities.v)
        wmax = maximum(ocean.model.velocities.w)

        step_time = 1e-9 * (time_ns() - wall_time[])

        msg1 = @sprintf("time: %s, iteration: %d, Δt: %s, ",
                        prettytime(sim), iteration(sim), prettytime(sim.Δt))
        msg2 = @sprintf("max(h): %.2e m, max(ℵ): %.2e ", hmax, ℵmax)
        msg3 = @sprintf("extrema(T, S): (%.2f, %.2f) ᵒC, (%.2f, %.2f) psu ",
                        Tmin, Tmax, Smin, Smax)
        msg4 = @sprintf("maximum(u): (%.2e, %.2e, %.2e) m/s, ", umax, vmax, wmax)
        msg5 = @sprintf("wall time: %s", prettytime(step_time))

        @info msg1 * msg2 * msg3 * msg4 * msg5

        if resource_probe_enabled
            gc = Base.gc_num()
            gc_time = Int64(Base.gc_time_ns())
            safepoint_time = safepoint_time_ns(gc)
            allocated = Int64(Base.gc_bytes())

            probe = @sprintf("PROBE iter=%d step=%.3fs gc=%.3fs pauses=%d full=%d safepoint=%.3fs allocd=%.3fGiB live=%.2fGiB host_free=%.2fGiB",
                             iteration(sim), step_time,
                             1e-9 * (gc_time - previous_gc_time[]),
                             gc.pause - previous_collections[],
                             gc.full_sweep - previous_full_sweeps[],
                             1e-9 * (safepoint_time - previous_safepoint_time[]),
                             (allocated - previous_allocated[]) / 2^30,
                             Base.gc_live_bytes() / 2^30,
                             Sys.free_memory() / 2^30)

            device = device_memory_report(architecture(ocean.model.grid))

            if !isnothing(device)
                probe *= @sprintf(" gpu_free=%.2fGiB gpu_total=%.2fGiB gpu_pool=%.2fGiB",
                                  device.free / 2^30, device.total / 2^30, device.reserved / 2^30)
            end

            @info probe

            previous_gc_time[]        = gc_time
            previous_collections[]    = gc.pause
            previous_full_sweeps[]    = gc.full_sweep
            previous_safepoint_time[] = safepoint_time
            previous_allocated[]      = allocated
        end

        # Determinism probe: hash the prognostic state at a few iterations.
        # Compare these hashes between two pickup-from-same-checkpoint runs:
        # the first iteration whose hashes differ pinpoints when divergence
        # is introduced. Cheap (~ms; one host copy of each field).
        iter = iteration(sim)
        if iter in (1, 5, 100, 1000)
            T  = Array(parent(interior(ocean.model.tracers.T)))
            S  = Array(parent(interior(ocean.model.tracers.S)))
            u  = Array(parent(interior(ocean.model.velocities.u)))
            h  = Array(parent(interior(sea_ice.model.ice_thickness)))
            @info @sprintf("STATE_HASH iter=%d  T=%016x  S=%016x  u=%016x  h=%016x",
                           iter, hash(T), hash(S), hash(u), hash(h))
        end

        wall_time[] = time_ns()

        return nothing
    end

    return progress
end
