using StaticArrays
using Thermodynamics
using SurfaceFluxes

using ..OceanSeaIceModels: reference_density,
                           heat_capacity,
                           sea_ice_thickness,
                           downwelling_radiation,
                           freshwater_flux

using ClimaSeaIce: SlabSeaIceModel

using Oceananigans: HydrostaticFreeSurfaceModel, architecture
using Oceananigans.Grids: inactive_node
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: ConstantField
using Oceananigans.Utils: launch!, Time

using Oceananigans.Operators: ℑxᶜᵃᵃ, ℑyᵃᶜᵃ

using KernelAbstractions: @kernel, @index

#####
##### Container for organizing information related to fluxes
#####

struct OceanSeaIceSurfaceFluxes{T, P, C, R, PI, PC, FT}
    turbulent :: T
    prescribed :: P
    total :: C
    radiation :: R
    previous_ice_thickness :: PI
    previous_ice_concentration :: PC
    # The ocean is Boussinesq, so these are _only_ coupled properties:
    ocean_reference_density :: FT
    ocean_heat_capacity :: FT
end

Base.summary(crf::OceanSeaIceSurfaceFluxes) = "OceanSeaIceSurfaceFluxes"
Base.show(io::IO, crf::OceanSeaIceSurfaceFluxes) = print(io, summary(crf))

function OceanSeaIceSurfaceFluxes(ocean, sea_ice=nothing;
                                  atmosphere = nothing,
                                  radiation = nothing,
                                  ocean_reference_density = reference_density(ocean),
                                  ocean_heat_capacity = heat_capacity(ocean))

    FT = eltype(ocean.model.grid)

    ocean_reference_density = convert(FT, ocean_reference_density)
    ocean_heat_capacity = convert(FT, ocean_heat_capacity)

    # It's the "thermodynamics gravitational acceleration"
    # (as opposed to the one used for the free surface)
    g = ocean.model.buoyancy.model.gravitational_acceleration
    turbulent_fluxes = SimilarityTheoryTurbulentFluxes(FT, gravitational_acceleration=g)

    prescribed_fluxes = nothing

    if isnothing(sea_ice)
        previous_ice_thickness = nothing
        previous_ice_concentration = nothing
    else
        previous_ice_thickness = deepcopy(sea_ice.model.ice_thickness)
        previous_ice_concentration = deepcopy(sea_ice.model.ice_concentration)
    end

    ocean_grid = ocean.model.grid
    ρₒ = ocean_reference_density
    Jᵘ = surface_flux(ocean.model.velocities.u)
    Jᵛ = surface_flux(ocean.model.velocities.v)
    Jᵘᶜᶜᶜ = Field{Center, Center, Nothing}(ocean_grid)
    Jᵛᶜᶜᶜ = Field{Center, Center, Nothing}(ocean_grid)

    ocean_momentum_fluxes = (u = Jᵘ,       # fluxes used in the model
                             v = Jᵛ,       #
                             τˣ = ρₒ * Jᵘ, # momentum fluxes multiplied by reference density
                             τʸ = ρₒ * Jᵛ, # 
                             uᶜᶜᶜ = Jᵘᶜᶜᶜ, # fluxes computed by bulk formula at cell centers
                             vᶜᶜᶜ = Jᵛᶜᶜᶜ)

    tracers = ocean.model.tracers
    ocean_tracer_fluxes = NamedTuple(name => surface_flux(tracers[name]) for name in keys(tracers))

    cₚ = ocean_heat_capacity
    ocean_heat_flux = ρₒ * cₚ * ocean_tracer_fluxes.T

    total_ocean_fluxes = (momentum = ocean_momentum_fluxes,
                          tracers = ocean_tracer_fluxes,
                          heat = ocean_heat_flux)

    total_fluxes = (; ocean=total_ocean_fluxes)

    return OceanSeaIceSurfaceFluxes(turbulent_fluxes,
                                    prescribed_fluxes,
                                    total_fluxes,
                                    radiation,
                                    previous_ice_thickness,
                                    previous_ice_concentration,
                                    ocean_reference_density,
                                    ocean_heat_capacity)
end

#####
##### Surface flux computation
#####

const c = Center()
const f = Face()

function compute_atmosphere_ocean_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice
    atmosphere = coupled_model.atmosphere

    # Basic model properties
    grid = ocean.model.grid
    arch = architecture(grid)
    clock = ocean.model.clock

    # Ocean, atmosphere, and sea ice state
    ocean_velocities = surface_velocities(ocean)
    ocean_tracers    = surface_tracers(ocean)

    atmosphere_velocities            = atmosphere.velocities
    atmosphere_tracers               = atmosphere.tracers
    atmosphere_pressure              = atmosphere.pressure
    atmosphere_downwelling_radiation = atmosphere.downwelling_radiation
    atmosphere_freshwater_flux       = atmosphere.freshwater_flux

    ice_thickness = sea_ice_thickness(sea_ice)

    # Fluxes, and flux contributors
    centered_velocity_fluxes = (u = coupled_model.fluxes.total.ocean.momentum.uᶜᶜᶜ,
                                v = coupled_model.fluxes.total.ocean.momentum.vᶜᶜᶜ)

    staggered_velocity_fluxes = (u = coupled_model.fluxes.total.ocean.momentum.u,
                                 v = coupled_model.fluxes.total.ocean.momentum.v)

    net_tracer_fluxes    = coupled_model.fluxes.total.ocean.tracers
    turbulent_fluxes     = coupled_model.fluxes.turbulent
    prescribed_fluxes    = coupled_model.fluxes.prescribed
    radiation_properties = coupled_model.fluxes.radiation

    ocean_state = merge(ocean_velocities, ocean_tracers)
    atmosphere_state = merge(atmosphere_velocities, atmosphere_tracers, (; p=atmosphere_pressure))

    launch!(arch, grid, :xy, compute_atmosphere_ocean_turbulent_fluxes!,
            grid, clock,
            centered_velocity_fluxes,
            net_tracer_fluxes,
            turbulent_fluxes,
            atmosphere_freshwater_flux,
            atmosphere_downwelling_radiation,
            radiation_properties,
            ocean_state,
            atmosphere_state,
            atmosphere.reference_height, # height at which the state is known
            atmosphere.thermodynamics_parameters,
            coupled_model.fluxes.ocean_reference_density,
            coupled_model.fluxes.ocean_heat_capacity,
            ice_thickness)

    # Note: I think this can be avoided if we modify the preceding kernel
    # to compute from 0:Nx+1, ie in halo regions
    fill_halo_regions!(centered_velocity_fluxes)

    launch!(arch, grid, :xy, accumulate_atmosphere_ocean_fluxes!,
            grid, clock,
            staggered_velocity_fluxes,
            net_tracer_fluxes,
            centered_velocity_fluxes,
            prescribed_fluxes,
            ocean_state,
            atmosphere_state,
            atmosphere_downwelling_radiation,
            radiation_properties,
            atmosphere_freshwater_flux,
            coupled_model.fluxes.ocean_reference_density,
            coupled_model.fluxes.ocean_heat_capacity,
            ice_thickness)

    return nothing
end

@kernel function compute_atmosphere_ocean_turbulent_fluxes!(grid,
                                                            clock,
                                                            centered_velocity_fluxes,
                                                            net_tracer_fluxes,
                                                            turbulent_fluxes,
                                                            prescribed_freshwater_flux,
                                                            downwelling_radiation,
                                                            radiation_properties,
                                                            ocean_state,
                                                            atmos_state,
                                                            atmosphere_reference_height,
                                                            atmosphere_thermodynamics_parameters,
                                                            ocean_reference_density,
                                                            ocean_heat_capacity,
                                                            ice_thickness)

    i, j = @index(Global, NTuple)

    time = Time(clock.time)

    # Extract state variables at cell centers
    @inbounds begin
        # Ocean state
        uₒ = ℑxᶜᵃᵃ(i, j, 1, grid, ocean_state.u)
        vₒ = ℑyᵃᶜᵃ(i, j, 1, grid, ocean_state.v)
        Uₒ = SVector(uₒ, vₒ)
        Tₒ = ocean_state.T[i, j, 1]
        Sₒ = ocean_state.S[i, j, 1]

        # Atmos state
        uₐ = atmos_state.u[i, j, 1, time]
        vₐ = atmos_state.v[i, j, 1, time]
        Uₐ = SVector(uₐ, vₐ)

        Tₐ = atmos_state.T[i, j, 1, time]
        pₐ = atmos_state.p[i, j, 1, time]
        qᵗₐ = atmos_state.q[i, j, 1] # total specific humidity
    end

    # Build atmospheric state
    ℂ = atmosphere_thermodynamics_parameters
    ψₐ = thermodynamic_atmospheric_state = AtmosphericThermodynamics.PhaseEquil_pTq(ℂ, pₐ, Tₐ, qᵗₐ)

    # Build surface state with saturated specific humidity
    surface_type = AtmosphericThermodynamics.Liquid()
    q★ = surface_saturation_specific_humidity(ℂ, Tₒ, ψₐ, surface_type)
    
    # Thermodynamic and dynamic state at the surface
    ψ₀ = thermodynamic_surface_state = AtmosphericThermodynamics.PhaseEquil_pTq(ℂ, pₐ, Tₒ, q★)
    Ψ₀ = dynamic_surface_state = SurfaceFluxes.StateValues(zero(grid), Uₒ, ψ₀)

    # Thermodynamic and dynamic state at reference level h above the surface
    h = atmosphere_reference_height # elevation of atmos variables relative to surface
    Ψₐ = dynamic_atmos_state = SurfaceFluxes.StateValues(h, Uₐ, ψₐ)

    # Roughness lengths...
    zᵐ = convert(eltype(grid), 1e-2)
    zʰ = convert(eltype(grid), 1e-2)

    values = SurfaceFluxes.ValuesOnly(Ψₐ, Ψ₀, zᵐ, zʰ)
    conditions = SurfaceFluxes.surface_conditions(turbulent_fluxes, values)
    
    # Compute heat fluxes, bulk flux first
    Qd = net_downwelling_radiation(i, j, grid, time, downwelling_radiation, radiation_properties)
    Qu = net_upwelling_radiation(i, j, grid, time, radiation_properties, ocean_state)
    Qs = conditions.shf
    Qℓ = conditions.lhf
    ΣQ = Qu + Qs + Qℓ

    E  = conditions.evaporation
    F  = cross_realm_flux(i, j, grid, time, prescribed_freshwater_flux)

    @show conditions

    Jᵘ = centered_velocity_fluxes.u
    Jᵛ = centered_velocity_fluxes.v
    Jᵀ = net_tracer_fluxes.T
    Jˢ = net_tracer_fluxes.S

    ρₒ = ocean_reference_density
    cₒ = ocean_heat_capacity

    atmos_ocean_Jᵘ = conditions.ρτxz / ρₒ
    atmos_ocean_Jᵛ = conditions.ρτyz / ρₒ
    atmos_ocean_Jᵀ = ΣQ / (ρₒ * cₒ)
    atmos_ocean_Jˢ = Sₒ * (E + F)

    kᴺ = size(grid, 3) # index of the top ocean cell
    inactive = inactive_node(i, j, kᴺ, grid, c, c, c)

    @inbounds begin
        Jᵘ[i, j, 1] = ifelse(inactive, zero(grid), atmos_ocean_Jᵘ)
        Jᵛ[i, j, 1] = ifelse(inactive, zero(grid), atmos_ocean_Jᵛ)
        Jᵀ[i, j, 1] = ifelse(inactive, zero(grid), atmos_ocean_Jᵀ)
        Jˢ[i, j, 1] = ifelse(inactive, zero(grid), atmos_ocean_Jˢ)
    end
end

@kernel function accumulate_atmosphere_ocean_fluxes!(grid,
                                                     clock,
                                                     staggered_velocity_fluxes,
                                                     centered_velocity_fluxes)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3) # index of the top ocean cell

    time = Time(clock.time)

    Jᵘ = staggered_velocity_fluxes.u
    Jᵛ = staggered_velocity_fluxes.v

    @inbounds begin
        Jᵘ[i, j, 1] = ℑxᶠᵃᵃ(i, j, k, grid, centered_velocity_fluxes.u)
        Jᵛ[i, j, 1] = ℑyᵃᶠᵃ(i, j, k, grid, centered_velocity_fluxes.v)
    end

    #=
    # Note: there could one or more formula(e)
    τˣ_formula = bulk_momentum_flux_formulae.u
    τʸ_formula = bulk_momentum_flux_formulae.v
    Q_formula = bulk_heat_flux_formulae
    F_formula = bulk_tracer_flux_formulae.S

    atmos_state_names = keys(atmos_state)
    ocean_state_names = keys(atmos_state)

    atmos_state_ij = stateindex(atmos_state, i, j, 1, time)
    ocean_state_ij = stateindex(ocean_state, i, j, 1, time)

    # Compute transfer velocity scale
    ΔUᶠᶜᶜ = bulk_velocity_scaleᶠᶜᶜ(i, j, grid, time, bulk_velocity, atmos_state, ocean_state)
    ΔUᶜᶠᶜ = bulk_velocity_scaleᶜᶠᶜ(i, j, grid, time, bulk_velocity, atmos_state, ocean_state)
    ΔUᶜᶜᶜ = bulk_velocity_scaleᶜᶜᶜ(i, j, grid, time, bulk_velocity, atmos_state, ocean_state)

    # Compute momentum fluxes
    τˣ = cross_realm_flux(i, j, grid, time, τˣ_formula, ΔUᶠᶜᶜ, atmos_state, ocean_state)
    τʸ = cross_realm_flux(i, j, grid, time, τʸ_formula, ΔUᶜᶠᶜ, atmos_state, ocean_state)

    # Compute heat fluxes, bulk flux first
    Qd = net_downwelling_radiation(i, j, grid, time, downwelling_radiation, radiation)
     
    Qu = net_upwelling_radiation(i, j, grid, time, radiation, ocean_state)
    Q★ = cross_realm_flux(i, j, grid, time, Q_formula, ΔUᶜᶜᶜ, atmos_state_ij, ocean_state_ij)
    Q = Q★ + Qd + Qu

    # Compute salinity fluxes, bulk flux first
    Fp = cross_realm_flux(i, j, grid, time, prescribed_freshwater_flux)
    F★ = cross_realm_flux(i, j, grid, time, F_formula, ΔUᶜᶜᶜ, atmos_state_ij, ocean_state_ij)
    F = F★ + Fp

    # Then the rest of the heat fluxes
    ρₒ = ocean_reference_density
    cₚ = ocean_heat_capacity

    atmos_ocean_Jᵘ = τˣ / ρₒ
    atmos_ocean_Jᵛ = τʸ / ρₒ
    atmos_ocean_Jᵀ = Q / (ρₒ * cₚ)

    S = ocean_state_ij.S
    atmos_ocean_Jˢ = S * F

    @inbounds begin
        # Set fluxes
        # TODO: should this be peripheral_node?
        Jᵘ[i, j, 1] = ifelse(inactive_node(i, j, kᴺ, grid, f, c, c), zero(grid), atmos_ocean_Jᵘ)
        Jᵛ[i, j, 1] = ifelse(inactive_node(i, j, kᴺ, grid, c, f, c), zero(grid), atmos_ocean_Jᵛ)
        Jᵀ[i, j, 1] = ifelse(inactive_node(i, j, kᴺ, grid, c, c, c), zero(grid), atmos_ocean_Jᵀ)
        Jˢ[i, j, 1] = ifelse(inactive_node(i, j, kᴺ, grid, c, c, c), zero(grid), atmos_ocean_Jˢ)
    end
    =#
end

@inline function net_downwelling_radiation(i, j, grid, time, downwelling_radiation, radiation)
    Qˢʷ = downwelling_radiation.shortwave
    Qˡʷ = downwelling_radiation.longwave
    α = stateindex(radiation.reflection.ocean, i, j, 1, time)

    return @inbounds - (1 - α) * Qˢʷ[i, j, 1, time] - Qˡʷ[i, j, 1, time]
end

@inline function net_upwelling_radiation(i, j, grid, time, radiation, ocean_state)
    σ = radiation.stefan_boltzmann_constant
    Tᵣ = radiation.reference_temperature
    ϵ = stateindex(radiation.emission.ocean, i, j, 1, time)

    # Ocean surface temperature (departure from reference, typically in ᵒC)
    Tₒ = @inbounds ocean_state.T[i, j, 1]

    # Note: positive implies _upward_ heat flux, and therefore cooling.
    return σ * ϵ * (Tₒ + Tᵣ)^4
end

@inline cross_realm_flux(i, j, grid, time, ::Nothing,        args...) = zero(grid)
@inline cross_realm_flux(i, j, grid, time, a::AbstractArray, args...) = stateindex(a, i, j, 1, time)
@inline cross_realm_flux(i, j, grid, time, nt::NamedTuple,   args...) = cross_realm_flux(i, j, grid, time, values(nt), args...)

@inline cross_realm_flux(i, j, grid, time, flux_tuple::Tuple{<:Any, <:Any}, args...) =
    cross_realm_flux(i, j, grid, time, flux_tuple[1], args...) +
    cross_realm_flux(i, j, grid, time, flux_tuple[2], args...)

@inline cross_realm_flux(i, j, grid, time, flux_tuple::Tuple{<:Any, <:Any, <:Any}, args...) =
    cross_realm_flux(i, j, grid, time, flux_tuple[1], args...) +
    cross_realm_flux(i, j, grid, time, flux_tuple[2], args...) +
    cross_realm_flux(i, j, grid, time, flux_tuple[3], args...)

@inline cross_realm_flux(i, j, grid, time, flux_tuple::Tuple{<:Any, <:Any, <:Any, <:Any}, args...) =
    cross_realm_flux(i, j, grid, time, flux_tuple[1], args...) +
    cross_realm_flux(i, j, grid, time, flux_tuple[2], args...) +
    cross_realm_flux(i, j, grid, time, flux_tuple[3], args...) +
    cross_realm_flux(i, j, grid, time, flux_tuple[4], args...)

#=
function default_atmosphere_ocean_fluxes(FT=Float64, tracers=tuple(:S))
    # Note: we are constantly coping with the fact that the ocean is ᵒC.
    ocean_reference_temperature = 273.15
    momentum_transfer_coefficient = 5e-3
    evaporation_transfer_coefficient = 1e-3
    sensible_heat_transfer_coefficient = 2e-3
    vaporization_enthalpy  = 2.5e-3

    momentum_transfer_coefficient      = convert(FT, momentum_transfer_coefficient)
    evaporation_transfer_coefficient   = convert(FT, evaporation_transfer_coefficient)
    sensible_heat_transfer_coefficient = convert(FT, sensible_heat_transfer_coefficient)
    vaporization_enthalpy              = convert(FT, vaporization_enthalpy)

    τˣ = BulkFormula(RelativeUVelocity(), momentum_transfer_coefficient)
    τʸ = BulkFormula(RelativeVVelocity(), momentum_transfer_coefficient)
    momentum_flux_formulae = (u=τˣ, v=τʸ)

    # Note: reference temperature comes in here
    water_specific_humidity_difference = SpecificHumidity(FT)
    evaporation = nothing #BulkFormula(SpecificHumidity(FT), evaporation_transfer_coefficient)
    tracer_flux_formulae = (; S = evaporation)

    latent_heat_difference = LatentHeat(specific_humidity_difference = water_specific_humidity_difference; vaporization_enthalpy)
    latent_heat_formula    = nothing #BulkFormula(latent_heat_difference,  evaporation_transfer_coefficient)

    sensible_heat_difference = SensibleHeat(FT; ocean_reference_temperature)
    sensible_heat_formula = BulkFormula(sensible_heat_difference, sensible_heat_transfer_coefficient)

    heat_flux_formulae = (sensible_heat_formula, latent_heat_formula)

    return CrossRealmFluxes(momentum = momentum_flux_formulae,
                            heat = heat_flux_formulae,
                            tracers = tracer_flux_formulae)
end
=#

#=
#####
##### Bulk formula
#####

"""
    BulkFormula(air_sea_difference, transfer_coefficient)

The basic structure of a flux `J` computed by a bulk formula is:

```math
J = - ρₐ * C * Δc * ΔU
```

where `ρₐ` is the density of air, `C` is the `transfer_coefficient`,
`Δc` is the air_sea_difference, and `ΔU` is the bulk velocity scale.
"""
struct BulkFormula{F, CD}
    air_sea_difference :: F
    transfer_coefficient :: CD
end

@inline function cross_realm_flux(i, j, grid, time, formula::BulkFormula, ΔU, atmos_state, ocean_state)
    ρₐ = stateindex(atmos_state.ρ, i, j, 1, time)
    C = formula.transfer_coefficient
    Δc = air_sea_difference(i, j, grid, time, formula.air_sea_difference, atmos_state, ocean_state)

    # Note the sign convention, which corresponds to positive upward fluxes:
    return - ρₐ * C * Δc * ΔU
end

#####
##### Air-sea differences
#####

@inline air_sea_difference(i, j, grid, time, air, sea) = stateindex(air, i, j, 1, time) -
                                                         stateindex(sea, i, j, 1, time)

struct RelativeUVelocity end
struct RelativeVVelocity end

@inline function air_sea_difference(i, j, grid, time, ::RelativeUVelocity, atmos_state, ocean_state)
    uₐ = atmos_state.u
    uₒ = ocean_state.u
    return air_sea_difference(i, j, grid, time, uₐ, uₒ)
end

@inline function air_sea_difference(i, j, grid, time, ::RelativeVVelocity, atmos_state, ocean_state)
    vₐ = atmos_state.v
    vₒ = ocean_state.v
    return air_sea_difference(i, j, grid, time, vₐ, vₒ)
end

struct SensibleHeat{FT}
    ocean_reference_temperature :: FT
end

SensibleHeat(FT::DataType=Float64; ocean_reference_temperature=273.15) =
    SensibleHeat(convert(FT, ocean_reference_temperature))

@inline function air_sea_difference(i, j, grid, time, Qs::SensibleHeat, atmos_state, ocean_state)
    cₚ = stateindex(atmos_state.cₚ, i, j, 1, time)
    Tₐ = atmos_state.T

    # Compute ocean temperature in degrees K
    Tᵣ = Qs.ocean_reference_temperature 
    Tₒᵢ = stateindex(ocean_state.T, i, j, 1, time)
    Tₒ = Tₒᵢ + Tᵣ
    
    ΔT = air_sea_difference(i, j, grid, time, Tₐ, Tₒ)

    return @inbounds cₚ[i, j, 1] * ΔT
end

struct SpecificHumidity{S}
    saturation_specific_humidity :: S

    @doc """
        SpecificHumidity(FT = Float64;
                               saturation_specific_humidity = LargeYeagerSaturationVaporFraction(FT))

    """
    function SpecificHumidity(FT = Float64;
                              saturation_specific_humidity = LargeYeagerSaturationVaporFraction(FT))
        S = typeof(saturation_specific_humidity)
        return new{S}(saturation_specific_humidity)
    end
end

struct LargeYeagerSaturationVaporFraction{FT}
    q₀ :: FT
    c₁ :: FT
    c₂ :: FT
    reference_temperature :: FT
end

"""
    LargeYeagerSaturationVaporFraction(FT = Float64;
                                       q₀ = 0.98,
                                       c₁ = 640380,
                                       c₂ = -5107.4,
                                       reference_temperature = 273.15)

"""
function LargeYeagerSaturationVaporFraction(FT = Float64;
                                            q₀ = 0.98,
                                            c₁ = 640380,
                                            c₂ = -5107.4,
                                            reference_temperature = 273.15)

    return LargeYeagerSaturationVaporFraction(convert(FT, q₀),
                                              convert(FT, c₁),
                                              convert(FT, c₂),
                                              convert(FT, reference_temperature))
end

@inline function saturation_specific_humidity(i, j, grid, time,
                                           ratio::LargeYeagerSaturationVaporFraction,
                                           atmos_state, ocean_state)

    Tₒ = stateindex(ocean_state.T, i, j, 1, time)
    ρₐ = stateindex(atmos_state.ρ, i, j, 1, time)
    Tᵣ = ratio.reference_temperature
    q₀ = ratio.q₀
    c₁ = ratio.c₁
    c₂ = ratio.c₂

    return q₀ * c₁ * exp(-c₂ / (Tₒ + Tᵣ))
end

@inline function air_sea_difference(i, j, grid, time, diff::SpecificHumidity, atmos_state, ocean_state)
    vapor_fraction = diff.saturation_specific_humidity 
    qₐ = stateindex(atmos_state.q, i, j, 1, time)
    qₛ = saturation_specific_humidity(i, j, grid, time, vapor_fraction, atmos_state, ocean_state)
    return qₐ - qₛ
end

struct LatentHeat{Q, FT}
    specific_humidity_difference :: Q
    vaporization_enthalpy :: FT
end

"""
    LatentHeat(FT = Float64;
               vaporization_enthalpy = 2.5e3 # J / g
               specific_humidity_difference = SpecificHumidity(FT))

"""
function LatentHeat(FT = Float64;
                    vaporization_enthalpy = 2.5e3, # J / g
                    specific_humidity_difference = SpecificHumidity(FT))

    vaporization_enthalpy = convert(FT, vaporization_enthalpy)
    return LatentHeat(specific_humidity_difference, vaporization_enthalpy)
end

@inline function air_sea_difference(i, j, grid, time, diff::LatentHeat, atmos, ocean)
    Δq = air_sea_difference(i, j, grid, time, diff.specific_humidity_difference, atmos, ocean)
    Λᵥ = diff.vaporization_enthalpy
    return Λᵥ * Δq
end

#####
##### Bulk velocity scales
#####

#####
##### Convenience containers for surface fluxes
##### 
##### "Cross realm fluxes" can refer to the flux _data_ (ie, fields representing
##### the total flux for a given variable), or to the flux _components_ / formula.
#####

struct CrossRealmFluxes{M, H, T}
    momentum :: M
    heat :: H
    tracers :: T
end

CrossRealmFluxes(; momentum=nothing, heat=nothing, tracers=nothing) =
    CrossRealmFluxes(momentum, heat, tracers)

Base.summary(osf::CrossRealmFluxes) = "CrossRealmFluxes"
Base.show(io::IO, osf::CrossRealmFluxes) = print(io, summary(osf))

# struct AtmosphereOnlyVelocityScale end
struct RelativeVelocityScale end

@inline function bulk_velocity_scaleᶠᶜᶜ(i, j, grid, time, ::RelativeVelocityScale, atmos_state, ocean_state)
    uₐ = atmos_state.u
    vₐ = atmos_state.v
    uₒ = ocean_state.u
    vₒ = ocean_state.v
    Δu = stateindex(uₐ, i, j, 1, time) - stateindex(uₒ, i, j, 1, time)
    Δv² = ℑxyᶠᶜᵃ(i, j, 1, grid, Δϕt², vₐ, vₒ, time)
    return sqrt(Δu^2 + Δv²)
end

@inline function bulk_velocity_scaleᶜᶠᶜ(i, j, grid, time, ::RelativeVelocityScale, atmos_state, ocean_state)
    uₐ = atmos_state.u
    vₐ = atmos_state.v
    uₒ = ocean_state.u
    vₒ = ocean_state.v
    Δu² = ℑxyᶜᶠᵃ(i, j, 1, grid, Δϕt², uₐ, uₒ, time)
    Δv = stateindex(vₐ, i, j, 1, time) - stateindex(vₒ, i, j, 1, time)
    return sqrt(Δu² + Δv^2)
end

@inline function bulk_velocity_scaleᶜᶜᶜ(i, j, grid, time, ::RelativeVelocityScale, atmos_state, ocean_state)
    uₐ = atmos_state.u
    vₐ = atmos_state.v
    uₒ = ocean_state.u
    vₒ = ocean_state.v
    Δu² = ℑxᶜᵃᵃ(i, j, 1, grid, Δϕt², uₐ, uₒ, time)
    Δv² = ℑyᵃᶜᵃ(i, j, 1, grid, Δϕt², vₐ, vₒ, time)
    return sqrt(Δu² + Δv²)
end


=#

