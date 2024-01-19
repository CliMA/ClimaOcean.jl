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

using Oceananigans.Operators: ℑxᶜᵃᵃ, ℑyᵃᶜᵃ, ℑxᶠᵃᵃ, ℑyᵃᶠᵃ

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
    # ocean_temperature_units
    # freshwater_density ?
end

# Possible units for temperature and salinity
struct DegreesCelsius end
struct DegreesKelvin end

const celsius_to_kelvin = 273.15
@inline convert_to_kelvin(::DegreesCelsius, T::FT) where FT = T + convert(FT, celsius_to_kelvin)
@inline convert_to_kelvin(::DegreesKelvin, T) = T

Base.summary(crf::OceanSeaIceSurfaceFluxes) = "OceanSeaIceSurfaceFluxes"
Base.show(io::IO, crf::OceanSeaIceSurfaceFluxes) = print(io, summary(crf))

function OceanSeaIceSurfaceFluxes(ocean, sea_ice=nothing;
                                  atmosphere = nothing,
                                  radiation = nothing,
                                  ocean_reference_density = reference_density(ocean),
                                  ocean_heat_capacity = heat_capacity(ocean))

    grid = ocean.model.grid
    FT = eltype(grid)

    ocean_reference_density = convert(FT, ocean_reference_density)
    ocean_heat_capacity = convert(FT, ocean_heat_capacity)

    # It's the "thermodynamics gravitational acceleration"
    # (as opposed to the one used for the free surface)
    gravitational_acceleration = ocean.model.buoyancy.model.gravitational_acceleration
    turbulent_fluxes = SimilarityTheoryTurbulentFluxes(grid; gravitational_acceleration)

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

    launch!(arch, grid, :xy, reconstruct_momentum_fluxes!,
            grid, staggered_velocity_fluxes, centered_velocity_fluxes)

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
        Tₒ = ocean_state.T[i, j, 1] + 273.15 # K
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
    water_mole_fraction = turbulent_fluxes.water_mole_fraction
    water_vapor_saturation = turbulent_fluxes.water_vapor_saturation
    q★ = seawater_saturation_specific_humidity(ℂ, Tₒ, Sₒ, ψₐ,
                                               water_mole_fraction,
                                               water_vapor_saturation,
                                               surface_type)
    
    # Thermodynamic and dynamic state at the surface
    ψ₀ = thermodynamic_surface_state = AtmosphericThermodynamics.PhaseEquil_pTq(ℂ, pₐ, Tₒ, q★)
    Ψ₀ = dynamic_surface_state = SurfaceFluxes.StateValues(zero(grid), Uₒ, ψ₀)

    # Thermodynamic and dynamic state at reference level h above the surface
    h = atmosphere_reference_height # elevation of atmos variables relative to surface
    Ψₐ = dynamic_atmos_state = SurfaceFluxes.StateValues(h, Uₐ, ψₐ)

    # Roughness lengths...
    zᵐ = convert(eltype(grid), 1e-2)
    zʰ = convert(eltype(grid), 1e-2)

    Uᵍ = zero(grid) # gustiness
    β = one(grid)   # surface "resistance"
    values = SurfaceFluxes.ValuesOnly(Ψₐ, Ψ₀, zᵐ, zʰ, Uᵍ, β)
    conditions = SurfaceFluxes.surface_conditions(turbulent_fluxes, values)
    update_turbulent_flux_fields!(turbulent_fluxes.fields, i, j, grid, conditions)
    
    # Compute heat fluxes, bulk flux first
    Qd = net_downwelling_radiation(i, j, grid, time, downwelling_radiation, radiation_properties)
    Qu = net_upwelling_radiation(i, j, grid, time, radiation_properties, ocean_state)
    Qs = conditions.shf # sensible heat flux
    Qℓ = conditions.lhf # latent heat flux
    ΣQ = Qd + Qu + Qs + Qℓ

    # Accumulate freshwater fluxes. Rain, snow, runoff -- all freshwater.
    M = cross_realm_flux(i, j, grid, time, prescribed_freshwater_flux) # mass fluxes apparently?
    ρᶠ = 1000 # density of freshwater?
    ΣF = M / ρᶠ # convert from a mass flux to a volume flux / velocity

    # Apparently, conditions.evaporation is a mass flux of water.
    # So, we divide by the density of freshwater.
    E = - conditions.evaporation / ρᶠ # ?

    # Clip evaporation. TODO: figure out why this is needed.
    E = min(zero(E), E) # why does this happen?
    ΣF += E

    # Compute fluxes for u, v, T, S from momentum, heat, and freshwater fluxes
    Jᵘ = centered_velocity_fluxes.u
    Jᵛ = centered_velocity_fluxes.v
    Jᵀ = net_tracer_fluxes.T
    Jˢ = net_tracer_fluxes.S

    ρₒ = ocean_reference_density
    cₒ = ocean_heat_capacity

    atmos_ocean_Jᵘ = conditions.ρτxz / ρₒ
    atmos_ocean_Jᵛ = conditions.ρτyz / ρₒ
    atmos_ocean_Jᵀ = ΣQ / (ρₒ * cₒ)
    atmos_ocean_Jˢ = Sₒ * ΣF

    # Mask fluxes over land for convenience
    kᴺ = size(grid, 3) # index of the top ocean cell
    inactive = inactive_node(i, j, kᴺ, grid, c, c, c)

    @inbounds begin
        Jᵘ[i, j, 1] = ifelse(inactive, zero(grid), atmos_ocean_Jᵘ)
        Jᵛ[i, j, 1] = ifelse(inactive, zero(grid), atmos_ocean_Jᵛ)
        Jᵀ[i, j, 1] = ifelse(inactive, zero(grid), atmos_ocean_Jᵀ)
        Jˢ[i, j, 1] = ifelse(inactive, zero(grid), atmos_ocean_Jˢ)
    end
end

@kernel function reconstruct_momentum_fluxes!(grid, J, Jᶜᶜᶜ)
    i, j = @index(Global, NTuple)

    @inbounds begin
        J.u[i, j, 1] = ℑxᶠᵃᵃ(i, j, 1, grid, Jᶜᶜᶜ.u)
        J.v[i, j, 1] = ℑyᵃᶠᵃ(i, j, 1, grid, Jᶜᶜᶜ.v)
    end
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
    Tₒ = @inbounds ocean_state.T[i, j, 1] + 273.15 # K

    # Note: positive implies _upward_ heat flux, and therefore cooling.
    return σ * ϵ * Tₒ^4
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

