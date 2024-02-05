using StaticArrays
using Thermodynamics
using SurfaceFluxes

using ..OceanSeaIceModels: reference_density,
                           heat_capacity,
                           sea_ice_concentration,
                           downwelling_radiation,
                           freshwater_flux

using ClimaSeaIce: SlabSeaIceModel

using Oceananigans: HydrostaticFreeSurfaceModel, architecture
using Oceananigans.Grids: inactive_node, node
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: ConstantField, interpolate
using Oceananigans.Utils: launch!, Time

using Oceananigans.Operators: ℑxᶜᵃᵃ, ℑyᵃᶜᵃ, ℑxᶠᵃᵃ, ℑyᵃᶠᵃ

using KernelAbstractions: @kernel, @index

#####
##### Container for organizing information related to fluxes
#####

struct OceanSeaIceSurfaceFluxes{T, P, C, R, PI, PC, FT, UN}
    turbulent :: T
    prescribed :: P
    total :: C
    radiation :: R
    previous_ice_thickness :: PI
    previous_ice_concentration :: PC
    # The ocean is Boussinesq, so these are _only_ coupled properties:
    ocean_reference_density :: FT
    ocean_heat_capacity :: FT
    freshwater_density :: FT
    ocean_temperature_units :: UN
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
                                  freshwater_density = 1000,
                                  ocean_temperature_units = DegreesCelsius(),
                                  ocean_reference_density = reference_density(ocean),
                                  ocean_heat_capacity = heat_capacity(ocean))

    grid = ocean.model.grid
    FT = eltype(grid)

    ocean_reference_density = convert(FT, ocean_reference_density)
    ocean_heat_capacity = convert(FT, ocean_heat_capacity)
    freshwater_density = convert(FT, freshwater_density)

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
                                    ocean_heat_capacity,
                                    freshwater_density,
                                    ocean_temperature_units)
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
    atmosphere_grid = atmosphere.grid

    # Basic model properties
    grid = ocean.model.grid
    arch = architecture(grid)
    clock = ocean.model.clock

    # Ocean, atmosphere, and sea ice state
    ocean_velocities  = surface_velocities(ocean)
    ocean_tracers     = surface_tracers(ocean)
    ice_concentration = sea_ice_concentration(sea_ice)

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
    atmosphere_state = merge(atmosphere.velocities, atmosphere.tracers, (; p=atmosphere.pressure))

    launch!(arch, grid, :xy, compute_atmosphere_ocean_turbulent_fluxes!,
            grid, clock,
            centered_velocity_fluxes,
            net_tracer_fluxes,
            turbulent_fluxes,
            atmosphere.freshwater_flux,
            atmosphere.downwelling_radiation,
            radiation_properties,
            ocean_state,
            atmosphere_state,
            atmosphere_grid,
            atmosphere.reference_height, # height at which the state is known
            atmosphere.thermodynamics_parameters,
            coupled_model.fluxes.ocean_reference_density,
            coupled_model.fluxes.ocean_heat_capacity,
            coupled_model.fluxes.freshwater_density,
            coupled_model.fluxes.ocean_temperature_units,
            ice_concentration)

    # Note: I think this can be avoided if we modify the preceding kernel
    # to compute from 0:Nx+1, ie in halo regions
    fill_halo_regions!(centered_velocity_fluxes)

    launch!(arch, grid, :xy, reconstruct_momentum_fluxes!,
            grid, staggered_velocity_fluxes, centered_velocity_fluxes)

    return nothing
end

const c = Center()

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
                                                            atmos_grid,
                                                            atmosphere_reference_height,
                                                            atmosphere_thermodynamics_parameters,
                                                            ocean_reference_density,
                                                            ocean_heat_capacity,
                                                            freshwater_density,
                                                            ocean_temperature_units,
                                                            ice_concentration)

    i, j = @index(Global, NTuple)

    time = Time(clock.time)

    # Extract state variables at cell centers
    @inbounds begin
        # Ocean state
        uₒ = ℑxᶜᵃᵃ(i, j, 1, grid, ocean_state.u)
        vₒ = ℑyᵃᶜᵃ(i, j, 1, grid, ocean_state.v)
        Tₒ = ocean_state.T[i, j, 1]
        Tₒ = convert_to_kelvin(ocean_temperature_units, Tₒ)
        Sₒ = ocean_state.S[i, j, 1]

        # Atmos state
        X = node(i, j, 1, grid, c, c, c)

        uₐ = interpolate_atmos_field_time_series(atmos_state.u, X, time, atmos_grid)
        vₐ = interpolate_atmos_field_time_series(atmos_state.v, X, time, atmos_grid)

        Tₐ = interpolate_atmos_field_time_series(atmos_state.T, X, time, atmos_grid)
        pₐ = interpolate_atmos_field_time_series(atmos_state.p, X, time, atmos_grid)
        qₐ = interpolate_atmos_field_time_series(atmos_state.q, X, time, atmos_grid)

        Qs = interpolate_atmos_field_time_series(downwelling_radiation.shortwave, X, time, atmos_grid)
        Qℓ = interpolate_atmos_field_time_series(downwelling_radiation.longwave,  X, time, atmos_grid)

        # Accumulate freshwater mass fluxes. Rain, snow, runoff -- all freshwater.
        M = interpolate_atmos_field_time_series(prescribed_freshwater_flux, X, time, atmos_grid)
    end

    # Build thermodynamic and dynamic states in the atmosphere and surface.
    # Notation:
    #   ⋅ ϕ ≡ thermodynamic state vector
    #   ⋅ Φ ≡ "dynamic" state vector (thermodynamics + reference height + velocity)
    ℂₐ = atmosphere_thermodynamics_parameters
    ϕₐ = thermodynamic_atmospheric_state = AtmosphericThermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₐ, qₐ)

    hₐ = atmosphere_reference_height # elevation of atmos variables relative to surface
    Uₐ = SVector(uₐ, vₐ)
    Φₐ = dynamic_atmos_state = SurfaceFluxes.StateValues(hₐ, Uₐ, ϕₐ)

    # Build surface state with saturated specific humidity
    surface_type = AtmosphericThermodynamics.Liquid()
    qₒ = seawater_saturation_specific_humidity(ℂₐ, Tₒ, Sₒ, ϕₐ,
                                               turbulent_fluxes.water_mole_fraction,
                                               turbulent_fluxes.water_vapor_saturation,
                                               surface_type)
    
    # Thermodynamic and dynamic surface state
    ϕ₀ = thermodynamic_surface_state = AtmosphericThermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₒ, qₒ)

    h₀ = zero(grid) # surface height
    Uₒ = SVector(uₒ, vₒ)
    Φ₀ = dynamic_surface_state = SurfaceFluxes.StateValues(h₀, Uₒ, ϕ₀)

    f = turbulent_fluxes.similarity_function
    ℓ = turbulent_fluxes.roughness_lengths
    fluxes = compute_turbulent_surface_fluxes(f, ℓ, turbulent_fluxes, Φₐ, Φ₀)

    #=
    # Initial guess for the roughness length.
    FT = eltype(grid)
    zᵐ = zʰ = convert(FT, 5e-4) # τ = 0.3 => u★ = sqrt(τ / ρₐ) ~ z₀ ~ 5e-4

    # Solve for the surface fluxes with initial roughness length guess
    Uᵍ = zero(grid) # gustiness
    β = one(grid)   # surface "resistance"
    values = SurfaceFluxes.ValuesOnly(Φₐ, Φ₀, zᵐ, zʰ, Uᵍ, β)
    conditions = SurfaceFluxes.surface_conditions(turbulent_fluxes, values)

    # It's like a fixed point iteration
    g = turbulent_fluxes.gravitational_acceleration
    α = convert(FT, 0.011) # Charnock parameter
    u★ = conditions.ustar
    zᵐ = zʰ = α * u★^2 / g
    values = SurfaceFluxes.ValuesOnly(Φₐ, Φ₀, zᵐ, zʰ, Uᵍ, β)
    conditions = SurfaceFluxes.surface_conditions(turbulent_fluxes, values)
    Qc = conditions.shf # sensible or "conductive" heat flux
    Qe = conditions.lhf # latent or "evaporative" heat flux
    =#

    Qc = fluxes.sensible_heat
    Qe = fluxes.latent_heat
    Mt = fluxes.freshwater
    τˣ = fluxes.zonal_momentum
    τʸ = fluxes.meridional_momentum
    
    # Compute heat fluxes, bulk flux first
    Qd = net_downwelling_radiation(i, j, grid, time, Qs, Qℓ, radiation_properties)
    Qu = net_upwelling_radiation(i, j, grid, time, radiation_properties, ocean_state, ocean_temperature_units)
    ΣQ = Qd + Qu + Qc + Qe

    # Convert from a mass flux to a volume flux (aka velocity).
    # Also switch the sign, for some reason we are given freshwater flux as positive down.
    ρᶠ = freshwater_density
    ΣF = - M / ρᶠ

    # Mt is the turbulent freshwater (moisture) flux, which we convert here
    # to a volume flux by dividing by the density of freshwater.
    Ft = Mt / ρᶠ
    ΣF += Ft

    update_turbulent_flux_fields!(turbulent_fluxes.fields, i, j, grid, fluxes, ρᶠ)

    # Compute fluxes for u, v, T, S from momentum, heat, and freshwater fluxes
    Jᵘ = centered_velocity_fluxes.u
    Jᵛ = centered_velocity_fluxes.v
    Jᵀ = net_tracer_fluxes.T
    Jˢ = net_tracer_fluxes.S

    ρₒ = ocean_reference_density
    cₒ = ocean_heat_capacity

    #atmos_ocean_Jᵘ = conditions.ρτxz / ρₒ
    #atmos_ocean_Jᵛ = conditions.ρτyz / ρₒ
    atmos_ocean_Jᵘ = τˣ / ρₒ
    atmos_ocean_Jᵛ = τʸ / ρₒ
    atmos_ocean_Jᵀ = ΣQ / (ρₒ * cₒ)
    atmos_ocean_Jˢ = - Sₒ * ΣF

    # Mask fluxes over land for convenience
    kᴺ = size(grid, 3) # index of the top ocean cell
    inactive = inactive_node(i, j, kᴺ, grid, c, c, c)

    @inbounds begin
        Jᵘ[i, j, 1] = ifelse(inactive, 0, atmos_ocean_Jᵘ)
        Jᵛ[i, j, 1] = ifelse(inactive, 0, atmos_ocean_Jᵛ)
        Jᵀ[i, j, 1] = ifelse(inactive, 0, atmos_ocean_Jᵀ)
        Jˢ[i, j, 1] = ifelse(inactive, 0, atmos_ocean_Jˢ)
    end
end

@kernel function reconstruct_momentum_fluxes!(grid, J, Jᶜᶜᶜ)
    i, j = @index(Global, NTuple)

    @inbounds begin
        J.u[i, j, 1] = ℑxᶠᵃᵃ(i, j, 1, grid, Jᶜᶜᶜ.u)
        J.v[i, j, 1] = ℑyᵃᶠᵃ(i, j, 1, grid, Jᶜᶜᶜ.v)
    end
end

@inline function net_downwelling_radiation(i, j, grid, time, Qs, Qℓ, radiation)
    α = stateindex(radiation.reflection.ocean, i, j, 1, time)

    return @inbounds - (1 - α) * Qs - Qℓ
end

@inline function net_upwelling_radiation(i, j, grid, time, radiation, ocean_state, ocean_temperature_units)
    σ = radiation.stefan_boltzmann_constant
    ϵ = stateindex(radiation.emission.ocean, i, j, 1, time)

    # Ocean surface temperature (departure from reference, typically in ᵒC)
    Tₒ = @inbounds ocean_state.T[i, j, 1]
    Tₒ = convert_to_kelvin(ocean_temperature_units, Tₒ)

    # Note: positive implies _upward_ heat flux, and therefore cooling.
    return ϵ * σ * Tₒ^4
end

#####
##### Utility for interpolating tuples of fields
#####

# Note: assumes loc = (c, c, c)
@inline interpolate_atmos_field_time_series(J, x, t, grid) =
    interpolate(x, t, J, (c, c, c), grid)

@inline interpolate_atmos_field_time_series(ΣJ::NamedTuple, args...) =
    interpolate_atmos_field_time_series(values(ΣJ), args...)

@inline interpolate_atmos_field_time_series(ΣJ::Tuple{<:Any}, args...) =
    interpolate_atmos_field_time_series(ΣJ[1], args...) +
    interpolate_atmos_field_time_series(ΣJ[2], args...)

@inline interpolate_atmos_field_time_series(ΣJ::Tuple{<:Any, <:Any}, args...) =
    interpolate_atmos_field_time_series(ΣJ[1], args...) +
    interpolate_atmos_field_time_series(ΣJ[2], args...)

@inline interpolate_atmos_field_time_series(ΣJ::Tuple{<:Any, <:Any, <:Any}, args...) =
    interpolate_atmos_field_time_series(ΣJ[1], args...) +
    interpolate_atmos_field_time_series(ΣJ[2], args...) +
    interpolate_atmos_field_time_series(ΣJ[3], args...)

@inline interpolate_atmos_field_time_series(ΣJ::Tuple{<:Any, <:Any, <:Any, <:Any}, args...) =
    interpolate_atmos_field_time_series(ΣJ[1], args...) +
    interpolate_atmos_field_time_series(ΣJ[2], args...) +
    interpolate_atmos_field_time_series(ΣJ[3], args...) +
    interpolate_atmos_field_time_series(ΣJ[4], args...)

