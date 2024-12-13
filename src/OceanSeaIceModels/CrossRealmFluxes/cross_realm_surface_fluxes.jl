using StaticArrays
using Thermodynamics
using SurfaceFluxes

using ..OceanSeaIceModels: reference_density,
                           heat_capacity,
                           sea_ice_concentration,
                           sea_ice_thickness,
                           downwelling_radiation,
                           freshwater_flux,
                           SeaIceSimulation


using ClimaSeaIce: SeaIceModel

using Oceananigans: HydrostaticFreeSurfaceModel, architecture
using Oceananigans.Grids: inactive_node, node
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: ConstantField, interpolate
using Oceananigans.Utils: launch!, Time, KernelParameters

# using Oceananigans.OutputReaders: extract_field_time_series, update_field_time_series!

using Oceananigans.Operators: ℑxᶜᵃᵃ, ℑyᵃᶜᵃ, ℑxᶠᵃᵃ, ℑyᵃᶠᵃ

using KernelAbstractions: @kernel, @index

#####
##### Container for organizing information related to fluxes
#####

struct CrossRealmSurfaceFluxes{T, P, C, R, PI, PC, FT, UN, ATM}
    turbulent :: T # the turbulent fluxes
    prescribed :: P
    # Add `components` which will also store components of the total fluxes
    # (eg latent, sensible heat flux)
    total :: C
    radiation :: R
    previous_ice_thickness :: PI
    previous_ice_concentration :: PC
    # The ocean is Boussinesq, so these are _only_ coupled properties:
    ocean_reference_density :: FT
    ocean_heat_capacity :: FT
    sea_ice_reference_density :: FT
    sea_ice_heat_capacity :: FT
    freshwater_density :: FT
    ocean_temperature_units :: UN
    sea_ice_temperature_units :: UN
    # Scratch space to store the atmosphere state at the surface 
    # interpolated to the ocean grid
    surface_atmosphere_state :: ATM
end

struct TurbulentFluxes{T, FT, C, W, I, M, F}
    thermodynamic_parameters :: T
    gravitational_acceleration :: FT
    coefficients :: C
    water_vapor_saturation :: W    # model for computing the saturation water vapor mass over ocean
    ice_vapor_saturation :: I      # model for computing the saturation water vapor mass over ice
    water_mole_fraction :: M       # mole fraction of H₂O in seawater
    fields :: F                    # fields that store turbulent fluxes
end

struct ClasiusClapyeronSaturation end
 
@inline function water_saturation_specific_humidity(::ClasiusClapyeronSaturation, ℂₐ, ρₛ, Tₛ)
    FT = eltype(ℂₐ)
    p★ = AtmosphericThermodynamics.saturation_vapor_pressure(ℂₐ, convert(FT, Tₛ), Liquid())
    q★ = AtmosphericThermodynamics.q_vap_saturation_from_density(ℂₐ, convert(FT, Tₛ), ρₛ, p★)
    return q★
end

const PATP = PrescribedAtmosphereThermodynamicsParameters

# Possible units for temperature and salinity
struct DegreesCelsius end
struct DegreesKelvin end

const celsius_to_kelvin = 273.15
@inline convert_to_kelvin(::DegreesCelsius, T::FT) where FT = T + convert(FT, celsius_to_kelvin)
@inline convert_to_kelvin(::DegreesKelvin, T) = T

Base.summary(crf::CrossRealmSurfaceFluxes) = "CrossRealmSurfaceFluxes"
Base.show(io::IO, crf::CrossRealmSurfaceFluxes) = print(io, summary(crf))

"""
    We need a docstring...

- `thermodynamics_parameters`: The thermodynamics parameters used to calculate atmospheric stability and
                               saturation pressure. Default: `PATP(FT)`, alias for `PrescribedAtmosphereThermodynamicsParameters`.
- `water_vapor_saturation`: The water vapor saturation law. Default: `ClasiusClapyeronSaturation()` that follows the 
                            Clasius-Clapyeron pressure formulation.
- `water_mole_fraction`: The water mole fraction used to calculate the `seawater_saturation_specific_humidity`. 
                         Default: 0.98, the rest is assumed to be other substances such as chlorine, sodium sulfide, and magnesium.
"""
function CrossRealmSurfaceFluxes(ocean, sea_ice=nothing;
                                 atmosphere = nothing,
                                 radiation = nothing,
                                 freshwater_density = 1000,
                                 prescribed_fluxes = nothing, # ?? Is this ever used ??
                                 ocean_temperature_units = DegreesCelsius(),
                                 turbulent_fluxes = nothing,
                                 water_vapor_saturation = ClasiusClapyeronSaturation(),
                                 ice_vapor_saturation = ClasiusClapyeronSaturation(),
                                 water_mole_fraction = 0.98,
                                 thermodynamics_parameters = nothing,
                                 ocean_reference_density = reference_density(ocean),
                                 ocean_heat_capacity = heat_capacity(ocean),
                                 sea_ice_reference_density = reference_density(sea_ice),
                                 sea_ice_heat_capacity = reference_density(sea_ice))

    ocean_grid = ocean.model.grid
    FT = eltype(ocean_grid)

    ocean_reference_density = convert(FT, ocean_reference_density)
    ocean_heat_capacity = convert(FT, ocean_heat_capacity)
    sea_ice_reference_density = isnothing(sea_ice_reference_density) ? nothing : convert(FT, sea_ice_reference_density)
    sea_ice_heat_capacity = isnothing(sea_ice_heat_capacity) ? nothing : convert(FT, sea_ice_heat_capacity)
    freshwater_density = convert(FT, freshwater_density)
    water_mole_fraction = convert(FT, water_mole_fraction)
    
    if !isnothing(atmosphere)
        if isnothing(thermodynamics_parameters)
            thermodynamics_parameters = PATP(FT)
        end

        # It's the "thermodynamics gravitational acceleration"
        # (as opposed to the one used for the free surface)
        gravitational_acceleration = ocean.model.buoyancy.formulation.gravitational_acceleration

        # Build turbulent fluxes if they do not exist
        if isnothing(turbulent_fluxes)
            ocean_fluxes = SimilarityTheoryFluxes()
            sea_ice_fluxes = if sea_ice isa SeaIceSimulation
                SimilarityTheoryFluxes()
            else
                nothing
            end
            coefficients = (ocean=ocean_fluxes, sea_ice=sea_ice_fluxes)
            ocean_fields = turbulent_fluxes_fields(ocean_grid)
            sea_ice_fields = if sea_ice isa SeaIceSimulation
                turbulent_fluxes_fields(sea_ice.model.grid)
            else
                nothing
            end

            fluxes_fields = (ocean=ocean_fields, sea_ice=sea_ice_fields)
            
            turbulent_fluxes = TurbulentFluxes(thermodynamics_parameters,
                                               gravitational_acceleration,
                                               coefficients,
                                               water_vapor_saturation,
                                               ice_vapor_saturation,
                                               water_mole_fraction,
                                               fluxes_fields)
        end
    end

    if sea_ice isa SeaIceSimulation
        previous_ice_thickness = deepcopy(sea_ice.model.ice_thickness)
        previous_ice_concentration = deepcopy(sea_ice.model.ice_concentration)
    else
        previous_ice_thickness = nothing
        previous_ice_concentration = nothing
    end

    total_ocean_fluxes = surface_model_fluxes(ocean.model,
                                              ocean_reference_density,
                                              ocean_heat_capacity)

                                        
    total_sea_ice_fluxes = if sea_ice isa SeaIceSimulation
        surface_model_fluxes(sea_ice.model,
                             sea_ice_reference_density,
                             sea_ice_heat_capacity)
    else
        nothing
    end

    # The actual fields in the boundary conditions of 
    # model-specific velocities and tracers
    total_fluxes = (; ocean=total_ocean_fluxes,
                    sea_ice=total_sea_ice_fluxes)

    surface_atmosphere_state = interpolated_surface_atmosphere_state(ocean_grid)

    return CrossRealmSurfaceFluxes(turbulent_fluxes,
                                   prescribed_fluxes,
                                   total_fluxes,
                                   radiation,
                                   previous_ice_thickness,
                                   previous_ice_concentration,
                                   ocean_reference_density,
                                   ocean_heat_capacity,
                                   sea_ice_reference_density,
                                   sea_ice_heat_capacity,
                                   freshwater_density,
                                   ocean_temperature_units,
                                   ocean_temperature_units,
                                   surface_atmosphere_state)
end

function surface_model_fluxes(model, ρₛ, cₛ)
    grid = model.grid
    τx = surface_flux(model.velocities.u)
    τy = surface_flux(model.velocities.v)
    τxᶜᶜᶜ = Field{Center, Center, Nothing}(grid)
    τyᶜᶜᶜ = Field{Center, Center, Nothing}(grid)

   surface_momentum_fluxes = (u    = τx,      # fluxes used in the model
                              v    = τy,      #
                              # Including these (which are only a user convenience, not needed for
                              # time-stepping) incurs about 100s in construction
                              # time for CrossRealmSurfaceFluxes:
                              # ρτx  = ρₒ * τx, # momentum fluxes multiplied by reference density
                              # ρτy  = ρₒ * τy, # user convenience 
                              uᶜᶜᶜ = τxᶜᶜᶜ,   # fluxes computed by bulk formula at cell centers
                              vᶜᶜᶜ = τyᶜᶜᶜ)

    tracers = model.tracers
    surface_tracer_fluxes = NamedTuple(name => surface_flux(tracers[name])
                                     for name in keys(tracers))

    surface_heat_flux = ρₛ * cₛ * surface_tracer_fluxes.T

    fluxes = (momentum = surface_momentum_fluxes,
              tracers = surface_tracer_fluxes,
              heat = surface_heat_flux)

    return fluxes
end

function turbulent_fluxes_fields(grid)
    water_vapor   = Field{Center, Center, Nothing}(grid)
    latent_heat   = Field{Center, Center, Nothing}(grid)
    sensible_heat = Field{Center, Center, Nothing}(grid)
    x_momentum    = Field{Center, Center, Nothing}(grid)
    y_momentum    = Field{Center, Center, Nothing}(grid)
    T_surface     = Field{Center, Center, Nothing}(grid)

    return (; latent_heat, sensible_heat, water_vapor, x_momentum, y_momentum, T_surface)
end

function interpolated_surface_atmosphere_state(ocean_grid)
    surface_atmosphere_state = (u  = Field{Center, Center, Nothing}(ocean_grid),
                                v  = Field{Center, Center, Nothing}(ocean_grid),
                                T  = Field{Center, Center, Nothing}(ocean_grid),
                                q  = Field{Center, Center, Nothing}(ocean_grid),
                                p  = Field{Center, Center, Nothing}(ocean_grid),
                                Qs = Field{Center, Center, Nothing}(ocean_grid),
                                Qℓ = Field{Center, Center, Nothing}(ocean_grid),
                                Mp = Field{Center, Center, Nothing}(ocean_grid))

    return surface_atmosphere_state
end
    
#####
##### Utility for interpolating tuples of fields
#####

# Note: assumes loc = (c, c, nothing) (and the third location should
# not matter.)
@inline interp_atmos_time_series(J, X, time, grid, args...) =
    interpolate(X, time, J, (c, c, nothing), grid, args...)

@inline interp_atmos_time_series(ΣJ::NamedTuple, args...) =
    interp_atmos_time_series(values(ΣJ), args...)

@inline interp_atmos_time_series(ΣJ::Tuple{<:Any}, args...) =
    interp_atmos_time_series(ΣJ[1], args...) +
    interp_atmos_time_series(ΣJ[2], args...)

@inline interp_atmos_time_series(ΣJ::Tuple{<:Any, <:Any}, args...) =
    interp_atmos_time_series(ΣJ[1], args...) +
    interp_atmos_time_series(ΣJ[2], args...)

@inline interp_atmos_time_series(ΣJ::Tuple{<:Any, <:Any, <:Any}, args...) =
    interp_atmos_time_series(ΣJ[1], args...) +
    interp_atmos_time_series(ΣJ[2], args...) +
    interp_atmos_time_series(ΣJ[3], args...)

@inline interp_atmos_time_series(ΣJ::Tuple{<:Any, <:Any, <:Any, <:Any}, args...) =
    interp_atmos_time_series(ΣJ[1], args...) +
    interp_atmos_time_series(ΣJ[2], args...) +
    interp_atmos_time_series(ΣJ[3], args...) +
    interp_atmos_time_series(ΣJ[4], args...)

