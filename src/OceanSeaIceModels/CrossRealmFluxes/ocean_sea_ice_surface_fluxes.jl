using StaticArrays
using Thermodynamics
using SurfaceFluxes

using ..OceanSeaIceModels: reference_density,
                           heat_capacity,
                           sea_ice_concentration,
                           downwelling_radiation,
                           freshwater_flux

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

struct OceanSeaIceSurfaceFluxes{T, P, C, R, PI, PC, FT, UN, ATM}
    turbulent :: T
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
    freshwater_density :: FT
    ocean_temperature_units :: UN
    # Scratch space to store the atmosphere state at the surface 
    # interpolated to the ocean grid
    surface_atmosphere_state :: ATM
end

# Possible units for temperature and salinity
struct DegreesCelsius end
struct DegreesKelvin end

const celsius_to_kelvin = 273.15
@inline convert_to_kelvin(::DegreesCelsius, T::FT) where FT = T + convert(FT, celsius_to_kelvin)
@inline convert_to_kelvin(::DegreesKelvin, T) = T

Base.summary(crf::OceanSeaIceSurfaceFluxes) = "OceanSeaIceSurfaceFluxes"

function Base.show(io::IO, crf::OceanSeaIceSurfaceFluxes)
    print(io, summary(crf), "\n")
    print(io, "├── radiation: ", summary(crf.radiation), "\n")
    print(io, "└── turbulent coefficients: ", summary(crf.turbulent), "\n")
    return nothing
end

const SeaIceSimulation = Simulation{<:SeaIceModel}

function OceanSeaIceSurfaceFluxes(ocean, sea_ice=nothing;
                                  atmosphere = nothing,
                                  radiation = nothing,
                                  freshwater_density = 1000,
                                  ocean_temperature_units = DegreesCelsius(),
                                  similarity_theory = nothing,
                                  ocean_reference_density = reference_density(ocean),
                                  ocean_heat_capacity = heat_capacity(ocean))

    ocean_grid = ocean.model.grid
    FT = eltype(ocean_grid)

    ocean_reference_density = convert(FT, ocean_reference_density)
    ocean_heat_capacity = convert(FT, ocean_heat_capacity)
    freshwater_density = convert(FT, freshwater_density)

    if !isnothing(atmosphere)
        # It's the "thermodynamics gravitational acceleration"
        # (as opposed to the one used for the free surface)
        gravitational_acceleration = ocean.model.buoyancy.formulation.gravitational_acceleration

        if isnothing(similarity_theory)
            similarity_theory = SimilarityTheoryTurbulentFluxes(ocean_grid; gravitational_acceleration)
        end
    end

    prescribed_fluxes = nothing

    if sea_ice isa SeaIceSimulation
        previous_ice_thickness = deepcopy(sea_ice.model.ice_thickness)
        previous_ice_concentration = deepcopy(sea_ice.model.ice_concentration)
    else
        previous_ice_thickness = nothing
        previous_ice_concentration = nothing
    end

    total_ocean_fluxes = ocean_model_fluxes(ocean.model,
                                            ocean_reference_density,
                                            ocean_heat_capacity)

    total_fluxes = (; ocean=total_ocean_fluxes)

    surface_atmosphere_state = interpolated_surface_atmosphere_state(ocean_grid)

    return OceanSeaIceSurfaceFluxes(similarity_theory,
                                    prescribed_fluxes,
                                    total_fluxes,
                                    radiation,
                                    previous_ice_thickness,
                                    previous_ice_concentration,
                                    ocean_reference_density,
                                    ocean_heat_capacity,
                                    freshwater_density,
                                    ocean_temperature_units,
                                    surface_atmosphere_state)
end

function ocean_model_fluxes(model, ρₒ, cₚ)
    grid = model.grid
    τx = surface_flux(model.velocities.u)
    τy = surface_flux(model.velocities.v)
    τxᶜᶜᶜ = Field{Center, Center, Nothing}(grid)
    τyᶜᶜᶜ = Field{Center, Center, Nothing}(grid)

    ocean_momentum_fluxes = (u    = τx,      # fluxes used in the model
                             v    = τy,      #
                             # Including these (which are only a user convenience, not needed for
                             # time-stepping) incurs about 100s in construction
                             # time for OceanSeaIceSurfaceFluxes:
                             # ρτx  = ρₒ * τx, # momentum fluxes multiplied by reference density
                             # ρτy  = ρₒ * τy, # user convenience 
                             uᶜᶜᶜ = τxᶜᶜᶜ,   # fluxes computed by bulk formula at cell centers
                             vᶜᶜᶜ = τyᶜᶜᶜ)

    tracers = model.tracers
    ocean_tracer_fluxes = NamedTuple(name => surface_flux(tracers[name])
                                     for name in keys(tracers))

    ocean_heat_flux = ρₒ * cₚ * ocean_tracer_fluxes.T
    upwelling_radiation = Field{Center, Center, Nothing}(grid)

    fluxes = (; momentum = ocean_momentum_fluxes,
                tracers = ocean_tracer_fluxes,
                heat = (; total = ocean_heat_flux, upwelling_radiation))
    
    return fluxes
end

function interpolated_surface_atmosphere_state(ocean_grid)
    surface_atmosphere_state = (u  = Field{Center, Center, Nothing}(ocean_grid),
                                v  = Field{Center, Center, Nothing}(ocean_grid),
                                T  = Field{Center, Center, Nothing}(ocean_grid),
                                q  = Field{Center, Center, Nothing}(ocean_grid),
                                p  = Field{Center, Center, Nothing}(ocean_grid),
                                Qs = Field{Center, Center, Nothing}(ocean_grid),
                                Qℓ = Field{Center, Center, Nothing}(ocean_grid),
                                Qu = Field{Center, Center, Nothing}(ocean_grid),
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

