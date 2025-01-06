using ClimaOcean
using Oceananigans
using SpeedyWeather
using ClimaOcean.OceanSeaIceModels.Atmospheres: HeatCapacityParameters, ConstitutiveParameters, PrognosticAtmosphere

import ClimaOcean.OceanSeaIceModels: time_step!

#####
##### Ocean model (idealized quarter degree model)
#####

Nx = 360
Ny = 180
Nz = 40

depth = 5000
φmax  = 75

arch = Oceananigans.CPU()
r_faces = ClimaOcean.exponential_z_faces(; Nz, depth)

# A quarter degree ocean model (idealized)
grid = LatitudeLongitudeGrid(arch, 
                             size=(Nx, Ny, Nz), 
                             latitude=(-φmax, φmax), 
                             longitude=(0, 360),
                             halo = (6, 6, 2),
                             z = r_faces)

ocean = ClimaOcean.ocean_simulation(grid)

# Initial conditions 
# - parabolic temperature profile with a stratification
# - constant salinity

Tᵢ(λ, φ, z) = 30.0 * cosd(φ)^2 * (1 + z / depth)
Sᵢ(λ, φ, z) = 35.0 

Oceananigans.set!(ocean.model, T=Tᵢ, S=Sᵢ)

#####
##### Atmospheric model
#####

spectral_grid = SpectralGrid(trunc=121, nlayers=8, Grid=FullClenshawGrid)
model      = PrimitiveWetModel(spectral_grid)
simulation = initialize!(model)
atmosphere = PrognosticAtmosphere(; simulation)

# Initialize the speedyweather timestepping correctly
#=

=# 

# Extend the PrognosticAtmosphere functions to work with SpeedyWeather
import ClimaOcean.OceanSeaIceModels.Atmospheres: 
                    regrid_fluxes_to_atmospheric_model!, 
                    interpolate_atmospheric_state!

# Out-source the time_step! to the prognostic atmosphere model
function time_step!(atmos::PrognosticAtmosphere{<:Any, <:SpeedyWeather.Simulation}) 
    progn = atmos.simulation.prognostic_variables
    diagn = atmos.simulation.diagnostic_variables
    model = atmos.simulation.model
    
    (; clock) = progn
    (; Δt, Δt_millisec) = model.time_stepping

    (; output, feedback) = model
    
    SpeedyWeather.timestep!(progn, diagn, 2Δt, model) # calculate tendencies and leapfrog forward
    SpeedyWeather.timestep!(clock, Δt_millisec)       # time of lf=2 and diagn after timestep!

    SpeedyWeather.progress!(feedback, progn)          # updates the progress meter bar
    SpeedyWeather.output!(output, progn, diagn, model)
    SpeedyWeather.callback!(model.callbacks, progn, diagn, model)

    return nothing
end

# Interpolate the atmospheric surface fields to the ocean/sea-ice model grid
function interpolate_atmospheric_state!(surface_atmosphere_state, 
                                        interpolated_prescribed_freshwater_flux, 
                                        atmos::PrognosticAtmosphere{<:Any, <:SpeedyWeather.Simulation}, 
                                        grid, clock)
    nothing
end

# Regrid the fluxes from the ocean/sea-ice grid to the atmospheric model grid
function regrid_fluxes_to_atmospheric_model!(atmos::PrognosticAtmosphere{<:Any, <:SpeedyWeather.Simulation}, net_tracer_fluxes, centered_velocity_fluxes)
    nothing
end

#####
##### Coupled model
##### 

radiation   = Radiation() # Need to fix this to be consistent with SpeedyWeather? 
earth_model = OceanSeaIceModel(ocean; atmosphere, radiation)
earth       = ClimaOcean.Simulation(earth_model, Δt=10minutes)

#####
##### Progress function
#####

function progress(earth)

    # Some outputs from speedyweather?


    # Some outputs from oceananigans?


end

#####
##### Run the coupled model
##### 

time_step!(earth) # If no issue we can think about the details