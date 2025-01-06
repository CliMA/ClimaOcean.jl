using ClimaOcean
using Oceananigans
using SpeedyWeather

#####
##### Ocean model (idealized quarter degree model)
#####

Nx = 1440
Ny = 600
Nz = 40

depth = 5000
φmax  = 75

arch = CPU()
r_faces = ClimaOcean.exponential_z_faces(; Nz, depth)

# A quarter degree ocean model (idealized)
grid = LatitudeLongitudeGrid(arch, 
                             size=(Nx, Ny, Nz), 
                             latitude=(-φmax, φmax), 
                             longitude=(0, 360),
                             z = r_faces)

ocean = ClimaOcean.ocean_simulation(grid)

# Initial conditions 
# - parabolic temperature profile with a stratification
# - constant salinity

Tᵢ(λ, φ, z) = 30.0 * (1 - 1 / φmax^2 * φ^2) * (1 + z / depth)
Sᵢ(λ, φ, z) = 35.0 

set!(ocean.model, T=Tᵢ, S=Sᵢ)

#####
##### Atmospheric model
#####

spectral_grid = SpectralGrid(trunc=119, nlev=8, Grid=FullClenshawGrid)
model = PrimitiveWetModel(spectral_grid)
simulation = initialize!(model)

# Extend the PrognosticAtmosphere functions to work with SpeedyWeather
import ClimaOcean.OceanSeaIceModels.Atmospheres: 
                    regrid_fluxes_to_atmospheric_model!, 
                    interpolate_atmospheric_state!

# Interpolate the atmospheric surface fields to the ocean/sea-ice model grid
function interpolate_atmospheric_state!(atmos::PrognosticAtmosphere{<:Any, <:Simulation}, surface_atmosphere_state, grid, clock)
    nothing
end

# Regrid the fluxes from the ocean/sea-ice grid to the atmospheric model grid
function regrid_fluxes_to_atmospheric_model!(atmos::PrognosticAtmosphere{<:Any, <:Simulation}, net_tracer_fluxes, centered_velocity_fluxes)
    nothing
end

#####
##### Coupled model
##### 

radiation = Radiation() # ? 

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)