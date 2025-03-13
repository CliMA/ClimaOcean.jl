using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using ClimaOcean.OceanSimulations
using ClimaOcean.ECCO
using ClimaOcean.DataWrangling
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium
using Printf

# A single column grid with one point in the Arctic ocean
grid = LatitudeLongitudeGrid(size = 1, 
                             latitude = 87.0,
                             longitude = 33.0,
                             z = (-10, 0),
                             topology = (Flat, Flat, Bounded))

#####
##### A Prescribed Ocean model
#####

dataset = ECCO4Monthly()
dates   = all_dates(dataset)[1:24]

T = ECCOFieldTimeSeries(:temperature, grid; dates, inpainting=nothing, time_indices_in_memory=24)
S = ECCOFieldTimeSeries(:salinity, grid;    dates, inpainting=nothing, time_indices_in_memory=24)

ocean_model = PrescribedOceanModel((; T, S); grid)
ocean = Simulation(ocean_model, Δt=30minutes, stop_time=730days)

##### 
##### A Prescribed Atmosphere model
#####

atmosphere = JRA55PrescribedAtmosphere(latitude = 87.0,
                                       longitude = 33.0,
                                       backend = InMemory())

#####
##### A Prognostic Sea-ice model with no dynamics
#####

# Remember to pass the SSS as a bottom bc to the sea ice!
SSS = view(ocean.model.tracers.S, :, :, 1)
bottom_heat_boundary_condition = IceWaterThermalEquilibrium(SSS)

sea_ice = sea_ice_simulation(grid; bottom_heat_boundary_condition) 

set!(sea_ice.model, h=Metadata(:sea_ice_thickness;     dataset, dates=first(dates)), 
                    ℵ=Metadata(:sea_ice_concentration; dataset, dates=first(dates)))

#####
##### Radiation of the Ocean and the Sea ice
#####

# to mimick snow, as in Semtner (1975) we use a sea ice albedo that is 0.64 for top temperatures 
# above -0.1 ᵒC and 0.75 for top temperatures below -0.1 ᵒC. We also disable radiation for the ocean.
radiation = Radiation(sea_ice_albedo=0.7, ocean_albedo=1, ocean_emissivity=0)

#####
##### Arctic coupled model
#####

arctic = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
arctic = Simulation(arctic, Δt=30minutes, stop_time=730days)

# Sea-ice variables
h = sea_ice.model.ice_thickness
ℵ = sea_ice.model.ice_concentration

# Fluxes
Tu = arctic.model.interfaces.atmosphere_sea_ice_interface.temperature
Qˡ = arctic.model.interfaces.atmosphere_sea_ice_interface.fluxes.latent_heat
Qˢ = arctic.model.interfaces.atmosphere_sea_ice_interface.fluxes.sensible_heat
Qⁱ = arctic.model.interfaces.sea_ice_ocean_interface.fluxes.interface_heat
Qᶠ = arctic.model.interfaces.sea_ice_ocean_interface.fluxes.frazil_heat
Qᵗ = arctic.model.interfaces.net_fluxes.sea_ice_top.heat
Qᴮ = arctic.model.interfaces.net_fluxes.sea_ice_bottom.heat

# Output writers
arctic.output_writers[:vars] = JLD2OutputWriter(sea_ice.model, (; h, ℵ, Tu, Qˡ, Qˢ, Qⁱ, Qᶠ, Qᵗ, Qᴮ),
                                                 filename = "sea_ice_quantities.jld2",
                                                 schedule = IterationInterval(12),
                                                 overwrite_existing=true)

arctic.output_writers[:averages] = JLD2OutputWriter(sea_ice.model, (; h, ℵ, Tu, Qˡ, Qˢ, Qⁱ, Qᶠ, Qᵗ, Qᴮ),
                                                    filename = "averaged_sea_ice_quantities.jld2",
                                                    schedule = AveragedTimeInterval(1days),
                                                    overwrite_existing=true)

wall_time = Ref(time_ns())

using Statistics

function progress(sim)
    sea_ice = sim.model.sea_ice
    h = first(sea_ice.model.ice_thickness)
    ℵ = first(sea_ice.model.ice_concentration)
    T = first(sim.model.interfaces.atmosphere_sea_ice_interface.temperature)

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg1 = @sprintf("time: %s, iteration: %d, Δt: %s, ", prettytime(sim), iteration(sim), prettytime(sim.Δt))
    msg2 = @sprintf("h: %.2e m, ℵ: %.2e , T: %.2e ᴼC", h, ℵ, T)
    msg3 = @sprintf("wall time: %s \n", prettytime(step_time))

    @info msg1 * msg2 * msg3

    wall_time[] = time_ns()
    return nothing
end

# And add it as a callback to the simulation.
add_callback!(arctic, progress, IterationInterval(10))

run!(arctic)