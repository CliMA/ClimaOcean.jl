# # One-degree tripolar ocean–sea ice simulation
#
# Coupled ocean–sea ice simulation on a 1° `TripolarGrid` (360×180) with CATKE,
# Gent-McWilliams, and biharmonic viscosity. Forced by repeat-year JRA55
# atmospheric reanalysis and run for two years.

using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Dates
using Printf
using CUDA
using CairoMakie

# ### Build ocean and sea ice

arch = GPU()
ocean   = one_degree_tripolar_ocean(arch; zstar=true)
sea_ice = one_degree_tripolar_sea_ice(ocean)

# ### Initial conditions from ECCO

date = DateTime(1993, 1, 1)
set!(ocean.model, T = Metadatum(:temperature; date, dataset = ECCO4Monthly()),
                  S = Metadatum(:salinity;    date, dataset = ECCO4Monthly()))

set!(sea_ice.model, h = Metadatum(:sea_ice_thickness;    date, dataset = ECCO4Monthly()),
                    ℵ = Metadatum(:sea_ice_concentration; date, dataset = ECCO4Monthly()))

# ### Atmospheric forcing

atmosphere = JRA55PrescribedAtmosphere(arch; backend = JRA55NetCDFBackend(41),
                                       include_rivers_and_icebergs = false)
radiation = Radiation(arch)

# ### Simulation

coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
simulation = Simulation(coupled_model; Δt = 20minutes, stop_time = 2 * 365days)

# ### Progress messenger

add_callback!(simulation, Progress(), IterationInterval(100))

# ### Output writers

ocean_outputs = merge(ocean.model.tracers, ocean.model.velocities)

ocean.output_writers[:surface] = JLD2Writer(ocean.model, ocean_outputs;
                                            schedule = TimeInterval(5days),
                                            filename = "one_degree_coupled_ocean_surface",
                                            indices = (:, :, ocean.model.grid.Nz),
                                            including = [:grid],
                                            overwrite_existing = true)

sea_ice_outputs = (; h = sea_ice.model.ice_thickness,
                     ℵ = sea_ice.model.ice_concentration)

sea_ice.output_writers[:fields] = JLD2Writer(sea_ice.model, sea_ice_outputs;
                                             schedule = TimeInterval(5days),
                                             filename = "one_degree_coupled_sea_ice",
                                             overwrite_existing = true)

# ### Run!

run!(simulation)

# ### Diagnostic report

simulation_report(ocean, filename = "one_degree_coupled_report.png")

# ![](one_degree_coupled_report.png)
