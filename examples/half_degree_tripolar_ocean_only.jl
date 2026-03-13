# # Half-degree tripolar ocean simulation
#
# Ocean-only simulation on a 0.5° `TripolarGrid` (720×360) with CATKE,
# Gent-McWilliams, and biharmonic viscosity. Forced by repeat-year JRA55
# atmospheric reanalysis and run for two years.

using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Dates
using Printf
using CUDA
using CairoMakie

# ### Build the ocean

arch = GPU()
ocean = half_degree_tripolar_ocean(arch; zstar=true)

# ### Initial conditions from ECCO

date = DateTime(1993, 1, 1)
set!(ocean.model, T = Metadatum(:temperature; date, dataset = ECCO4Monthly()),
                  S = Metadatum(:salinity;    date, dataset = ECCO4Monthly()))

# ### Atmospheric forcing

atmosphere = JRA55PrescribedAtmosphere(arch; backend = JRA55NetCDFBackend(41),
                                       include_rivers_and_icebergs = false)
radiation = Radiation(arch)

# ### Simulation

coupled_model = OceanOnlyModel(ocean; atmosphere, radiation)
simulation = Simulation(coupled_model; Δt = 10minutes, stop_time = 2 * 365days)

# ### Progress messenger

add_callback!(simulation, Progress(), IterationInterval(100))

# ### Output writers

outputs = merge(ocean.model.tracers, ocean.model.velocities)

ocean.output_writers[:surface] = JLD2Writer(ocean.model, outputs;
                                            schedule = TimeInterval(5days),
                                            filename = "half_degree_ocean_only_surface",
                                            indices = (:, :, ocean.model.grid.Nz),
                                            including = [:grid],
                                            overwrite_existing = true)

# ### Run!

run!(simulation)

# ### Diagnostic report

simulation_report(ocean, filename = "half_degree_ocean_only_report.png")

# ![](half_degree_ocean_only_report.png)
