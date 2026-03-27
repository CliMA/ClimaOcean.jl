# Sixth-degree tripolar ocean–sea ice simulation (distributed, 4 GPUs)
#
# Run with: srun -n 4 julia --project examples/sixth_degree_tripolar_ocean_sea_ice.jl
#
# Coupled ocean–sea ice simulation on a 1/6° TripolarGrid (2160×1080)
# distributed across 4 GPUs with Partition(2, 2). Forced by repeat-year
# JRA55 atmospheric reanalysis and run for two years. Output is saved
# to JLD2 files for post-processing.

using MPI
MPI.Init()

using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Oceananigans.DistributedComputations
using Dates
using Printf

# ### Distributed architecture

arch = Distributed(GPU(), partition = Partition(2, 2))

# ### Build ocean and sea ice

ocean   = sixth_degree_tripolar_ocean(arch; zstar=true)
sea_ice = sixth_degree_tripolar_sea_ice(ocean)

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
simulation = Simulation(coupled_model; Δt = 5minutes, stop_time = 2 * 365days)

# ### Progress messenger

add_callback!(simulation, Progress(), IterationInterval(100))

# ### Output writers

ocean_outputs = merge(ocean.model.tracers, ocean.model.velocities)

ocean.output_writers[:surface] = JLD2Writer(ocean.model, ocean_outputs;
                                            schedule = TimeInterval(5days),
                                            filename = "sixth_degree_coupled_ocean_surface",
                                            indices = (:, :, ocean.model.grid.Nz),
                                            including = [:grid],
                                            overwrite_existing = true)

sea_ice_outputs = (; h = sea_ice.model.ice_thickness,
                     ℵ = sea_ice.model.ice_concentration)

sea_ice.output_writers[:fields] = JLD2Writer(sea_ice.model, sea_ice_outputs;
                                             schedule = TimeInterval(5days),
                                             filename = "sixth_degree_coupled_sea_ice",
                                             overwrite_existing = true)

# ### Run!

run!(simulation)

# ### Diagnostic report (rank 0 only)

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    using CairoMakie
    simulation_report(ocean, filename = "sixth_degree_coupled_report.png")
end

# ![](sixth_degree_coupled_report.png)
