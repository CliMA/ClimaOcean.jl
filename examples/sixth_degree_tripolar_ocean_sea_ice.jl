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
T_meta = Metadatum(:temperature;          date, dataset = ECCO4Monthly())
S_meta = Metadatum(:salinity;             date, dataset = ECCO4Monthly())
h_meta = Metadatum(:sea_ice_thickness;    date, dataset = ECCO4Monthly())
ℵ_meta = Metadatum(:sea_ice_concentration; date, dataset = ECCO4Monthly())

# Pre-download via the NumericalEarthArtifacts mirror so the build survives
# transient outages of the ECCO drive.
foreach(download_with_fallback, (T_meta, S_meta, h_meta, ℵ_meta))

set!(ocean.model,   T = T_meta, S = S_meta)
set!(sea_ice.model, h = h_meta, ℵ = ℵ_meta)

# ### Atmospheric forcing

atmosphere = JRA55PrescribedAtmosphere(arch; time_indices_in_memory = 41)
radiation  = JRA55PrescribedRadiation(arch; time_indices_in_memory = 41)

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
