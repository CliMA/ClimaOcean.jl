# # Twelfth-degree OMIP simulation
#
# The OMIP-2 simulation on the eddy-resolving 1/12ᵒ NEMO eORCA12 mesh, distributed over four GPUs:
#
# ```
# mpiexec -n 4 julia --project=examples examples/twelfth_degree_omip.jl
# ```

using MPI
MPI.Init()

using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Oceananigans.DistributedComputations
using CUDA

configuration = TwelfthDegree()
arch = Distributed(GPU(), partition = Partition(1, 4))

simulation = omip_simulation(configuration; 
                             arch, 
                             output_directory = "twelfth_degree_omip",
                             file_splitting_interval = 180days)

simulation.stop_time = 300 * 365days
run!(simulation)
