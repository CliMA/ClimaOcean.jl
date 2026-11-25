# # Quarter-degree OMIP simulation
#
# The OMIP-2 simulation of `one_degree_omip.jl` on the eddy-permitting 1/4ᵒ NEMO eORCA025 mesh. At this resolution
# the eddies are resolved, so the configuration `QuarterDegree()` carries no Gent-McWilliams or Redi diffusivity
# and no explicit horizontal viscosity, and starts with a short time step that lengthens as the ocean adjusts.

using ClimaOcean
using Oceananigans
using Oceananigans.Units
using CUDA

configuration = QuarterDegree()

simulation = omip_simulation(configuration; 
                             arch = GPU(), 
                             output_directory = "quarter_degree_omip")

simulation.stop_time = 300 * 365days
run!(simulation)
