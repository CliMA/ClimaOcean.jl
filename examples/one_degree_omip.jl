# # One-degree OMIP simulation
#
# A global ocean--sea ice simulation on the 1ᵒ NEMO eORCA1 mesh, in the configuration of the Ocean Model
# Intercomparison Project (OMIP-2): forced by the multi-year JRA55-do reanalysis from 1958, with the surface
# salinity restored to the World Ocean Atlas climatology. The configuration and its defaults are those of
# `OneDegree()`; see `?omip_simulation` for the physics common to every resolution.

using ClimaOcean
using Oceananigans
using Oceananigans.Units
using CUDA

configuration = OneDegree()

simulation = omip_simulation(configuration; 
                             arch = GPU(), 
                             output_directory = "one_degree_omip")

# Five cycles of the 1958--2017 forcing. A run that was stopped restarts from its latest checkpoint
# with `run!(simulation; pickup = true)`.

simulation.stop_time = 300 * 365days
run!(simulation)
