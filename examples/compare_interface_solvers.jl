# # Comparing Interface State Solvers
#
# This example compares the performance and accuracy of two iterative solvers
# for computing surface fluxes in ClimaOcean:
#
# 1. **Fixed-point (Picard) iteration** - The traditional approach
# 2. **Good Broyden's method** - A quasi-Newton method with faster convergence
#
# Both solvers compute the characteristic scales (u*, theta*, q*) used in
# Monin-Obukhov similarity theory for turbulent flux calculations.

using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.JRA55
using ClimaOcean.Oceans
using BenchmarkTools
using ClimaOcean.OceanSeaIceModels.InterfaceComputations:
    SimilarityTheoryFluxes,
    FixedPointSolver,
    BroydenSolver,
    ConvergenceStopCriteria

using Oceananigans
using Dates
using Printf
using Statistics

# # Setup: Build the grid and load data
#
# We use a subset of the ECCO grid for faster computation.
# The full global grid can be used for more comprehensive benchmarking.
# Load atmosphere data (first two time indices from JRA55)

atmosphere = JRA55PrescribedAtmosphere(; backend = JRA55NetCDFBackend(2))

# # Define the solvers to compare
#
# We create two flux formulations with different solvers:
# - Fixed-point solver with default settings
# - Broyden solver (Val(3) for 3 variables: u*, theta*, q*)

tolerance = 1e-8
maxiter = 100

## Fixed-point solver
fp_solver = FixedPointSolver(Float64; tolerance, maxiter)
fp_fluxes = SimilarityTheoryFluxes(Float64; solver = fp_solver)

## Broyden solver (3 variables for BulkTemperature)
br_solver = BroydenSolver(Val(3), Float64; tolerance, maxiter)
br_fluxes = SimilarityTheoryFluxes(Float64; solver = br_solver)

println("Solver configurations:")
println("  Fixed-point: ", fp_fluxes.solver)
println("  Broyden:     ", br_fluxes.solver)

# # Run comparison: Fixed-point solver
#
# First, we create an ocean simulation and coupled model using the fixed-point solver.

println("\n" * "="^60)
println("Running with Fixed-Point Solver")
println("="^60)

## Set ocean state from ECCO data
T_metadata = ECCOMetadatum(:temperature; date = DateTime(1993, 1, 1))
S_metadata = ECCOMetadatum(:salinity;    date = DateTime(1993, 1, 1))

grid  = Field(T_metadata).grid
ocean = ocean_simulation(grid; closure = nothing, momentum_advection = nothing, tracer_advection = nothing)

set!(ocean.model; T = T_metadata, S = S_metadata)

fp_interfaces = ComponentInterfaces(atmosphere, ocean; atmosphere_ocean_fluxes = fp_fluxes)
br_interfaces = ComponentInterfaces(atmosphere, ocean; atmosphere_ocean_fluxes = br_fluxes)

## Create coupled model with fixed-point solver
## Time the flux computation
coupled_model_fp = OceanSeaIceModel(ocean;
                                    atmosphere,
                                    interfaces = fp_interfaces)

t_fp = @benchmark begin
    OceanSeaIceModel(ocean;
                     atmosphere,
                     interfaces = fp_interfaces)
end

t_fp = t_fp.times[1]

println(@sprintf("  Flux computation time: %.3f seconds", t_fp))

## Extract computed fluxes
fluxes_fp = coupled_model_fp.interfaces.atmosphere_ocean_interface.fluxes

# # Run comparison: Broyden solver
#
# Now we repeat with the Broyden solver.

println("\n" * "="^60)
println("Running with Broyden Solver")
println("="^60)

coupled_model_br = OceanSeaIceModel(ocean;
                                    atmosphere,
                                    interfaces = br_interfaces)

## Create coupled model with Broyden solver
## Time the flux computation
t_br = @benchmark begin
    OceanSeaIceModel(ocean;
                     atmosphere,
                     interfaces = br_interfaces)
end

t_br = t_br.times[1]

println(@sprintf("  Flux computation time: %.3f seconds", t_br))

## Extract computed fluxes
fluxes_br = coupled_model_br.interfaces.atmosphere_ocean_interface.fluxes

# # Compare results
#
# We compare the flux fields computed by both solvers to verify they produce
# consistent results.

println("\n" * "="^60)
println("Comparing Results")
println("="^60)

## Helper function to compute statistics
function compare_fields(field_fp, field_br, name)
    fp_data = interior(field_fp, :, :, 1)
    br_data = interior(field_br, :, :, 1)

    ## Filter out NaN values (land points)
    valid = .!isnan.(fp_data) .& .!isnan.(br_data)
    fp_valid = fp_data[valid]
    br_valid = br_data[valid]

    if length(fp_valid) > 0
        diff = abs.(fp_valid .- br_valid)
        rel_diff = diff ./ (abs.(fp_valid) .+ 1e-10)

        max_abs_diff = maximum(diff)
        mean_abs_diff = mean(diff)
        max_rel_diff = maximum(rel_diff)
        mean_rel_diff = mean(rel_diff)

        println(@sprintf("  %s:", name))
        println(@sprintf("    Max absolute difference:  %.2e", max_abs_diff))
        println(@sprintf("    Mean absolute difference: %.2e", mean_abs_diff))
        println(@sprintf("    Max relative difference:  %.2e", max_rel_diff))
        println(@sprintf("    Mean relative difference: %.2e", mean_rel_diff))
    else
        println("  $name: No valid data points")
    end
end

compare_fields(fluxes_fp.sensible_heat, fluxes_br.sensible_heat, "Sensible heat flux")
compare_fields(fluxes_fp.latent_heat, fluxes_br.latent_heat, "Latent heat flux")
compare_fields(fluxes_fp.x_momentum, fluxes_br.x_momentum, "Zonal momentum flux")
compare_fields(fluxes_fp.y_momentum, fluxes_br.y_momentum, "Meridional momentum flux")
compare_fields(fluxes_fp.water_vapor, fluxes_br.water_vapor, "Water vapor flux")

# # Performance summary

println("\n" * "="^60)
println("Performance Summary")
println("="^60)
println(@sprintf("  Fixed-point solver time: %.3f seconds", t_fp))
println(@sprintf("  Broyden solver time:     %.3f seconds", t_br))
speedup = t_fp / t_br
if speedup > 1
    println(@sprintf("  Broyden speedup: %.2fx faster", speedup))
else
    println(@sprintf("  Fixed-point speedup: %.2fx faster", 1/speedup))
end

# # Visualize the comparison
#
# Create side-by-side plots of the sensible heat flux from both solvers,
# plus a difference plot.

using CairoMakie

fig = Figure(size = (1200, 800), fontsize = 12)

lambda, phi, z = nodes(fluxes_fp.sensible_heat)

## Fixed-point result
ax1 = Axis(fig[1, 1], title = "Sensible Heat Flux (Fixed-Point)", ylabel = "Latitude")
hm1 = heatmap!(ax1, lambda, phi, interior(fluxes_fp.sensible_heat, :, :, 1);
               colormap = :bwr, colorrange = (-200, 200))
Colorbar(fig[1, 2], hm1, label = "W/m^2")

## Broyden result
ax2 = Axis(fig[1, 3], title = "Sensible Heat Flux (Broyden)", ylabel = "Latitude")
hm2 = heatmap!(ax2, lambda, phi, interior(fluxes_br.sensible_heat, :, :, 1);
               colormap = :bwr, colorrange = (-200, 200))
Colorbar(fig[1, 4], hm2, label = "W/m^2")

## Difference
diff_sensible = interior(fluxes_fp.sensible_heat, :, :, 1) .- interior(fluxes_br.sensible_heat, :, :, 1)
ax3 = Axis(fig[2, 1], title = "Difference (Fixed-Point - Broyden)",
           xlabel = "Longitude", ylabel = "Latitude")
hm3 = heatmap!(ax3, lambda, phi, diff_sensible;
               colormap = :bwr, colorrange = (-1e-6, 1e-6))
Colorbar(fig[2, 2], hm3, label = "W/m^2")

## Latent heat comparison
ax4 = Axis(fig[2, 3], title = "Latent Heat Flux (Fixed-Point)",
           xlabel = "Longitude", ylabel = "Latitude")
hm4 = heatmap!(ax4, lambda, phi, interior(fluxes_fp.latent_heat, :, :, 1);
               colormap = :bwr, colorrange = (-400, 400))
Colorbar(fig[2, 4], hm4, label = "W/m^2")

save("solver_comparison.png", fig)
println("\nFigure saved to solver_comparison.png")

# ![](solver_comparison.png)

# # Conclusion
#
# Both solvers should produce nearly identical results (within numerical tolerance).
# The Broyden solver typically converges in fewer iterations than fixed-point iteration,
# which can lead to performance improvements, especially for difficult convergence cases.
#
# Key observations:
# - Results should match to within the solver tolerance (1e-8)
# - Broyden's method uses quasi-Newton updates which can accelerate convergence
# - For well-behaved problems, both solvers perform similarly
# - For stiff or slowly converging problems, Broyden may offer significant speedup
