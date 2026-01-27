using ClimaOcean, Oceananigans, JLD2

arch = CPU()

Nx = 720 # longitudinal direction
Ny = 360 # meridional direction
Nz = 100

z_faces = ExponentialDiscretization(Nz, -6000, 0; scale=1800, mutable=true)

grid = TripolarGrid(arch;
                    size = (Nx, Ny, Nz),
                    z = z_faces,
                    halo = (7, 7, 7))

bottom_height = Oceananigans.on_architecture(arch, jldopen("bottom_height.jld2")["bottom"])
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

file = jldopen("halfdegree_iteration957000.jld2")
u = set!(XFaceField(grid), file["uo"])
v = set!(YFaceField(grid), file["vo"])
η = set!(ZFaceField(grid, indices=(:, :, 100)), file["ηo"])

using Oceananigans.Utils: launch!
using Oceananigans.Models.HydrostaticFreeSurfaceModels: _update_zstar_scaling!, surface_kernel_parameters

Oceananigans.Utils.launch!(CPU(), grid, surface_kernel_parameters(grid), _update_zstar_scaling!, η, grid)
parent(grid.z.σᶜᶜ⁻) .= parent(grid.z.σᶜᶜⁿ)

Oceananigans.BoundaryConditions.fill_halo_regions!((u, v))

using ClimaOcean.Diagnostics: BrokenLineSet, atlantic_ocean_mask, compute_streamfunction

# Create broken lines for Atlantic AMOC
alb = ClimaOcean.Diagnostics.BrokenLineSet(grid, 0:0.5:65; basin_mask = atlantic_ocean_mask(grid))

# Compute the streamfunction
ψ = ClimaOcean.Diagnostics.compute_streamfunction(u, v, alb)

# Convert to Sverdrups
ψ_Sv = ψ ./ 1e6

# Get the latitudes for plotting
lats = ClimaOcean.Diagnostics.band_latitudes(alb)

# Print some diagnostics
println("AMOC streamfunction computed")
println("  Latitude range: ", minimum(lats), "° to ", maximum(lats), "°")
println("  Number of latitudes: ", length(lats))
println("  Streamfunction range: ", minimum(ψ_Sv), " to ", maximum(ψ_Sv), " Sv")
println("  Max AMOC strength: ", maximum(ψ_Sv), " Sv")