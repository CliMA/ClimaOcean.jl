# # [One-degree global ocean--sea ice simulation](@id one-degree-ocean-seaice)
#
# This example configures a global ocean--sea ice simulation at 1ᵒ horizontal resolution with
# realistic bathymetry and a few closures including the "Gent-McWilliams" `IsopycnalSkewSymmetricDiffusivity`.
# The simulation is forced by repeat-year JRA55 atmospheric reanalysis
# and initialized by temperature, salinity, sea ice concentration, and sea ice thickness
# from the ECCO state estimate.
#
# For this example, we need Oceananigans, ClimaOcean, Dates, CUDA, and
# CairoMakie to visualize the simulation.

using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Dates
using Printf
using Statistics
using CUDA

# ### Grid and Bathymetry

# We start by constructing an underlying TripolarGrid at ~1 degree resolution,

arch = GPU()
Nx = 360
Ny = 180
Nz = 50

depth = 5000meters
z = ExponentialDiscretization(Nz, -depth, 0; scale = depth/4)

underlying_grid = TripolarGrid(arch; size = (Nx, Ny, Nz), halo = (5, 5, 4), z)

# Next, we build bathymetry on this grid, using interpolation passes to smooth the bathymetry.
# With 2 major basins, we keep the Mediterranean (though we need to manually open the Gibraltar
# Strait to connect it to the Atlantic):

bottom_height = regrid_bathymetry(underlying_grid;
                                  minimum_depth = 10,
                                  interpolation_passes = 10,
                                  major_basins = 2)

# We then incorporate the bathymetry into an ImmersedBoundaryGrid,

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height);
                            active_cells_map=true)

# ### Closures
#
# We include a Gent-McWilliams isopycnal diffusivity as a parameterization for the mesoscale
# eddy fluxes. For vertical mixing at the upper-ocean boundary layer we include the CATKE
# parameterization.

using Oceananigans.TurbulenceClosures: IsopycnalSkewSymmetricDiffusivity, AdvectiveFormulation

eddy_closure = IsopycnalSkewSymmetricDiffusivity(κ_skew=1e3, κ_symmetric=1e3, skew_flux_formulation=AdvectiveFormulation()) 
@inline νhb(i, j, k, grid, ℓx, ℓy, ℓz, clock, fields, λ) = Oceananigans.Operators.Az(i, j, k, grid, ℓx, ℓy, ℓz)^2 / λ
horizontal_viscosity = HorizontalScalarBiharmonicDiffusivity(ν=νhb, discrete_form=true, parameters=15days)
vertical_mixing = ClimaOcean.Oceans.default_ocean_closure()

# ### Ocean simulation
# Now we bring everything together to construct the ocean simulation.
# We use a split-explicit timestepping with 70 substeps for the barotropic mode.

free_surface       = SplitExplicitFreeSurface(grid; substeps=70)
momentum_advection = WENOVectorInvariant(order=5)
tracer_advection   = WENO(order=5)

ocean = ocean_simulation(grid; momentum_advection, tracer_advection, free_surface,
                         closure=(eddy_closure, horizontal_viscosity, vertical_mixing))

@info "We've built an ocean simulation with model:"
@show ocean.model

# ### Sea Ice simulation
#
# We also build a sea ice simulation. We use the default configuration:
# EVP rheology and a zero-layer thermodynamic model that advances thickness
# and concentration.

sea_ice = sea_ice_simulation(grid, ocean; advection=tracer_advection)

# ### Initial condition

# We initialize the ocean and sea ice models with data from the ECCO state estimate.

date = DateTime(1993, 1, 1)
dataset = ECCO4Monthly()
ecco_temperature           = Metadatum(:temperature; date, dataset)
ecco_salinity              = Metadatum(:salinity; date, dataset)
ecco_sea_ice_thickness     = Metadatum(:sea_ice_thickness; date, dataset)
ecco_sea_ice_concentration = Metadatum(:sea_ice_concentration; date, dataset)

set!(ocean.model, T=ecco_temperature, S=ecco_salinity)
set!(sea_ice.model, h=ecco_sea_ice_thickness, ℵ=ecco_sea_ice_concentration)

# ### Atmospheric forcing

# We force the simulation with a JRA55-do atmospheric reanalysis.
radiation  = Radiation(arch)
atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(80),
                                       include_rivers_and_icebergs = false)

# ### Coupled simulation

# Now we are ready to build the coupled ocean--sea ice model and bring everything
# together into a `simulation`.

# With Runge-Kutta 3rd order time-stepping we can safely use a timestep of 20 minutes.

coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
simulation = Simulation(coupled_model; Δt=20minutes, stop_time=365days)

# ### A progress messenger
#
# We write a function that prints out a helpful progress message while the simulation runs.

wall_time = Ref(time_ns())

function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T
    e = ocean.model.tracers.e
    Tmin, Tmax, Tavg = minimum(T), maximum(T), mean(view(T, :, :, ocean.model.grid.Nz))
    emax = maximum(e)
    umax = (maximum(abs, u), maximum(abs, v), maximum(abs, w))

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg1 = @sprintf("time: %s, iter: %d", prettytime(sim), iteration(sim))
    msg2 = @sprintf(", max|uo|: (%.1e, %.1e, %.1e) m s⁻¹", umax...)
    msg3 = @sprintf(", extrema(To): (%.1f, %.1f) ᵒC, mean(To(z=0)): %.1f ᵒC", Tmin, Tmax, Tavg)
    msg4 = @sprintf(", max(e): %.2f m² s⁻²", emax)
    msg5 = @sprintf(", wall time: %s \n", prettytime(step_time))

    @info msg1 * msg2 * msg3 * msg4 * msg5

    wall_time[] = time_ns()

     return nothing
end

# And add it as a callback to the simulation.
add_callback!(simulation, progress, TimeInterval(5days))

# ### Output
#
# We are almost there! We need to save some output. Below we choose to save _only surface_
# fields using the `indices` keyword argument. We save all the velocity and tracer components.
# Note, that besides temperature and salinity, the CATKE vertical mixing parameterization
# also uses a prognostic turbulent kinetic energy, `e`, to diagnose the vertical mixing length.

ocean_outputs = merge(ocean.model.tracers, ocean.model.velocities)
sea_ice_outputs = merge((h = sea_ice.model.ice_thickness,
                         ℵ = sea_ice.model.ice_concentration,
                         T = sea_ice.model.ice_thermodynamics.top_surface_temperature),
                         sea_ice.model.velocities)

ocean.output_writers[:surface] = JLD2Writer(ocean.model, ocean_outputs;
                                            schedule = TimeInterval(1days),
                                            filename = "ocean_one_degree_surface_fields",
                                            indices = (:, :, grid.Nz),
                                            overwrite_existing = true)

sea_ice.output_writers[:surface] = JLD2Writer(sea_ice.model, sea_ice_outputs;
                                              schedule = TimeInterval(1days),
                                              filename = "sea_ice_one_degree_surface_fields",
                                              overwrite_existing = true)

# ### Ready to run

# We are ready to press the big red button and run the simulation.
run!(simulation)

# ### A movie
#
# We load the saved output and make a movie of the simulation. First we plot a snapshot:
using CairoMakie

# We suffix the ocean fields with "o":
uo = FieldTimeSeries("ocean_one_degree_surface_fields.jld2",  "u"; backend = OnDisk())
vo = FieldTimeSeries("ocean_one_degree_surface_fields.jld2",  "v"; backend = OnDisk())
To = FieldTimeSeries("ocean_one_degree_surface_fields.jld2",  "T"; backend = OnDisk())
eo = FieldTimeSeries("ocean_one_degree_surface_fields.jld2",  "e"; backend = OnDisk())

# and sea ice fields with "i":
ui = FieldTimeSeries("sea_ice_one_degree_surface_fields.jld2", "u"; backend = OnDisk())
vi = FieldTimeSeries("sea_ice_one_degree_surface_fields.jld2", "v"; backend = OnDisk())
hi = FieldTimeSeries("sea_ice_one_degree_surface_fields.jld2", "h"; backend = OnDisk())
ℵi = FieldTimeSeries("sea_ice_one_degree_surface_fields.jld2", "ℵ"; backend = OnDisk())
Ti = FieldTimeSeries("sea_ice_one_degree_surface_fields.jld2", "T"; backend = OnDisk())

times = uo.times
Nt = length(times)
n = Observable(Nt)

# We create a land mask and use it to fill land points with `NaN`s.
land = interior(To.grid.immersed_boundary.bottom_height) .≥ 0

Toₙ = @lift begin
    Tₙ = interior(To[$n])
    Tₙ[land] .= NaN
    view(Tₙ, :, :, 1)
end

eoₙ = @lift begin
    eₙ = interior(eo[$n])
    eₙ[land] .= NaN
    view(eₙ, :, :, 1)
end

heₙ = @lift begin
    hₙ = interior(hi[$n])
    ℵₙ = interior(ℵi[$n])
    hₙ[land] .= NaN
    view(hₙ, :, :, 1) .* view(ℵₙ, :, :, 1)
end

# We compute the surface speeds for the ocean and the sea ice.
uoₙ = Field{Face, Center, Nothing}(uo.grid)
voₙ = Field{Center, Face, Nothing}(vo.grid)

uiₙ = Field{Face, Center, Nothing}(ui.grid)
viₙ = Field{Center, Face, Nothing}(vi.grid)

so = Field(sqrt(uoₙ^2 + voₙ^2))
si = Field(sqrt(uiₙ^2 + viₙ^2))

soₙ = @lift begin
    parent(uoₙ) .= parent(uo[$n])
    parent(voₙ) .= parent(vo[$n])
    compute!(so)
    soₙ = interior(so)
    soₙ[land] .= NaN
    view(soₙ, :, :, 1)
end

siₙ = @lift begin
    parent(uiₙ) .= parent(ui[$n])
    parent(viₙ) .= parent(vi[$n])
    compute!(si)
    siₙ = interior(si)
    hₙ = interior(hi[$n])
    ℵₙ = interior(ℵi[$n])
    he = hₙ .* ℵₙ
    siₙ[he .< 1e-7] .= 0
    siₙ[land] .= NaN
    view(siₙ, :, :, 1)
end

# Finally, we plot a snapshot of the surface speed, temperature, and the turbulent
# eddy kinetic energy from the CATKE vertical mixing parameterization as well as the
# sea ice speed and the effective sea ice thickness.
fig = Figure(size=(1200, 1000))

title = @lift string("Global 1ᵒ ocean simulation after ", prettytime(times[$n] - times[1]))

axso = Axis(fig[1, 1])
axsi = Axis(fig[1, 3])
axTo = Axis(fig[2, 1])
axhi = Axis(fig[2, 3])
axeo = Axis(fig[3, 1])

hmo = heatmap!(axso, soₙ, colorrange = (0, 0.5), colormap = :deep,  nan_color=:lightgray)
hmi = heatmap!(axsi, siₙ, colorrange = (0, 0.5), colormap = :greys, nan_color=:lightgray)
Colorbar(fig[1, 2], hmo, label = "Ocean Surface speed (m s⁻¹)")
Colorbar(fig[1, 4], hmi, label = "Sea ice speed (m s⁻¹)")

hmo = heatmap!(axTo, Toₙ, colorrange = (-1, 32), colormap = :magma, nan_color=:lightgray)
hmi = heatmap!(axhi, heₙ, colorrange =  (0, 4),  colormap = :blues, nan_color=:lightgray)
Colorbar(fig[2, 2], hmo, label = "Surface Temperature (ᵒC)")
Colorbar(fig[2, 4], hmi, label = "Effective ice thickness (m)")

hm = heatmap!(axeo, eoₙ, colorrange = (0, 1e-3), colormap = :solar, nan_color=:lightgray)
Colorbar(fig[3, 2], hm, label = "Turbulent Kinetic Energy (m² s⁻²)")

for ax in (axso, axsi, axTo, axhi, axeo)
    hidedecorations!(ax)
end

Label(fig[0, :], title)

save("global_snapshot.png", fig)
nothing #hide

# ![](global_snapshot.png)

# And now a movie:

CairoMakie.record(fig, "one_degree_global_ocean_surface.mp4", 1:Nt, framerate = 8) do nn
    n[] = nn
end
nothing #hide

# ![](one_degree_global_ocean_surface.mp4)
