using ClimaOcean
using Oceananigans
using Oceananigans.Units
using CFTime
using Dates
using Printf
using Glob

# ### Grid configuration 
#
# We define a global grid with a horizontal resolution of 1/4 degree and 40 vertical levels.
# The grid is a `LatitudeLongitudeGrid` spanning latitudes from 75°S to 75°N.
# We use an exponential vertical spacing to better resolve the upper-ocean layers.
# The total depth of the domain is set to 6000 meters.
# Finally, we specify the architecture for the simulation, which in this case is a GPU.

arch = CPU()

Nx = 144
Ny = 60
Nz = 40

depth = 6000meters
z_faces = exponential_z_faces(; Nz, depth)

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             z = z_faces,
                             latitude  = (-75, 75),
                             longitude = (0, 360))

# ### Ocean model configuration
#
# We build our ocean model using `ocean_simulation`,

ocean = ocean_simulation(grid)

# which uses the default `ocean.model`,

ocean.model

# We initialize the ocean model with ECCO2 temperature and salinity for January 1, 1993.

date = DateTimeProlepticGregorian(1993, 1, 1)
set!(ocean.model, T=ECCOMetadata(:temperature; dates=date),
                  S=ECCOMetadata(:salinity; dates=date))

# ### Prescribed atmosphere and radiation
#
# Next we build a prescribed atmosphere state and radiation model,
# which will drive the ocean simulation. We use the default `Radiation` model,

# The radiation model specifies an ocean albedo emissivity to compute the net radiative
# fluxes. The default ocean albedo is based on Payne (1982) and depends on cloud cover
# (calculated from the ratio of maximum possible incident solar radiation to actual
# incident solar radiation) and latitude. The ocean emissivity is set to 0.97.

radiation = Radiation(arch)

# The atmospheric data is prescribed using the JRA55 dataset.
# The JRA55 dataset provides atmospheric data such as temperature, humidity, and winds
# to calculate turbulent fluxes using bulk formulae, see [`CrossRealmFluxes`](@ref).
# The number of snapshots that are loaded into memory is determined by
# the `backend`. Here, we load 41 snapshots at a time into memory.

atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(41))

# ## The coupled simulation

# Next we assemble the ocean, atmosphere, and radiation
# into a coupled model,

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)

# We then create a coupled simulation. We start with a small-ish time step of 90 seconds.
# We run the simulation for 10 days with this small-ish time step.

simulation = Simulation(coupled_model; Δt=10, stop_iteration=50)

# We define a callback function to monitor the simulation's progress,

wall_time = Ref(time_ns())

function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T

    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))

    umax = (maximum(abs, interior(u)),
            maximum(abs, interior(v)),
            maximum(abs, interior(w)))

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg = @sprintf("Iter: %d, time: %s, Δt: %s", iteration(sim), prettytime(sim), prettytime(sim.Δt))
    msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e) m s⁻¹, extrema(T): (%.2f, %.2f) ᵒC, wall time: %s",
                    umax..., Tmax, Tmin, prettytime(step_time))

    @info msg 

    wall_time[] = time_ns()
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

# ### Set up output writers
#
# We define output writers to save the simulation data at regular intervals.
# In this case, we save the surface fluxes and surface fields at a relatively high frequency (every day).
# The `indices` keyword argument allows us to save only a slice of the three dimensional variable.
# Below, we use `indices` to save only the values of the variables at the surface, which corresponds to `k = grid.Nz`



outputs = merge(ocean.model.tracers, ocean.model.velocities)
ocean.output_writers[:surface] = JLD2OutputWriter(ocean.model, outputs;
                                                  schedule = IterationInterval(10),
                                                  filename = "checkpointer_mwe_surface",
                                                  indices = (:, :, grid.Nz),
                                                  with_halos = true,
                                                  overwrite_existing = true,
                                                  array_type = Array{Float32})
output_dir = "/g/data/v46/txs156/ClimaOcean.jl-checkpointer/examples/"
prefix = "checkpointer_mwe"

ocean.output_writers[:checkpoint] = Checkpointer(ocean.model;
                                                  schedule = IterationInterval(20),
                                                  prefix = prefix,
                                                  cleanup = true,
                                                  dir = output_dir,
                                                  verbose = true,
                                                  overwrite_existing = true)

# ### Spinning up the simulation
#
# We spin up the simulation with a small-ish time-step to resolve the "initialization shock"
# associated with starting from ECCO2 initial conditions that are both interpolated and also
# satisfy a different dynamical balance than our simulation.

# We check if a checkpointer already exists - if not, we can run the initial start up
pattern = prefix * "*"
checkpoint_file = glob(pattern, output_dir)
if !isempty(checkpoint_file)
    # If checkpoint exists, load the simulation state
    @info ("Checkpoint found, resuming the simulation from the checkpoint.")

    set!(simulation, checkpoint_file[1])
    
    coupled_model = OceanSeaIceModel(simulation.model.ocean; atmosphere, radiation)
    
    simulation = Simulation(coupled_model; Δt=10, stop_iteration=100)

    run!(simulation)

else
   @info ("Checkpoint not found, spinning up simulation from scratch.") 
   run!(simulation)
end
