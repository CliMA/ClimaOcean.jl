using ClimaOcean
using Oceananigans
using Oceananigans.Units

Ny = 41
Nx = 80

grid = RectilinearGrid(size=(Nx, Ny, 1), 
                          x=(0, 400kilometers), 
                          y=(0, 200kilometers), 
                          z=(-10, 0),
                   topology=(Periodic, Bounded, Bounded))

####
#### Ocean simulation
####

μf = 5.4e-2
Tf = - μf * 30

#  parabolic profile in Y, max @ j=4, min @ j=ny, amplitude=1.K
function T_restoring(i, j, k, grid, clock, fields, p)  
    Ny = size(grid, 2)
    j′ = (j - 4) / (Ny - 4) 
    Tr = p.Tf + 0.5 * j′^2
    Ti = @inbounds fields.T[i, j, k]
    return p.rate * (Tr - Ti)
end

FT = Forcing(T_restoring, discrete_form=true, parameters=(; rate=1/2days, Tf))

ocean = ocean_simulation(grid;
                         momentum_advection = nothing,
                         tracer_advection = nothing,
                         free_surface = nothing,
                         closure = nothing,
                         bottom_drag_coefficient = 0,
                         equation_of_state = LinearEquationOfState(thermal_expansion=2e-4, haline_contraction=0),
                         forcing = (; T=FT,)
)

# Make sure we use FE!
ocean.model.timestepper.χ = - 0.5

# Ocean initial conditions

set!(ocean.model, T=Tf, u=0.2)

####
#### Sea ice simulation
####

sea_ice = sea_ice_simulation(grid)

ℵi = ones(size(grid))
ℵi[:, 1]    .= 0.0
ℵi[:, 2]    .= 0.1
ℵi[:, Ny]   .= 0.0
ℵi[:, Ny-1] .= 0.1

set!(sea_ice.model, h=0.2, ℵ=ℵi)

####
#### Atmosphere simulation
####

atmos_times = range(0, 360Oceananigans.Units.days, length=10)
atmosphere  = PrescribedAtmosphere(grid, atmos_times)

Cf = 640380.0
Ce = 5107.4
ρa = 1.2
rh = 0.7

Ta(x, y) = 273.15 + 4 * sin(π * (1 + 2 * x / grid.Lx))
Tb(x, y) = Cf * exp(-Ce / Ta(x, y))
qa(x, y) = rh * Tb(x, y) / ρa

for t in eachindex(atmos_times)
    set!(atmosphere.tracers.T[t],    Ta)
    set!(atmosphere.tracers.q[t],    qa)
    set!(atmosphere.velocities.u[t], 10)

    Oceananigans.BoundaryConditions.fill_halo_regions!(atmosphere.tracers.T[t])
    Oceananigans.BoundaryConditions.fill_halo_regions!(atmosphere.tracers.q[t])
    Oceananigans.BoundaryConditions.fill_halo_regions!(atmosphere.velocities.u[t])
end

####
#### Coupling
####

radiation = Radiation(sea_ice_albedo=0.6, ocean_albedo=0.1)
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
coupled_simulation = Simulation(coupled_model, Δt=1hour, stop_time=30days)

progress(sim) = 
    @info "Time: " * prettytime(sim)

add_callback!(coupled_simulation, progress, IterationInterval(100))

run!(coupled_simulation)