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
τR = 2592000 

#  parabolic profile in Y, max @ j=4, min @ j=ny, amplitude=1.K
function T_restoring(i, j, k, grid, clock, fields, p)  
    Ny = size(grid, 2)
    j′ = (Ny - j + 4) / (Ny - 4) 
    Tr = p.Tf + 0.5 * j′^2
    Ti = @inbounds fields.T[i, j, k]
    return p.rate * (Tr - Ti)
end

FT = Forcing(T_restoring, discrete_form=true, parameters=(; rate=1/τR, Tf))

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

set!(ocean.model, T=Tf, u=0.2, S=30)

fix_salinity!(sim) = 
    fill!(sim.model.tracers.S, 30)
    
add_callback!(ocean, fix_salinity!, IterationInterval(1))

####
#### Sea ice simulation
####

sea_ice = sea_ice_simulation(grid; ice_salinity=4)

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

Ta(x, y) = 273.15 - 4 * sin(π * (1 + 2 * x / grid.Lx))
Tb(x, y) = Cf * exp(-Ce / Ta(x, y))
qa(x, y) = rh * Tb(x, y) / ρa

for t in eachindex(atmos_times)
    set!(atmosphere.tracers.T[t],    Ta)
    set!(atmosphere.tracers.q[t],    qa)
    set!(atmosphere.velocities.u[t], 10)
    set!(atmosphere.downwelling_radiation.longwave[t],  250)
    set!(atmosphere.downwelling_radiation.shortwave[t], 100)

    Oceananigans.BoundaryConditions.fill_halo_regions!(atmosphere.tracers.T[t])
    Oceananigans.BoundaryConditions.fill_halo_regions!(atmosphere.tracers.q[t])
    Oceananigans.BoundaryConditions.fill_halo_regions!(atmosphere.velocities.u[t])
end

####
#### Coupling
####

radiation = Radiation(sea_ice_albedo=0.6, ocean_albedo=0.05)
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
coupled_simulation = Simulation(coupled_model, Δt=1hour, stop_time=30days)

progress(sim) = 
    @info "Time: " * prettytime(sim)

add_callback!(coupled_simulation, progress, IterationInterval(100))

run!(coupled_simulation)

#####
##### Reading MITgcm results
#####

# Qas * SIarea + Qao * (1 - SIarea)

function load_mitgcm_results(path, iter)
    iter_string = string(iter, pad=10)
    file = joinpath(path, "snapshot.$(iter_string).data")

    # 'SIarea  ' 'SIheff  ' 'THETA   ' 'SItices ' 'SIqnet  ' 'SIqsw   ' 'SIempmr ' 'oceSflux' 'SIatmQnt' 'SIatmFW '
    variables = zeros(80 * 42 * 10)
    read!(file, variables)

    variables = bswap.(variables)
    variables = reshape(variables, 80, 42, 10)

    ℵ   = variables[:, :, 1]
    h   = variables[:, :, 2] 
    T   = variables[:, :, 3]
    Ts  = variables[:, :, 4]
    Q   = variables[:, :, 5]
    Qsw = variables[:, :, 6]
    emp = variables[:, :, 7]
    Jˢ  = variables[:, :, 8]
    Qˢ  = variables[:, :, 9]
    Fˢ  = variables[:, :, 10]

    Ts[Ts .== -999] .= NaN

    return (ℵ=ℵ, h=h, T=T, Ts=Ts, Q=Q, Qsw=Qsw, emp=emp, Jˢ=Jˢ, Qˢ=Qˢ, Fˢ=Fˢ)
end

mitgmc_thermo = load_mitgcm_results("mitgcm_results/run_thermo/res_30d", 719)

indices = (:, 35, 1)
fig = Figure()
ax  = Axis(fig[1, 1])
lines!(ax, mitgmc_thermo.T[indices...], color=:blue, label="MITgcm")
lines!(ax, ocean.model.tracers.T[indices...], color=:red, label="Oceananigans")

ax  = Axis(fig[1, 2])
lines!(ax, mitgmc_thermo.h[indices...], color=:blue, label="MITgcm")
lines!(ax, interior(sea_ice.model.ice_thickness, indices...) .* interior(sea_ice.model.ice_concentration, indices...), color=:red, label="Oceananigans")

ax  = Axis(fig[2, 1])
lines!(ax, mitgmc_thermo.Ts[indices...] .- 273.15, color=:blue, label="MITgcm")
lines!(ax, interior(coupled_model.interfaces.atmosphere_sea_ice_interface.temperature, indices...) .* interior(sea_ice.model.ice_concentration, indices...), color=:red, label="Oceananigans")

ax  = Axis(fig[2, 2])
lines!(ax, mitgmc_thermo.Qsw[indices...], color=:blue, label="MITgcm")
ℵ = interior(sea_ice.model.ice_concentration, indices...)
α = ℵ .* 0.4 + 0.9 * (1 .- ℵ)
lines!(ax, - (1 .-  ℵ).* 100 .* α, color=:red, label="Oceananigans")

ax  = Axis(fig[1, 3])
lines!(ax, mitgmc_thermo.ℵ[indices...], color=:blue, label="MITgcm")
lines!(ax, ℵ, color=:red, label="Oceananigans")
