using Printf
using Oceananigans.Units
using ClimaOcean
using Oceananigans
using Oceananigans.Operators
using ClimaOcean.DataWrangling
using ClimaOcean.OceanSimulations
using Oceananigans.Units
using ClimaOcean.DataWrangling.ECCO
using ClimaOcean.DataWrangling.ECCO: inpainted_ecco_field, ECCO2Daily, ECCO4Monthly
using ClimaOcean.DataWrangling.JRA55: JRA55_prescribed_atmosphere
using ClimaOcean.OceanSeaIceModels: Radiation
using CFTime

# Upload ECCO fields
T = inpainted_ecco_field(:temperature; maxiter = Inf)
S = inpainted_ecco_field(:salinity;    maxiter = Inf)
u = inpainted_ecco_field(:u_velocity;  maxiter = Inf)
v = inpainted_ecco_field(:v_velocity;  maxiter = Inf)

include("ecco_immersed_grid.jl")
grid = ecco_immersed_grid(metadata)

atmosphere  = JRA55_prescribed_atmosphere(1:2; grid = grid.underlying_grid)

ocean = ocean_simulation(grid; momentum_advection = nothing,
                                 tracer_advection = nothing)

ocean_model = ocean.model

# setting ecco variables in the model
set!(ocean_model, T = T, S = S, u = u, v = v)

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation = Radiation())

function centered_surface_u_velocity(u)
    𝒰ᶜᶜᶜ = KernelFunctionOperation{Center, Center, Center}(ℑxᶜᵃᵃ, grid, u)
    uᶜᶜᶜ = Field(𝒰ᶜᶜᶜ)
    compute!(uᶜᶜᶜ)

    Nz = size(uᶜᶜᶜ.grid, 3)

    return Array(interior(uᶜᶜᶜ, :, :, Nz))
end

function centered_surface_v_velocity(v)
    𝒱ᶜᶜᶜ = KernelFunctionOperation{Center, Center, Center}(ℑyᵃᶜᵃ, grid, v)
    vᶜᶜᶜ = Field(𝒱ᶜᶜᶜ)
    compute!(vᶜᶜᶜ)

    Nz = size(vᶜᶜᶜ.grid, 3)

    return Array(interior(vᶜᶜᶜ, :, :, Nz))
end

# Write down variables to be read by MATLAB!
using MAT

matfile = matopen("surface_fluxes.mat", "w")

# Ocean variables
uo = centered_surface_u_velocity(ocean_model.velocities.u)
vo = centered_surface_v_velocity(ocean_model.velocities.v)
To = Array(interior(ocean_model.tracers.T, :, :, grid.Nz))
So = Array(interior(ocean_model.tracers.S, :, :, grid.Nz))
write(matfile, "uo", uo)
write(matfile, "vo", vo)
write(matfile, "To", To)
write(matfile, "So", So)

# Atmospheric variables
ua = Array(interior(atmosphere.velocities.u[1], :, :, 1))
va = Array(interior(atmosphere.velocities.v[1], :, :, 1))
Ta = Array(interior(atmosphere.tracers.T[1], :, :, 1))
qa = Array(interior(atmosphere.tracers.q[1], :, :, 1))
ra = Array(interior(atmosphere.tracers.r[1], :, :, 1))
pa = Array(interior(atmosphere.pressure[1], :, :, 1))
Rs = Array(interior(atmosphere.downwelling_radiation.shortwave[1], :, :, 1))
Rl = Array(interior(atmosphere.downwelling_radiation.longwave[1], :, :, 1))
write(matfile, "ua", ua)
write(matfile, "va", va)
write(matfile, "Ta", Ta)
write(matfile, "qa", qa)
write(matfile, "ra", ra)
write(matfile, "pa", pa)
write(matfile, "Rs", Rs)
write(matfile, "Rl", Rl)

import ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: net_downwelling_radiation

@inline net_downwelling_radiation(i::Int, j::Int, k::Int, grid, time, Qs::Field, Ql::Field, args...) = 
        net_downwelling_radiation(i, j,  grid, time, Qs[i, j, k], Ql[i, j, k], args...)

Rad = KernelFunctionOperation{Center, Center, Center}(net_downwelling_radiation, 
                                                      grid, 
                                                      Time(0),
                                                      atmosphere.downwelling_radiation.shortwave[1], 
                                                      atmosphere.downwelling_radiation.longwave[1],
                                                      coupled_model.fluxes.radiation)

Rin = compute!(Field(Rad))
Rin = Array(interior(Rin, :, :, 1))
Rou = 0.97 * 5.67e-8 .* (Array(interior(ocean.model.tracers.T, :, :, 50)) .+ 273.15) .^ 4

# Turbulent fluxes
Ql = Array(interior(coupled_model.fluxes.turbulent.fields.latent_heat,       :, :, 1))
Qs = Array(interior(coupled_model.fluxes.turbulent.fields.sensible_heat,     :, :, 1))
Mv = Array(interior(coupled_model.fluxes.turbulent.fields.water_vapor,       :, :, 1))
tx = Array(interior(coupled_model.fluxes.turbulent.fields.x_momentum,        :, :, 1))
ty = Array(interior(coupled_model.fluxes.turbulent.fields.y_momentum,        :, :, 1))
Tf = Array(interior(ocean.model.tracers.T.boundary_conditions.top.condition, :, :, 1))
write(matfile, "Ql", Ql)
write(matfile, "Qs", Qs)
write(matfile, "Mv", Mv)
write(matfile, "tx", tx)
write(matfile, "ty", ty)
write(matfile, "Tf", Tf)
write(matfile, "Ri", Rin)
write(matfile, "Ro", Rou)

close(matfile)
