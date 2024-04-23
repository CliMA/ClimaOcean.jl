using Printf
using Oceananigans.Units
using ClimaOcean
using Oceananigans
using Oceananigans.Operators
using ClimaOcean.ECCO2
using ClimaOcean.OceanSimulations
using Oceananigans.Units
using ClimaOcean.JRA55: JRA55_prescribed_atmosphere
using KernelAbstractions: @index, @kernel
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

import Oceananigans.OutputReaders: interpolate!, FieldTimeSeries
using Oceananigans.Grids: architecture, location, node, with_halo
using Oceananigans.BoundaryConditions
using Oceananigans.Utils: launch!, Time
using Oceananigans.Fields: interpolate, instantiate

@kernel function _interpolate_field_time_series_new!(target_fts, target_grid, target_location, target_times,
                                                     source_fts, source_grid, source_location)

    # 4D index, cool!
    i, j, k, n = @index(Global, NTuple)

    target_node = node(i, j, k, target_grid, target_location...)
    at_time     = @inbounds target_times[n]

    @inbounds target_fts[i, j, k, n] = interpolate(target_node, at_time,
                                                   source_fts, source_location, source_grid)
end

function interpolate!(target_fts::FieldTimeSeries, source_fts::FieldTimeSeries)

    target_grid = target_fts.grid
    source_grid = source_fts.grid

    @assert architecture(target_grid) == architecture(source_grid)
    arch = architecture(target_grid)

    # Make locations
    source_location = map(instantiate, location(source_fts))
    target_location = map(instantiate, location(target_fts))

    target_times = target_fts.times

    launch!(arch, target_grid, size(target_fts),
            _interpolate_field_time_series_new!,
            target_fts.data, target_grid, target_location, Time.(target_times),
            source_fts, source_grid, source_location)

    fill_halo_regions!(target_fts)

    return nothing
end

# Upload ECCO2 fields
T = ECCO2.ecco2_field(:temperature)
S = ECCO2.ecco2_field(:salinity)
u = ECCO2.ecco2_field(:u_velocity)
v = ECCO2.ecco2_field(:v_velocity)

include("ecco2_immersed_grid.jl")
grid = ecco2_immersed_grid()

# Let's leave out the radiation for the moment (too simple to test)
atmosphere  = JRA55_prescribed_atmosphere(1:2; backend = InMemory(), grid = grid.underlying_grid)

ocean = ocean_simulation(grid; momentum_advection = nothing,
                                 tracer_advection = nothing)

ocean_model = ocean.model

# setting ecco variables in the model
set!(ocean_model, T = T, S = S, u = u, v = v)

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation = nothing)

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
So = Array(interior(ocean_model.tracers.T, :, :, grid.Nz))
write(matfile, "uo", uo)
write(matfile, "vo", vo)
write(matfile, "To", To)
write(matfile, "So", So)

# Atmospheric variables
ua = Array(interior(atmosphere.velocities.u[1], :, :, 1))
va = Array(interior(atmosphere.velocities.v[1], :, :, 1))
Ta = Array(interior(atmosphere.tracers.T[1], :, :, 1))
qa = Array(interior(atmosphere.tracers.q[1], :, :, 1))
pa = Array(interior(atmosphere.pressure[1], :, :, 1))
write(matfile, "ua", ua)
write(matfile, "va", va)
write(matfile, "Ta", Ta)
write(matfile, "qa", qa)
write(matfile, "pa", pa)

# Turbulent fluxes
Ql = Array(interior(coupled_model.fluxes.turbulent.fields.latent_heat,   :, :, 1))
Qs = Array(interior(coupled_model.fluxes.turbulent.fields.sensible_heat, :, :, 1))
Mv = Array(interior(coupled_model.fluxes.turbulent.fields.water_vapor,   :, :, 1))
tx = Array(interior(coupled_model.fluxes.turbulent.fields.x_momentum,    :, :, 1))
ty = Array(interior(coupled_model.fluxes.turbulent.fields.y_momentum,    :, :, 1))
write(matfile, "Ql", Ql)
write(matfile, "Qs", Qs)
write(matfile, "Mv", Mv)
write(matfile, "tx", tx)
write(matfile, "ty", ty)

close(matfile)
