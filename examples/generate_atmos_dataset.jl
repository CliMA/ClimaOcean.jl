using ClimaOcean
using Oceananigans
using JLD2

time_indices = 1:1

qt = ClimaOcean.JRA55.jra55_field_time_series(:specific_humidity; time_indices)
Tt = ClimaOcean.JRA55.jra55_field_time_series(:temperature; time_indices)
pt = ClimaOcean.JRA55.jra55_field_time_series(:sea_level_pressure; time_indices)

Nx, Ny, Nz = size(qt[1])

q = Array(interior(qt[1], 1:4:Nx, 1:4:Ny, 1))
T = Array(interior(Tt[1], 1:4:Nx, 1:4:Ny, 1))
p = Array(interior(pt[1], 1:4:Nx, 1:4:Ny, 1))

@save "atmospheric_state.jld2" q T p

