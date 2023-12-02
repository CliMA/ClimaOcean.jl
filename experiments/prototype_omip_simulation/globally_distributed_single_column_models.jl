using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.Fields: ConstantField, ZeroField, interpolate!
using Oceananigans.Models.HydrostaticFreeSurfaceModels: ColumnEnsembleSize

using ClimaOcean
using ClimaOcean.OceanSeaIceModels: TwoStreamDownwellingRadiation,
                                    PrescribedAtmosphere,
                                    SurfaceRadiation

using ClimaOcean.JRA55: jra55_field_time_series

using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

using NCDatasets
using GLMakie
using Printf

using Downloads: download

temperature_filename = "THETA.1440x720x50.19920102.nc"
salinity_filename = "SALT.1440x720x50.19920102.nc"
ice_thickness_filename = "SIheff.1440x720.19920102.nc"

# Downloaded from https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N

temperature_url = "https://www.dropbox.com/scl/fi/01h96yo2fhnnvt2zkmu0d/" *
                  "THETA.1440x720x50.19920102.nc?rlkey=ycso2v09gc6v2qb5j0lff0tjs&dl=0"

salinity_url = "https://www.dropbox.com/scl/fi/t068we10j5skphd461zg8/" *
               "SALT.1440x720x50.19920102.nc?rlkey=r5each0ytdtzh5icedvzpe7bw&dl=0"

ice_thickness_url = "https://www.dropbox.com/scl/fi/x0v9gjrfebwsef4tv1dvn/" *
                    "SIheff.1440x720.19920102.nc?rlkey=2uel3jtzbsplr28ejcnx3u6am&dl=0"

isfile(temperature_filename)   || download(temperature_url,   temperature_filename)
isfile(salinity_filename)      || download(salinity_url,      salinity_filename)
isfile(ice_thickness_filename) || download(ice_thickness_url, ice_thickness_filename)

temperature_ds = Dataset(temperature_filename)
salinity_ds = Dataset(salinity_filename)
ice_thickness_ds = Dataset(ice_thickness_filename)

# Construct vertical coordinate
depth = temperature_ds["DEPTH_T"][:]
zc = -reverse(depth)

# Interface depths from cell center depths
zf = (zc[1:end-1] .+ zc[2:end]) ./ 2
push!(zf, 0)

Δz = zc[2] - zc[1]
pushfirst!(zf, zf[1] - Δz)

Tᵢ = temperature_ds["THETA"][:, :, :, 1]
Sᵢ = salinity_ds["SALT"][:, :, :, 1]
ℋᵢ = ice_thickness_ds["SIheff"][:, :, 1]

Nx′, Ny′, Nz = size(Tᵢ)

#####
##### Construct the grid
#####

arch = CPU()

Nx = 1
target_longitude = 180
western_limit = target_longitude
eastern_limit = target_longitude + 1/4
i₁ = 4 * western_limit
i₂ = i₁ + 1

southern_limit = -70
northern_limit = +55
j₁ = 4 * (90 + southern_limit)
j₂ = 720 - 4 * (90 - northern_limit) + 1

Ny = j₂ - j₁ + 1

Tᵢ = Tᵢ[i₁:i₂, j₁:j₂, :]
Sᵢ = Sᵢ[i₁:i₂, j₁:j₂, :]

Tᵢ = convert(Array{Float32, 3}, Tᵢ)
Sᵢ = convert(Array{Float32, 3}, Sᵢ)
ℋᵢ = convert(Array{Float32, 2}, ℋᵢ)

Tᵢ = reverse(Tᵢ, dims=3)
Sᵢ = reverse(Sᵢ, dims=3)

# missing_value = Float32(-9.9e22)

grid = LatitudeLongitudeGrid(arch,
                             size = (Nx, Ny, Nz),
                             halo = (1, 1, 1),
                             longitude = (western_limit, eastern_limit),
                             latitude = (southern_limit, northern_limit),
                             z = zf,
                             topology = (Periodic, Bounded, Bounded))

include("omip_atmosphere.jl")

column_ensemble_size = ColumnEnsembleSize(Nz=Nz, ensemble=(Nx, Ny))
column_ensemble_halo_size = ColumnEnsembleSize(Nz=0, Hz=1)

single_column_grid = RectilinearGrid(arch,
                                     size = column_ensemble_size,
                                     halo = column_ensemble_halo_size,
                                     z = zf,
                                     topology = (Flat, Flat, Bounded))

Ω = Oceananigans.Coriolis.Ω_Earth
f(i, j) = 2Ω * sind(φnode(i, j, 1, grid))
coriolis = [FPlane(f=f(i, j)) for i=1:Nx, j=1:Ny]

top_ocean_heat_flux          = Qᵀ = Field{Center, Center, Nothing}(single_column_grid)
top_salt_flux                = Fˢ = Field{Center, Center, Nothing}(single_column_grid)
top_zonal_momentum_flux      = τˣ = Field{Face, Center, Nothing}(single_column_grid)
top_meridional_momentum_flux = τʸ = Field{Center, Face, Nothing}(single_column_grid)

ocean_boundary_conditions = (u = FieldBoundaryConditions(top=FluxBoundaryCondition(τˣ)),
                             v = FieldBoundaryConditions(top=FluxBoundaryCondition(τʸ)),
                             T = FieldBoundaryConditions(top=FluxBoundaryCondition(Qᵀ)),
                             S = FieldBoundaryConditions(top=FluxBoundaryCondition(Fˢ)))

# Model construction
teos10 = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(equation_of_state=teos10)
closure = CATKEVerticalDiffusivity()

ocean_model = HydrostaticFreeSurfaceModel(; buoyancy, closure, coriolis,
                                          grid = single_column_grid,
                                          tracers = (:T, :S, :e),
                                          boundary_conditions = ocean_boundary_conditions)

#=
#####
##### Setup JRA55 atmosphere
#####

time_indices = 1:10

ua_jra55  = jra55_field_time_series(:eastward_velocity;               time_indices, architecture=arch)
va_jra55  = jra55_field_time_series(:northward_velocity;              time_indices, architecture=arch)
Ta_jra55  = jra55_field_time_series(:temperature;                     time_indices, architecture=arch)
qa_jra55  = jra55_field_time_series(:relative_humidity;               time_indices, architecture=arch)
Fr_jra55  = jra55_field_time_series(:freshwater_rain_flux;            time_indices, architecture=arch)
Fs_jra55  = jra55_field_time_series(:freshwater_snow_flux;            time_indices, architecture=arch)
Qlw_jra55 = jra55_field_time_series(:downwelling_longwave_radiation;  time_indices, architecture=arch)
Qsw_jra55 = jra55_field_time_series(:downwelling_shortwave_radiation; time_indices, architecture=arch)

times = ua_jra55.times

u_bcs = FieldBoundaryConditions(grid, (Face, Center, Nothing))
v_bcs = FieldBoundaryConditions(grid, (Center, Face, Nothing))
c_bcs = FieldBoundaryConditions(grid, (Center, Center, Nothing))

ua_jra55  = FieldTimeSeries{Face,   Center, Nothing}(grid, times; boundary_conditions=u_bcs)
va_jra55  = FieldTimeSeries{Center,   Face, Nothing}(grid, times; boundary_conditions=v_bcs)
Ta_jra55  = FieldTimeSeries{Center, Center, Nothing}(grid, times; boundary_conditions=c_bcs)
qa_jra55  = FieldTimeSeries{Center, Center, Nothing}(grid, times; boundary_conditions=c_bcs)
Fr_jra55  = FieldTimeSeries{Center, Center, Nothing}(grid, times; boundary_conditions=c_bcs)
Fs_jra55  = FieldTimeSeries{Center, Center, Nothing}(grid, times; boundary_conditions=c_bcs)
Qlw_jra55 = FieldTimeSeries{Center, Center, Nothing}(grid, times; boundary_conditions=c_bcs)
Qsw_jra55 = FieldTimeSeries{Center, Center, Nothing}(grid, times; boundary_conditions=c_bcs)

interpolate!(ua_jra55,   u_jra55_native)
interpolate!(va_jra55,   v_jra55_native)
interpolate!(Ta_jra55,   T_jra55_native)
interpolate!(qa_jra55,   q_jra55_native)
interpolate!(Fr_jra55,  Fr_jra55_native)
interpolate!(Fs_jra55,  Fs_jra55_native)
interpolate!(Qlw_jra55, Qlw_jra55_native)
interpolate!(Qsw_jra55, Qsw_jra55_native)

 ua_scm = FieldTimeSeries{Face,   Center, Nothing}(single_column_grid, times)
 va_scm = FieldTimeSeries{Center,   Face, Nothing}(single_column_grid, times)
 Ta_scm = FieldTimeSeries{Center, Center, Nothing}(single_column_grid, times)
 qa_scm = FieldTimeSeries{Center, Center, Nothing}(single_column_grid, times)
 Fr_scm = FieldTimeSeries{Center, Center, Nothing}(single_column_grid, times)
 Fs_scm = FieldTimeSeries{Center, Center, Nothing}(single_column_grid, times)
Qlw_scm = FieldTimeSeries{Center, Center, Nothing}(single_column_grid, times)
Qsw_scm = FieldTimeSeries{Center, Center, Nothing}(single_column_grid, times)

interior(  u_scm, :, :, :) .= interior(  u_jra55, :, :, :)
interior(  v_scm, :, :, :) .= interior(  v_jra55, :, :, :)
interior(  T_scm, :, :, :) .= interior(  T_jra55, :, :, :)
interior(  q_scm, :, :, :) .= interior(  q_jra55, :, :, :)
interior( Fr_scm, :, :, :) .= interior( Fr_jra55, :, :, :)
interior( Fs_scm, :, :, :) .= interior( Fs_jra55, :, :, :)
interior(Qlw_scm, :, :, :) .= interior(Qlw_jra55, :, :, :)
interior(Qsw_scm, :, :, :) .= interior(Qsw_jra55, :, :, :)

velocities = (u = u_scm, v = v_scm)
tracers = (T = T_scm, q = q_scm)
freshwater_flux = (rain = Fr_scm, snow = Fs_scm)
downwelling_radiation = TwoStreamDownwellingRadiation(shortwave=Qsw_scm, longwave=Qsw_scm)
atmosphere = PrescribedAtmosphere(times; velocities, freshwater_flux, tracers, downwelling_radiation)
=#

set!(ocean_model, T=Tᵢ, S=Sᵢ)

surface_radiation = SurfaceRadiation()

coupled_model = OceanSeaIceModel(ocean;
                                 atmosphere = atmosphere
                                 surface_radiation = surface_radiation)

coupled_simulation = Simulation(coupled_model, Δt=5minutes, stop_iteration=2) #stop_time=30days)

run!(coupled_simulation)
