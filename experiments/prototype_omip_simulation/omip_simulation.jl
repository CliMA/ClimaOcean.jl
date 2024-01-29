using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.Fields: ConstantField, ZeroField, interpolate!

using ClimaOcean
using ClimaOcean.OceanSeaIceModels:
    adjust_ice_covered_ocean_temperature!,
    TwoStreamDownwellingRadiation,
    PrescribedAtmosphere,
    SurfaceRadiation

using ClimaOcean.JRA55: jra55_field_time_series

using ClimaSeaIce
using ClimaSeaIce: IceWaterThermalEquilibrium

using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

using NCDatasets
using GLMakie
using Printf

using Downloads: download

start_time = time_ns()

temperature_filename = "THETA.1440x720x50.19920102.nc"
salinity_filename = "SALT.1440x720x50.19920102.nc"
ice_thickness_filename = "SIheff.1440x720.19920102.nc"
#=

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

Nx, Ny′, Nz = size(Tᵢ)

elapsed = time_ns() - start_time
@info "Initial condition built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

#####
##### Construct the grid
#####

arch = CPU()
southern_limit = -79
northern_limit = -30
j₁ = 4 * (90 + southern_limit)
j₂ = 720 - 4 * (90 - northern_limit) + 1
Ny = j₂ - j₁ + 1

Tᵢ = Tᵢ[:, j₁:j₂, :]
Sᵢ = Sᵢ[:, j₁:j₂, :]
ℋᵢ = ℋᵢ[:, j₁:j₂]

Tᵢ = convert(Array{Float32, 3}, Tᵢ)
Sᵢ = convert(Array{Float32, 3}, Sᵢ)
ℋᵢ = convert(Array{Float32, 2}, ℋᵢ)

Tᵢ = reverse(Tᵢ, dims=3)
Sᵢ = reverse(Sᵢ, dims=3)

missing_value = Float32(-9.9e22)

# Construct bottom_height depth by analyzing T
Nx, Ny, Nz = size(Tᵢ)
bottom_height = ones(Nx, Ny) .* (zf[1] - Δz)

for i = 1:Nx, j = 1:Ny
    @inbounds for k = Nz:-1:1
        if Tᵢ[i, j, k] < -10
            bottom_height[i, j] = zf[k+1]
            break
        end
    end
end

grid = LatitudeLongitudeGrid(arch,
                             size = (Nx, Ny, Nz),
                             longitude = (0, 360),
                             halo = (7, 7, 7),
                             latitude = (southern_limit, northern_limit),
                             z = zf,
                             topology = (Periodic, Bounded, Bounded))

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

elapsed = time_ns() - start_time
@info "Grid constructed including evaluation of bottom height from initial condition data. " *
      prettytime(elapsed * 1e-9)
start_time = time_ns()

# Defines `ocean`, an `Oceananigans.Simulation`
include("omip_ocean_component.jl")

elapsed = time_ns() - start_time
@info "Ocean component built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

ocean_model.clock.time = 0
ocean_model.clock.iteration = 0
set!(ocean_model, T=Tᵢ, S=Sᵢ, e=1e-6)

# Defines `sea_ice`, an `Oceananigans.Simulation`
include("omip_sea_ice_component.jl")

elapsed = time_ns() - start_time
@info "Sea ice component built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

# Defines `atmosphere`, a `ClimaOcean.OceanSeaIceModels.PrescribedAtmosphere`
# also defines `radiation`, a `ClimaOcean.OceanSeaIceModels.Radiation`
include("omip_atmosphere.jl")

elapsed = time_ns() - start_time
@info "Atmosphere built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

surface_radiation = SurfaceRadiation()
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, surface_radiation)
coupled_simulation = Simulation(coupled_model, Δt=5minutes, stop_iteration=2)

elapsed = time_ns() - start_time
@info "Coupled simulation built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

adjust_ice_covered_ocean_temperature!(coupled_model)

wall_clock = Ref(time_ns())

function progress(sim)

    msg1 = string("Iter: ", iteration(sim), ", time: ", prettytime(sim))

    elapsed = 1e-9 * (time_ns() - wall_clock[])
    msg2 = string(", wall time: ", prettytime(elapsed))
    wall_clock[] = time_ns()

    u, v, w = sim.model.ocean.model.velocities
    msg3 = @sprintf(", max|u|: (%.2e, %.2e, %.2e)", 
                    maximum(abs, u),
                    maximum(abs, v),
                    maximum(abs, w))

    @info msg1 * msg2 * msg3
end

coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))

using Oceananigans.Operators: ζ₃ᶠᶠᶜ
u, v, w = ocean_model.velocities
ζ = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, u, v)
a = @at (Center, Center, Center) (u^2 + v^2) / 2
ℋ = sea_ice_model.ice_thickness

# Build flux outputs
Jᵘ = coupled_model.surfaces.ocean.momentum.u
Jᵛ = coupled_model.surfaces.ocean.momentum.v
Jᵀ = coupled_model.surfaces.ocean.tracers.T
F  = coupled_model.surfaces.ocean.tracers.S
ρₒ = coupled_model.ocean_reference_density
cₚ = coupled_model.ocean_heat_capacity

Q = ρₒ * cₚ * Jᵀ
τˣ = ρₒ * Jᵘ
τʸ = ρₒ * Jᵛ

fluxes = (; τˣ, τʸ, Q, F)

fields = merge(ocean_model.velocities, ocean_model.tracers,
               (; ζ = Field(ζ), ℋ ))

# Slice fields at the surface
#fields = NamedTuple(name => view(fields[name], :, :, Nz) for name in keys(fields))
outputs = merge(fields, fluxes)

filename = "omip_surface_fields.jld2"

coupled_simulation.output_writers[:surface] = JLD2OutputWriter(ocean_model, outputs; filename,
                                                               #schedule = TimeInterval(1day),
                                                               schedule = IterationInterval(1),
                                                               indices = (:, :, Nz),
                                                               overwrite_existing = true)

run!(coupled_simulation)
=#

filename = "omip_surface_fields.jld2"

τˣt = FieldTimeSeries(filename, "τˣ")
τʸt = FieldTimeSeries(filename, "τʸ")
Qt  = FieldTimeSeries(filename, "Q")
Ft  = FieldTimeSeries(filename, "F")

Tt = FieldTimeSeries(filename, "T")
St = FieldTimeSeries(filename, "S")
et = FieldTimeSeries(filename, "e")
ζt = FieldTimeSeries(filename, "ζ")

function nan_land!(ψt)
    Nt = length(ψt.times)
    land = interior(ψt[1], :, :, 1) .== 0
    for n = 2:Nt
        ψn = interior(ψt[n], :, :, 1)
        ψn[land] .= NaN
    end
    return nothing
end

for ψt in (τˣt, τʸt, Qt, Ft, Tt, St, et, ζt)
    nan_land!(ψt)
end

λf, φc, zc = nodes(τˣt)
λc, φf, zc = nodes(τʸt)
λc, φc, zc = nodes(Qt)
λf, φf, zc = nodes(ζt)

fig = Figure(size=(2400, 1200))

Nt = length(Tt.times)
slider = Slider(fig[5, 2:3], range=1:Nt, startvalue=1)
n = slider.value #Observable(1)

τˣn = @lift interior(τˣt[$n], :, :, 1)
τʸn = @lift interior(τʸt[$n], :, :, 1)
Qn  = @lift interior(Qt[$n], :, :, 1)
Fn  = @lift interior(Ft[$n], :, :, 1)

Tn  = @lift interior(Tt[$n], :, :, 1)
Sn  = @lift interior(St[$n], :, :, 1)
en  = @lift interior(et[$n], :, :, 1)
ζn  = @lift interior(ζt[$n], :, :, 1)

Qmax = 1000
τmax = 1.0
Fmax = 1e-4

Tmax = 32
Tmin = -2
Smax = 35
Smin = 20
elim = 1e-2
ζlim = 1e-4

τlim = 3τmax / 4
Qlim = 3Qmax / 4
Flim = 3Fmax / 4

axx = Axis(fig[1, 2])
axy = Axis(fig[2, 2])
axQ = Axis(fig[3, 2])
axF = Axis(fig[4, 2])

axT = Axis(fig[1, 3])
axS = Axis(fig[2, 3])
axe = Axis(fig[3, 3])
axz = Axis(fig[4, 3])

hmx = heatmap!(axx, λf, φc, τˣn, colorrange=(-τlim, τlim), colormap=:balance, nan_color=:gray)
hmy = heatmap!(axy, λc, φf, τʸn, colorrange=(-τlim, τlim), colormap=:balance, nan_color=:gray)
hmQ = heatmap!(axQ, λc, φc, Qn,  colorrange=(-Qlim, Qlim), colormap=:balance, nan_color=:gray)
hmF = heatmap!(axF, λc, φc, Fn,  colorrange=(-Flim, Flim), colormap=:balance, nan_color=:gray)

Colorbar(fig[1, 1], hmx, flipaxis=false, label="Eastward momentum flux (N m⁻²)")
Colorbar(fig[2, 1], hmy, flipaxis=false, label="Northward momentum flux (N m⁻²)")
Colorbar(fig[3, 1], hmQ, flipaxis=false, label="Heat flux (W m⁻²)")
Colorbar(fig[4, 1], hmF, flipaxis=false, label="Salt flux (m s⁻¹ psu)")

hmT = heatmap!(axT, λf, φc, Tn, colorrange=(Tmin, Tmax), colormap=:thermal,  nan_color=:gray)
hmS = heatmap!(axS, λc, φf, Sn, colorrange=(Smin, Smax), colormap=:haline,   nan_color=:gray)
hme = heatmap!(axe, λc, φc, en, colorrange=(0, elim),    colormap=:solar,    nan_color=:gray)
hmz = heatmap!(axz, λf, φf, ζn, colorrange=(-ζlim, ζlim), colormap=:balance, nan_color=:gray)

Colorbar(fig[1, 4], hmx, label="Temperature (ᵒC)")
Colorbar(fig[2, 4], hmy, label="Salinity (psu)")
Colorbar(fig[3, 4], hmQ, label="Turbulent kinetic energy (m² s⁻²)")
Colorbar(fig[4, 4], hmF, label="Vorticity (s⁻¹)")

display(fig)

