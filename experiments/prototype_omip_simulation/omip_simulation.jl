using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using ClimaOcean
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
hᵢ = ice_thickness_ds["SIheff"][:, :, 1]

Tᵢ = convert(Array{Float32, 3}, Tᵢ)
Sᵢ = convert(Array{Float32, 3}, Sᵢ)

Tᵢ = reverse(Tᵢ, dims=3)
Sᵢ = reverse(Sᵢ, dims=3)

missing_value = Float32(-9.9e22)

# Construct bottom_height depth by analyzing T
Nx, Ny′, Nz = size(Tᵢ)
bottom_height = ones(Nx, Ny′) .* (zf[1] - Δz)

for i = 1:Nx, j = 1:Ny′
    @inbounds for k = Nz:-1:1
        if Tᵢ[i, j, k] < -10
            bottom_height[i, j] = zf[k+1]
            break
        end
    end
end

#####
##### Construct the grid
#####

arch = GPU()
southern_limit = -79
northern_limit = 85
j₁ = 4 * (90 + southern_limit)
j₂ = 720 - 4 * (90 - northern_limit) + 1
Ny = j₂ - j₁ + 1

Tᵢ = Tᵢ[:, j₁:j₂, :]
Sᵢ = Sᵢ[:, j₁:j₂, :]
bottom_height = bottom_height[:, j₁:j₂]


grid = LatitudeLongitudeGrid(arch,
                             size = (Nx, Ny, Nz),
                             longitude = (0, 360),
                             halo = (7, 7, 7),
                             latitude = (southern_limit, northern_limit),
                             z = zf,
                             topology = (Periodic, Bounded, Bounded))

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

#####
##### Setup ice model
#####

ice_ocean_heat_flux      = Field{Center, Center, Nothing}(grid)
top_ocean_heat_flux = Qᵀ = Field{Center, Center, Nothing}(grid)
top_salt_flux       = Qˢ = Field{Center, Center, Nothing}(grid)

ocean_boundary_conditions = (T = FieldBoundaryConditions(top=FluxBoundaryCondition(Qᵀ)),
                             S = FieldBoundaryConditions(top=FluxBoundaryCondition(Qˢ)))

#=
using Oceananigans.Grids: φnodes, λnodes

λ = λnodes(grid, Center())
φ = φnodes(grid, Center())

fig = Figure()
ax = Axis(fig[1, 1])
heatmap!(ax, λ, φ, bottom_height)
=#

# Model construction
teos10 = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(equation_of_state=teos10)

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: MixingLength
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: TurbulentKineticEnergyEquation

mixing_length = MixingLength(Cᵇ=0.01)
turbulent_kinetic_energy_equation = TurbulentKineticEnergyEquation(Cᵂϵ=1.0)
closure = CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation)

tracer_advection = WENO()
momentum_advection = VectorInvariant(vorticity_scheme = WENO(),
                                     divergence_scheme = WENO(),
                                     vertical_scheme = WENO())

ocean_model = HydrostaticFreeSurfaceModel(; grid, buoyancy, closure,
                                          tracer_advection, momentum_advection,
                                          tracers = (:T, :S, :e),
                                          free_surface = SplitExplicitFreeSurface(cfl=0.2; grid),
                                          boundary_conditions = ocean_boundary_conditions,
                                          coriolis = HydrostaticSphericalCoriolis())

set!(ocean_model, T=Tᵢ, S=Sᵢ)

#####
##### Setup ice model
#####

Nz = size(grid, 3)
So = ocean_model.tracers.S
ocean_surface_salinity = view(So, :, :, Nz)
bottom_bc = IceWaterThermalEquilibrium(ConstantField(30)) #ocean_surface_salinity)

u, v, w = ocean_model.velocities
ocean_surface_velocities = (u = view(u, :, :, Nz), #interior(u, :, :, Nz),
                            v = view(v, :, :, Nz), #interior(v, :, :, Nz),
                            w = ZeroField())

ice_model = SlabSeaIceModel(ice_grid;
                            velocities = ocean_surface_velocities,
                            advection = WENO(),
                            ice_consolidation_thickness = 0.05,
                            ice_salinity = 4,
                            internal_thermal_flux = ConductiveFlux(conductivity=2),
                            top_thermal_flux = ConstantField(0), # W m⁻²
                            top_thermal_boundary_condition = PrescribedTemperature(-2),
                            bottom_thermal_boundary_condition = bottom_bc,
                            bottom_thermal_flux = ice_ocean_heat_flux)

set!(ice_model, h=hᵢ) 

#simulation = Simulation(model, Δt=5minutes, stop_time=360days)
ocean_simulation = Simulation(ocean_model; Δt=5minutes, verbose=false)
ice_simulation = Simulation(ice_model, Δt=5minutes, verbose=false)
coupled_model = IceOceanModel(ice_simulation, ocean_simulation)
coupled_simulation = Simulation(coupled_model, Δt=5minutes, stop_time=5days)

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

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

using Oceananigans.Operators: ζ₃ᶠᶠᶜ
u, v, w = ocean_model.velocities
ζ = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, u, v)

outputs = merge(ocean_model.velocities, ocean_model.tracers, (; ζ))
filename = "omip_surface_fields.jld2"

ocean_simulation.output_writers[:surface] = JLD2OutputWriter(ocean_model, outputs; filename,
                                                             schedule = TimeInterval(12hour),
                                                             indices = (:, :, Nz),
                                                             overwrite_existing = true)

run!(coupled_simulation)

#####
##### Visualize
#####

et = FieldTimeSeries(filename, "e")
ζt = FieldTimeSeries(filename, "ζ")
times = et.times
Nt = length(times)

fig = Figure(resolution=(1800, 1200))

axe = Axis(fig[1, 1])
axZ = Axis(fig[2, 1])
slider = Slider(fig[3, 1], range=1:Nt, startvalue=1)
n = slider.value

en = @lift interior(et[$n], :, :, 1)
ζn = @lift interior(ζt[$n], :, :, 1)

#heatmap!(axT, λ, φ, Tn, colorrange=(-1, 25), colormap=:thermal)
#heatmap!(axS, λ, φ, Sn, colorrange=(28, 35), colormap=:haline)
heatmap!(axe, λ, φ, en, colorrange=(0, 1e-4), colormap=:solar)
heatmap!(axZ, λ, φ, ζn, colorrange=(-2e-5, 2e-5), colormap=:redblue)

display(fig)

