using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Oceananigans.BuoyancyFormulations: buoyancy_frequency
using Oceananigans.Units: Time
using Dates
using Printf
using GLMakie

location_name = "equatorial_shallow_regions"
z_faces = ExponentialCoordinate(60, -6200)
z_faces = Oceananigans.Grids.MutableVerticalDiscretization(z_faces)
λ★ = 123.21302547843129
φ★ = 9.894610398995553

grid = LatitudeLongitudeGrid(size = (3, 3, 60),
                             longitude = (λ★-0.5, λ★+0.5),
                             latitude  = (φ★-0.5, φ★+0.5),
                             z = z_faces,
                             topology = (Bounded, Bounded, Bounded))

bottom_height = [0.0       0.0    0.0;
                 0.0       0.0    0.0;
              -384.522  -161.595  0.0]

grid  = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

@info "Building the simulation"

ocean = ocean_simulation(grid; 
                         Δt=20minutes, 
                         momentum_advection=nothing, 
                         tracer_advection=nothing)

set!(ocean.model, T=Metadatum(:temperature, dataset=ECCO4Monthly()),
                  S=Metadatum(:salinity,    dataset=ECCO4Monthly()))

atmosphere = JRA55PrescribedAtmosphere(end_date = DateTime(1958, 12, 30), # Last day of the simulation
                                       backend  = JRA55NetCDFBackend(100),
                                       dataset  = MultiYearJRA55(),
                                       dir      = "./")
radiation = Radiation()
coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)
simulation = Simulation(coupled_model, Δt=ocean.Δt, stop_time=7200days)

wall_clock = Ref(time_ns())

function progress(sim)
    msg = "Ocean Station Papa"
    msg *= string(", iter: ", iteration(sim), ", time: ", prettytime(sim))

    elapsed = 1e-9 * (time_ns() - wall_clock[])
    msg *= string(", wall time: ", prettytime(elapsed))
    wall_clock[] = time_ns()

    u, v, w = sim.model.ocean.model.velocities
    msg *= @sprintf(", max|u|: (%.2e, %.2e)", maximum(abs, u), maximum(abs, v))

    T = sim.model.ocean.model.tracers.T
    S = sim.model.ocean.model.tracers.S
    e = if sim.model.ocean.model.closure isa CATKEVerticalDiffusivity
        sim.model.ocean.model.tracers.e
    else
        sim.model.ocean.model.tracers.S
    end

    τx = first(sim.model.interfaces.net_fluxes.ocean_surface.u)
    τy = first(sim.model.interfaces.net_fluxes.ocean_surface.v)
    Q  = first(sim.model.interfaces.net_fluxes.ocean_surface.Q)

    u★ = sqrt(sqrt(τx^2 + τy^2))

    Nz = size(T, 3)
    msg *= @sprintf(", u★: %.2f m s⁻¹", u★)
    msg *= @sprintf(", Q: %.2f W m⁻²",  Q)
    msg *= @sprintf(", T₀: %.2f ᵒC", first(interior(T, 1, 1, Nz)))
    msg *= @sprintf(", extrema(T): (%.2f, %.2f) ᵒC", minimum(T), maximum(T))
    msg *= @sprintf(", S₀: %.2f g/kg", first(interior(S, 1, 1, Nz)))
    msg *= @sprintf(", e₀: %.2e m² s⁻²", first(interior(e, 1, 1, Nz)))

    @info msg

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# Build flux outputs
τx  = simulation.model.interfaces.net_fluxes.ocean_surface.u
τy  = simulation.model.interfaces.net_fluxes.ocean_surface.v
JT  = simulation.model.interfaces.net_fluxes.ocean_surface.T
Js  = simulation.model.interfaces.net_fluxes.ocean_surface.S
E   = simulation.model.interfaces.atmosphere_ocean_interface.fluxes.water_vapor
Qc  = simulation.model.interfaces.atmosphere_ocean_interface.fluxes.sensible_heat
Qv  = simulation.model.interfaces.atmosphere_ocean_interface.fluxes.latent_heat
Qu  = simulation.model.interfaces.atmosphere_ocean_interface.fluxes.upwelling_longwave
Ql  = simulation.model.interfaces.atmosphere_ocean_interface.fluxes.downwelling_longwave
Qs  = simulation.model.interfaces.atmosphere_ocean_interface.fluxes.downwelling_shortwave
ρₒ  = simulation.model.interfaces.ocean_properties.reference_density
cₚ  = simulation.model.interfaces.ocean_properties.heat_capacity
ua  = simulation.model.interfaces.exchanger.exchange_atmosphere_state.u
va  = simulation.model.interfaces.exchanger.exchange_atmosphere_state.v
Ta  = simulation.model.interfaces.exchanger.exchange_atmosphere_state.T
qa  = simulation.model.interfaces.exchanger.exchange_atmosphere_state.q
Fwf = simulation.model.interfaces.exchanger.exchange_atmosphere_state.Mp
u★  = simulation.model.interfaces.atmosphere_ocean_interface.fluxes.friction_velocity
θ★  = simulation.model.interfaces.atmosphere_ocean_interface.fluxes.temperature_scale
q★  = simulation.model.interfaces.atmosphere_ocean_interface.fluxes.water_vapor_scale
Ts  = simulation.model.interfaces.atmosphere_ocean_interface.temperature
qs  = simulation.model.interfaces.atmosphere_ocean_interface.humidity

Q = ρₒ * cₚ * JT
ρτx = ρₒ * τx
ρτy = ρₒ * τy
N² = buoyancy_frequency(ocean.model)
κc = ocean.model.diffusivity_fields.κc

fluxes = (; ρτx, ρτy, E, Js, JT, Q, Qv, Qc, ua, va, Ta, qa, Qs, Ql, Qu, Fwf, u★, θ★, q★, Ts, qs)
auxiliary_fields = (; N², κc)
fields = merge(ocean.model.velocities, ocean.model.tracers, auxiliary_fields)

# Slice fields at the surface
outputs = merge(fields, fluxes)

filename = "single_column_omip_$(location_name)_take2"

simulation.output_writers[:jld2] = JLD2Writer(ocean.model, outputs; filename,
                                              schedule = TimeInterval(8hours),
                                              overwrite_existing = true)

run!(simulation)

filename *= ".jld2"

ua  = FieldTimeSeries(filename, "ua")
va  = FieldTimeSeries(filename, "va")
Ta  = FieldTimeSeries(filename, "Ta")
qa  = FieldTimeSeries(filename, "qa")
Ql  = FieldTimeSeries(filename, "Ql")
Qs  = FieldTimeSeries(filename, "Qs")
Qu  = FieldTimeSeries(filename, "Qu")
Pt  = FieldTimeSeries(filename, "Fwf")

u  = FieldTimeSeries(filename, "u")
v  = FieldTimeSeries(filename, "v")
T  = FieldTimeSeries(filename, "T")
S  = FieldTimeSeries(filename, "S")
e  = if ocean.model.closure isa CATKEVerticalDiffusivity
    FieldTimeSeries(filename, "e")
else
    FieldTimeSeries(filename, "S") # Use salinity as a proxy for e
end

N² = FieldTimeSeries(filename, "N²")
κ  = FieldTimeSeries(filename, "κc")

Qv  = FieldTimeSeries(filename, "Qv")
Qc  = FieldTimeSeries(filename, "Qc")
Q   = FieldTimeSeries(filename, "Q")
Js  = FieldTimeSeries(filename, "Js")
Ev  = FieldTimeSeries(filename, "E")
ρτx = FieldTimeSeries(filename, "ρτx")
ρτy = FieldTimeSeries(filename, "ρτy")
u★  = FieldTimeSeries(filename, "u★")
θ★  = FieldTimeSeries(filename, "θ★")
q★  = FieldTimeSeries(filename, "q★") 
Ts  = FieldTimeSeries(filename, "Ts") 
qs  = FieldTimeSeries(filename, "qs") 

Nz = size(T, 3)
times = Qc.times
ρₒ = coupled_model.interfaces.ocean_properties.reference_density

fig = Figure(size=(1800, 1800))

axτ = Axis(fig[1, 1:3], xlabel="Days since Oct 1 1992", ylabel="Wind stress (N m⁻²)")
axQ = Axis(fig[1, 4:6], xlabel="Days since Oct 1 1992", ylabel="Heat flux (W m⁻²)")
axu = Axis(fig[2, 1:3], xlabel="Days since Oct 1 1992", ylabel="Velocities (m s⁻¹)")
axT = Axis(fig[2, 4:6], xlabel="Days since Oct 1 1992", ylabel="Surface temperature (ᵒC)")
axF = Axis(fig[3, 1:3], xlabel="Days since Oct 1 1992", ylabel="Freshwater volume flux (m s⁻¹)")
axS = Axis(fig[3, 4:6], xlabel="Days since Oct 1 1992", ylabel="Surface salinity (g kg⁻¹)")

axuz = Axis(fig[4:5, 1:2], xlabel="Velocities (m s⁻¹)",                ylabel="z (m)")
axTz = Axis(fig[4:5, 3:4], xlabel="Temperature (ᵒC)",                  ylabel="z (m)")
axSz = Axis(fig[4:5, 5:6], xlabel="Salinity (g kg⁻¹)",                 ylabel="z (m)")
axNz = Axis(fig[6:7, 1:2], xlabel="Buoyancy frequency (s⁻²)",          ylabel="z (m)")
axκz = Axis(fig[6:7, 3:4], xlabel="Eddy diffusivity (m² s⁻¹)",         ylabel="z (m)") #, xscale=log10)
axez = Axis(fig[6:7, 5:6], xlabel="Turbulent kinetic energy (m² s⁻²)", ylabel="z (m)") #, xscale=log10)

title = @sprintf("Single-column simulation at %.2f, %.2f", φ★, λ★)
Label(fig[0, 1:6], title)

n = Observable(1)

times = (times .- times[1]) ./days
Nt = length(times)
tn = @lift times[$n]

colors = Makie.wong_colors()

τx = interior(ρτx, 3, 2, 1, :) ./ ρₒ
τy = interior(ρτy, 3, 2, 1, :) ./ ρₒ

lines!(axu, times, interior(u, 3, 2, Nz, :), color=colors[1], label="Zonal")
lines!(axu, times, interior(v, 3, 2, Nz, :), color=colors[2], label="Meridional")
lines!(axu, times, interior(u★, 3, 2, 1, :), color=colors[3], label="Ocean-side u★")
vlines!(axu, tn, linewidth=4, color=(:black, 0.5))
axislegend(axu)

lines!(axτ, times, interior(ρτx, 3, 2, 1, :), label="Zonal")
lines!(axτ, times, interior(ρτy, 3, 2, 1, :), label="Meridional")
vlines!(axτ, tn, linewidth=4, color=(:black, 0.5))
axislegend(axτ)

lines!(axT, times, interior(Ta, 3, 2, 1, 1:Nt) .- 273.15, color=colors[1], linewidth=2, linestyle=:dash, label="Atmosphere temperature")
lines!(axT, times, interior(T, 3, 2, Nz, :), color=colors[2], linewidth=4, label="Ocean surface temperature")
vlines!(axT, tn, linewidth=4, color=(:black, 0.5))
axislegend(axT)

lines!(axQ, times, interior(Qc, 3, 2, 1, 1:Nt),   color=colors[2], label="Sensible",        linewidth=2)
lines!(axQ, times, interior(Qv, 3, 2, 1, 1:Nt),   color=colors[3], label="Latent",          linewidth=2)
lines!(axQ, times, interior(Qs, 3, 2, 1, 1:Nt),   color=colors[4], label="Shortwave",       linewidth=2)
lines!(axQ, times, interior(Ql, 3, 2, 1, 1:Nt),   color=colors[5], label="Longwave",        linewidth=2)
lines!(axQ, times, interior(Qu, 3, 2, 1, 1:Nt),   color=colors[6], label="Upwelling",       linewidth=2)
lines!(axQ, times, interior(Q, 3, 2, 1, 1:Nt),    color=colors[7], label="Total heat flux", linewidth=4)
vlines!(axQ, tn, linewidth=4, color=(:black, 0.5))
axislegend(axQ)

lines!(axF, times, interior(Pt, 3, 2, 1, 1:Nt),   label="Prescribed freshwater flux")
lines!(axF, times, - interior(Ev, 3, 2, 1, 1:Nt), label="Evaporation")
vlines!(axF, tn, linewidth=4, color=(:black, 0.5))
axislegend(axF)

lines!(axS, times, interior(S, 3, 2, Nz, :))
vlines!(axS, tn, linewidth=4, color=(:black, 0.5))

zc = znodes(T)
zf = znodes(κ)
un  = @lift interior(u[$n],  3, 2, :)
vn  = @lift interior(v[$n],  3, 2, :)
Tn  = @lift interior(T[$n],  3, 2, :)
Sn  = @lift interior(S[$n],  3, 2, :)
κn  = @lift log10.(interior(κ[$n], 3, 2, :))
en  = @lift interior(e[$n],  3, 2, :)
N²n = @lift interior(N²[$n], 3, 2, :)

scatterlines!(axuz, un,  zc, label="u")
scatterlines!(axuz, vn,  zc, label="v")
scatterlines!(axTz, Tn,  zc)
scatterlines!(axSz, Sn,  zc)
scatterlines!(axez, en,  zc)
scatterlines!(axNz, N²n, zf)
scatterlines!(axκz, κn,  zf)
ylims!(axuz, (-400, 0))
ylims!(axTz, (-400, 0))
ylims!(axSz, (-400, 0))
ylims!(axez, (-400, 0))
ylims!(axNz, (-400, 0))
ylims!(axκz, (-400, 0))

xlims!(axTz, (15, 45))
xlims!(axSz, (32, 37))
xlims!(axez, (0, 0.00035))
xlims!(axκz, (-8, 1))
xlims!(axSz, (32, 37))

axislegend(axuz)

record(fig, "single_column_profiles.mp4", 1:Nt, framerate=24) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end
nothing #hide

# ![](single_column_profiles.mp4)
