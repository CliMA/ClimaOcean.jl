# This simulation tests the difference between using a
# `BulkTemperature` and a `SkinTemperature`

using ClimaOcean
using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: SimilarityTheoryTurbulentFluxes, SkinTemperature
using Oceananigans
using Oceananigans.Units
using Oceananigans.BuoyancyModels: buoyancy_frequency
using Oceananigans.Units: Time
using Printf

# Ocean station papa location
location_name = "ocean_station_papa"
λ★, φ★ = 35.1, 50.1

grid = RectilinearGrid(size = 200,
                       x = λ★,
                       y = φ★,
                       z = (-400, 0),
                       topology = (Flat, Flat, Bounded))

ocean1 = ocean_simulation(grid; Δt=10minutes, coriolis=FPlane(latitude = φ★))
ocean2 = ocean_simulation(grid; Δt=10minutes, coriolis=FPlane(latitude = φ★))

# We set initial conditions from ECCO:
set!(ocean1.model.tracers.T, ECCOMetadata(:temperature))
set!(ocean1.model.tracers.S, ECCOMetadata(:salinity))
set!(ocean2.model.tracers.T, ECCOMetadata(:temperature)) 
set!(ocean2.model.tracers.S, ECCOMetadata(:salinity)) 

simulation_days = 31
snapshots_per_day = 8 # corresponding to JRA55's 3-hour frequency
last_time = simulation_days * snapshots_per_day
atmosphere = JRA55_prescribed_atmosphere(1:last_time;
                                         longitude = λ★,
                                         latitude = φ★,
                                         backend = InMemory())
                                         
# We continue constructing a simulation.

radiation = Radiation()
coupled_model_prescribed = OceanSeaIceModel(ocean1; atmosphere, radiation)

similarity_theory = SimilarityTheoryTurbulentFluxes(grid; surface_temperature_type = SkinTemperature(; κ = 0.5))
coupled_model_diagnostic = OceanSeaIceModel(ocean2; atmosphere, radiation, similarity_theory)

for (coupled_model, suffix) in zip([coupled_model_diagnostic, coupled_model_prescribed],
                                   ["diagnostic", "prescribed"])

    simulation = Simulation(coupled_model, Δt=ocean1.Δt, stop_time=10days)

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
        e = sim.model.ocean.model.tracers.e

        τx = first(sim.model.fluxes.total.ocean.momentum.u)
        τy = first(sim.model.fluxes.total.ocean.momentum.v)
        Q = first(sim.model.fluxes.total.ocean.heat)

        u★ = sqrt(sqrt(τx^2 + τy^2))
        Ts = first(interior(sim.model.fluxes.turbulent.fields.T_surface, 1, 1, 1))

        Nz = size(T, 3)
        msg *= @sprintf(", u★: %.2f m s⁻¹", u★)
        msg *= @sprintf(", Q: %.2f W m⁻²",  Q)
        msg *= @sprintf(", T₀: %.2f ᵒC", Ts)
        msg *= @sprintf(", S₀: %.2f g/kg", first(interior(S, 1, 1, Nz)))
        msg *= @sprintf(", e₀: %.2e m² s⁻²", first(interior(e, 1, 1, Nz)))

        @info msg
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

    # Build flux outputs
    τx = coupled_model.fluxes.total.ocean.momentum.u
    τy = coupled_model.fluxes.total.ocean.momentum.v
    JT = coupled_model.fluxes.total.ocean.tracers.T
    Js = coupled_model.fluxes.total.ocean.tracers.S
    E  = coupled_model.fluxes.turbulent.fields.water_vapor
    Qc = coupled_model.fluxes.turbulent.fields.sensible_heat
    Qv = coupled_model.fluxes.turbulent.fields.latent_heat
    Ts = coupled_model.fluxes.turbulent.fields.T_surface
    ρₒ = coupled_model.fluxes.ocean_reference_density
    cₚ = coupled_model.fluxes.ocean_heat_capacity

    Q = ρₒ * cₚ * JT
    ρτx = ρₒ * τx
    ρτy = ρₒ * τy
    N² = buoyancy_frequency(coupled_model.ocean.model)
    κc = coupled_model.ocean.model.diffusivity_fields.κc

    fluxes = (; ρτx, ρτy, E, Js, Qv, Qc, Ts)
    auxiliary_fields = (; N², κc)
    fields = merge(coupled_model.ocean.model.velocities, coupled_model.ocean.model.tracers, auxiliary_fields)

    # Slice fields at the surface
    outputs = merge(fields, fluxes)

    filename = "single_column_omip_$(location_name)_$(suffix)"

    simulation.output_writers[:jld2] = JLD2OutputWriter(coupled_model.ocean.model, outputs; filename,
                                                        schedule = TimeInterval(3hours),
                                                        overwrite_existing = true)

    run!(simulation)
end

#####
##### Visualization
#####

filename_prescribed = "single_column_omip_$(location_name)_prescribed.jld2"
filename_diagnostic = "single_column_omip_$(location_name)_diagnostic.jld2"

# Diagnosed
ud  = FieldTimeSeries(filename_diagnostic, "u")
vd  = FieldTimeSeries(filename_diagnostic, "v")
Td  = FieldTimeSeries(filename_diagnostic, "v")
ed  = FieldTimeSeries(filename_diagnostic, "T")
Sd  = FieldTimeSeries(filename_diagnostic, "S")
Nd² = FieldTimeSeries(filename_diagnostic, "N²")
κd  = FieldTimeSeries(filename_diagnostic, "κc")

Tsd = FieldTimeSeries(filename_diagnostic, "Ts")
Qvd  = FieldTimeSeries(filename_diagnostic, "Qv")
Qcd  = FieldTimeSeries(filename_diagnostic, "Qc")
Jsd  = FieldTimeSeries(filename_diagnostic, "Js")
Evd  = FieldTimeSeries(filename_diagnostic, "E")
ρτxd = FieldTimeSeries(filename_diagnostic, "ρτx")
ρτyd = FieldTimeSeries(filename_diagnostic, "ρτy")

# Prescribed
up  = FieldTimeSeries(filename_prescribed, "u")
vp  = FieldTimeSeries(filename_prescribed, "v")
Tp  = FieldTimeSeries(filename_prescribed, "T")
Sp  = FieldTimeSeries(filename_prescribed, "S")
ep  = FieldTimeSeries(filename_prescribed, "e")
Np² = FieldTimeSeries(filename_prescribed, "N²")
κp  = FieldTimeSeries(filename_prescribed, "κc")

Tsp  = FieldTimeSeries(filename_prescribed, "Ts")
Qvp  = FieldTimeSeries(filename_prescribed, "Qv")
Qcp  = FieldTimeSeries(filename_prescribed, "Qc")
Jsp  = FieldTimeSeries(filename_prescribed, "Js")
Evp  = FieldTimeSeries(filename_prescribed, "E")
ρτxp = FieldTimeSeries(filename_prescribed, "ρτx")
ρτyp = FieldTimeSeries(filename_prescribed, "ρτy")

Nz = size(ud, 3)
times = Qcd.times

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
axκz = Axis(fig[6:7, 3:4], xlabel="Eddy diffusivity (m² s⁻¹)",         ylabel="z (m)", xscale=log10)
axez = Axis(fig[6:7, 5:6], xlabel="Turbulent kinetic energy (m² s⁻²)", ylabel="z (m)", xscale=log10)

title = @sprintf("Single-column simulation at %.2f, %.2f", φ★, λ★)
Label(fig[0, 1:6], title)

n = Observable(1)

times = (times .- times[1]) ./days
Nt = length(times)
tn = @lift times[$n]

colors = Makie.wong_colors()

ρₒ = coupled_model_prescribed.fluxes.ocean_reference_density
τxp = interior(ρτxp, 1, 1, 1, :) ./ ρₒ
τyp = interior(ρτyp, 1, 1, 1, :) ./ ρₒ
up★ = @. (τxp^2 + τyp^2)^(1/4)
τxd = interior(ρτxd, 1, 1, 1, :) ./ ρₒ
τyd = interior(ρτyd, 1, 1, 1, :) ./ ρₒ
ud★ = @. (τxd^2 + τyd^2)^(1/4)

lines!(axu, times, interior(ud, 1, 1, Nz, :), color=colors[1], label="Zonal")
lines!(axu, times, interior(vd, 1, 1, Nz, :), color=colors[2], label="Meridional")
lines!(axu, times, interior(up, 1, 1, Nz, :), color=colors[1], linestyle = :dash, label="Zonal")
lines!(axu, times, interior(vp, 1, 1, Nz, :), color=colors[2], linestyle = :dash, label="Meridional")
lines!(axu, times, ud★, color=colors[3], label="Ocean-side u★") 
lines!(axu, times, up★, color=colors[3], label="Ocean-side u★", linestyle = :dash) 
vlines!(axu, tn, linewidth=4, color=(:black, 0.5))
axislegend(axu)

lines!(axτ, times, interior(ρτxd, 1, 1, 1, :), label="Zonal")
lines!(axτ, times, interior(ρτyd, 1, 1, 1, :), label="Meridional")
lines!(axτ, times, interior(ρτxp, 1, 1, 1, :), linestyle=:dash, label="Zonal")
lines!(axτ, times, interior(ρτyp, 1, 1, 1, :), linestyle=:dash, label="Meridional")
vlines!(axτ, tn, linewidth=4, color=(:black, 0.5))
axislegend(axτ)

lines!(axT, times, interior(Tsd, 1, 1, 1, :), color=colors[2], linewidth=4, label="Ocean surface temperature")
lines!(axT, times, interior(Tsp, 1, 1, 1, :), color=colors[2], linewidth=4, linestyle = :dash, label="Ocean surface temperature")
vlines!(axT, tn, linewidth=4, color=(:black, 0.5))
axislegend(axT)

lines!(axQ, times, interior(Qvd, 1, 1, 1, 1:Nt), color=colors[2], label="Sensible",  linewidth=2)
lines!(axQ, times, interior(Qcd, 1, 1, 1, 1:Nt), color=colors[3], label="Latent",    linewidth=2)
lines!(axQ, times, interior(Qvp, 1, 1, 1, 1:Nt), color=colors[2], linestyle=:dash, label="Sensible",  linewidth=2)
lines!(axQ, times, interior(Qcp, 1, 1, 1, 1:Nt), color=colors[3], linestyle=:dash, label="Latent",    linewidth=2)
vlines!(axQ, tn, linewidth=4, color=(:black, 0.5))
axislegend(axQ)

lines!(axF, times, - interior(Evd, 1, 1, 1, 1:Nt), label="Evaporation")
lines!(axF, times, - interior(Evd, 1, 1, 1, 1:Nt), linestyle=:dash, label="Evaporation")
vlines!(axF, tn, linewidth=4, color=(:black, 0.5))
axislegend(axF)

lines!(axS, times, interior(Sd, 1, 1, Nz, :))
lines!(axS, times, interior(Sp, 1, 1, Nz, :), linestyle=:dash)
vlines!(axS, tn, linewidth=4, color=(:black, 0.5))

zc = znodes(Tp)
zf = znodes(κp)
upn  = @lift interior(up[$n],  1, 1, :)
vpn  = @lift interior(vp[$n],  1, 1, :)
Tpn  = @lift interior(Tp[$n],  1, 1, :)
Spn  = @lift interior(Sp[$n],  1, 1, :)
κpn  = @lift interior(κp[$n],  1, 1, :)
epn  = @lift interior(ep[$n],  1, 1, :)
Np²n = @lift interior(Np²[$n], 1, 1, :)

udn  = @lift interior(ud[$n],  1, 1, :)
vdn  = @lift interior(vd[$n],  1, 1, :)
Tdn  = @lift interior(Td[$n],  1, 1, :)
Sdn  = @lift interior(Sd[$n],  1, 1, :)
κdn  = @lift interior(κd[$n],  1, 1, :)
edn  = @lift interior(ed[$n],  1, 1, :)
Nd²n = @lift interior(Nd²[$n], 1, 1, :)

scatterlines!(axuz, udn,  zc, label="u") 
scatterlines!(axuz, vdn,  zc, label="v") 
scatterlines!(axTz, Tdn,  zc) 
scatterlines!(axSz, Sdn,  zc) 
scatterlines!(axez, edn,  zc) 
scatterlines!(axNz, Nd²n, zf) 
scatterlines!(axκz, κdn,  zf) 
scatterlines!(axuz, upn,  zc, linestyle=:dash, label="u") 
scatterlines!(axuz, vpn,  zc, linestyle=:dash, label="v") 
scatterlines!(axTz, Tpn,  zc, linestyle=:dash) 
scatterlines!(axSz, Spn,  zc, linestyle=:dash) 
scatterlines!(axez, epn,  zc, linestyle=:dash) 
scatterlines!(axNz, Np²n, zf, linestyle=:dash) 
scatterlines!(axκz, κpn,  zf, linestyle=:dash) 

axislegend(axuz)

ulim = max(maximum(abs, up), maximum(abs, vp))
xlims!(axuz, -ulim, ulim)

Tmax = maximum(interior(Tp))
Tmin = minimum(interior(Tp))
xlims!(axTz, Tmin - 0.1, Tmax + 0.1)

Nmax = maximum(interior(Np²))
xlims!(axNz, -Nmax/10, Nmax * 1.05)

κmax = maximum(interior(κp))
xlims!(axκz, 1e-9, κmax * 1.1)

emax = maximum(interior(ep))
xlims!(axez, 1e-11, emax * 1.1)

Smax = maximum(interior(Sp))
Smin = minimum(interior(Sp))
xlims!(axSz, Smin - 0.2, Smax + 0.2)

record(fig, "single_column_profiles.mp4", 1:Nt, framerate=24) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end
nothing #hide
# ![](single_column_profiles.mp4)
