using Oceananigans
using Oceananigans.Units
using Oceananigans.BuoyancyModels: buoyancy_frequency
using Oceananigans.Units: Time

using ClimaOcean
using ClimaOcean.OceanSeaIceModels: Radiation
using ClimaOcean.DataWrangling.JRA55: jra55_prescribed_atmosphere
using ClimaOcean.DataWrangling.ECCO2: ecco2_field

using GLMakie
using Printf
using Dates

locations = (
    #eastern_mediterranean = (λ =  30, φ = 32), 
    ocean_station_papa = (λ = 215, φ = 50), 
    north_atlantic = (λ = 325, φ = 50), 
    drake_passage = (λ = 300, φ = -60), 
    weddell_sea = (λ = 325, φ = -70), 
    tasman_southern_ocean = (λ = 145, φ = -55), 
)

for location in keys(locations)

#location = :ocean_station_papa
    start_time = time_ns()

    include("omip_ocean_component.jl")

    epoch = Date(1992, 1, 1)
    date = Date(1992, 10, 1)
    start_seconds = Second(date - epoch).value
    uᵢ = ecco2_field(:u_velocity, date)
    vᵢ = ecco2_field(:v_velocity, date)
    Tᵢ = ecco2_field(:temperature, date)
    Sᵢ = ecco2_field(:salinity, date)

    land = interior(Tᵢ) .< -10
    interior(Tᵢ)[land] .= NaN
    interior(Sᵢ)[land] .= NaN

    teos10 = TEOS10EquationOfState()
    buoyancy = SeawaterBuoyancy(equation_of_state=teos10)
    tracers = (T=Tᵢ, S=Sᵢ)
    N²_op = buoyancy_frequency(buoyancy, Tᵢ.grid, tracers)
    N² = Field(N²_op)
    compute!(N²)

    elapsed = time_ns() - start_time
    @info "Initial condition built. " * prettytime(elapsed * 1e-9)
    start_time = time_ns()

    #####
    ##### Construct the grid
    #####

    zc = znodes(Tᵢ)
    zf = znodes(N²)

    arch = CPU()

    Δ = 1/4 # resolution in degrees
    φ₁ = -90 + Δ/2
    φ₂ = +90 - Δ/2
    λ₁ = 0   + Δ/2
    λ₂ = 360 - Δ/2
    φe = φ₁:Δ:φ₂
    λe = λ₁:Δ:λ₂

    λ★, φ★ = locations[location]

    i★ = searchsortedfirst(λe, λ★)
    j★ = searchsortedfirst(φe, φ★)

    longitude = (λe[i★] - Δ/2, λe[i★] + Δ/2)
    latitude  = (φe[j★] - Δ/2, φe[j★] + Δ/2)

    # Column
    uc = interior(uᵢ, i★:i★, j★:j★, :)
    vc = interior(vᵢ, i★:i★, j★:j★, :)
    Tc = interior(Tᵢ, i★:i★, j★:j★, :)
    Sc = interior(Sᵢ, i★:i★, j★:j★, :)

    # Find bottom
    zm = -400
    zf = znodes(Tᵢ.grid, Face())
    kb = findlast(T -> T < -20, Tc[1, 1, :])
    km = findlast(z -> z < zm, zf)
    k★ = isnothing(kb) ? km : max(kb + 1, km)

    Nz = size(Tc, 3)
    kf = k★:Nz+1
    kc = k★:Nz
    zf = zf[kf]
    uc = uc[:, :, kc]
    vc = vc[:, :, kc]
    Tc = Tc[:, :, kc]
    Sc = Sc[:, :, kc]
    Nz′ = length(kc)

    grid = LatitudeLongitudeGrid(arch; longitude, latitude,
                                 size = (1, 1, Nz′),
                                 z = zf,
                                 topology = (Periodic, Periodic, Bounded))

    elapsed = time_ns() - start_time
    @info "Grid constructed. " * prettytime(elapsed * 1e-9)
    start_time = time_ns()

    ocean = omip_ocean_component(grid)
    elapsed = time_ns() - start_time
    @info "Ocean component built. " * prettytime(elapsed * 1e-9)
    start_time = time_ns()

    Ndays = 365
    Nt = 8 * Ndays
    atmosphere = jra55_prescribed_atmosphere(grid, 1:Nt) #, 1:21)
    elapsed = time_ns() - start_time
    @info "Atmosphere built. " * prettytime(elapsed * 1e-9)
    start_time = time_ns()

    ocean.model.clock.time = start_seconds
    ocean.model.clock.iteration = 0
    set!(ocean.model, T=Tc, S=Sc, e=1e-6)

    ua = atmosphere.velocities.u
    va = atmosphere.velocities.v
    Ta = atmosphere.tracers.T
    qa = atmosphere.tracers.q
    times = ua.times

    #=
    fig = Figure(resolution=(1200, 1800))
    axu = Axis(fig[1, 1])
    axT = Axis(fig[2, 1])
    axq = Axis(fig[3, 1])

    lines!(axu, times ./ days, interior(ua, 1, 1, 1, :))
    lines!(axu, times ./ days, interior(va, 1, 1, 1, :))
    lines!(axT, times ./ days, interior(Ta, 1, 1, 1, :))
    lines!(axq, times ./ days, interior(qa, 1, 1, 1, :))

    display(fig)
    =#

    sea_ice = nothing
    radiation = Radiation()
    coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

    coupled_simulation = Simulation(coupled_model, Δt=10minutes, stop_time=start_seconds + 60days)

    elapsed = time_ns() - start_time
    @info "Coupled simulation built. " * prettytime(elapsed * 1e-9)
    start_time = time_ns()

    wall_clock = Ref(time_ns())

    function progress(sim)
        msg = string("(", location, ")")
        msg *= string(", iter: ", iteration(sim), ", time: ", prettytime(sim))

        elapsed = 1e-9 * (time_ns() - wall_clock[])
        msg *= string(", wall time: ", prettytime(elapsed))
        wall_clock[] = time_ns()

        u, v, w = sim.model.ocean.model.velocities
        msg *= @sprintf(", max|u|: (%.2e, %.2e)", maximum(abs, u), maximum(abs, v))

        T = sim.model.ocean.model.tracers.T
        S = sim.model.ocean.model.tracers.S
        e = sim.model.ocean.model.tracers.e

        τˣ = first(sim.model.fluxes.total.ocean.momentum.τˣ)
        τʸ = first(sim.model.fluxes.total.ocean.momentum.τʸ)
        u★ = (τˣ^2 + τʸ^2)^(1/4)
        Q = first(sim.model.fluxes.total.ocean.heat)

        Nz = size(T, 3)
        msg *= @sprintf(", u★: %.2f m s⁻¹", u★)
        msg *= @sprintf(", Q: %.2f W m⁻²", Q)
        msg *= @sprintf(", T₀: %.2f ᵒC",     first(interior(T, 1, 1, Nz)))
        msg *= @sprintf(", extrema(T): (%.2f, %.2f) ᵒC", minimum(T), maximum(T))
        msg *= @sprintf(", S₀: %.2f g/kg",   first(interior(S, 1, 1, Nz)))
        msg *= @sprintf(", e₀: %.2e m² s⁻²", first(interior(e, 1, 1, Nz)))

        @info msg
    end

    coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

    # Build flux outputs
    Jᵘ  = coupled_model.fluxes.total.ocean.momentum.u
    Jᵛ  = coupled_model.fluxes.total.ocean.momentum.v
    Jᵀ  = coupled_model.fluxes.total.ocean.tracers.T
    Jˢ  = coupled_model.fluxes.total.ocean.tracers.S
    E   = coupled_model.fluxes.turbulent.fields.freshwater
    Qse = coupled_model.fluxes.turbulent.fields.sensible_heat
    Qla = coupled_model.fluxes.turbulent.fields.latent_heat
    ρₒ  = coupled_model.fluxes.ocean_reference_density
    cₚ  = coupled_model.fluxes.ocean_heat_capacity

    Q = ρₒ * cₚ * Jᵀ
    τˣ = ρₒ * Jᵘ
    τʸ = ρₒ * Jᵛ
    N² = buoyancy_frequency(ocean.model)
    κᶜ = ocean.model.diffusivity_fields.κᶜ

    fluxes = (; τˣ, τʸ, E, F, Q, Qse, Qla)

    auxiliary_fields = (; N², κᶜ)
    fields = merge(ocean.model.velocities, ocean.model.tracers, auxiliary_fields)

    # Slice fields at the surface
    outputs = merge(fields, fluxes)

    output_attributes = Dict(
        "κᶜ"  => Dict("long_name" => "Tracer diffusivity",          "units" => "m^2 / s"),
        "Q"   => Dict("long_name" => "Net heat flux",               "units" => "W / m^2", "convention" => "positive upwards"),
        "Qla" => Dict("long_name" => "Latent heat flux",            "units" => "W / m^2", "convention" => "positive upwards"),
        "Qse" => Dict("long_name" => "Sensible heat flux",          "units" => "W / m^2", "convention" => "positive upwards"),
        "Jˢ"  => Dict("long_name" => "Salt flux",                   "units" => "g kg⁻¹ m s⁻¹", "convention" => "positive upwards"),
        "E"   => Dict("long_name" => "Freshwater evaporation flux", "units" => "m s⁻¹", "convention" => "positive upwards"),
        "e"   => Dict("long_name" => "Turbulent kinetic energy",    "units" => "m^2 / s^2"),
        "τˣ"  => Dict("long_name" => "Zonal momentum flux",         "units" => "m^2 / s^2"),
        "τʸ"  => Dict("long_name" => "Meridional momentum flux",    "units" => "m^2 / s^2"),
    )

    filename = "single_column_omip_$location"

    coupled_simulation.output_writers[:jld2] = JLD2OutputWriter(ocean.model, outputs; filename,
                                                                schedule = TimeInterval(3hours),
                                                                overwrite_existing = true)

    coupled_simulation.output_writers[:nc] = NetCDFOutputWriter(ocean.model, outputs; filename,
                                                                schedule = AveragedTimeInterval(1days),
                                                                output_attributes,
                                                                overwrite_existing = true)

    run!(coupled_simulation)

    filename *= ".jld2"

    ut = FieldTimeSeries(filename, "u")
    vt = FieldTimeSeries(filename, "v")
    Tt = FieldTimeSeries(filename, "T")
    St = FieldTimeSeries(filename, "S")
    et = FieldTimeSeries(filename, "e")
    N²t = FieldTimeSeries(filename, "N²")
    κt = FieldTimeSeries(filename, "κᶜ")

    Qt = FieldTimeSeries(filename, "Q")
    Qset = FieldTimeSeries(filename, "Qse")
    Qlat = FieldTimeSeries(filename, "Qla")
    Ft = FieldTimeSeries(filename, "Jˢ")
    Et = FieldTimeSeries(filename, "E")
    τˣt = FieldTimeSeries(filename, "τˣ")
    τʸt = FieldTimeSeries(filename, "τʸ")

    Nz = size(Tt, 3)
    times = Qt.times

    ua = atmosphere.velocities.u
    va = atmosphere.velocities.v
    Ta = atmosphere.tracers.T
    qa = atmosphere.tracers.q
    Qlw = atmosphere.downwelling_radiation.longwave
    Qsw = atmosphere.downwelling_radiation.shortwave
    Pr = atmosphere.freshwater_flux.rain
    Ps = atmosphere.freshwater_flux.snow

    Nt = length(times)
    uat = zeros(Nt)
    vat = zeros(Nt)
    Tat = zeros(Nt)
    qat = zeros(Nt)
    Qswt = zeros(Nt)
    Qlwt = zeros(Nt)
    Pt = zeros(Nt)

    for n = 1:Nt
        t = times[n]
        uat[n] = ua[1, 1, 1, Time(t)]
        vat[n] = va[1, 1, 1, Time(t)]
        Tat[n] = Ta[1, 1, 1, Time(t)]
        qat[n] = qa[1, 1, 1, Time(t)]
        Qswt[n] = Qsw[1, 1, 1, Time(t)]
        Qlwt[n] = Qlw[1, 1, 1, Time(t)]
        Pt[n] = Pr[1, 1, 1, Time(t)] + Ps[1, 1, 1, Time(t)]
    end

    set_theme!(Theme(linewidth=3))

    fig = Figure(resolution=(2400, 1800))

    axτ = Axis(fig[1, 1:2], xlabel="Days since Oct 1 1992", ylabel="Wind stress (N m⁻²)")
    axu = Axis(fig[2, 1:2], xlabel="Days since Oct 1 1992", ylabel="Velocities (m s⁻¹)")
    axQ = Axis(fig[1, 3:4], xlabel="Days since Oct 1 1992", ylabel="Heat flux (W m⁻²)")
    axT = Axis(fig[2, 3:4], xlabel="Days since Oct 1 1992", ylabel="Surface temperature (ᵒC)")
    axF = Axis(fig[1, 5:6], xlabel="Days since Oct 1 1992", ylabel="Freshwater volume flux (m s⁻¹)")
    axS = Axis(fig[2, 5:6], xlabel="Days since Oct 1 1992", ylabel="Surface salinity (g kg⁻¹)")

    axuz = Axis(fig[3, 1], xlabel="Velocities (m s⁻¹)",                ylabel="z (m)")
    axTz = Axis(fig[3, 2], xlabel="Temperature (ᵒC)",                  ylabel="z (m)")
    axSz = Axis(fig[3, 3], xlabel="Salinity (g kg⁻¹)",                 ylabel="z (m)")
    axNz = Axis(fig[3, 4], xlabel="Buoyancy frequency (s⁻²)",          ylabel="z (m)")
    axκz = Axis(fig[3, 5], xlabel="Eddy diffusivity (m² s⁻¹)",         ylabel="z (m)", xscale=log10)
    axez = Axis(fig[3, 6], xlabel="Turbulent kinetic energy (m² s⁻²)", ylabel="z (m)", xscale=log10)

    title = @sprintf("Single column simulation at %.2f, %.2f", φ★, λ★)
    Label(fig[0, 1:6], title)

    slider = Slider(fig[4, 1:6], range=1:Nt, startvalue=1)
    n = slider.value

    times = (times .- times[1]) ./days
    tn = @lift times[$n]

    colors = Makie.wong_colors()

    #lines!(axu, times, uat, color=colors[1])
    #lines!(axu, times, vat, color=colors[2])

    ρₒ = coupled_model.fluxes.ocean_reference_density
    Jᵘt = interior(τˣt, 1, 1, 1, :) ./ ρₒ
    Jᵛt = interior(τʸt, 1, 1, 1, :) ./ ρₒ
    u★ = @. (Jᵘt^2 + Jᵛt^2)^(1/4)

    lines!(axu, times, interior(ut, 1, 1, Nz, :), color=colors[1], label="Zonal")
    lines!(axu, times, interior(vt, 1, 1, Nz, :), color=colors[2], label="Meridional")
    lines!(axu, times, u★, color=colors[3], label="Ocean-side u★") 
    vlines!(axu, tn, linewidth=4, color=(:black, 0.5))
    axislegend(axu)

    lines!(axτ, times, interior(τˣt, 1, 1, 1, :), label="Zonal")
    lines!(axτ, times, interior(τʸt, 1, 1, 1, :), label="Meridional")
    vlines!(axτ, tn, linewidth=4, color=(:black, 0.5))
    axislegend(axτ)

    lines!(axT, times, Tat .- 273.15,             color=colors[1], linewidth=2, linestyle=:dash, label="Atmosphere temperature")
    lines!(axT, times, interior(Tt, 1, 1, Nz, :), color=colors[2], linewidth=4, label="Ocean surface temperature")
    vlines!(axT, tn, linewidth=4, color=(:black, 0.5))
    axislegend(axT)

    lines!(axQ, times, interior(Qt, 1, 1, 1, :),   color=colors[1], label="Total",     linewidth=6)
    lines!(axQ, times, interior(Qset, 1, 1, 1, :), color=colors[2], label="Sensible",  linewidth=2)
    lines!(axQ, times, interior(Qlat, 1, 1, 1, :), color=colors[3], label="Latent",    linewidth=2)
    lines!(axQ, times, - Qswt,                     color=colors[4], label="Shortwave", linewidth=2)
    lines!(axQ, times, - Qlwt,                     color=colors[5], label="Longwave",  linewidth=2)
    vlines!(axQ, tn, linewidth=4, color=(:black, 0.5))
    axislegend(axQ)

    #lines!(axF, times, interior(Ft, 1, 1, 1, :), label="Net freshwater flux")
    lines!(axF, times, Pt, label="Prescribed freshwater flux")
    lines!(axF, times, - interior(Et, 1, 1, 1, :), label="Evaporation")
    vlines!(axF, tn, linewidth=4, color=(:black, 0.5))
    axislegend(axF)

    lines!(axS, times, interior(St, 1, 1, Nz, :))
    vlines!(axS, tn, linewidth=4, color=(:black, 0.5))

    zc = znodes(Tt)
    zf = znodes(κt)
    un = @lift interior(ut[$n], 1, 1, :)
    vn = @lift interior(vt[$n], 1, 1, :)
    Tn = @lift interior(Tt[$n], 1, 1, :)
    Sn = @lift interior(St[$n], 1, 1, :)
    κn = @lift interior(κt[$n], 1, 1, :)
    en = @lift max.(1e-6, interior(et[$n], 1, 1, :))
    N²n = @lift interior(N²t[$n], 1, 1, :)

    scatterlines!(axuz, un,  zc, label="u") 
    scatterlines!(axuz, vn,  zc, label="v") 
    scatterlines!(axTz, Tn,  zc) 
    scatterlines!(axSz, Sn,  zc) 
    scatterlines!(axez, en,  zc) 
    scatterlines!(axNz, N²n, zf) 
    scatterlines!(axκz, κn,  zf) 

    axislegend(axuz)

    Tmax = maximum(interior(Tt))
    Tmin = minimum(interior(Tt))
    xlims!(axTz, Tmin - 0.1, Tmax + 0.1)

    Nmax = maximum(interior(N²t))
    Nmin = minimum(interior(N²t))
    xlims!(axNz, Nmin / 2, Nmin * 1.1)

    emax = maximum(interior(et))
    xlims!(axez, 8e-7, emax * 1.1)
    xlims!(axκz, 1e-7, 10)

    Smax = maximum(interior(St))
    Smin = minimum(interior(St))
    xlims!(axSz, Smin - 0.2, Smax + 0.2)

    display(fig)

    record(fig, "$(location)_single_column_simulation.mp4", 1:Nt, framerate=24) do nn
        @info "Drawing frame $nn of $Nt..."
        n[] = nn
    end
end
