using JLD2
using Oceananigans
using OrthogonalSphericalShellGrids
using SeawaterPolynomials
using SeawaterPolynomials: TEOS10EquationOfState

using Oceananigans.Fields: ConstantField
using Oceananigans.Models: seawater_density

function calculate_climatology(filename; start_day = 0)
    geopotential_height = ConstantField(Float32(-200.0))
    equation_of_state   = TEOS10EquationOfState(Float32)

    T_series = FieldTimeSeries(filename, "T"; backend = OnDisk())
    S_series = FieldTimeSeries(filename, "S"; backend = OnDisk())
    u_series = FieldTimeSeries(filename, "u"; backend = OnDisk())
    v_series = FieldTimeSeries(filename, "v"; backend = OnDisk())
    w_series = FieldTimeSeries(filename, "w"; backend = OnDisk())
    e_series = FieldTimeSeries(filename, "e"; backend = OnDisk())

    grid = T_series.grid

    Na = size(T_series[1])
    Nu = size(u_series[1])
    Nv = size(v_series[1])
    Nw = size(w_series[1])

    times = T.times 
    tidx  = findfirst(t -> t > start_day * 86400, times)

    T  = zeros(Float32, Na)
    T² = zeros(Float32, Na)
    S  = zeros(Float32, Na)
    S² = zeros(Float32, Na)
    u  = zeros(Float32, Nu)
    u² = zeros(Float32, Nu)
    v  = zeros(Float32, Nv)
    v² = zeros(Float32, Nv)
    w  = zeros(Float32, Nw)
    w² = zeros(Float32, Nw)
    e  = zeros(Float32, Na)
    ρ  = zeros(Float32, Na)
    ρ² = zeros(Float32, Na)

    mdata = zeros(12)

    for month in 1:12
        for t in tidx:length(times)
            mnth = findfirst(x -> x >= day, total_days)
            @show day, mnth, month
            if mnth == month
                @info "adding day $day to month $mnth"
                Ttmp = interior(T_series[t])
                Stmp = interior(S_series[t])
                utmp = interior(u_series[t])
                vtmp = interior(v_series[t])
                wtmp = interior(w_series[t])
                T  .+= Ttmp
                S  .+= Stmp
                T² .+= Ttmp .* Ttmp
                S² .+= Stmp .* Stmp
                ρtmp = seawater_density(grid, equation_of_state, Ttmp, Stmp, geopotential_height)
                ρ  .+= ρtmp
                ρ² .+= ρtmp .* ρtmp
                u  .+= utmp
                v  .+= vtmp
                w  .+= wtmp
                u² .+= utmp .* utmp
                v² .+= vtmp .* vtmp
                w² .+= wtmp .* wtmp
                e  .+= e_series[t]
      
                mdata[mnth] += 1
            end
            GC.gc()
        end

        T  ./= mdata[month]
        S  ./= mdata[month]
        T² ./= mdata[month]
        S² ./= mdata[month]
        ρ  ./= mdata[month]
        ρ² ./= mdata[month]
        u  ./= mdata[month]
        v  ./= mdata[month]
        w  ./= mdata[month]
        u² ./= mdata[month]
        v² ./= mdata[month]
        w² ./= mdata[month]
        e  ./= mdata[month]

        jldopen("climatology_month_$(month).jl"; T, S, T², S², ρ, ρ², u, v, w, u², v², w², e)

        for v in [T, S, T², S², ρ, ρ², u, v, w, u², v², w², e]
            fill!(v, 0)
        end
    end
end