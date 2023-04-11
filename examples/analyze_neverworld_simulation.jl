using Oceananigans
using Oceananigans.Units
using Printf
using JLD2

dir = "archive"
backend = OnDisk()
bt = FieldTimeSeries("$dir/neverworld_xyz.jld2", "b"; backend)
t = bt.times
Nt = length(t)
grid = bt.grid

Nyears = 20
t_end = 200 * 360 * day # 200 "years"
t_start = (200 - Nyears) * 360 * day # 200 "years"

t = bt.times
n_stop  = searchsortedfirst(t, t_stop)
n_start = searchsortedfirst(t, t_start) - 1

for season = 1:4
    n_average_start = n_start + season - 1
    n_average = n_average_start:4:n_stop

    for name in ("b", "u", "v", "w", "e", "κᶜ")
        ψt = FieldTimeSeries("$dir/neverworld_xyz.jld2", name; backend)
        LX, LY, LZ = location(ψt)
        ψ_avg = Field{LX, LY, LZ}(grid)

        for n in n_average
            @info name season n
            ψn = ψt[n]
            ψ_avg .+= ψn / length(n_average)
        end
        
        filename = string("climatology_", name, "_season_", season, ".jld2")
        file = jldopen(filename, "a+")
        file[name] = ψ_avg
        file["season"] = season
        close(file)
    end
end

