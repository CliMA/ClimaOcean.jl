using Oceananigans
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using Oceananigans.Units
using Printf
using JLD2

backend = OnDisk()
bt = FieldTimeSeries("neverworld_xyz.jld2", "b"; backend)
ut = FieldTimeSeries("neverworld_xyz.jld2", "u"; backend)
ζt = FieldTimeSeries("neverworld_xyz.jld2", "ζ"; backend)
et = FieldTimeSeries("neverworld_xyz.jld2", "e"; backend)
κt = FieldTimeSeries("neverworld_xyz.jld2", "κᶜ"; backend)

t = bt.times
Nt = length(t)
grid = bt.grid

averaging_years = 5

for month in 1:12
    n_start = Nt - 12 * averaging_years - month + 1
    n_final = Nt - month + 1
    n_average = n_start:12:n_final

    for name in ("b", "u", "ζ", "e", "κᶜ")
        ψt = FieldTimeSeries("neverworld_xyz.jld2", name; backend)
        LX, LY, LZ = location(ψt)
        ψ_avg = Field{LX, LY, LZ}(grid)

        for n in n_average
            @info name month n
            ψn = ψt[n]
            ψ_avg .+= ψn / length(n_average)
        end
        
        filename = string("climatology_", name, "_month_", month, ".jld2")
        file = jldopen(filename, "a+")
        file[name] = ψ_avg
        file["month"] = month
        close(file)
    end
end

