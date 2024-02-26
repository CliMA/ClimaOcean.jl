using Oceananigans
using Oceananigans.OutputReaders: Cyclical, update_field_time_series!
using Oceananigans.Utils: Time
using GLMakie

grid = RectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1))
times = 0:6
path = "test_field_time_series.jld2"
rm(path, force=true)
name = "c"
c = CenterField(grid)
odct = FieldTimeSeries{Center, Center, Center}(grid; backend=OnDisk(), path, name)
for n in times
    set!(c, n)
    set!(odct, c, n, n)
end

timct = FieldTimeSeries(path, name; backend=InMemory(), time_indexing=Cyclical())
pimct = FieldTimeSeries(path, name; backend=InMemory(3), time_indexing=Cyclical())

ts = -2.1:0.1:17.1
Nt = length(ts)
pci = zeros(Nt)

for (n, t) in enumerate(ts)
    update_field_time_series!(pimct, Time(t))
    pci[n] = pimct[1, 1, 1, Time(t)]
end

tci = [timct[1, 1, 1, Time(t)] for t in ts]

fig = Figure()
ax = Axis(fig[1, 1])
scatterlines!(ax, ts, tci, marker='s', color=:blue)
scatterlines!(ax, ts, pci, marker=:square, color=:pink)
scatter!(ax, timct.times, timct[1, 1, 1, :], marker='o', markersize=20)
display(fig)
