using Test
using ClimaOcean
using Oceananigans
using Oceananigans.Fields: CenterField, interior
using ClimaOcean.OMIPConfigurations: overflow_restoring_mask, overflow_restoring_tendency, overflow_restoring_forcing

@testset "Overflow restoring" begin

    grid = ImmersedBoundaryGrid(
        LatitudeLongitudeGrid(size = (16, 16, 12), longitude = (-40, -20), latitude = (58, 70),
                              z = (-3000, 0), halo = (4, 4, 4)),
        GridFittedBottom((λ, φ) -> λ < -34 ? -400.0 : -3000.0))     # a shelf to the west, deep to the east

    longitude, latitude, minimum_depth = (-32.0, -24.0), (62.0, 66.5), 1500

    mask = overflow_restoring_mask(grid; longitude, latitude, minimum_depth)
    m = Array(interior(mask))

    @test all(x -> x == 0 || x == 1, m)
    @test any(x -> x == 1, m)                       # the region is not empty on this grid

    # Every selected cell must satisfy all three criteria, and every cell satisfying them be selected.
    λs = Array(λnodes(grid, Center())); φs = Array(φnodes(grid, Center())); zs = Array(znodes(grid, Center()))
    for i in axes(m, 1), j in axes(m, 2), k in axes(m, 3)
        wet = !Oceananigans.ImmersedBoundaries.inactive_node(i, j, k, grid, Center(), Center(), Center())
        want = wet && longitude[1] <= λs[i] <= longitude[2] &&
                      latitude[1]  <= φs[j] <= latitude[2]  && -zs[k] > minimum_depth
        @test (m[i, j, k] == 1) == want
    end

    # The tendency relaxes toward the target inside the mask and vanishes outside it.
    T = CenterField(grid); set!(T, 5.0)
    S = CenterField(grid); set!(S, 35.0)
    fields = (; T, S)
    clock = Clock(time = 0.0)
    target, timescale = 2.0, 10 * 86400.0
    parameters = (; mask, rate = 1 / timescale, target, tracer_name = Val(:T))

    G = [overflow_restoring_tendency(i, j, k, grid, clock, fields, parameters)
         for i in axes(m,1), j in axes(m,2), k in axes(m,3)]

    @test all(G[m .== 0] .== 0)                                  # inert outside
    @test all(G[m .== 1] .≈ (target - 5.0) / timescale)          # relaxes toward the target
    @test all(G[m .== 1] .< 0)                                   # 5 °C is warmer than the 2 °C target

    # Halving the timescale doubles the rate.
    fast = (; mask, rate = 2 / timescale, target, tracer_name = Val(:T))
    Gf = [overflow_restoring_tendency(i, j, k, grid, clock, fields, fast)
          for i in axes(m,1), j in axes(m,2), k in axes(m,3)]
    @test all(Gf[m .== 1] .≈ 2 .* G[m .== 1])

    # The convenience constructor is off by default and returns both tracers when on.
    @test overflow_restoring_forcing(grid, nothing) == NamedTuple()
    f = overflow_restoring_forcing(grid, 10 * 86400.0; longitude, latitude, minimum_depth)
    @test haskey(f, :T) && haskey(f, :S)
end
