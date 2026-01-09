using Test

using Oceananigans
using ClimaOcean

@testset "PrescribedAtmosphere set!" begin
    grid = RectilinearGrid(CPU();
                           size = (4, 4, 1),
                           x = (0, 1),
                           y = (0, 1),
                           z = (0, 1),
                           topology = (Periodic, Periodic, Bounded))

    times = [0.0, 3600.0, 7200.0]
    atmosphere = PrescribedAtmosphere(grid, times)

    set!(atmosphere;
         T = t -> 300 + t / 3600,
         u = t -> 5.0,
         v = t -> -2.0,
         p = t -> 101000.0)

    Tparent = parent(atmosphere.tracers.T)
    uparent = parent(atmosphere.velocities.u)
    vparent = parent(atmosphere.velocities.v)
    pparent = parent(atmosphere.pressure)

    @test size(Tparent)[end] == length(times)
    @test size(uparent)[end] == length(times)
    @test size(vparent)[end] == length(times)
    @test size(pparent)[end] == length(times)

    @test Tparent[1, 1, 1, 1] == 300.0
    @test Tparent[1, 1, 1, 2] == 301.0
    @test Tparent[1, 1, 1, 3] == 302.0

    @test uparent[1, 1, 1, 1] == 5.0
    @test vparent[1, 1, 1, 1] == -2.0
    @test pparent[1, 1, 1, 1] == 101000.0
end


