using Test
using ClimaOcean
using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: CenterField, interior
using Oceananigans.Operators: Azᶜᶜᶜ, Δzᶜᶜᶜ
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

using ClimaOcean.OMIPConfigurations: BottomBoundaryLayer, bottom_boundary_layer_tendency, bottom_level,
                       update_bottom_boundary_layer!

const Nx, Ny, Nz = 16, 4, 20

deepening(x, y) = -200 - 800 * (x / 160e3)   # shallow at x = 0, deep at x = 160 km
level(x, y)     = -600.0

function slope_grid(bottom_height)
    grid = RectilinearGrid(size = (Nx, Ny, Nz), x = (0, 160e3), y = (0, 40e3), z = (-1000, 0),
                           topology = (Bounded, Periodic, Bounded))
    return ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))
end

# Total tendency-weighted volume, the peak tendency, and the number of cells the scheme touches.
function tendency_summary(bottom_height, temperature; diffusivity = 5000)
    grid = slope_grid(bottom_height)

    T = CenterField(grid)
    S = CenterField(grid)
    set!(T, temperature)
    set!(S, (x, y, z) -> 34.8)
    fill_halo_regions!(T)
    fill_halo_regions!(S)
    fields = (; T, S)

    bbl = BottomBoundaryLayer(grid, TEOS10EquationOfState(); diffusivity)
    update_bottom_boundary_layer!(bbl, grid, fields)
    parameters = (bottom_boundary_layer = bbl, tracer_name = Val(:T))
    clock = Clock(time = 0.0)

    integral = 0.0
    peak = 0.0
    active = 0
    tendency = zeros(Nx, Ny, Nz)

    for i in 1:Nx, j in 1:Ny, k in 1:Nz
        G = bottom_boundary_layer_tendency(i, j, k, grid, clock, fields, parameters)
        tendency[i, j, k] = G
        integral += G * Azᶜᶜᶜ(i, j, k, grid) * Δzᶜᶜᶜ(i, j, k, grid)
        peak = max(peak, abs(G))
        G != 0 && (active += 1)
    end

    return (; integral, peak, active, tendency, bbl)
end

dense_upslope(x, y, z)   = x < 55e3 ? -1.0 : 4.0
dense_downslope(x, y, z) = x > 55e3 ? -1.0 : 4.0
uniform(x, y, z)         = 4.0

@testset "Bottom boundary layer" begin

    @testset "activation" begin
        # Dense water upslope of a deeper neighbour is the one configuration that drives a flux.
        upslope = tendency_summary(deepening, dense_upslope)
        @test upslope.peak > 0
        @test upslope.active > 0

        # Reversing the density gradient must switch it off: this guards the sign of the slope
        # criterion, the one error that would still produce a plausible-looking run.
        @test tendency_summary(deepening, dense_downslope).peak == 0

        # A level bottom has no slope sign, so the criterion can never be met.
        @test tendency_summary(level, dense_upslope).peak == 0

        # No buoyancy contrast, no flux.
        @test tendency_summary(deepening, uniform).peak == 0
    end

    @testset "conservation" begin
        upslope = tendency_summary(deepening, dense_upslope)
        @test upslope.integral == 0
    end

    @testset "direction" begin
        upslope = tendency_summary(deepening, dense_upslope)
        bbl = upslope.bbl
        j = 2
        # The front sits between i = 5 and i = 6; the dense anomaly must propagate downslope, so the
        # cold cell warms and its deeper neighbour cools by an equal and opposite amount.
        cold = upslope.tendency[5, j, bottom_level(bbl, 5, j)]
        warm = upslope.tendency[6, j, bottom_level(bbl, 6, j)]
        @test cold > 0
        @test warm < 0
        @test cold ≈ -warm
    end

    @testset "transport magnitude" begin
        # Transport is linear in the coefficient, so halving κ halves the tendency everywhere.
        upslope = tendency_summary(deepening, dense_upslope; diffusivity = 5000)
        halved  = tendency_summary(deepening, dense_upslope; diffusivity = 2500)
        @test halved.peak ≈ upslope.peak / 2
    end
end
