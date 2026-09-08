using Test
using ClimaOcean
using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: CenterField, interior
using Oceananigans.Operators: Azᶜᶜᶜ, Δzᶜᶜᶜ
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

using ClimaOcean.OMIPConfigurations: AdvectiveBottomBoundaryLayer, BottomBoundaryLayer, bottom_level,
                       advective_bottom_boundary_layer_tendency, bottom_boundary_layer_tendency,
                       update_advective_bottom_boundary_layer!, update_bottom_boundary_layer!,
                       merge_tracer_forcings

const Nx, Ny, Nz = 8, 4, 20

# A single step: the left half of the domain is a shelf at 300 m, the right half is 900 m deep.
step_bottom(x, y) = x < 40e3 ? -300.0 : -900.0
level_bottom(x, y) = -600.0

function step_grid(bottom_height)
    grid = RectilinearGrid(size = (Nx, Ny, Nz), x = (0, 80e3), y = (0, 40e3), z = (-1000, 0),
                           topology = (Bounded, Periodic, Bounded))
    return ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))
end

function state(grid, temperature)
    T = CenterField(grid)
    S = CenterField(grid)
    set!(T, temperature)
    set!(S, (x, y, z) -> 34.8)
    fill_halo_regions!(T)
    fill_halo_regions!(S)
    return (; T, S)
end

# Total tendency-weighted volume and the peak tendency, for the activation and conservation checks.
function tendency_summary(bottom_height, temperature; transport_coefficient = 10)
    grid = step_grid(bottom_height)
    fields = state(grid, temperature)

    bbl = AdvectiveBottomBoundaryLayer(grid, TEOS10EquationOfState(); transport_coefficient)
    update_advective_bottom_boundary_layer!(bbl, grid, fields)
    parameters = (bottom_boundary_layer = bbl, tracer_name = Val(:T))
    clock = Clock(time = 0.0)

    integral = 0.0
    peak = 0.0
    for i in 1:Nx, j in 1:Ny, k in 1:Nz
        G = advective_bottom_boundary_layer_tendency(i, j, k, grid, clock, fields, parameters)
        integral += G * Azᶜᶜᶜ(i, j, k, grid) * Δzᶜᶜᶜ(i, j, k, grid)
        peak = max(peak, abs(G))
    end

    return (; integral, peak, bbl, grid, fields)
end

dense_upslope(x, y, z)   = x < 40e3 ? -1.0 : 4.0
dense_downslope(x, y, z) = x > 40e3 ? -1.0 : 4.0
uniform(x, y, z)         = 4.0

@testset "Advective bottom boundary layer" begin

    @testset "activation" begin
        # Dense water upslope of a deeper neighbour is the one configuration that opens the circuit.
        @test tendency_summary(step_bottom, dense_upslope).peak > 0

        # Reversing the density contrast must switch it off: `max(0, δρ)` allows downslope flow only.
        @test tendency_summary(step_bottom, dense_downslope).peak == 0

        # A level bottom has no downslope direction, so there is no step to descend.
        @test tendency_summary(level_bottom, dense_upslope).peak == 0

        # No buoyancy contrast, no transport.
        @test tendency_summary(step_bottom, uniform).peak == 0
    end

    @testset "conservation" begin
        # The three limbs of the circuit must sum to zero over the cells they touch, or the
        # overturning would be a tracer source.
        @test tendency_summary(step_bottom, dense_upslope).integral ≈ 0 atol = 1e-8
    end

    @testset "transport scales with the coefficient" begin
        full   = tendency_summary(step_bottom, dense_upslope; transport_coefficient = 10)
        halved = tendency_summary(step_bottom, dense_upslope; transport_coefficient = 5)
        @test halved.peak ≈ full.peak / 2
    end

    # What actually separates this scheme from the diffusive one is the *topology* of the circuit,
    # not its strength. Both drive the deep bottom cell towards the shelf value; they differ in which
    # cells the displaced water passes through, and in where the water returning to the shelf comes
    # from. Those are the properties worth pinning down, and they are the ones a refactor could break.
    @testset "overturning topology" begin
        grid = step_grid(step_bottom)
        clock = Clock(time = 0.0)

        # A deep column stratified continuously, so that "the deep bottom cell" and "the deep column at
        # the shelf's own level" are distinguishable, and every level of the upwelling limb carries a
        # nonzero difference. A piecewise-constant profile would leave the limb present but valued zero,
        # since advecting a uniform tracer does nothing.
        stratified(x, y, z) = x < 40e3 ? -1.0 : 10.0 + 0.01z
        fields = state(grid, stratified)

        advective = AdvectiveBottomBoundaryLayer(grid, TEOS10EquationOfState(); transport_coefficient = 10)
        diffusive = BottomBoundaryLayer(grid, TEOS10EquationOfState(); diffusivity = 5000)
        update_advective_bottom_boundary_layer!(advective, grid, fields)
        update_bottom_boundary_layer!(diffusive, grid, fields)

        Gadv(i, j, k) = advective_bottom_boundary_layer_tendency(i, j, k, grid, clock, fields,
                            (bottom_boundary_layer = advective, tracer_name = Val(:T)))
        Gdif(i, j, k) = bottom_boundary_layer_tendency(i, j, k, grid, clock, fields,
                            (bottom_boundary_layer = diffusive, tracer_name = Val(:T)))

        shelf, deep = 4, 5                     # the two columns either side of the step
        kˢ = bottom_level(advective, shelf, 2)
        kᵈ = bottom_level(advective, deep, 2)
        @test kˢ > kᵈ                          # the shelf is the shallower column

        # The diffusive scheme touches only the two bottom cells of the pair; the overturning
        # circuit also carries the displaced water up through every level between them.
        touched(G) = count(!=(0), [G(i, 2, k) for i in 1:Nx, k in 1:Nz])
        @test touched(Gadv) > touched(Gdif)
        @test all(Gadv(deep, 2, k) != 0 for k in (kᵈ+1):kˢ)
        @test all(Gdif(deep, 2, k) == 0 for k in (kᵈ+1):kˢ)

        # Both drive the deep bottom cell towards the shelf value: that part is common.
        @test Gadv(deep, 2, kᵈ) < 0            # cooled by the dense shelf water arriving
        @test Gdif(deep, 2, kᵈ) < 0

        # The difference is what returns to the shelf. The diffusive exchange sends the deep *bottom*
        # water back, so the shelf is contaminated by the abyssal water it just made; the overturning
        # circuit returns ambient from the shelf's own level instead.
        Vˢ = Azᶜᶜᶜ(shelf, 2, kˢ, grid) * Δzᶜᶜᶜ(shelf, 2, kˢ, grid)
        Qadv = advective.transport_x[shelf, 2, 1]
        Qdif = diffusive.transport_x[shelf, 2, 1]
        @test Qadv > 0 && Qdif > 0

        Tˢ  = fields.T[shelf, 2, kˢ]
        Tᵈˢ = fields.T[deep,  2, kˢ]           # ambient at the shelf's level
        Tᵈ  = fields.T[deep,  2, kᵈ]           # abyssal water
        @test Tᵈˢ != Tᵈ                        # the stratification makes the two targets differ

        @test Gadv(shelf, 2, kˢ) ≈ Qadv * (Tᵈˢ - Tˢ) / Vˢ
        @test Gdif(shelf, 2, kˢ) ≈ Qdif * (Tᵈ  - Tˢ) / Vˢ
    end

    @testset "merge_tracer_forcings" begin
        @test merge_tracer_forcings(NamedTuple(), NamedTuple()) == NamedTuple()
        @test merge_tracer_forcings((T = 1,), NamedTuple()) == (T = 1,)

        both = merge_tracer_forcings((T = 1, S = 2), (T = 3, S = 4))
        @test both.T == (1, 3)
        @test both.S == (2, 4)
    end
end
