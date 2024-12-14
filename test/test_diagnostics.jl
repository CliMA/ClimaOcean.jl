include("runtests_setup.jl")

using SeawaterPolynomials: TEOS10EquationOfState
using Oceananigans.BuoyancyFormulations: buoyancy
using Oceananigans: location
using ClimaOcean.DataWrangling.ECCO: ECCO_field, ECCOFieldTimeSeries
using ClimaOcean.Diagnostics: MixedLayerDepthField, MixedLayerDepthOperand

@testset "MixedLayerDepthField" begin
    for arch in test_architectures
        grid = LatitudeLongitudeGrid(arch;
                                     size = (3, 3, 100),
                                     latitude  = (0, 30),
                                     longitude = (150, 180),
                                     z = (-1000, 0))

        bottom_height = regrid_bathymetry(grid; 
                                          minimum_depth = 10,
                                          interpolation_passes = 5,
                                          major_basins = 1)

        grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

        start = DateTimeProlepticGregorian(1993, 1, 1)
        stop  = DateTimeProlepticGregorian(1993, 2, 1)
        dates = range(start; stop, step=Month(1))

        Tmeta = ECCOMetadata(:temperature; dates)
        Smeta = ECCOMetadata(:salinity; dates)

        Tt = ECCOFieldTimeSeries(Tmeta, grid; time_indices_in_memory=2)
        St = ECCOFieldTimeSeries(Smeta, grid; time_indices_in_memory=2)

        equation_of_state = TEOS10EquationOfState()
        sb = SeawaterBuoyancy(; equation_of_state)
        tracers = (T=Tt[1], S=St[1])
        h = MixedLayerDepthField(sb, grid, tracers)

        @test h isa Field
        @test location(h) == (Center, Center, Nothing)
        @test h.operand isa MixedLayerDepthOperand
        @test h.operand.buoyancy_perturbation isa KernelFunctionOperation

        compute!(h)
        @test @allowscalar h[1, 1, 1] ≈ 16.255836 # m

        tracers = (T=Tt[2], S=St[2])
        h.operand.buoyancy_perturbation = buoyancy(sb, grid, tracers)
        compute!(h)
        @test @allowscalar h[1, 1, 1] ≈ 9.295890287 # m
    end
end
