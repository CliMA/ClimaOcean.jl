include("runtests_setup.jl")

using SeawaterPolynomials: TEOS10EquationOfState
using Oceananigans: location
using Oceananigans.Models: buoyancy_operation
using ClimaOcean.Diagnostics: MixedLayerDepthField, MixedLayerDepthOperand
using ClimaOcean.Diagnostics: MeridionalStreamfunction, compute_streamfunction
using ConservativeRegridding

for arch in test_architectures, dataset in (ECCO4Monthly(),)
    A = typeof(arch)
    @info "Testing MixedLayerDepthField with $(typeof(dataset)) on $A"

    @testset "MixedLayerDepthField" begin
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

        Tmeta = Metadata(:temperature; dataset, dates)
        Smeta = Metadata(:salinity; dataset, dates)

        Tt = FieldTimeSeries(Tmeta, grid; time_indices_in_memory=2)
        St = FieldTimeSeries(Smeta, grid; time_indices_in_memory=2)

        equation_of_state = TEOS10EquationOfState()
        sb = SeawaterBuoyancy(; equation_of_state)
        tracers = (T=Tt[1], S=St[1])
        h = MixedLayerDepthField(sb, grid, tracers)

        @test h isa Field
        @test location(h) == (Center, Center, Nothing)
        @test h.operand isa MixedLayerDepthOperand
        @test h.operand.buoyancy_perturbation isa KernelFunctionOperation

        compute!(h)
        if dataset isa ECCO4Monthly
            @test @allowscalar h[1, 1, 1] ≈ 16.2558363 # m
        end

        tracers = (T=Tt[2], S=St[2])
        h.operand.buoyancy_perturbation = buoyancy_operation(sb, grid, tracers)
        compute!(h)
        if dataset isa ECCO4Monthly
            @test @allowscalar h[1, 1, 1] ≈ 9.2957298 # m
        end
    end
end

for arch in test_architectures
    A = typeof(arch)
    @info "Testing MeridionalStreamfunction on $A"

    @testset "MeridionalStreamfunction" begin
        # Create a simple LatitudeLongitudeGrid
        grid = LatitudeLongitudeGrid(arch;
                                     size = (36, 18, 3),
                                     longitude = (0, 360),
                                     latitude = (-60, 60),
                                     z = (-500, 0))

        # Create a w field with a simple sinusoidal pattern
        w = CenterField(grid)
        set!(w, (λ, φ, z) -> sin(deg2rad(φ)))

        # Test constructor
        moc = MeridionalStreamfunction(w; resolution=5, latitude=(-50, 50))

        @test moc.latitude_longitude_grid isa LatitudeLongitudeGrid
        @test moc.regridded_w isa Field

        # Test compute_streamfunction
        result = compute_streamfunction(moc, w)

        @test haskey(result, :ψ)
        @test haskey(result, :latitude)
        @test haskey(result, :depth)
        @test size(result.ψ, 1) == length(result.latitude)
        @test size(result.ψ, 2) == length(result.depth)

        # Check that latitude and depth arrays have reasonable values
        @test minimum(result.latitude) >= -50
        @test maximum(result.latitude) <= 50
        @test minimum(result.depth) >= -500
        @test maximum(result.depth) <= 0
    end
end
