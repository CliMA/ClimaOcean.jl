using ClimaOcean
using ClimaOcean: SAARAtlas, absolute_salinity, time_step
using Oceananigans.Units
using Test

@testset "Absolute salinity reproduces the GSW check values" begin
    # The example of `gsw_SA_from_SP` in the GSW-Matlab documentation
    Sᴾ = [34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600]
    p = [10, 50, 125, 250, 600, 1000]
    Sᴬ = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]

    atlas = SAARAtlas()
    @test [absolute_salinity(atlas, Sᴾ[n], p[n], 188, 4) for n in eachindex(p)] ≈ Sᴬ atol = 1e-4
end

@testset "Time-step schedules" begin
    schedule = TimeStepSchedule(4minutes, 60days => 30minutes, 365days => 45minutes)

    @test time_step(schedule, 0) == 4minutes
    @test time_step(schedule, 60days - 1) == 4minutes
    @test time_step(schedule, 60days) == 30minutes
    @test time_step(schedule, 400days) == 45minutes
    @test time_step(TimeStepSchedule(90minutes), 1000days) == 90minutes
end

@testset "OMIP configurations" begin
    for configuration in (OneDegree(), QuarterDegree(), TwelfthDegree())
        @test configuration isa OMIPConfiguration
    end

    @test isnothing(QuarterDegree().ocean.eddy_closure)
    @test OneDegree(biharmonic_timescale = 10days).ocean.biharmonic_timescale == 10days
    @test QuarterDegree(sea_ice_thickness_categories = 1).sea_ice.thickness_categories == 1
    @test TwelfthDegree().ocean.background_diffusivity === henyey_diffusivity
    @test henyey_diffusivity(0, 0, 0, 0) == 2e-6
    @test henyey_diffusivity(0, 90, 0, 0) ≈ 1e-5
end
