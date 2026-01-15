using Test
using Reactant
using Oceananigans.Models: initialization_update_state!
using Oceananigans.Architectures: ReactantState
using ClimaOcean

gpu_test = get(ENV, "GPU_TEST", "false") == "true"

if gpu_test
    Reactant.set_default_backend("gpu")
else
    Reactant.set_default_backend("cpu")
end

@testset "Reactant extension tests" begin
    arch = ReactantState()
    grid = LatitudeLongitudeGrid(arch;
                                 size = (256, 128, 10),
                                 longitude = (0, 360),
                                 latitude = (-80, 80),
                                 halo = (7, 7, 7),
                                 z = (-6000, 0))

    free_surface = SplitExplicitFreeSurface(substeps=10)
    ocean = ocean_simulation(grid; Δt=300, free_surface)

    # We use an idealized atmosphere to avoid downloading the whole JRA55 data
    atmos_grid  = LatitudeLongitudeGrid(arch, Float32; size=(320, 200), 
                                                       latitude=(-90, 90), 
                                                       longitude=(0, 360), 
                                                       topology=(Periodic, Bounded, Flat))

    atmos_times = range(0, 360Oceananigans.Units.days, length=10)
    atmosphere  = PrescribedAtmosphere(atmos_grid, atmos_times)

    radiation = Radiation(arch)
    coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)

    # Test that Reactant does _not_ initialize in the constructor for OceanSeaIceModel
    exchanger = coupled_model.interfaces.exchanger.atmosphere
    state     = exchanger.state
    regridder = exchanger.regridder
    @test all(regridder.i .== 0)
    @test all(regridder.j .== 0)

    # This tests that update_state! is not called
    ue = state.u
    @test all(ue .== 0) # not initialized with Reactant
end
