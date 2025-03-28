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
    ocean = ocean_simulation(grid; free_surface)
    backend = JRA55NetCDFBackend(4)
    atmosphere = JRA55PrescribedAtmosphere(arch; backend)
    radiation = Radiation(arch)
    coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)

    # Test that Reactant does _not_ initialize in the constructor for OceanSeaIceModel
    exchange_state = coupled_model.interfaces.exchanger.exchange_atmosphere_state
    atmos_exchanger = coupled_model.interfaces.exchanger.atmosphere_exchanger
    @test all(atmos_exchanger.i .== 0)
    @test all(atmos_exchanger.j .== 0)

    # This tests that update_state! is not called
    ue = exchange_state.u
    @test all(ue .== 0) # not initialized with Reactant
end

