using ClimaOcean
using ClimaOcean: compute_total_water!
using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: MutableVerticalDiscretization
using Oceananigans.TurbulenceClosures: IsopycnalSkewSymmetricDiffusivity
using JLD2
using Test

@testset "Total water normalization and OMIP diagnostics" begin
    grid = LatitudeLongitudeGrid(CPU(); size = (12, 10, 4), halo = (7, 7, 7), longitude = (0, 360), latitude = (-60, 60),
                                 z = MutableVerticalDiscretization(collect(range(-400, 0, length = 5))))

    closure = (NumericalEarth.Oceans.default_ocean_closure(),
               IsopycnalSkewSymmetricDiffusivity(κ_skew = 500, κ_symmetric = 800))

    ocean = ocean_simulation(grid; closure, free_surface = SplitExplicitFreeSurface(grid; substeps = 20), warn = false)
    set!(ocean.model, T = 10, S = 35)
    sea_ice = sea_ice_simulation(grid, ocean)

    atmosphere_grid = LatitudeLongitudeGrid(CPU(); size = (36, 18), longitude = (0, 360), latitude = (-90, 90),
                                            topology = (Periodic, Bounded, Flat))
    atmosphere = PrescribedAtmosphere(atmosphere_grid, [0, 1days])
    parent(atmosphere.precipitation_flux.rain) .= 1e-4
    radiation = PrescribedRadiation(atmosphere_grid, [0, 1days])

    coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
    simulation = Simulation(coupled_model; Δt = 10minutes, stop_iteration = 12)

    normalization = TotalWaterNormalization(coupled_model)
    add_callback!(simulation, normalization, IterationInterval(1))

    output_directory = mktempdir()
    add_omip_diagnostics!(simulation; output_directory, prefix = "test",
                          surface_averaging_interval = 1hour,
                          field_averaging_interval = 1hour,
                          averages_interval = 1hour)

    run!(simulation)

    # Without the normalization the rain would add about 10¹¹ m³ of water over the run
    @test only(Array(compute_total_water!(normalization, coupled_model))) ≈
          only(Array(normalization.reference_water)) atol = 1

    for (file, names) in (("test_surface.jld2", ("tos", "sos", "zos", "mlotst", "tauuo", "hfds", "siconc", "sisnthick")),
                          ("test_fields.jld2", ("to", "so", "uvol", "vvol", "vvolgm", "tke")),
                          ("test_averages.jld2", ("zosga", "to_h", "voco", "soco", "sivol", "snvol")))
        timeseries = jldopen(joinpath(output_directory, file)) do io
            keys(io["timeseries"])
        end
        @test all(name -> name ∈ timeseries, names)
    end
end
