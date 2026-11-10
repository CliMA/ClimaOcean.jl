using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Dates: DateTime
using TOML

regression_architecture = parse(Bool, get(ENV, "GPU_TEST", "false")) ? GPU() : CPU()
regression_reference_path = joinpath(@__DIR__, "one_degree_regression_reference.toml")

function one_degree_regression_simulation(arch)
    simulation = omip_simulation(OneDegree(); arch, diagnostics = false,
                                 dates = (DateTime(1958, 1, 1), DateTime(1958, 2, 1)))
    simulation.stop_iteration = 100
    return simulation
end

global_mean(c) = only(Array(interior(compute!(Field(Average(c))))))
global_variance(c) = global_mean((c - global_mean(c))^2)
global_integral(c) = only(Array(interior(compute!(Field(Integral(c))))))

function regression_summary(simulation)
    ocean = simulation.model.ocean.model
    sea_ice = simulation.model.sea_ice.model
    T, S = ocean.tracers.T, ocean.tracers.S
    u, v, w = ocean.velocities
    Nz = size(ocean.grid, 3)

    # The free surface lives on a grid with extended halos, so it is averaged on the ocean grid
    η = Field{Center, Center, Nothing}(ocean.grid)
    interior(η) .= interior(ocean.free_surface.displacement)

    return Dict("mean_temperature" => global_mean(T),
                "mean_salinity" => global_mean(S),
                "mean_surface_temperature" => global_mean(view(T, :, :, Nz)),
                "mean_sea_level" => global_mean(η),
                "temperature_variance" => global_variance(T),
                "salinity_variance" => global_variance(S),
                "sea_level_variance" => global_variance(η),
                "zonal_velocity_variance" => global_variance(u),
                "meridional_velocity_variance" => global_variance(v),
                "vertical_velocity_variance" => global_variance(w),
                "maximum_zonal_velocity" => maximum(abs, u),
                "maximum_meridional_velocity" => maximum(abs, v),
                "maximum_vertical_velocity" => maximum(abs, w),
                "maximum_turbulent_kinetic_energy" => maximum(ocean.tracers.e),
                "sea_ice_volume" => global_integral(sea_ice.ice_thickness * sea_ice.ice_concentration))
end
