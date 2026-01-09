include("runtests_setup.jl")

using ClimaOcean.OceanSeaIceModels: IceBathHeatFlux,
                                     TwoEquationHeatFlux,
                                     ThreeEquationHeatFlux,
                                     MomentumBasedFrictionVelocity

using ClimaOcean.OceanSeaIceModels.InterfaceComputations: compute_interface_heat_flux,
                                                          get_friction_velocity,
                                                          solve_interface_conditions,
                                                          SeaIceOceanInterface,
                                                          ComponentInterfaces

using ClimaSeaIce.SeaIceThermodynamics: LinearLiquidus, melting_temperature

@testset "Sea ice-ocean heat flux formulations" begin

    @testset "IceBathHeatFlux construction" begin
        flux = IceBathHeatFlux()
        @test flux.characteristic_melting_speed == 1e-5

        flux2 = IceBathHeatFlux(characteristic_melting_speed = 1e-4)
        @test flux2.characteristic_melting_speed == 1e-4

        flux3 = IceBathHeatFlux(Float32)
        @test flux3.characteristic_melting_speed isa Float32
    end

    @testset "TwoEquationHeatFlux construction" begin
        flux = TwoEquationHeatFlux()
        @test flux.heat_transfer_coefficient == 0.006

        flux2 = TwoEquationHeatFlux(heat_transfer_coefficient = 0.01,
                                    friction_velocity = MomentumBasedFrictionVelocity())

        @test flux2.heat_transfer_coefficient == 0.01
        @test flux2.friction_velocity isa MomentumBasedFrictionVelocity

        flux3 = TwoEquationHeatFlux(Float32)
        @test flux3.heat_transfer_coefficient isa Float32
    end

    @testset "ThreeEquationHeatFlux construction" begin
        flux = ThreeEquationHeatFlux()
        @test flux.heat_transfer_coefficient == 0.0095  # Default from Hieronymus et al. (2021)
        @test flux.salt_transfer_coefficient ≈ 0.0095 / 35  # R = 35

        flux2 = ThreeEquationHeatFlux(heat_transfer_coefficient = 0.01,
                                      salt_transfer_coefficient = 0.001,
                                      friction_velocity = MomentumBasedFrictionVelocity())
        @test flux2.heat_transfer_coefficient == 0.01
        @test flux2.salt_transfer_coefficient == 0.001
        @test flux2.friction_velocity isa MomentumBasedFrictionVelocity

        flux3 = ThreeEquationHeatFlux(Float32)
        @test flux3.heat_transfer_coefficient isa Float32
        @test flux3.salt_transfer_coefficient isa Float32
    end

    @testset "solve_interface_conditions" begin
        # Test with a linear liquidus
        liquidus = LinearLiquidus(Float64)

        Tₒ = 1.0    # Ocean temperature (above freezing)
        Sₒ = 35.0   # Ocean salinity
        Sᵢ = 5.0    # Ice salinity
        Gₕ = 1e-7   # Small growth rate
        αₕ = 0.006  # Heat transfer coefficient
        αₛ = 0.0001 # Salt transfer coefficient
        u★ = 0.01   # Friction velocity

        Tᵦ, Sᵦ = solve_interface_conditions(Tₒ, Sₒ, Sᵢ, Gₕ, αₕ, αₛ, u★, liquidus)

        # Interface salinity should be between ice and ocean salinity
        @test Sᵦ >= Sᵢ
        @test Sᵦ <= Sₒ

        # Interface temperature should be at freezing point of interface salinity
        Tₘ = melting_temperature(liquidus, Sᵦ)
        @test Tᵦ ≈ Tₘ
    end
end

@testset "Coupled model with different heat flux formulations" begin
    for arch in test_architectures
        A = typeof(arch)
        @info "Testing heat flux formulations on $A"

        grid = LatitudeLongitudeGrid(arch,
                                     size = (10, 10, 1),
                                     latitude = (-80, 80),
                                     longitude = (0, 360),
                                     z = (-100, 0))

        ocean = ocean_simulation(grid, momentum_advection=nothing, closure=nothing, tracer_advection=nothing)
        sea_ice = sea_ice_simulation(grid, ocean)

        backend = JRA55NetCDFBackend(4)
        atmosphere = JRA55PrescribedAtmosphere(arch; backend)
        radiation = Radiation(arch)

        # Test with ThreeEquationHeatFlux (default)
        @test begin
            coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
            flux_form = coupled_model.interfaces.sea_ice_ocean_interface.flux_formulation
            flux_form isa ThreeEquationHeatFlux
        end

        # Test with IceBathHeatFlux via ComponentInterfaces
        @test begin
            flux = IceBathHeatFlux()
            interfaces = ComponentInterfaces(atmosphere, ocean, sea_ice;
                                              radiation,
                                              sea_ice_ocean_heat_flux = flux)
            coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation, interfaces)
            flux_form = coupled_model.interfaces.sea_ice_ocean_interface.flux_formulation
            flux_form isa IceBathHeatFlux
        end

        # Test with TwoEquationHeatFlux via ComponentInterfaces
        @test begin
            flux = TwoEquationHeatFlux()
            interfaces = ComponentInterfaces(atmosphere, ocean, sea_ice;
                                             radiation,
                                             sea_ice_ocean_heat_flux = flux)
            coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation, interfaces)
            flux_form = coupled_model.interfaces.sea_ice_ocean_interface.flux_formulation
            flux_form isa TwoEquationHeatFlux
        end

        # Test time stepping with each formulation
        for sea_ice_ocean_heat_flux in [IceBathHeatFlux(),
                                        TwoEquationHeatFlux(),
                                        ThreeEquationHeatFlux()]

            @testset "Time stepping with $(nameof(typeof(sea_ice_ocean_heat_flux)))" begin
                interfaces = ComponentInterfaces(atmosphere, ocean, sea_ice;
                                                 radiation,
                                                 sea_ice_ocean_heat_flux)
                coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation, interfaces)
                @test begin
                    time_step!(coupled_model, 60)
                    true
                end
            end
        end
    end
end
