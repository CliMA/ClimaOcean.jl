include("runtests_setup.jl")

using ClimaOcean.OceanSeaIceModels: IceBathHeatFlux,
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
        @test flux.heat_transfer_coefficient == 0.006

        flux2 = IceBathHeatFlux(heat_transfer_coefficient = 0.01,
                                friction_velocity = MomentumBasedFrictionVelocity())

        @test flux2.heat_transfer_coefficient == 0.01
        @test flux2.friction_velocity isa MomentumBasedFrictionVelocity

        flux3 = IceBathHeatFlux(Float32)
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
        # Test parameters
        liquidus = LinearLiquidus(Float64)
        αₕ = 0.0095  # Heat transfer coefficient
        αₛ = αₕ / 35  # Salt transfer coefficient (R = 35)
        u★ = 0.002   # Friction velocity
        L  = 334e3   # Latent heat of fusion (J/kg)
        ρₒ = 1025.0  # Ocean reference density (kg/m³)
        cₒ = 3991.0  # Ocean heat capacity (J/kg/K)

        @testset "Warm ocean (melting conditions)" begin
            Tₒ = 2.0    # Ocean temperature well above freezing
            Sₒ = 35.0   # Ocean salinity
            Sᵢ = 5.0    # Ice salinity

            Tᵦ, Sᵦ, q = solve_interface_conditions(Tₒ, Sₒ, Sᵢ, αₕ, αₛ, u★, L, ρₒ, cₒ, liquidus)

            # Interface salinity should be between ice and ocean salinity
            @test Sᵦ >= Sᵢ
            @test Sᵦ <= Sₒ

            # Interface temperature should be at freezing point of interface salinity
            Tₘ = melting_temperature(liquidus, Sᵦ)
            @test Tᵦ ≈ Tₘ

            # Warm ocean should cause melting (q > 0)
            @test q > 0
        end

        @testset "Cool ocean (weak melting)" begin
            # Ocean just above freezing - weak melting conditions
            Sₒ = 35.0
            Tₘ_ocean = melting_temperature(liquidus, Sₒ)
            Tₒ = Tₘ_ocean + 0.5  # Ocean 0.5°C above freezing
            Sᵢ = 5.0

            Tᵦ, Sᵦ, q = solve_interface_conditions(Tₒ, Sₒ, Sᵢ, αₕ, αₛ, u★, L, ρₒ, cₒ, liquidus)

            # Interface salinity should be between ice and ocean salinity
            @test Sᵦ >= Sᵢ
            @test Sᵦ <= Sₒ

            # Interface temperature should be at freezing point
            Tₘ = melting_temperature(liquidus, Sᵦ)
            @test Tᵦ ≈ Tₘ

            # Ocean above freezing should still cause melting (q > 0)
            @test q > 0
        end

        @testset "Ocean near freezing point" begin
            Sₒ = 35.0
            Tₘ_ocean = melting_temperature(liquidus, Sₒ)
            Tₒ = Tₘ_ocean + 0.01  # Just slightly above freezing
            Sᵢ = 5.0

            Tᵦ, Sᵦ, q = solve_interface_conditions(Tₒ, Sₒ, Sᵢ, αₕ, αₛ, u★, L, ρₒ, cₒ, liquidus)

            # Interface salinity should still be bounded
            @test Sᵦ >= Sᵢ
            @test Sᵦ <= Sₒ

            # Interface temperature should be at freezing point
            Tₘ = melting_temperature(liquidus, Sᵦ)
            @test Tᵦ ≈ Tₘ

            # Very small melt rate expected
            @test abs(q) < 1e-6
        end

        @testset "Various salinity conditions" begin
            Tₒ = 1.0  # Warm ocean
            Sᵢ = 5.0  # Ice salinity

            for Sₒ in [30.0, 33.0, 35.0, 37.0, 40.0]
                Tᵦ, Sᵦ, q = solve_interface_conditions(Tₒ, Sₒ, Sᵢ, αₕ, αₛ, u★, L, ρₒ, cₒ, liquidus)

                # Interface salinity must always be bounded
                @test Sᵦ >= Sᵢ
                @test Sᵦ <= Sₒ

                # Interface temperature at freezing point
                Tₘ = melting_temperature(liquidus, Sᵦ)
                @test Tᵦ ≈ Tₘ
            end
        end

        @testset "Zero ice salinity" begin
            Tₒ = 1.0
            Sₒ = 35.0
            Sᵢ = 0.0  # Fresh ice

            Tᵦ, Sᵦ, q = solve_interface_conditions(Tₒ, Sₒ, Sᵢ, αₕ, αₛ, u★, L, ρₒ, cₒ, liquidus)

            @test Sᵦ >= Sᵢ
            @test Sᵦ <= Sₒ
            @test Tᵦ ≈ melting_temperature(liquidus, Sᵦ)
        end

        @testset "High friction velocity" begin
            Tₒ = 1.0
            Sₒ = 35.0
            Sᵢ = 5.0
            u★_high = 0.1  # High turbulence

            Tᵦ, Sᵦ, q = solve_interface_conditions(Tₒ, Sₒ, Sᵢ, αₕ, αₛ, u★_high, L, ρₒ, cₒ, liquidus)

            @test Sᵦ >= Sᵢ
            @test Sᵦ <= Sₒ
            @test Tᵦ ≈ melting_temperature(liquidus, Sᵦ)
        end

        @testset "Low friction velocity" begin
            Tₒ = 1.0
            Sₒ = 35.0
            Sᵢ = 5.0
            u★_low = 0.0001  # Very low turbulence

            Tᵦ, Sᵦ, q = solve_interface_conditions(Tₒ, Sₒ, Sᵢ, αₕ, αₛ, u★_low, L, ρₒ, cₒ, liquidus)

            @test Sᵦ >= Sᵢ
            @test Sᵦ <= Sₒ
            @test Tᵦ ≈ melting_temperature(liquidus, Sᵦ)
        end
    end
end

@testset "Salt flux sign conventions in coupled model" begin
    # Test that computed salt fluxes have the correct sign based on ocean temperature:
    # - Warm ocean (T > Tₘ) → melting → q > 0 → Jˢ > 0 (fresh meltwater dilutes ocean)
    # - Cold ocean (T < Tₘ) → freezing → q < 0 → Jˢ < 0 (brine rejection adds salt)
    # Sign convention: Jˢ > 0 means salinity is extracted from ocean

    for arch in test_architectures
        A = typeof(arch)
        @info "Testing salt flux sign conventions on $A"

        grid = LatitudeLongitudeGrid(arch,
                                     size = (4, 4, 1),
                                     latitude = (-80, 80),
                                     longitude = (0, 360),
                                     z = (-100, 0))

        ocean = ocean_simulation(grid, momentum_advection=nothing, closure=nothing, tracer_advection=nothing)
        sea_ice = sea_ice_simulation(grid, ocean)

        backend = JRA55NetCDFBackend(4)
        atmosphere = JRA55PrescribedAtmosphere(arch; backend)
        radiation = Radiation(arch)

        for sea_ice_ocean_heat_flux in [IceBathHeatFlux(), ThreeEquationHeatFlux()]
            @testset "Salt flux with $(nameof(typeof(sea_ice_ocean_heat_flux)))" begin
                interfaces = ComponentInterfaces(atmosphere, ocean, sea_ice;
                                                 radiation,
                                                 sea_ice_ocean_heat_flux)
                coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation, interfaces)

                # Test melting conditions: warm ocean above freezing
                # Freezing point at S=35 is about -1.9°C
                set!(ocean.model, T=2.0, S=35.0)  # Warm ocean
                set!(sea_ice.model, h=1.0, ℵ=1.0, S=5.0)

                time_step!(coupled_model, 60)

                # Get the computed fluxes
                Jˢ = coupled_model.interfaces.sea_ice_ocean_interface.fluxes.salt
                Qᵢ = coupled_model.interfaces.sea_ice_ocean_interface.fluxes.interface_heat

                # Warm ocean should cause melting → Qᵢ > 0 (heat into ice)
                Qᵢ_cpu = Array(interior(Qᵢ, :, :, 1))
                @test all(Qᵢ_cpu .> 0)

                # During melting, fresh meltwater dilutes ocean → Jˢ > 0
                Jˢ_cpu = Array(interior(Jˢ, :, :, 1))
                @test all(Jˢ_cpu .> 0)
            end
        end
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

        # Test time stepping with each formulation
        for sea_ice_ocean_heat_flux in [IceBathHeatFlux(),
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
