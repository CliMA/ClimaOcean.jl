include("runtests_setup.jl")

using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes:    
                                    celsius_to_kelvin, 
                                    convert_to_kelvin, 
                                    SimilarityScales,
                                    seawater_saturation_specific_humidity,
                                    surface_flux

using Thermodynamics
using CUDA

import ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: water_saturation_specific_humidity

struct FixedSpecificHumidity{FT}
    qₒ :: FT
end

@inline water_saturation_specific_humidity(h::FixedSpecificHumidity, args...) = h.qₒ

@testset "Test surface fluxes" begin
    @info " Testing zero fluxes..."
    for arch in test_architectures
        grid = LatitudeLongitudeGrid(arch;
                                    size = (1, 1), 
                                latitude = 0, 
                               longitude = (-0.5, 0.5), 
                                       z = (-1, 0),
                                topology = (Periodic, Flat, Bounded))
        
        ocean = ocean_simulation(grid; momentum_advection = nothing, 
                                        tracer_advection = nothing, 
                                                closure = nothing,
                                bottom_drag_coefficient = 0.0)

        atmosphere = JRA55_prescribed_atmosphere(1:2; grid, backend = InMemory())
        
        CUDA.@allowscalar begin
            h  = atmosphere.reference_height
            pₐ = atmosphere.pressure[1][1, 1, 1]

            Tₐ = atmosphere.tracers.T[1][1, 1, 1]
            qₐ = atmosphere.tracers.q[1][1, 1, 1]
            
            uₐ = atmosphere.velocities.u[1][1, 1, 1]
            vₐ = atmosphere.velocities.v[1][1, 1, 1]

            ℂₐ = atmosphere.thermodynamics_parameters

            # Force the saturation humidity of the ocean to be 
            # equal to the atmospheric saturation humidity
            water_vapor_saturation = FixedSpecificHumidity(qₐ)
            water_mole_fraction = 1

            # turbulent fluxes that force a specific humidity at the ocean's surface
            similarity_theory = SimilarityTheoryTurbulentFluxes(grid; water_vapor_saturation, water_mole_fraction)

            # Thermodynamic parameters of the atmosphere
            g  = similarity_theory.gravitational_acceleration
            𝒬ₐ = Thermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₐ, qₐ)
            cp = Thermodynamics.cp_m(ℂₐ, 𝒬ₐ)
            ρₐ = Thermodynamics.air_density(ℂₐ, 𝒬ₐ)
            ℰv = Thermodynamics.latent_heat_vapor(ℂₐ, 𝒬ₐ)

            # Ensure that the ΔT between atmosphere and ocean is zero 
            # Note that the Δθ accounts for the "lapse rate" at height h
            Tₒ = Tₐ - celsius_to_kelvin + h / cp * g
            
            set!(ocean.model, u = uₐ, v = vₐ, T = Tₒ)

            # Compute the turbulent fluxes (neglecting radiation)
            coupled_model    = OceanSeaIceModel(ocean; atmosphere, similarity_theory)
            turbulent_fluxes = coupled_model.fluxes.turbulent.fields

            # Make sure all fluxes are (almost) zero!
            @test turbulent_fluxes.x_momentum[1, 1, 1]    < eps(eltype(grid))
            @test turbulent_fluxes.y_momentum[1, 1, 1]    < eps(eltype(grid))
            @test turbulent_fluxes.sensible_heat[1, 1, 1] < eps(eltype(grid))
            @test turbulent_fluxes.latent_heat[1, 1, 1]   < eps(eltype(grid))
            @test turbulent_fluxes.water_vapor[1, 1, 1]   < eps(eltype(grid))

            @info " Testing neutral fluxes..."
            
            # Constructing very special fluxes that do not account for stability of
            # the atmosphere, have zero gustiness and a constant roughness length of
            # `1e-4` for momentum, water vapor and temperature
            # For this case we can compute the fluxes by hand.
            ℓ = 1e-4
            
            @inline zero_stability_function(ζ) = zero(ζ)

            stability_functions = SimilarityScales(zero_stability_function, 
                                                zero_stability_function, 
                                                zero_stability_function)

            roughness_lengths = SimilarityScales(ℓ, ℓ, ℓ)
            similarity_theory = SimilarityTheoryTurbulentFluxes(grid; 
                                                                roughness_lengths, 
                                                                gustiness_parameter = 0,
                                                                stability_functions)

            # mid-latitude ocean conditions
            set!(ocean.model, u = 0, v = 0, T = 15, S = 30)
            
            coupled_model = OceanSeaIceModel(ocean; atmosphere, similarity_theory)

            # Now manually compute the fluxes:
            Tₒ = ocean.model.tracers.T[1, 1, 1] + celsius_to_kelvin
            Sₒ = ocean.model.tracers.S[1, 1, 1]
            qₒ = seawater_saturation_specific_humidity(ℂₐ, Tₒ, Sₒ, 𝒬ₐ,
                                                    similarity_theory.water_mole_fraction,
                                                    similarity_theory.water_vapor_saturation,
                                                    Thermodynamics.Liquid())
            
            𝒬ₒ = Thermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₒ, qₒ)
            qₒ = Thermodynamics.vapor_specific_humidity(ℂₐ, 𝒬ₒ)

            # Differences!
            Δu = uₐ
            Δv = vₐ
            ΔU = sqrt(Δu^2 + Δv^2)
            Δθ = Tₐ - Tₒ + h / cp * g
            Δq = qₐ - qₒ
            ϰ  = similarity_theory.von_karman_constant
            
            # Characteristic scales
            u★ = ϰ / log(h / ℓ) * ΔU
            θ★ = ϰ / log(h / ℓ) * Δθ
            q★ = ϰ / log(h / ℓ) * Δq

            τx = - ρₐ * u★^2 * Δu / sqrt(Δu^2 + Δv^2)
            τy = - ρₐ * u★^2 * Δv / sqrt(Δu^2 + Δv^2)
            Qs = - ρₐ * cp * u★ * θ★
            Mv = - ρₐ * u★ * q★
            Ql = - ρₐ * u★ * q★ * ℰv

            turbulent_fluxes = coupled_model.fluxes.turbulent.fields

            # Make sure fluxes agree with the hand-calculated ones
            @test turbulent_fluxes.x_momentum[1, 1, 1]    ≈ τx
            @test turbulent_fluxes.y_momentum[1, 1, 1]    ≈ τy
            @test turbulent_fluxes.sensible_heat[1, 1, 1] ≈ Qs
            @test turbulent_fluxes.latent_heat[1, 1, 1]   ≈ Ql
            @test turbulent_fluxes.water_vapor[1, 1, 1]   ≈ Mv

        end

        @info " Testing FreezingLimitedOceanTemperature..." 

        grid = LatitudeLongitudeGrid(arch;
                                    size = (2, 2, 10), 
                                latitude = (-0.5, 0.5), 
                               longitude = (-0.5, 0.5), 
                                       z = (-1, 0),
                                topology = (Bounded, Bounded, Bounded))

        ocean = ocean_simulation(grid; momentum_advection = nothing, 
                                        tracer_advection = nothing, 
                                                closure = nothing,
                                bottom_drag_coefficient = 0.0)

        atmosphere = JRA55_prescribed_atmosphere(1:2; grid, backend = InMemory())

        fill!(ocean.model.tracers.T, -2.0)

        CUDA.@allowscalar begin
            ocean.model.tracers.T[1, 2, 10] = 1.0
            ocean.model.tracers.T[2, 1, 10] = 1.0

            # Cap all fluxes exept for heating ones where T < 0
            sea_ice = FreezingLimitedOceanTemperature()

            # Always cooling!
            fill!(atmosphere.tracers.T, 273.15 - 20)

            coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation=nothing)

            turbulent_fluxes = coupled_model.fluxes.turbulent.fields

            # Make sure that the fluxes are zero when the temperature is below the minimum
            # but not zero when it is above
            u, v, _ = ocean.model.velocities
            T, S    = ocean.model.tracers

            for field in (u, v, T, S)
                flux = surface_flux(field)
                @test flux[1, 2, 1] != 0.0 # below freezing and cooling, no flux
                @test flux[2, 1, 1] != 0.0 # below freezing and cooling, no flux
                @test flux[1, 1, 1] == 0.0 # above freezing and cooling
                @test flux[2, 2, 1] == 0.0 # above freezing and cooling
            end

            # Test that the temperature has snapped up to freezing
            @test minimum(ocean.model.tracers.T) == 0
        end
    end
end


