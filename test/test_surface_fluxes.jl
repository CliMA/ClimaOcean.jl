include("runtests_setup.jl")

using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes:    
                                    celsius_to_kelvin, 
                                    convert_to_kelvin, 
                                    SimilarityScales,
                                    seawater_saturation_specific_humidity,
                                    surface_flux,
                                    SkinTemperature, 
                                    BulkTemperature

using Thermodynamics
using CUDA
using KernelAbstractions: @kernel, @index
using Oceananigans.TimeSteppers: update_state!
using Oceananigans.Units: hours, days

import ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: water_saturation_specific_humidity

using Statistics: mean, std

struct FixedSpecificHumidity{FT}
    qₒ :: FT
end

@inline water_saturation_specific_humidity(h::FixedSpecificHumidity, args...) = h.qₒ

# TODO: Remove this when https://github.com/CliMA/Oceananigans.jl/pull/3923 is merged
import Oceananigans.Fields: _fractional_indices
_fractional_indices(at_node, grid, ::Nothing, ::Nothing, ::Nothing) = (nothing, nothing, nothing)

@testset "Test surface fluxes" begin
    for arch in test_architectures
        grid = LatitudeLongitudeGrid(arch;
                                     size = 1, 
                                     latitude = 0, 
                                     longitude = 0,
                                     z = (-1, 0),
                                     topology = (Flat, Flat, Bounded))
        
        ocean = ocean_simulation(grid;
                                 momentum_advection = nothing, 
                                 tracer_advection = nothing, 
                                 closure = nothing,
                                 bottom_drag_coefficient = 0.0)

        atmosphere = JRA55PrescribedAtmosphere(1:2; grid, architecture = arch, backend = InMemory()) 
        
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

            # Thermodynamic parameters of the atmosphere
            𝒬ₐ = Thermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₐ, qₐ)
            cp = Thermodynamics.cp_m(ℂₐ, 𝒬ₐ)
            ρₐ = Thermodynamics.air_density(ℂₐ, 𝒬ₐ)
            ℰv = Thermodynamics.latent_heat_vapor(ℂₐ, 𝒬ₐ)

            # turbulent fluxes that force a specific humidity at the ocean's surface
            for Tmode in (BulkTemperature, SkinTemperature)
                @info " Testing zero fluxes with $(Tmode)..."

                turbulent_coefficients = SimilarityTheoryFluxes(surface_temperature_type = Tmode())

                g = ocean.model.buoyancy.formulation.gravitational_acceleration

                # Ensure that the ΔT between atmosphere and ocean is zero 
                # Note that the Δθ accounts for the "lapse rate" at height h
                Tₒ = Tₐ - celsius_to_kelvin + h / cp * g
                
                set!(ocean.model, u = uₐ, v = vₐ, T = Tₒ)

                # Compute the turbulent fluxes (neglecting radiation)
                coupled_model    = OceanSeaIceModel(ocean; atmosphere, turbulent_coefficients, water_vapor_saturation, water_mole_fraction)
                turbulent_fluxes = coupled_model.fluxes.turbulent.fields.ocean

                # Make sure all fluxes are (almost) zero!
                @test turbulent_fluxes.x_momentum[1, 1, 1]    < eps(eltype(grid))
                @test turbulent_fluxes.y_momentum[1, 1, 1]    < eps(eltype(grid))
                @test turbulent_fluxes.sensible_heat[1, 1, 1] < eps(eltype(grid))
                @test turbulent_fluxes.latent_heat[1, 1, 1]   < eps(eltype(grid))
                @test turbulent_fluxes.water_vapor[1, 1, 1]   < eps(eltype(grid))
            end

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
            similarity_theory = SimilarityTheoryFluxes(; roughness_lengths, 
                                                         gustiness_parameter = 0,
                                                         stability_functions)

            # mid-latitude ocean conditions
            set!(ocean.model, u = 0, v = 0, T = 15, S = 30)
            
            coupled_model = OceanSeaIceModel(ocean; atmosphere, turbulent_coefficients=similarity_theory)

            # Now manually compute the fluxes:
            Tₒ = ocean.model.tracers.T[1, 1, 1] + celsius_to_kelvin
            Sₒ = ocean.model.tracers.S[1, 1, 1]
            qₒ = seawater_saturation_specific_humidity(ℂₐ, Tₒ, Sₒ, 𝒬ₐ,
                                                       coupled_model.fluxes.turbulent.water_mole_fraction,
                                                       coupled_model.fluxes.turbulent.water_vapor_saturation,
                                                       Thermodynamics.Liquid())
            
            𝒬ₒ = Thermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₒ, qₒ)
            qₒ = Thermodynamics.vapor_specific_humidity(ℂₐ, 𝒬ₒ)
            g  = ocean.model.buoyancy.formulation.gravitational_acceleration

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

            turbulent_fluxes = coupled_model.fluxes.turbulent.fields.ocean

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

        atmosphere = JRA55PrescribedAtmosphere(1:2; grid, architecture = arch, backend = InMemory())

        fill!(ocean.model.tracers.T, -2.0)

        CUDA.@allowscalar begin
            ocean.model.tracers.T[1, 2, 10] = 1.0
            ocean.model.tracers.T[2, 1, 10] = 1.0

            # Cap all fluxes exept for heating ones where T < 0
            sea_ice = FreezingLimitedOceanTemperature()

            # Always cooling!
            fill!(atmosphere.tracers.T, 273.15 - 20)

            coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation=nothing)

            # Make sure that temperature fluxes are zero when the temperature 
            # is below the minimum but not zero when it is above
            Jᵀ = surface_flux(ocean.model.tracers.T)

            @test Jᵀ[1, 2, 1] != 0.0 # below freezing and cooling, no flux
            @test Jᵀ[2, 1, 1] != 0.0 # below freezing and cooling, no flux
            @test Jᵀ[1, 1, 1] == 0.0 # above freezing and cooling
            @test Jᵀ[2, 2, 1] == 0.0 # above freezing and cooling

            # Test that the temperature has snapped up to freezing
            @test minimum(ocean.model.tracers.T) == 0
        end
    end
end

@testset "Fluxes regression" begin
    for arch in test_architectures
        @info "Testing fluxes regression..."

        grid = LatitudeLongitudeGrid(arch; 
                                     size = (20, 20, 20), 
                                 latitude = (-60, 60), 
                                longitude = (0, 360),
                                        z = (-5000, 0))

        # Speed up compilation by removing all the unnecessary stuff
        momentum_advection = nothing
        tracer_advection   = nothing
        tracers  = (:T, :S)
        buoyancy = nothing
        closure  = nothing
        coriolis = nothing

        ocean = ocean_simulation(grid; momentum_advection, tracer_advection, closure, tracers, coriolis)

        T_metadata = ECCOMetadata(:temperature)
        S_metadata = ECCOMetadata(:salinity)

        set!(ocean.model; T=T_metadata, S=S_metadata)

        atmosphere = JRA55PrescribedAtmosphere(1:10; grid, architecture = arch, backend = InMemory())
        radiation  = Radiation(ocean_albedo=0.1, ocean_emissivity=1.0)
        sea_ice    = nothing

        coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

        @show coupled_model.sea_ice
        
        times = 0:1hours:1days
        Ntimes = length(times)

        # average the fluxes over one day
        Jᵀ = interior(ocean.model.tracers.T.boundary_conditions.top.condition, :, :, 1) ./ Ntimes
        Jˢ = interior(ocean.model.tracers.S.boundary_conditions.top.condition, :, :, 1) ./ Ntimes
        τˣ = interior(ocean.model.velocities.u.boundary_conditions.top.condition, :, :, 1) ./ Ntimes
        τʸ = interior(ocean.model.velocities.v.boundary_conditions.top.condition, :, :, 1) ./ Ntimes

        for time in times[2:end]
            coupled_model.clock.time = time
            update_state!(coupled_model)
            Jᵀ .+= interior(ocean.model.tracers.T.boundary_conditions.top.condition, :, :, 1) ./ Ntimes
            Jˢ .+= interior(ocean.model.tracers.S.boundary_conditions.top.condition, :, :, 1) ./ Ntimes
            τˣ .+= interior(ocean.model.velocities.u.boundary_conditions.top.condition, :, :, 1) ./ Ntimes
            τʸ .+= interior(ocean.model.velocities.v.boundary_conditions.top.condition, :, :, 1) ./ Ntimes
        end

        Jᵀ_mean = mean(Jᵀ) 
        Jˢ_mean = mean(Jˢ) 
        τˣ_mean = mean(τˣ) 
        τʸ_mean = mean(τʸ) 

        Jᵀ_std = std(Jᵀ)
        Jˢ_std = std(Jˢ)
        τˣ_std = std(τˣ)
        τʸ_std = std(τʸ)

        # Regression test
        @test Jᵀ_mean ≈ -3.526464713488678e-5
        @test Jˢ_mean ≈ 1.1470078542716042e-6
        @test τˣ_mean ≈ -1.0881334225579832e-5
        @test τʸ_mean ≈ 5.653281786086694e-6
            
        @test Jᵀ_std ≈ 7.477575901188957e-5
        @test Jˢ_std ≈ 3.7416720607945508e-6
        @test τˣ_std ≈ 0.00011349625113971719
        @test τʸ_std ≈ 7.627885224680635e-5
    end
end

