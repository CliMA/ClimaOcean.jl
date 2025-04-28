include("runtests_setup.jl")

using ClimaOcean.OceanSeaIceModels.InterfaceComputations:
                                   ComponentInterfaces,
                                   celsius_to_kelvin,
                                   convert_to_kelvin,
                                   SimilarityScales,
                                   saturation_specific_humidity,
                                   surface_flux,
                                   SkinTemperature,
                                   BulkTemperature,
                                   DiffusiveFlux

using Thermodynamics
using CUDA
using KernelAbstractions: @kernel, @index
using Oceananigans.TimeSteppers: update_state!
using Oceananigans.Units: hours, days
using ClimaOcean.DataWrangling: all_dates

using ClimaSeaIce.SeaIceMomentumEquations
using ClimaSeaIce.Rheologies

import ClimaOcean.OceanSeaIceModels.InterfaceComputations: saturation_specific_humidity

using Statistics: mean, std

struct FixedSpecificHumidity{FT}
    qₒ :: FT
end

@inline saturation_specific_humidity(h::FixedSpecificHumidity, args...) = h.qₒ

@testset "Test surface fluxes" begin
    for arch in test_architectures
        grid = LatitudeLongitudeGrid(arch;
                                     size = 1,
                                     latitude = 10,
                                     longitude = 10,
                                     z = (-1, 0),
                                     topology = (Flat, Flat, Bounded))

        ocean = ocean_simulation(grid;
                                 momentum_advection = nothing,
                                 tracer_advection = nothing,
                                 closure = nothing,
                                 bottom_drag_coefficient = 0.0)

        dates = all_dates(RepeatYearJRA55(), :temperature)
        atmosphere = JRA55PrescribedAtmosphere(arch, Float64; end_date=dates[2], backend = InMemory())

        CUDA.@allowscalar begin
            h  = atmosphere.surface_layer_height
            pₐ = atmosphere.pressure[1][1, 1, 1]

            Tₐ = 15 + celsius_to_kelvin
            qₐ = 0.003

            uₐ = atmosphere.velocities.u[1][1, 1, 1]
            vₐ = atmosphere.velocities.v[1][1, 1, 1]

            ℂₐ = atmosphere.thermodynamics_parameters

            fill!(parent(atmosphere.tracers.T),    Tₐ)
            fill!(parent(atmosphere.tracers.q),    qₐ)
            fill!(parent(atmosphere.velocities.u), uₐ)
            fill!(parent(atmosphere.velocities.v), vₐ)
            fill!(parent(atmosphere.pressure),     pₐ)

            # Force the saturation humidity of the ocean to be
            # equal to the atmospheric saturation humidity
            atmosphere_ocean_interface_specific_humidity = FixedSpecificHumidity(qₐ)

            # Thermodynamic parameters of the atmosphere
            𝒬ₐ = Thermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₐ, qₐ)
            cp = Thermodynamics.cp_m(ℂₐ, 𝒬ₐ)
            ρₐ = Thermodynamics.air_density(ℂₐ, 𝒬ₐ)
            ℰv = Thermodynamics.latent_heat_vapor(ℂₐ, 𝒬ₐ)

            # No radiation equivalent
            radiation = Radiation(ocean_emissivity=0, ocean_albedo=1)

            # turbulent fluxes that force a specific humidity at the ocean's surface
            for atmosphere_ocean_interface_temperature in (BulkTemperature(), SkinTemperature(DiffusiveFlux(1, 1e-2)))
                @info " Testing zero fluxes with $(atmosphere_ocean_interface_temperature)..."

                interfaces = ComponentInterfaces(atmosphere, ocean;
                                                 radiation,
                                                 atmosphere_ocean_interface_specific_humidity,
                                                 atmosphere_ocean_interface_temperature)

                g = ocean.model.buoyancy.formulation.gravitational_acceleration

                # Ensure that the ΔT between atmosphere and ocean is zero
                # Note that the Δθ accounts for the "lapse rate" at height h
                Tₒ = Tₐ - celsius_to_kelvin + h / cp * g

                fill!(parent(ocean.model.velocities.u), uₐ)
                fill!(parent(ocean.model.velocities.v), vₐ)
                fill!(parent(ocean.model.tracers.T), Tₒ)

                # Compute the turbulent fluxes (neglecting radiation)
                coupled_model    = OceanSeaIceModel(ocean; atmosphere, interfaces)
                turbulent_fluxes = coupled_model.interfaces.atmosphere_ocean_interface.fluxes

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

            interfaces = ComponentInterfaces(atmosphere, ocean;
                                             atmosphere_ocean_flux_formulation=similarity_theory)

            # mid-latitude ocean conditions
            set!(ocean.model, u = 0, v = 0, T = 15, S = 30)

            coupled_model = OceanSeaIceModel(ocean; atmosphere, interfaces)

            # Now manually compute the fluxes:
            Tₒ = ocean.model.tracers.T[1, 1, 1] + celsius_to_kelvin
            Sₒ = ocean.model.tracers.S[1, 1, 1]

            interface_properties = interfaces.atmosphere_ocean_interface.properties
            q_formulation = interface_properties.specific_humidity_formulation
            qₒ = saturation_specific_humidity(q_formulation, ℂₐ, 𝒬ₐ.ρ, Tₒ, Sₒ)
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

            turbulent_fluxes = coupled_model.interfaces.atmosphere_ocean_interface.fluxes

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
                                                 coriolis = nothing,
                                                  closure = nothing,
                                  bottom_drag_coefficient = 0.0)

        dates = all_dates(RepeatYearJRA55(), :temperature)
        atmosphere = JRA55PrescribedAtmosphere(arch; end_date=dates[2], backend = InMemory())

        fill!(ocean.model.tracers.T, -2.0)

        CUDA.@allowscalar begin
            ocean.model.tracers.T[1, 2, 10] = 1.0
            ocean.model.tracers.T[2, 1, 10] = 1.0

            # Cap all fluxes exept for heating ones where T < 0
            sea_ice = FreezingLimitedOceanTemperature()

            # Always cooling!
            fill!(atmosphere.tracers.T, 273.15 - 20)

            coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

            # Test that the temperature has snapped up to freezing
            @test minimum(ocean.model.tracers.T) == 0
        end

        @info "Testing Surface Fluxes with sea ice..."

        grid = RectilinearGrid(arch;
                               size = (2, 2, 2),
                             extent = (1, 1, 1),
                           topology = (Periodic, Periodic, Bounded))

        ocean = ocean_simulation(grid; momentum_advection = nothing,
                                         tracer_advection = nothing,
                                                 coriolis = nothing,
                                                  closure = nothing,
                                  bottom_drag_coefficient = 0.0)

        SSU = view(ocean.model.velocities.u, :, :, grid.Nz)
        SSV = view(ocean.model.velocities.v, :, :, grid.Nz)

        τo  = SemiImplicitStress(uₑ=SSU, vₑ=SSV, Cᴰ=0.001, ρₑ=1000.0)
        τua = Field{Face, Center, Nothing}(grid)
        τva = Field{Center, Face, Nothing}(grid)

        dynamics = SeaIceMomentumEquation(grid;
                                          top_momentum_stress = (u=τua, v=τva),
                                          bottom_momentum_stress = τo,
                                          rheology = nothing,
                                          solver = ExplicitSolver())

        sea_ice = sea_ice_simulation(grid; dynamics, advection=Centered())

        # Set a velocity for the ocean
        fill!(ocean.model.velocities.u, 0.1)
        fill!(ocean.model.velocities.v, 0.2)
        fill!(ocean.model.tracers.T,   -2.0)

        # Test that we populate the sea-ice ocean stress
        earth = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation=Radiation())

        τx = earth.interfaces.sea_ice_ocean_interface.fluxes.x_momentum
        τy = earth.interfaces.sea_ice_ocean_interface.fluxes.y_momentum

        CUDA.@allowscalar begin
            @test τx[1, 1, 1] == sqrt(0.1^2 + 0.2^2) * 0.1
            @test τy[1, 1, 1] == sqrt(0.1^2 + 0.2^2) * 0.2
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

        T_metadata = Metadatum(:temperature, date=DateTimeProlepticGregorian(1993, 1, 1), dataset=ECCO4Monthly())
        S_metadata = Metadatum(:salinity,    date=DateTimeProlepticGregorian(1993, 1, 1), dataset=ECCO4Monthly())

        set!(ocean.model; T=T_metadata, S=S_metadata)

        end_date   = all_dates(RepeatYearJRA55(), :temperature)[10]
        atmosphere = JRA55PrescribedAtmosphere(arch; end_date, backend = InMemory())
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
        @test_broken Jᵀ_mean ≈ -3.526464713488678e-5
        @test_broken Jˢ_mean ≈ 1.1470078542716042e-6
        @test_broken τˣ_mean ≈ -1.0881334225579832e-5
        @test_broken τʸ_mean ≈ 5.653281786086694e-6

        @test_broken Jᵀ_std ≈ 7.477575901188957e-5
        @test_broken Jˢ_std ≈ 3.7416720607945508e-6
        @test_broken τˣ_std ≈ 0.00011349625113971719
        @test_broken τʸ_std ≈ 7.627885224680635e-5
    end
end

