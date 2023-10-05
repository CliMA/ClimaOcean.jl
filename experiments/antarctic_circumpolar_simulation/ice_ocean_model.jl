using Oceananigans.Operators

using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Models: AbstractModel
using Oceananigans.TimeSteppers: tick!
using Oceananigans.Utils: launch!

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

# Simulations interface
import Oceananigans: fields, prognostic_fields
import Oceananigans.Fields: set!
import Oceananigans.Models: timestepper, NaNChecker, default_nan_checker
import Oceananigans.OutputWriters: default_included_properties
import Oceananigans.Simulations: reset!, initialize!, iteration
import Oceananigans.TimeSteppers: time_step!, update_state!, time
import Oceananigans.Utils: prettytime

struct IceOceanModel{FT, C, G, I, O, S, PI, PC} <: AbstractModel{Nothing}
    clock :: C
    grid :: G # TODO: make it so simulation does not require this
    ice :: I
    previous_ice_thickness :: PI
    previous_ice_concentration :: PC
    ocean :: O
    solar_insolation :: S
    ocean_density :: FT
    ocean_heat_capacity :: FT
    ocean_emissivity :: FT
    stefan_boltzmann_constant :: FT
    reference_temperature :: FT
end

const IOM = IceOceanModel

Base.summary(::IOM) = "IceOceanModel"
prettytime(model::IOM) = prettytime(model.clock.time)
iteration(model::IOM) = model.clock.iteration
timestepper(::IOM) = nothing
reset!(::IOM) = nothing
initialize!(::IOM) = nothing
default_included_properties(::IOM) = tuple()
update_state!(::IOM) = nothing
prognostic_fields(cm::IOM) = nothing
fields(::IOM) = NamedTuple()

function IceOceanModel(ice, ocean; clock = Clock{Float64}(0, 0, 1))
    
    previous_ice_thickness = deepcopy(ice.model.ice_thickness)
    previous_ice_concentration = deepcopy(ice.model.ice_concentration)

    grid = ocean.model.grid
    ice_ocean_thermal_flux = Field{Center, Center, Nothing}(grid)
    ice_ocean_salt_flux = Field{Center, Center, Nothing}(grid)
    solar_insolation = Field{Center, Center, Nothing}(grid)

    ocean_density = 1024
    ocean_heat_capacity = 3991
    ocean_emissivity = 1
    reference_temperature = 273.15
    stefan_boltzmann_constant = 5.67e-8

    # How would we ensure consistency?
    try
        if ice.model.external_thermal_fluxes.top isa RadiativeEmission
            radiation = ice.model.external_thermal_fluxes.top
        else
            radiation = filter(flux isa RadiativeEmission, ice.model.external_thermal_fluxes.top) |> first
        end

        stefan_boltzmann_constant = radiation.stefan_boltzmann_constant
        reference_temperature = radiation.reference_temperature
    catch
    end

    FT = eltype(ocean.model.grid)

    return IceOceanModel(clock,
                         ocean.model.grid,
                         ice,
                         previous_ice_thickness,
                         previous_ice_concentration,
                         ocean,
                         solar_insolation,
                         convert(FT, ocean_density),
                         convert(FT, ocean_heat_capacity),
                         convert(FT, ocean_emissivity),
                         convert(FT, stefan_boltzmann_constant),
                         convert(FT, reference_temperature))
end

time(coupled_model::IceOceanModel) = coupled_model.clock.time

function compute_air_sea_flux!(coupled_model)
    ocean = coupled_model.ocean
    ice = coupled_model.ice

    T = ocean.model.tracers.T
    Nx, Ny, Nz = size(ocean.model.grid)

    grid = ocean.model.grid
    arch = architecture(grid)

    σ = coupled_model.stefan_boltzmann_constant
    ρₒ = coupled_model.ocean_density
    cₒ = coupled_model.ocean_heat_capacity
    Tᵣ = coupled_model.reference_temperature
    Qᵀ = T.boundary_conditions.top.condition
    Tₒ = ocean.model.tracers.T
    hᵢ = ice.model.ice_thickness
    ℵᵢ = ice.model.ice_concentration
    I₀ = coupled_model.solar_insolation

    launch!(arch, grid, :xy, _compute_air_sea_flux!,
            Qᵀ, grid, Tₒ, hᵢ, ℵᵢ, I₀,
            σ, ρₒ, cₒ, Tᵣ)
end

@kernel function _compute_air_sea_flux!(temperature_flux,
                                        grid,
                                        ocean_temperature,
                                        ice_thickness,
                                        ice_concentration,
                                        solar_insolation,
                                        σ, ρₒ, cₒ, Tᵣ)

    i, j = @index(Global, NTuple)

    Nz = size(grid, 3)

    @inbounds begin
        T₀ = ocean_temperature[i, j, Nz] # at the surface
        h = ice_thickness[i, j, 1]
        ℵ = ice_concentration[i, j, 1]
        I₀ = solar_insolation[i, j, 1]
    end

    # Radiation model
    ϵ = 1 # ocean emissivity
    ΣQᵀ = ϵ * σ * (T₀ + Tᵣ)^4 / (ρₒ * cₒ)

    # Also add solar insolation
    ΣQᵀ += I₀ / (ρₒ * cₒ)

    # Set the surface flux only if ice-free
    Qᵀ = temperature_flux

    # @inbounds Qᵀ[i, j, 1] = 0 #(1 - ℵ) * ΣQᵀ
end

function time_step!(coupled_model::IceOceanModel, Δt; callbacks=nothing)
    ocean = coupled_model.ocean
    ice = coupled_model.ice
    liquidus = ice.model.phase_transitions.liquidus
    grid = ocean.model.grid
    ice.Δt = Δt
    ocean.Δt = Δt

    fill_halo_regions!(h)

    # Initialization
    if coupled_model.clock.iteration == 0
        h⁻ = coupled_model.previous_ice_thickness
        hⁿ = coupled_model.ice.model.ice_thickness
        parent(h⁻) .= parent(hⁿ)
    end

    time_step!(ice)

    # TODO: put this in update_state!
    compute_ice_ocean_salinity_flux!(coupled_model)
    ice_ocean_latent_heat!(coupled_model)
    #compute_solar_insolation!(coupled_model)
    #compute_air_sea_flux!(coupled_model)

    time_step!(ocean)

    # TODO:
    # - Store fractional ice-free / ice-covered _time_ for more
    #   accurate flux computation?
    # - Or, input "excess heat flux" into ocean after the ice melts
    # - Currently, non-conservative for heat due bc we don't account for excess
        
    # TODO after ice time-step:
    #   - Adjust ocean temperature if the ice completely melts?
   
    tick!(coupled_model.clock, Δt)
    
    return nothing
end

function compute_ice_ocean_salinity_flux!(coupled_model)
    # Compute salinity increment due to changes in ice thickness

    ice = coupled_model.ice
    ocean = coupled_model.ocean
    grid = ocean.model.grid
    arch = architecture(grid)
    Qˢ = ocean.model.tracers.S.boundary_conditions.top.condition
    Sₒ = ocean.model.tracers.S
    Sᵢ = ice.model.ice_salinity
    Δt = ocean.Δt
    hⁿ = ice.model.ice_thickness
    h⁻ = coupled_model.previous_ice_thickness

    launch!(arch, grid, :xy, _compute_ice_ocean_salinity_flux!,
            Qˢ, grid, hⁿ, h⁻, Sᵢ, Sₒ, Δt)

    return nothing
end


@kernel function _compute_ice_ocean_salinity_flux!(ice_ocean_salinity_flux,
                                                   grid,
                                                   ice_thickness,
                                                   previous_ice_thickness,
                                                   ice_salinity,
                                                   ocean_salinity,
                                                   Δt)
    i, j = @index(Global, NTuple)

    Nz = size(grid, 3)

    hⁿ = ice_thickness
    h⁻ = previous_ice_thickness
    Qˢ = ice_ocean_salinity_flux
    Sᵢ = ice_salinity
    Sₒ = ocean_salinity

    @inbounds begin
        # Thickness of surface grid cell
        Δh = hⁿ[i, j, 1] - h⁻[i, j, 1]

        # Update surface salinity flux.
        # Note: the Δt below is the ocean time-step, eg.
        # ΔS = ⋯ - ∮ Qˢ dt ≈ ⋯ - Δtₒ * Qˢ 
        Qˢ[i, j, 1] = Δh / Δt * (Sᵢ[i, j, 1] - Sₒ[i, j, Nz])

        # Update previous ice thickness
        h⁻[i, j, 1] = hⁿ[i, j, 1]
    end
end

function ice_ocean_latent_heat!(coupled_model)
    ocean = coupled_model.ocean
    ice = coupled_model.ice
    ρₒ = coupled_model.ocean_density
    cₒ = coupled_model.ocean_heat_capacity
    Qₒ = ice.model.external_thermal_fluxes.bottom
    Tₒ = ocean.model.tracers.T
    Sₒ = ocean.model.tracers.S
    Δt = ocean.Δt
    hᵢ = ice.model.ice_thickness

    liquidus = ice.model.phase_transitions.liquidus
    grid = ocean.model.grid
    arch = architecture(grid)

    launch!(arch, grid, :xy, _compute_ice_ocean_latent_heat!,
            Qₒ, grid, hᵢ, Tₒ, Sₒ, liquidus, ρₒ, cₒ, Δt)

    return nothing
end

@kernel function _compute_ice_ocean_latent_heat!(latent_heat,
                                                 grid,
                                                 ice_thickness,
                                                 ocean_temperature,
                                                 ocean_salinity,
                                                 liquidus,
                                                 ρₒ, cₒ, Δt)

    i, j = @index(Global, NTuple)

    Nz = size(grid, 3)
    Qₒ = latent_heat
    hᵢ = ice_thickness
    Tₒ = ocean_temperature
    Sₒ = ocean_salinity

    δQ = zero(grid)
    icy_cell = @inbounds hᵢ[i, j, 1] > 0 # make ice bath approximation then

    @unroll for k = Nz:-1:1
        @inbounds begin
            # Various quantities
            Δz = Δzᶜᶜᶜ(i, j, k, grid)
            Tᴺ = Tₒ[i, j, k]
            Sᴺ = Sₒ[i, j, k]
        end

        # Melting / freezing temperature at the surface of the ocean
        Tₘ = melting_temperature(liquidus, Sᴺ)
                                 
        # Conditions for non-zero ice-ocean flux:
        #   - the ocean is below the freezing temperature, causing formation of ice.
        freezing = Tᴺ < Tₘ 

        #   - We are at the surface and the cell is covered by ice.
        icy_surface_cell = (k == Nz) & icy_cell

        # When there is a non-zero ice-ocean flux, we will instantaneously adjust the
        # temperature of the grid cells accordingly.
        adjust_temperature = freezing | icy_surface_cell

        # Compute change in ocean thermal energy.
        #
        #   - When Tᴺ < Tₘ, we heat the ocean back to melting temperature by extracting heat from the ice,
        #     assuming that the heat flux (which is carried by nascent ice crystals called frazil ice) floats
        #     instantaneously to the surface.
        #
        #   - When Tᴺ > Tₘ and we are in a surface cell covered by ice, we assume equilibrium
        #     and cool the ocean by injecting excess heat into the ice.
        # 
        δEₒ = adjust_temperature * ρₒ * cₒ * (Tₘ - Tᴺ)

        # Perform temperature adjustment
        @inline Tₒ[i, j, k] = ifelse(adjust_temperature, Tₘ, Tᴺ)

        # Compute the heat flux from ocean into ice.
        #
        # A positive value δQ > 0 implies that the ocean is cooled; ie heat
        # is fluxing upwards, into the ice. This occurs when applying the
        # ice bath equilibrium condition to cool down a warm ocean (δEₒ < 0).
        #
        # A negative value δQ < 0 implies that heat is fluxed from the ice into
        # the ocean, cooling the ice and heating the ocean (δEₒ > 0). This occurs when
        # frazil ice is formed within the ocean.
        
        δQ -= δEₒ * Δz / Δt
    end

    # Store ice-ocean flux
    @inbounds Qₒ[i, j, 1] = δQ
end

# Check for NaNs in the first prognostic field (generalizes to prescribed velocitries).
function default_nan_checker(model::IceOceanModel)
    u_ocean = model.ocean.model.velocities.u
    nan_checker = NaNChecker((; u_ocean))
    return nan_checker
end
