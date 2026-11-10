#####
##### Holding the water of the ocean, sea ice and snow at its initial value
#####

"""
    TotalWaterNormalization(coupled_model)

A callback that holds the total water of the ocean, sea ice and snow at its initial value. A forced configuration
cannot close its surface water budget, because precipitation and runoff are prescribed while evaporation follows
the model state, so most OMIP-2 models remove the net imbalance at every step [Tsujino et al. (2020)](@cite tsujino2020evaluation).
The free surface absorbs the excess, spread uniformly over the wet area, and the temperature and salinity of each
column are rescaled by the same factor so that their content is unchanged.
"""
struct TotalWaterNormalization{H, O, I, N, A, FT}
    surface_height :: H
    ocean_water :: O
    sea_ice_water :: I
    snow_water :: N
    reference_water :: A
    total_water :: A
    correction :: A
    sea_ice_density_ratio :: FT
    snow_density_ratio :: FT
    wet_area :: FT
end

@inline constant_value(ρ) = ρ
@inline constant_value(ρ::ConstantField) = ρ.constant

function TotalWaterNormalization(coupled_model)
    grid = coupled_model.ocean.model.grid
    FT = eltype(grid)
    sea_ice = coupled_model.sea_ice.model
    ρᵒᶜ = coupled_model.interfaces.ocean_properties.reference_density

    unit_area = Field{Center, Center, Nothing}(grid)
    set!(unit_area, 1)
    wet_area = only(Array(interior(compute!(Field(Integral(unit_area))))))

    surface_height = Field{Center, Center, Nothing}(grid)
    ocean_water = Field(Integral(surface_height))
    sea_ice_water = Field(Integral(sea_ice.ice_thickness * sea_ice.ice_concentration))
    snow_water = Field(Integral(sea_ice.snow_thickness * sea_ice.ice_concentration))

    reference_water, total_water, correction = (similar(interior(ocean_water)) for _ in 1:3)

    normalization = TotalWaterNormalization(surface_height, ocean_water, sea_ice_water, snow_water,
                                            reference_water, total_water, correction,
                                            convert(FT, constant_value(sea_ice.sea_ice_density) / ρᵒᶜ),
                                            convert(FT, constant_value(sea_ice.snow_density) / ρᵒᶜ),
                                            convert(FT, wet_area))

    compute_total_water!(normalization, coupled_model)
    reference_water .= total_water

    return normalization
end

# Freshwater-equivalent volume of the ocean, sea ice and snow, up to the constant volume below the free surface
function compute_total_water!(normalization, coupled_model)
    # The free surface displacement lives on a grid with extended halos, whose reductions count them as area
    interior(normalization.surface_height) .= interior(coupled_model.ocean.model.free_surface.displacement)

    normalization.total_water .= interior(compute!(normalization.ocean_water)) .+
        normalization.sea_ice_density_ratio .* interior(compute!(normalization.sea_ice_water)) .+
        normalization.snow_density_ratio .* interior(compute!(normalization.snow_water))

    return normalization.total_water
end

@kernel function _correct_sea_level!(η, T, S, grid, δ, kᴺ)
    i, j = @index(Global, NTuple)
    h = static_column_depthᶜᶜᵃ(i, j, grid)
    @inbounds ηⁿ = η[i, j, kᴺ]
    @inbounds δᵢ = δ[1, 1, 1]
    H = h + ηⁿ
    correctable = (h > 0) & (H - δᵢ > 0)
    f = ifelse(correctable, H / (H - δᵢ), one(grid))

    for k in 1:size(grid, 3)
        @inbounds T[i, j, k] *= f
        @inbounds S[i, j, k] *= f
    end

    @inbounds η[i, j, kᴺ] = ηⁿ - ifelse(correctable, δᵢ, zero(grid))
end

function (normalization::TotalWaterNormalization)(simulation)
    coupled_model = simulation.model
    ocean = coupled_model.ocean.model
    grid = ocean.grid
    Nx, Ny, Nz = size(grid)

    compute_total_water!(normalization, coupled_model)
    normalization.correction .= (normalization.total_water .- normalization.reference_water) ./ normalization.wet_area

    launch!(architecture(grid), grid, KernelParameters(1:Nx, 1:Ny), _correct_sea_level!,
            ocean.free_surface.displacement, ocean.tracers.T, ocean.tracers.S, grid, normalization.correction, Nz + 1)

    update_zstar_scaling!(grid, ocean.free_surface.displacement)

    return nothing
end

Oceananigans.prognostic_state(normalization::TotalWaterNormalization) =
    (; reference_water = only(Array(normalization.reference_water)))

function Oceananigans.restore_prognostic_state!(normalization::TotalWaterNormalization, state)
    normalization.reference_water .= state.reference_water
    return normalization
end
