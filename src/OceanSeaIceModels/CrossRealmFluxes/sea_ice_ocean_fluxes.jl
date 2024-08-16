using Oceananigans.Operators: Δzᶜᶜᶜ

function compute_sea_ice_ocean_fluxes!(coupled_model)
    return nothing
end

function compute_sea_ice_ocean_salinity_flux!(coupled_model)
    # Compute salinity increment due to changes in ice thickness

    sea_ice = coupled_model.se
    return nothing
end


function sea_ice_ocean_latent_heat_flux!(coupled_model)
    return nothing
end