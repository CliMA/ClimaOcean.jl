using ClimaSeaIce: SeaIceModel

sea_ice_thickness(sea_ice::Simulation{<:SeaIceModel}) = sea_ice.model.ice_thickness
sea_ice_thickness(::Nothing) = nothing

# Nothing yet...