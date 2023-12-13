using ClimaSeaIce: SlabSeaIceModel

sea_ice_thickness(sea_ice::Simulation{<:SlabSeaIceModel}) = sea_ice.model.ice_thickness
sea_ice_thickness(::Nothing) = nothing

# Nothing yet...
