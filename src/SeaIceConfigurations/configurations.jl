"""
    latitude_longitude_sea_ice(ocean; kwargs...)

Construct a sea ice `Simulation` for a latitude-longitude ocean.
The grid is extracted from the ocean simulation.
All keyword arguments are forwarded to `sea_ice_simulation`.
"""
latitude_longitude_sea_ice(ocean; kwargs...) =
    sea_ice_simulation(ocean.model.grid, ocean; kwargs...)

"""
    half_degree_tripolar_sea_ice(ocean; kwargs...)

Construct a sea ice `Simulation` for a half-degree tripolar ocean.
All keyword arguments are forwarded to `sea_ice_simulation`.
"""
half_degree_tripolar_sea_ice(ocean; kwargs...) =
    sea_ice_simulation(ocean.model.grid, ocean; kwargs...)

"""
    one_degree_tripolar_sea_ice(ocean; kwargs...)

Construct a sea ice `Simulation` for a one-degree tripolar ocean.
All keyword arguments are forwarded to `sea_ice_simulation`.
"""
one_degree_tripolar_sea_ice(ocean; kwargs...) =
    sea_ice_simulation(ocean.model.grid, ocean; kwargs...)

"""
    sixth_degree_tripolar_sea_ice(ocean; kwargs...)

Construct a sea ice `Simulation` for a sixth-degree tripolar ocean.
All keyword arguments are forwarded to `sea_ice_simulation`.
"""
sixth_degree_tripolar_sea_ice(ocean; kwargs...) =
    sea_ice_simulation(ocean.model.grid, ocean; kwargs...)

"""
    orca_sea_ice(ocean; kwargs...)

Construct a sea ice `Simulation` for an ORCA ocean.
All keyword arguments are forwarded to `sea_ice_simulation`.
"""
orca_sea_ice(ocean; kwargs...) =
    sea_ice_simulation(ocean.model.grid, ocean; kwargs...)
