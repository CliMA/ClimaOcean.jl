module SeaIceConfigurations

using Oceananigans
using Oceananigans.Units

using NumericalEarth.SeaIces: sea_ice_simulation

export latitude_longitude_sea_ice,
       half_degree_tripolar_sea_ice,
       one_degree_tripolar_sea_ice,
       sixth_degree_tripolar_sea_ice,
       orca_sea_ice

include("configurations.jl")

end # module
