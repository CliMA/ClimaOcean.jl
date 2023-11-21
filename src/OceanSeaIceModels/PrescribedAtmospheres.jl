module PrescribedAtmospheres

import ..OceanSeaIceModels: surface_velocities

struct PrescribedAtmosphere{U, F, T}
    velocities :: U
    freshwater_fluxes :: F
    tracers :: F
    times :: T
end

surface_velocities(pa::PrescribedAtmosphere) = pa.velocities

end # module
