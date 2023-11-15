module PrescribedAtmospheres

import ..OceanSeaIceModels: surface_velocities

struct PrescribedAtmosphere{U, T}
    velocities :: U
    times :: T
end

surface_velocities(pa::PrescribedAtmosphere) = pa.velocities

end # module
