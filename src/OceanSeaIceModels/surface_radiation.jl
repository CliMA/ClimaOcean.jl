struct SurfaceRadiation{FT, E, R}
    emission :: E
    reflection :: R
    stefan_boltzmann_constant :: FT
    reference_temperature :: FT
end

function SurfaceRadiation(FT=Float64;
                   ocean_emissivity = 0.97,
                   sea_ice_emissivity = 1.0,
                   ocean_albedo = 0.3,
                   sea_ice_albedo = 0.7,
                   reference_temperature = 273.15,
                   stefan_boltzmann_constant = 5.67e-8)

    ocean_emissivity isa Number && (ocean_emissivity = convert(FT, ocean_emissivity))
    ice_emissivity   isa Number && (ice_emissivity   = convert(FT, ice_emissivity))
    ocean_albedo     isa Number && (ocean_albedo     = convert(FT, ocean_albedo))
    ice_albedo       isa Number && (ice_albedo       = convert(FT, ice_albedo))

    emission = SurfaceProperties(ocean_emissivity, sea_ice_emissivity)
    reflection = SurfaceProperties(ocean_albedo, sea_ice_albedo)

    return SurfaceRadiation(emission,
                            reflection,
                            convert(FT, stefan_boltzmann_constant),
                            convert(FT, reference_temperature))
end

Base.summary(r::SurfaceRadiation) = "SurfaceRadiation"
Base.show(io::IO, r::SurfaceRadiation) = print(io, summary(osf))

struct SurfaceProperties{O, I}
    ocean :: O
    sea_ice :: I
end

