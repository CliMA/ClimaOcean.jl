struct Radiation{FT, E, R}
    emission :: E
    reflection :: R
    stefan_boltzmann_constant :: FT
end

function Radiation(arch = CPU(), FT=Float64;
                   ocean_emissivity = 0.97,
                   sea_ice_emissivity = 1.0,
                   ocean_albedo = TabulatedAlbedo(arch, FT),
                   sea_ice_albedo = 0.7,
                   stefan_boltzmann_constant = 5.67e-8)

    ocean_emissivity   isa Number && (ocean_emissivity   = convert(FT, ocean_emissivity))
    sea_ice_emissivity isa Number && (sea_ice_emissivity = convert(FT, sea_ice_emissivity))
    ocean_albedo       isa Number && (ocean_albedo       = convert(FT, ocean_albedo))
    sea_ice_albedo     isa Number && (sea_ice_albedo     = convert(FT, sea_ice_albedo))

    emission = SurfaceProperties(ocean_emissivity, sea_ice_emissivity)
    reflection = SurfaceProperties(ocean_albedo, sea_ice_albedo)

    return Radiation(emission,
                     reflection,
                     convert(FT, stefan_boltzmann_constant))
end

Base.summary(r::Radiation) = "Radiation"
Base.show(io::IO, r::Radiation) = print(io, summary(r))

struct SurfaceProperties{O, I}
    ocean :: O
    sea_ice :: I
end

