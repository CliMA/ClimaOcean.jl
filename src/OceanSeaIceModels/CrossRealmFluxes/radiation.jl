using Oceananigans.Grids: φnode

@inline hack_cosd(φ) = cos(π * φ / 180)
@inline hack_sind(φ) = sin(π * φ / 180)

struct Radiation{FT, E, R}
    emission :: E
    reflection :: R
    stefan_boltzmann_constant :: FT
end

Adapt.adapt_structure(to, r :: Radiation) = 
            Radiation(Adapt.adapt(to, r.emission),
                      Adapt.adapt(to, r.reflection),
                      Adapt.adapt(to, r.stefan_boltzmann_constant))

"""
    Radiation([arch = CPU(), FT=Float64];
              ocean_emissivity = 0.97,
              sea_ice_emissivity = 1.0,
              ocean_albedo = TabulatedAlbedo(arch, FT),
              sea_ice_albedo = 0.7,
              stefan_boltzmann_constant = 5.67e-8)

Constructs a `Radiation` object that represents the radiation properties of the ocean and sea ice.

# Arguments
===========

- `arch`: The architecture of the system (default: `CPU()`).
- `FT`: The floating-point type to use (default: `Float64`).

# Keyword Arguments
===================

- `ocean_emissivity`: The emissivity of the ocean surface (default: `0.97`).
- `sea_ice_emissivity`: The emissivity of the sea ice surface (default: `1.0`).
- `ocean_albedo`: The albedo of the ocean surface (default: `LatitudeDependentAlbedo(FT)`).
- `sea_ice_albedo`: The albedo of the sea ice surface (default: `0.7`).
- `stefan_boltzmann_constant`: The Stefan-Boltzmann constant (default: `5.67e-8`).
"""
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

Base.summary(r::Radiation{FT}) where FT = "Radiation{$FT}"

function Base.show(io::IO, r::Radiation)
    σ = r.stefan_boltzmann_constant

    print(io, summary(r), ":", '\n')
    print(io, "├── stefan_boltzmann_constant: ", prettysummary(σ), '\n')
    print(io, "├── emission: ", summary(r.emission), '\n')
    print(io, "│   ├── ocean: ", prettysummary(r.emission.ocean), '\n')
    print(io, "│   └── sea_ice: ", prettysummary(r.emission.ocean), '\n')
    print(io, "└── reflection: ", summary(r.reflection), '\n')
    print(io, "    ├── ocean: ", prettysummary(r.reflection.ocean), '\n')
    print(io, "    └── sea_ice: ", prettysummary(r.reflection.sea_ice))
end

struct LatitudeDependentAlbedo{FT}
    direct :: FT
    diffuse :: FT
end

"""
    LatitudeDependentAlbedo([FT::DataType=Float64]; diffuse = 0.069, direct = 0.011)

Constructs a `LatitudeDependentAlbedo` object. The albedo of the ocean surface is assumed to be a function of the latitude,
obeying the following formula (Large and Yeager, 2009):

    α(φ) = α.diffuse - α.direct * cos(2φ)

where `φ` is the latitude, `α.diffuse` is the diffuse albedo, and `α.direct` is the direct albedo.

# Arguments
===========

- `FT::DataType`: The data type of the albedo values. Default is `Float64`.

# Keyword Arguments
===================
- `diffuse`: The diffuse albedo value. Default is `0.069`.
- `direct`: The direct albedo value. Default is `0.011`.
"""
function LatitudeDependentAlbedo(FT::DataType=Float64; 
                                 diffuse = 0.069, 
                                 direct = 0.011) 
    
    return LatitudeDependentAlbedo(convert(FT, direct),
                                   convert(FT, diffuse))
end

Adapt.adapt_structure(to, α::LatitudeDependentAlbedo) = 
    LatitudeDependentAlbedo(Adapt.adapt(to, α.direct),                       
                            Adapt.adapt(to, α.diffuse))

function Base.summary(α::LatitudeDependentAlbedo{FT}) where FT
    formula = string(prettysummary(α.diffuse), - prettysummary(α.direct), " cos(2φ)")
    return string("LatitudeDepedendentAlbedo{$FT}: ", formula)
end

Base.show(io::IO, α::LatitudeDependentAlbedo) = print(io, summary(α))

@inline function stateindex(α::LatitudeDependentAlbedo, i, j, k, grid, time) 
    φ = φnode(i, j, k, grid, Center(), Center(), Center())
    α_diffuse = α.diffuse
    direct_correction = α.direct * hack_cosd(2φ)

    return α_diffuse - direct_correction
end

struct SurfaceProperties{O, I}
    ocean :: O
    sea_ice :: I
end

Adapt.adapt_structure(to, s :: SurfaceProperties) = 
    SurfaceProperties(Adapt.adapt(to, s.ocean),
                      Adapt.adapt(to, s.sea_ice))

Base.summary(properties::SurfaceProperties) = "SurfaceProperties"

function Base.show(io::IO, properties::SurfaceProperties)
    print(io, "SurfaceProperties:", '\n')
    print(io, "├── ocean: ", summary(properties.ocean), '\n')
    print(io, "└── sea_ice: ", summary(properties.sea_ice))
end

