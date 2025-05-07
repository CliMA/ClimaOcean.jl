using Oceananigans.Grids: φnode

@inline hack_cosd(φ) = cos(π * φ / 180)
@inline hack_sind(φ) = sin(π * φ / 180)

struct Radiation{FT, E, R}
    emission :: E
    reflection :: R
    stefan_boltzmann_constant :: FT
end

Adapt.adapt_structure(to, r :: Radiation) =  Radiation(Adapt.adapt(to, r.emission),
                                                       Adapt.adapt(to, r.reflection),
                                                       Adapt.adapt(to, r.stefan_boltzmann_constant))

"""
    Radiation([arch = CPU(), FT=Float64];
              ocean_emissivity = 0.97,
              sea_ice_emissivity = 1.0,
              ocean_albedo = LatitudeDependentAlbedo(FT),
              sea_ice_albedo = 0.7,
              stefan_boltzmann_constant = 5.67e-8)

Constructs a `Radiation` object that represents the radiation properties of the ocean and sea ice.

Arguments
=========

- `arch`: The architecture of the system. Default: `CPU()`.
- `FT`: The floating-point type to use. Default: `Float64`.

Keyword Arguments
=================

- `ocean_emissivity`: The emissivity of the ocean surface. Default: `0.97`.
- `sea_ice_emissivity`: The emissivity of the sea ice surface. Default: `1.0`.
- `ocean_albedo`: The albedo of the ocean surface. Default: `LatitudeDependentAlbedo(FT)`.
- `sea_ice_albedo`: The albedo of the sea ice surface. Default: `0.7`.
- `stefan_boltzmann_constant`: The Stefan-Boltzmann constant. Default: `5.67e-8`.
"""
function Radiation(arch = CPU(), FT=Oceananigans.defaults.FloatType;
                   ocean_emissivity = 0.97,
                   sea_ice_emissivity = 1.0,
                   ocean_albedo = LatitudeDependentAlbedo(FT),
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

const CCC = (Center, Center, Center)

@inline function upwelling_radiation(i, j, k, grid, time, T, σ, ϵ) 
    ϵi = stateindex(ϵ, i, j, k, grid, time, CCC)
    return σ * ϵi * T^4
end

# Split the individual bands
@inline function downwelling_longwave_radiation(i, j, k, grid, time, ϵ, Qℓ)  
    ϵi = stateindex(ϵ, i, j, k, grid, time, CCC)
    return - ϵi * Qℓ
end

@inline function downwelling_shortwave_radiation(i, j, k, grid, time, α, Qs)  
    αi = stateindex(α, i, j, k, grid, time, CCC, Qs) 
    return - (1 - αi) * Qs 
end

@inline net_downwelling_radiation(i, j, k, grid, time, α, ϵ, Qs, Qℓ) = 
    downwelling_shortwave_radiation(i, j, k, grid, time, α, Qs) + 
    downwelling_longwave_radiation(i, j, k, grid, time, ϵ, Qℓ)

# Inside the solver we lose the grid and time information 
@inline net_downwelling_radiation(Qs, Qℓ, α, ϵ) = - (1 - α) * Qs - ϵ * Qℓ

@inline upwelling_radiation(T, σ, ϵ) = σ * ϵ * T^4