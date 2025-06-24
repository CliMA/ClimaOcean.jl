"""
Implementation of several vertical grid options.
"""
module VerticalGrids

export ExponentialFaces,
       StretchedFaces,
       PowerLawStretching,
       LinearStretching

struct PowerLawStretching{T}
    power :: T
end

function (stretching::PowerLawStretching)(Δz, z)
    γ = stretching.power
    return Δz^γ
end

struct LinearStretching{T}
    coefficient :: T
end

function (stretching::LinearStretching)(Δz, z)
    c = stretching.coefficient
    return (1 + c) * Δz
end

struct StretchedFaces{FT, S, A} <: Function
    extent :: FT
    top_layer_minimum_spacing :: FT
    top_layer_height :: FT
    constant_bottom_spacing_depth :: FT
    maximum_spacing :: FT
    stretching :: S
    faces :: A

    function StretchedFaces(extent,
                            top_layer_minimum_spacing,
                            top_layer_height,
                            constant_bottom_spacing_depth,
                            maximum_spacing,
                            stretching;
                            rounding_digits=2)

        z_faces = calculate_stretched_faces(; extent,
                                            top_layer_minimum_spacing,
                                            top_layer_height,
                                            constant_bottom_spacing_depth,
                                            maximum_spacing,
                                            stretching,
                                            rounding_digits)

        FT = typeof(extent/2)
        S = typeof(stretching)
        A = typeof(z_faces)
        return new{FT, S, A}(extent, top_layer_minimum_spacing, top_layer_height,
                             constant_bottom_spacing_depth, maximum_spacing, stretching, z_faces)
    end
end

function calculate_stretched_faces(; extent = 5000,
                                   top_layer_minimum_spacing = 5.0,
                                   top_layer_height = 100.0,
                                   constant_bottom_spacing_depth = Inf,
                                   maximum_spacing = Inf,
                                   stretching = PowerLawStretching(1.02),
                                   rounding_digits = 2)

    Δz₀ = top_layer_minimum_spacing
    h₀  = top_layer_height

    # Generate surface layer grid
    z_faces = [-Δz₀ * (k-1) for k = 1:ceil(h₀ / Δz₀)]

    # Generate stretched interior grid
    Lz₀ = extent

    while z_faces[end] > - Lz₀
        Δz_above = z_faces[end-1] - z_faces[end]

        if z_faces[end] > - constant_bottom_spacing_depth
            Δz = stretching(Δz_above, z_faces[end])
            Δz = min(maximum_spacing, Δz)
        else
            Δz = Δz_above
        end

        push!(z_faces, round(z_faces[end] - Δz, digits=rounding_digits))
    end

    # Reverse grid to be right-side-up
    z_faces = reverse(z_faces)

    return z_faces
end

"""
    StretchedFaces(; depth = 5000,
                   surface_layer_Δz = 5.0,
                   surface_layer_height = 100.0,
                   constant_bottom_spacing_depth = Inf,
                   maximum_spacing = Inf,
                   stretching = PowerLawStretching(1.02),
                   rounding_digits = 2)

Return a type that describes a one-dimensional grid with `surface_layer_Δz` spacing
in a surface layer of extent `surface_layer_height`, and stretched according to
the `stretching` down to `depth`.
The interfaces extend from `depth = -z[1]` to `0 = z[end]`, where `Lz ≥ depth`.

The grid spacing `Δz` is limited to be less than `maximum_Δz`.
The grid is also uniformly-spaced below `constant_bottom_spacing_depth`.

`rounding_digits` controls the accuracy with which the grid face positions are saved.
"""
function StretchedFaces(; depth = 5000,
                        surface_layer_Δz = 5.0,
                        surface_layer_height = 100.0,
                        constant_bottom_spacing_depth = Inf,
                        maximum_Δz = Inf,
                        stretching = PowerLawStretching(1.02),
                        rounding_digits = 2)

    return StretchedFaces(depth,
                          surface_layer_Δz,
                          surface_layer_height,
                          constant_bottom_spacing_depth,
                          maximum_Δz,
                          stretching;
                          rounding_digits)
end

(g::StretchedFaces)(k) = @inbounds g.faces[k]

Base.length(g::StretchedFaces) = length(g.faces)-1

@inline exponential_profile(z, L, h) = @. expm1((z + L) / h) / expm1(L / h)

struct ExponentialFaces{FT} <: Function
    size :: Int
    extent :: FT
    scale :: FT

    @doc """
        ExponentialFaces(size::Int, extent; scale=extent/5)

    Return a type that describes a one-dimensional coordinate with `N+1` faces (i.e., `N` cells) that
    are exponentially spaced (or, equivalently, with spacings that grow linearly with depth)
    with `extent`. The coordinate spans `[-depth, 0]`. The exponential scaling is controlled by
    keyword argument `scale` (default: `extent/5`).
    """
    function ExponentialFaces(size::Int, extent; scale=extent/5)
        FT = typeof(scale)
        return new{FT}(size, extent, scale)
    end
end

Base.summary(g::ExponentialFaces) = "ExponentialFaces"

function Base.show(io::IO, g::ExponentialFaces)
    print(io, summary(g), '\n')
    print(io, "├── size: ", prettysummary(g.size), '\n')
    print(io, "├── extent: ", prettysummary(g.extent), '\n')
    print(io, "└── scale: ", prettysummary(g.scale), '\n')
end

function (g::ExponentialFaces)(k)
    Nz = g.size
    depth, scale = g.extent, g.scale

    scale_index = Nz * scale / depth

    ztop = exponential_profile(1, Nz, scale_index)
    zbottom = exponential_profile(Nz+1, Nz, scale_index)

    # use reverse index so that, e.g., k = Nz + 1 corresponds to the top
    z = exponential_profile(Nz+2-k, Nz, scale_index)

    # Normalize
    z -= ztop
    z *= - depth / (zbottom - ztop)

    if abs(z[1]) < 10eps(Float32)
        z = 0.0
    end

    if abs(z + depth) < 10eps(Float32)
        z = - depth
    end

    return z
end

# Vertical grid with 49 levels.
# Stretched from 10 meters spacing at surface
# to 400 meter at the bottom.
z_49_levels_10_to_400_meter_spacing = [
    -5244.5,
    -4834.0,
    -4446.5,
    -4082.0,
    -3740.5,
    -3422.0,
    -3126.5,
    -2854.0,
    -2604.5,
    -2378.0,
    -2174.45,
    -1993.62,
    -1834.68,
    -1695.59,
    -1572.76,
    -1461.43,
    -1356.87,
    -1255.54,
    -1155.53,
    -1056.28,
     -958.03,
     -861.45,
     -767.49,
     -677.31,
     -592.16,
     -513.26,
     -441.68,
     -378.18,
     -323.18,
     -276.68,
     -238.26,
     -207.16,
     -182.31,
     -162.49,
     -146.45,
     -133.04,
     -121.27,
     -110.47,
     -100.15,
      -90.04,
      -80.01,
      -70.0,
      -60.0,
      -50.0,
      -40.0,
      -30.0,
      -20.0,
      -10.0,
        0.0
]

end # module
