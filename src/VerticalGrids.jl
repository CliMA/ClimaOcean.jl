"""
Implementation of several vertical grid options.
"""
module VerticalGrids

export z_faces,
       z_centers,
       ExponentialFaces,
       exponential_vertical_faces,
       stretched_vertical_faces,
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

struct StretchedFaces{FT, S}
    extent :: FT
    top_layer_minimum_spacing :: FT
    top_layer_height :: FT
    constant_bottom_spacing_depth :: FT
    maximum_spacing :: FT
    stretching :: S

    function StretchedFaces(extent,
                            top_layer_minimum_spacing,
                            top_layer_height,
                            constant_bottom_spacing_depth,
                            maximum_spacing,
                            stretching)

        FT = typeof(extent/2)
        S = typeof(stretching)
        return new{FT, S}(extent, top_layer_minimum_spacing, top_layer_height, constant_bottom_spacing_depth, maximum_spacing, stretching)
    end
end

"""
    stretched_vertical_faces(; depth = 5000,
                             surface_layer_Δz = 5.0,
                             surface_layer_height = 100.0,
                             constant_bottom_spacing_depth = Inf,
                             maximum_spacing = Inf,
                             stretching = PowerLawStretching(1.02))

Return a type that describes a one-dimensional grid with `surface_layer_Δz` spacing
in a surface layer of extent `surface_layer_height`, and stretched according to
the `stretching` down to `depth`.
The interfaces extend from `depth = -z[1]` to `0 = z[end]`, where `Lz ≥ depth`.

The grid spacing `Δz` is limited to be less than `maximum_Δz`.
The grid is also uniformly-spaced below `constant_bottom_spacing_depth`.
"""
function stretched_vertical_faces(; depth = 5000,
                                  surface_layer_Δz = 5.0,
                                  surface_layer_height = 100.0,
                                  constant_bottom_spacing_depth = Inf,
                                  maximum_Δz = Inf,
                                  stretching = PowerLawStretching(1.02))

    return StretchedFaces(depth, surface_layer_Δz, surface_layer_height,
                          constant_bottom_spacing_depth, maximum_Δz, stretching)
end

"""
    z_centers(zgrid; rounding_digits = 2)

Return an array of ``z``-centers for a grid of `zgrid` type.
"""
function z_centers(zgrid; rounding_digits = 2)
    zf = z_faces(zgrid; rounding_digits)
    zc = [(zf[k] + zf[k+1])/2 for k in 1:length(zf)-1]
    return zc
end

"""
    z_faces(zgrid; rounding_digits = 2)

Return an array of ``z``-centers for a grid of `zgrid` type.
"""
function z_faces(zgrid::StretchedFaces; rounding_digits = 2)

    constant_bottom_spacing_depth = zgrid.constant_bottom_spacing_depth
    maximum_Δz = zgrid.maximum_spacing
    stretching = zgrid.stretching
    depth = zgrid.extent

    Δz₀ = zgrid.top_layer_minimum_spacing
    h₀  = zgrid.top_layer_height

    # Generate surface layer grid
    z_faces = [-Δz₀ * (k-1) for k = 1:ceil(h₀ / Δz₀)]

    # Generate stretched interior grid
    Lz₀ = depth

    while z_faces[end] > - Lz₀
        Δz_above = z_faces[end-1] - z_faces[end]

        if z_faces[end] > - constant_bottom_spacing_depth
            Δz = stretching(Δz_above, z_faces[end])
            Δz = min(maximum_Δz, Δz)
        else
            Δz = Δz_above
        end

        push!(z_faces, round(z_faces[end] - Δz, digits=rounding_digits))
    end

    # Reverse grid to be right-side-up
    z_faces = reverse(z_faces)

    return z_faces
end

@inline exponential_profile(z, L, h) = expm1((z + L) / h) / expm1(L / h)

struct ExponentialFaces{FT}
    size :: Int
    extent :: FT
    scale :: FT

    function ExponentialFaces(size::Int, extent, scale)
        FT = typeof(scale)
        return new{FT}(size, extent, scale)
    end
end

"""
    exponential_vertical_faces(; Nz, depth, scale=depth/5)

Return a type that describes a one-dimensional vertical grid with faces that are exponentially
spaced (or, equivalently, with spacings that grow linearly with depth) that has `Nz` cells,
goes down to `depth`, and the exponential scaling is controlled by `scale`.
"""
exponential_vertical_faces(; Nz, depth, scale=depth/5) = ExponentialFaces(Nz, depth, scale)

function z_faces(zgrid::ExponentialFaces; rounding_digits = 2)

    if rounding_digits ≥ 6
        @warn "rounding_digits = $(rounding_digits) seems excessive. It's beyond Float32 accuracy. Reducing rounding_digits to 5."
        rounding_digits=5
    end

    Nz = zgrid.size
    depth, scale = zgrid.extent, zgrid.scale

    k = collect(1:Nz+1)
    scale_index = Nz * scale / depth
    z_faces = exponential_profile.(k, Nz, scale_index)

    # Normalize
    z_faces .-= z_faces[1]
    z_faces .*= - depth / z_faces[end]

    if abs(z_faces[1]) < 10eps(Float32)
        z_faces[1] = 0.0
    end

    if abs(z_faces[end] + depth) < 10eps(Float32)
        z_faces[end] = - depth
    end

    @. z_faces = round(z_faces, digits=rounding_digits)

    return reverse(z_faces)
end

# function KDS_z_faces(; depth,
#                      surface_layer_Δz = 5.0,
#                      maximum_layer_Δz = 100.0)

#     ε = 1e-3 # m
#     s = 1
#     spacing(z) = maximum_layer_Δz * tanh( -z * π / (s * depth)) + ε

#     z_faces = Float64[]
#     Δz = Float64[]

#     push!(z_faces, 0.0)
#     push!(Δz, ε)

#     return z_faces, Δz
# end

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
