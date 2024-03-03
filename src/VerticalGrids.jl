module VerticalGrids

export stretched_vertical_faces, PowerLawStretching, LinearStretching

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

"""
    stretched_vertical_faces(; surface_layer_Δz = 5.0,
                               surface_layer_height = 100.0,
                               constant_bottom_spacing_depth = Inf,
                               maximum_Δz = Inf,
                               stretching = PowerLawStretching(1.02),
                               rounding_digits = 1,
                               depth = 5000)

Return an array of cell interfaces with `surface_layer_Δz` spacing in
a surface layer of height `surface_layer_height`, and stretched according to
the function `stretching(Δz_above, z_above)` down to `depth`.
The interfaces extends from `Lz = -z[1]` to `0 = z[end]`, where `Lz ≥ depth`.

The grid spacing `Δz` is limited to be less than `maximum_Δz`.
The grid is also uniformly-spaced below `constant_bottom_spacing_depth`.

`rounding_digits` controls the accuracy with which the grid face positions are saved.
"""
function stretched_vertical_faces(; surface_layer_Δz = 5.0,
                                    surface_layer_height = 100.0,
                                    constant_bottom_spacing_depth = Inf,
                                    maximum_Δz = Inf,
                                    stretching = PowerLawStretching(1.02),
                                    rounding_digits = 1,
                                    depth = 5000)

    Δz₀ = surface_layer_Δz
    h₀ = surface_layer_height

    # Generate surface layer grid
    z = [-Δz₀ * (k-1) for k = 1:ceil(h₀ / Δz₀)]

    # Generate stretched interior grid
    Lz₀ = depth

    while z[end] > - Lz₀
        Δz_above = z[end-1] - z[end]

        if z[end] > - constant_bottom_spacing_depth
            Δz = stretching(Δz_above, z[end])
            Δz = min(maximum_Δz, Δz)
        else
            Δz = Δz_above
        end

        push!(z, round(z[end] - Δz, digits=rounding_digits))
    end

    # Reverse grid to be right-side-up
    z = reverse(z)

    return z
end

@inline exponential_profile(z; Lz, h) = (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h)) 

function exponential_z_faces(; Nz, depth, h = Nz / 4.5)

    z_faces = exponential_profile.((1:Nz+1); Lz = Nz, h)

    # Normalize
    z_faces .-= z_faces[1]
    z_faces .*= - depth / z_faces[end]
    
    z_faces[1] = 0.0

    return reverse(z_faces)
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
