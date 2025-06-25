"""
Implementation of several vertical grid options.
"""
module VerticalGrids

export ExponentialInterfaces,
       StretchedInterfaces,
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

struct StretchedInterfaces{S, A} <: Function
    extent :: Float64
    top_layer_minimum_spacing :: Float64
    top_layer_height :: Float64
    constant_bottom_spacing_depth :: Float64
    maximum_spacing :: Float64
    stretching :: S
    faces :: A

    function StretchedInterfaces(extent,
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
        S = typeof(stretching)
        A = typeof(z_faces)
        return new{S, A}(extent, top_layer_minimum_spacing, top_layer_height,
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
    StretchedInterfaces(; depth = 5000,
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

`rounding_digits` controls the accuracy with which the grid interface positions are saved.

Example
=======

```jldoctest stretchedinterfaces
using ClimaOcean

z = StretchedInterfaces(depth = 200,
                        surface_layer_Δz = 20.0,
                        surface_layer_height = 100.0)

[z(k) for k in 1:length(z)+1]

# output

10-element Vector{Float64}:
 -200.74
 -173.42
 -147.82
 -123.8
 -101.23
  -80.0
  -60.0
  -40.0
  -20.0
   -0.0
```
"""
function StretchedInterfaces(; depth = 5000,
                             surface_layer_Δz = 5.0,
                             surface_layer_height = 100.0,
                             constant_bottom_spacing_depth = Inf,
                             maximum_Δz = Inf,
                             stretching = PowerLawStretching(1.02),
                             rounding_digits = 2)

    return StretchedInterfaces(depth,
                               surface_layer_Δz,
                               surface_layer_height,
                               constant_bottom_spacing_depth,
                               maximum_Δz,
                               stretching;
                               rounding_digits)
end

(g::StretchedInterfaces)(k) = g.faces[k]

Base.length(g::StretchedInterfaces) = length(g.faces)-1

@inline rightbiased_exponential_mapping(z, l, r, h) = @. r - (r - l) * expm1((r - z) / h) / expm1((r - l) / h)
@inline leftbiased_exponential_mapping(z, l, r, h)  = @. l + (r - l) * expm1((z - l) / h) / expm1((r - l) / h)

struct ExponentialInterfaces <: Function
    size :: Int
    left :: Float64
    right :: Float64
    scale :: Float64
    bias :: Symbol
end

"""
    ExponentialInterfaces(size::Int, left, right=0; scale=(right-left)/5, bias=:right)

Return a type that describes a one-dimensional coordinate with `N + 1` faces (i.e., `N` cells) that
are exponentially spaced (or, equivalently, with spacings that grow linearly).
The coordinate spans `[left, right]`. The exponential e-folding is controlled by `scale`.
The coordinate interfaces are stacked on the `bias`-side of the domain.

Arguments
=========
- `size`: The number of cells in the coordinate.
- `left`: The left-most interface of the coordinate.
- `right`: The right-most interface of the coordinate. Default: 0.

Keyword Arguments
=================
- `scale` :: The length scale of the exponential e-folding. Default: `(right - left) / 5`
- `bias :: Symbol`: Determine whether left or right biased. Default: `:right`.

Examples
========

```jldoctest exponentialinterfaces
using ClimaOcean

Nz = 10
left = -1000
right = 100

z = ExponentialInterfaces(Nz, left, right)

[z(k) for k in 1:Nz+1]

# output

11-element Vector{Float64}:
 -1000.0
  -564.247649441104
  -299.95048878528615
  -139.64615757253702
   -42.41666580727582
    16.55600197663209
    52.324733072619736
    74.0195651413529
    87.17814594835643
    95.15922864611028
   100.0
```

Above, the default `bias` is `:right`. We can get a left-biased grid via:

```jldoctest exponentialinterfaces
z = ExponentialInterfaces(Nz, left, right, bias=:left)

[z(k) for k in 1:Nz+1]

# output

11-element Vector{Float64}:
 -1000.0
  -995.1592286461103
  -987.1781459483565
  -974.0195651413529
  -952.3247330726198
  -916.556001976632
  -857.5833341927241
  -760.353842427463
  -600.0495112147139
  -335.75235055889596
   100.0
```
"""
ExponentialInterfaces(size::Int, left, right=0; scale=(right-left)/5, bias=:right) =
    ExponentialInterfaces(size, left, right, scale, bias)

function (g::ExponentialInterfaces)(k)
    Nz, left, right, scale = g.size, g.left, g.right, g.scale

    # uniform coordinate
    ξₖ = left + (k-1) * (right - left) / Nz

    # mapped coordinate
    if g.bias === :right
       zₖ = rightbiased_exponential_mapping(ξₖ, left, right, scale)
    elseif g.bias === :left
       zₖ = leftbiased_exponential_mapping(ξₖ, left, right, scale)
    end

    if abs(zₖ - left) < 10eps(Float32)
        zₖ = left
    elseif abs(zₖ - right) < 10eps(Float32)
        zₖ = right
    end

    return zₖ
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
