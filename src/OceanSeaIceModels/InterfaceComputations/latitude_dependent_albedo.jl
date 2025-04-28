using Oceananigans.Grids: ηnode

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

Arguments
=========

- `FT::DataType`: The data type of the albedo values. Default is `Float64`.

Keyword Arguments
=================

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
    formula = string(prettysummary(α.diffuse), " - ", prettysummary(α.direct), " ⋅ cos(2φ)")
    return string("LatitudeDepedendentAlbedo{$FT}: ", formula)
end

Base.show(io::IO, α::LatitudeDependentAlbedo) = print(io, summary(α))

@inline function stateindex(α::LatitudeDependentAlbedo, i, j, k, grid, time)
    φ = ηnode(i, j, k, grid, Center(), Center(), Center())
    α₀ = α.diffuse
    α₁ = α.direct
    return α₀ - α₁ * hack_cosd(2φ)
end

