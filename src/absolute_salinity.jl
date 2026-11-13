#####
##### TEOS-10 Absolute Salinity from Practical Salinity: a port of `gsw_sa_from_sp` of GSW-C
##### (https://github.com/TEOS-10/GSW-C) by way of SeawaterPolynomials.jl#48.
#####

# Standard Ocean Reference Salinity [g kg⁻¹] and the factor uₚₛ = Sₒ / 35 relating it to Practical Salinity
const Sₒ = 35.16504
const uₚₛ = Sₒ / 35

# The atlas flags land with ±9e90 and the lookups return 9e15 when they cannot resolve a point
const invalid_value = 9e15
const invalid_threshold = 1e10
const valid_ratio_threshold = 100

"""
    SAARAtlas

The TEOS-10 Absolute Salinity Anomaly Ratio atlas of [McDougall et al. (2012)](@cite mcdougall2012global):
`ratio[k, j, i]` on the pressures `p` [dbar], latitudes `φ` [°N] and longitudes `λ` [°E], and the deepest valid
pressure index `deepest_level[j, i]` of each column.
"""
struct SAARAtlas
    p :: Vector{Float64}
    φ :: Vector{Float64}
    λ :: Vector{Float64}
    ratio :: Array{Float64, 3}
    deepest_level :: Matrix{Float64}
end

function SAARAtlas(path = joinpath(artifact"gsw_saar_data", "gsw_saar_data.bin"))
    Np, Nφ, Nλ = 45, 45, 91
    p = Vector{Float64}(undef, Np)
    φ = Vector{Float64}(undef, Nφ)
    λ = Vector{Float64}(undef, Nλ)
    ratio = Array{Float64}(undef, Np, Nφ, Nλ)
    deepest_level = Matrix{Float64}(undef, Nφ, Nλ)

    open(path) do io
        read!(io, p)
        read!(io, φ)
        read!(io, λ)
        read!(io, ratio)
        read!(io, deepest_level)
    end

    return SAARAtlas(p, φ, λ, ratio, deepest_level)
end

# The Panama isthmus, below which the anomaly ratios of the Pacific and the Atlantic must not be mixed
const panama_longitudes = [260.00, 272.59, 276.50, 278.65, 280.73, 292.00]
const panama_latitudes  = [ 19.55,  13.97,   9.60,   8.10,   9.33,   3.40]

# Corners of the atlas cell, counterclockwise from the south-west one
const corner_offsets_i = (0, 1, 1, 0)
const corner_offsets_j = (0, 0, 1, 1)

function piecewise_linear(xs, ys, x)
    k = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
    return ys[k] + (x - xs[k]) * (ys[k+1] - ys[k]) / (xs[k+1] - xs[k])
end

isthmus_latitude(λ) = piecewise_linear(panama_longitudes, panama_latitudes, λ)

function mean_of_valid_corners(corners, valid)
    n = count(valid)
    mean = n == 0 ? 0.0 : sum(corners[c] for c in 1:4 if valid[c]) / n
    return ntuple(c -> valid[c] ? corners[c] : mean, 4)
end

function corner_ratios(atlas, i, j, k, λ, φ, inside_panama)
    corners = ntuple(c -> atlas.ratio[k, j + corner_offsets_j[c], i + corner_offsets_i[c]], 4)
    valid = ntuple(c -> abs(corners[c]) <= valid_ratio_threshold, 4)

    if inside_panama
        λʷ, λᵉ = atlas.λ[i], atlas.λ[i+1]
        φˢ, φⁿ = atlas.φ[j], atlas.φ[j+1]
        north_of_isthmus = isthmus_latitude(λ) ≤ φ
        corner_is_north = (isthmus_latitude(λʷ) ≤ φˢ, isthmus_latitude(λᵉ) ≤ φˢ,
                           isthmus_latitude(λᵉ) ≤ φⁿ, isthmus_latitude(λʷ) ≤ φⁿ)
        same_side = ntuple(c -> corner_is_north[c] == north_of_isthmus, 4)
        return mean_of_valid_corners(corners, valid .& same_side)
    end

    return abs(sum(corners)) >= invalid_threshold ? mean_of_valid_corners(corners, valid) : corners
end

function bilinear_ratio(atlas, i, j, k, λ, φ, ξ, η, inside_panama)
    r = corner_ratios(atlas, i, j, k, λ, φ, inside_panama)
    return (1 - η) * (r[1] + ξ * (r[2] - r[1])) + η * (r[4] + ξ * (r[3] - r[4]))
end

"""
    absolute_salinity_anomaly_ratio(atlas, p, λ, φ)

The Absolute Salinity Anomaly Ratio at sea pressure `p` [dbar], longitude `λ` [°E] and latitude `φ` [°N],
interpolated trilinearly in `atlas` (`gsw_saar`).
"""
function absolute_salinity_anomaly_ratio(atlas::SAARAtlas, p, λ, φ)
    -86 ≤ φ ≤ 90 || return invalid_value

    λ = mod(λ, 360)
    Nλ, Nφ = length(atlas.λ), length(atlas.φ)
    i = min(floor(Int, (Nλ - 1) * (λ - atlas.λ[1]) / (atlas.λ[end] - atlas.λ[1])) + 1, Nλ - 1)
    j = min(floor(Int, (Nφ - 1) * (φ - atlas.φ[1]) / (atlas.φ[end] - atlas.φ[1])) + 1, Nφ - 1)

    deepest_level = 0
    for c in 1:4
        level = atlas.deepest_level[j + corner_offsets_j[c], i + corner_offsets_i[c]]
        deepest_level = 0 < level < 1e90 ? max(deepest_level, Int(level)) : deepest_level
    end
    deepest_level == 0 && return 0.0

    p = min(p, atlas.p[deepest_level])
    k = min(searchsortedlast(atlas.p, p), length(atlas.p) - 1)

    ξ = (λ - atlas.λ[i]) / (atlas.λ[i+1] - atlas.λ[i])
    η = (φ - atlas.φ[j]) / (atlas.φ[j+1] - atlas.φ[j])
    ζ = (p - atlas.p[k]) / (atlas.p[k+1] - atlas.p[k])

    inside_panama = (panama_longitudes[1] ≤ λ ≤ panama_longitudes[end] - 0.001) &&
                    (panama_latitudes[end] ≤ φ ≤ panama_latitudes[1])

    r⁺ = bilinear_ratio(atlas, i, j, k,   λ, φ, ξ, η, inside_panama)
    r⁻ = bilinear_ratio(atlas, i, j, k+1, λ, φ, ξ, η, inside_panama)
    r⁻ = abs(r⁻) >= invalid_threshold ? r⁺ : r⁻

    r = r⁺ + ζ * (r⁻ - r⁺)
    return abs(r) >= invalid_threshold ? invalid_value : r
end

# The Baltic Sea polygon of `gsw_sa_from_sp_baltic`, after Feistel et al. (2010)
const baltic_west_longitudes = [12.6, 7.0, 26.0]
const baltic_west_latitudes  = [50.0, 59.0, 69.0]
const baltic_east_longitudes = [45.0, 26.0]
const baltic_east_latitudes  = [50.0, 69.0]

function baltic_absolute_salinity(Sᴾ, λ, φ)
    λ = mod(λ, 360)
    inside_box = (baltic_west_longitudes[2] < λ < baltic_east_longitudes[1]) &&
                 (baltic_west_latitudes[1] < φ < baltic_west_latitudes[3])
    inside_box || return invalid_value

    λʷ = piecewise_linear(baltic_west_latitudes, baltic_west_longitudes, φ)
    λᵉ = piecewise_linear(baltic_east_latitudes, baltic_east_longitudes, φ)

    return λʷ ≤ λ ≤ λᵉ ? (Sₒ - 0.087) / 35 * Sᴾ + 0.087 : invalid_value
end

"""
    absolute_salinity(atlas, Sᴾ, p, λ, φ)

TEOS-10 Absolute Salinity [g kg⁻¹] from Practical Salinity `Sᴾ` at sea pressure `p` [dbar], longitude `λ` [°E]
and latitude `φ` [°N], ``Sᴬ = uₚₛ Sᴾ (1 + r)`` with the anomaly ratio ``r`` of `atlas`, and the analytic
correction of [Feistel et al. (2010)](@cite feistel2010thermophysical) in the Baltic Sea (`gsw_sa_from_sp`).
"""
function absolute_salinity(atlas::SAARAtlas, Sᴾ, p, λ, φ)
    Sᴬ = baltic_absolute_salinity(Sᴾ, λ, φ)
    abs(Sᴬ) < invalid_threshold && return Sᴬ

    r = absolute_salinity_anomaly_ratio(atlas, p, λ, φ)
    abs(r) >= invalid_threshold && return NaN

    return uₚₛ * Sᴾ * (1 + r)
end
