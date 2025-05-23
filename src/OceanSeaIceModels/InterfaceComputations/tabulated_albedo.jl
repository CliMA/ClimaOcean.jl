using Oceananigans.Fields: interpolator
using Oceananigans.Grids: on_architecture
using Oceananigans.Utils: Time
using Base

using ClimaOcean.OceanSeaIceModels:
    PrescribedAtmosphere,
    TwoBandDownwellingRadiation

# Bilinear interpolation of the albedo α in α_table based on a
# transmissivity value (𝓉_values) and latitude (φ_values)
struct TabulatedAlbedo{FT, M, P, T}
    α_table :: M
    φ_values :: P
    𝓉_values :: T
    S₀ :: FT # Solar constant W / m^2
    day_to_radians :: FT
    noon_in_seconds :: Int
end

Adapt.adapt_structure(to, α::TabulatedAlbedo) =
    TabulatedAlbedo(Adapt.adapt(to, α.α_table),
                    Adapt.adapt(to, α.φ_values),
                    Adapt.adapt(to, α.𝓉_values),
                    Adapt.adapt(to, α.S₀),
                    Adapt.adapt(to, α.day_to_radians),
                    Adapt.adapt(to, α.noon_in_seconds))

# Tabulated from Payne (1972) https://doi.org/10.1175/1520-0469(1972)029<0959:AOTSS>2.0.CO;2
const α_payne = [ 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.06
                  0.062 0.062 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.06
                  0.072 0.070 0.068 0.065 0.065 0.063 0.062 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.061 0.060 0.061 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.06
                  0.087 0.083 0.079 0.073 0.070 0.068 0.066 0.065 0.064 0.063 0.062 0.061 0.061 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.06
                  0.115 0.108 0.098 0.086 0.082 0.077 0.072 0.071 0.067 0.067 0.065 0.063 0.062 0.061 0.061 0.060 0.060 0.060 0.060 0.061 0.061 0.061 0.061 0.060 0.059 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.060 0.059 0.059 0.05
                  0.163 0.145 0.130 0.110 0.101 0.092 0.084 0.079 0.072 0.072 0.068 0.067 0.064 0.063 0.062 0.061 0.061 0.061 0.060 0.060 0.060 0.060 0.060 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.05
                  0.235 0.198 0.174 0.150 0.131 0.114 0.103 0.094 0.083 0.080 0.074 0.074 0.070 0.067 0.065 0.064 0.063 0.062 0.061 0.060 0.060 0.060 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.058 0.058 0.05
                  0.318 0.263 0.228 0.192 0.168 0.143 0.127 0.113 0.099 0.092 0.084 0.082 0.076 0.072 0.070 0.067 0.065 0.064 0.062 0.062 0.060 0.060 0.060 0.059 0.059 0.059 0.059 0.059 0.059 0.059 0.058 0.058 0.058 0.058 0.058 0.058 0.058 0.058 0.057 0.058 0.058 0.058 0.058 0.057 0.057 0.05
                  0.395 0.336 0.290 0.248 0.208 0.176 0.151 0.134 0.117 0.107 0.097 0.091 0.085 0.079 0.075 0.071 0.068 0.067 0.065 0.063 0.062 0.061 0.060 0.060 0.060 0.059 0.059 0.058 0.058 0.058 0.057 0.057 0.057 0.057 0.057 0.057 0.057 0.056 0.056 0.056 0.056 0.056 0.056 0.056 0.056 0.05
                  0.472 0.415 0.357 0.306 0.252 0.210 0.176 0.154 0.135 0.125 0.111 0.102 0.094 0.086 0.081 0.076 0.072 0.071 0.068 0.066 0.065 0.063 0.062 0.061 0.060 0.059 0.058 0.057 0.057 0.057 0.056 0.055 0.055 0.055 0.055 0.055 0.055 0.054 0.053 0.054 0.053 0.053 0.054 0.054 0.053 0.05
                  0.542 0.487 0.424 0.360 0.295 0.242 0.198 0.173 0.150 0.136 0.121 0.110 0.101 0.093 0.086 0.081 0.076 0.073 0.069 0.067 0.065 0.064 0.062 0.060 0.059 0.058 0.057 0.056 0.055 0.055 0.054 0.053 0.053 0.052 0.052 0.052 0.051 0.051 0.050 0.050 0.050 0.050 0.051 0.050 0.050 0.05
                  0.604 0.547 0.498 0.407 0.331 0.272 0.219 0.185 0.160 0.141 0.127 0.116 0.105 0.097 0.089 0.083 0.077 0.074 0.069 0.066 0.063 0.061 0.059 0.057 0.056 0.055 0.054 0.053 0.053 0.052 0.051 0.050 0.050 0.049 0.049 0.049 0.048 0.047 0.047 0.047 0.046 0.046 0.047 0.047 0.046 0.04
                  0.655 0.595 0.556 0.444 0.358 0.288 0.236 0.190 0.164 0.145 0.130 0.119 0.107 0.098 0.090 0.084 0.076 0.073 0.068 0.064 0.060 0.058 0.056 0.054 0.053 0.051 0.050 0.049 0.048 0.048 0.047 0.046 0.046 0.045 0.045 0.045 0.044 0.043 0.043 0.043 0.042 0.042 0.043 0.042 0.042 0.04
                  0.693 0.631 0.588 0.469 0.375 0.296 0.245 0.193 0.165 0.145 0.131 0.118 0.106 0.097 0.088 0.081 0.074 0.069 0.065 0.061 0.057 0.055 0.052 0.050 0.049 0.047 0.046 0.046 0.044 0.044 0.043 0.042 0.042 0.041 0.041 0.040 0.040 0.039 0.039 0.039 0.038 0.038 0.038 0.038 0.038 0.03
                  0.719 0.656 0.603 0.480 0.385 0.300 0.250 0.193 0.164 0.145 0.131 0.116 0.103 0.092 0.084 0.076 0.071 0.065 0.061 0.057 0.054 0.051 0.049 0.047 0.045 0.043 0.043 0.042 0.041 0.040 0.039 0.039 0.038 0.038 0.037 0.036 0.036 0.035 0.035 0.034 0.034 0.034 0.034 0.034 0.034 0.03
                  0.732 0.670 0.592 0.474 0.377 0.291 0.246 0.190 0.162 0.144 0.130 0.114 0.100 0.088 0.080 0.072 0.067 0.062 0.058 0.054 0.050 0.047 0.045 0.043 0.041 0.039 0.039 0.038 0.037 0.036 0.036 0.035 0.035 0.034 0.033 0.032 0.032 0.032 0.031 0.031 0.031 0.030 0.030 0.030 0.030 0.03
                  0.730 0.652 0.556 0.444 0.356 0.273 0.235 0.188 0.160 0.143 0.129 0.113 0.097 0.086 0.077 0.069 0.064 0.060 0.055 0.051 0.047 0.044 0.042 0.039 0.037 0.035 0.035 0.035 0.034 0.033 0.033 0.032 0.032 0.032 0.029 0.029 0.029 0.029 0.028 0.028 0.028 0.028 0.027 0.027 0.028 0.02
                  0.681 0.602 0.488 0.386 0.320 0.252 0.222 0.185 0.159 0.142 0.127 0.111 0.096 0.084 0.075 0.067 0.062 0.058 0.054 0.050 0.046 0.042 0.040 0.036 0.035 0.033 0.032 0.032 0.031 0.030 0.030 0.030 0.030 0.029 0.027 0.027 0.027 0.027 0.026 0.026 0.026 0.026 0.026 0.026 0.026 0.02
                  0.581 0.494 0.393 0.333 0.288 0.237 0.211 0.182 0.158 0.141 0.126 0.110 0.095 0.083 0.074 0.066 0.061 0.057 0.053 0.049 0.045 0.041 0.039 0.034 0.033 0.032 0.031 0.030 0.029 0.028 0.028 0.028 0.028 0.027 0.026 0.026 0.026 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.02
                  0.453 0.398 0.342 0.301 0.266 0.226 0.205 0.180 0.157 0.140 0.125 0.109 0.095 0.083 0.074 0.065 0.061 0.057 0.052 0.048 0.044 0.040 0.038 0.033 0.032 0.031 0.030 0.029 0.028 0.027 0.027 0.026 0.026 0.026 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.02
                  0.425 0.370 0.325 0.290 0.255 0.220 0.200 0.178 0.157 0.140 0.122 0.108 0.095 0.083 0.074 0.065 0.061 0.056 0.052 0.048 0.044 0.040 0.038 0.033 0.032 0.031 0.030 0.029 0.028 0.027 0.026 0.026 0.026 0.026 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.02]

"""
    TabulatedAlbedo(arch = CPU(), FT = Float64;
                    S₀ = convert(FT, 1365),
                    α_table  = α_payne,
                    φ_values = (0:2:90) ./ 180 * π,
                    𝓉_values = 0:0.05:1)

Constructs a `TabulatedAlbedo` object that interpolated the albedo from a value table `α_table` that
is function of latitude `φ` and atmospheric transimissivity `𝓉`.

Note: `TabulatedAlbedo` assumes that the latitude and the transissivity in the table are uniformly spaced.

The transmissivity of the atmosphere is calculated as the ratio of the downwelling solar radiation to the
maximum possible downwelling solar radiation for a transparent atmosphere, function of hour of the day, latitude,
and day in the year.

Arguments
=========

- `arch`: The architecture to use. Default: `CPU()`.
- `FT`: The floating-point type to use. Default: `Float64`.

Keyword Arguments
=================

- `S₀`: The solar constant. Default: `convert(FT, 1365)`.
- `α_table`: The table of albedo values. Default: `α_payne`.
- `φ_values`: The latitude values for the table. Default: `(0:2:90) ./ 180 * π`.
- `𝓉_values`: The transmissivity values for the table. Default: `0:0.05:1`.
"""
function TabulatedAlbedo(arch = CPU(), FT = Oceananigans.defaults.FloatType;
                         S₀ = convert(FT, 1365),
                         α_table  = α_payne,
                         φ_values = (0:2:90) ./ 180 * π,
                         𝓉_values = 0:0.05:1,
                         day_to_radians  = convert(FT, 2π / 86400),
                         noon_in_seconds = 86400 ÷ 2) # assumes that midnight is at t = 0 seconds

    # Make everything GPU - ready
    α_table  = on_architecture(arch, convert.(FT, α_table))
    φ_values = on_architecture(arch, convert.(FT, φ_values))
    𝓉_values = on_architecture(arch, convert.(FT, 𝓉_values))

    return TabulatedAlbedo(α_table,
                           φ_values,
                           𝓉_values,
                           convert(FT, S₀),
                           convert(FT, day_to_radians),
                           noon_in_seconds)
end

Base.eltype(::TabulatedAlbedo{FT}) where FT = FT
Base.summary(::TabulatedAlbedo{FT}) where FT = "TabulatedAlbedo{$FT}"
Base.show(io::IO, α::TabulatedAlbedo) = print(io, summary(α))

@inline ϕ₁(ξ, η) = (1 - ξ) * (1 - η)
@inline ϕ₂(ξ, η) = (1 - ξ) *      η
@inline ϕ₃(ξ, η) =      ξ  * (1 - η)
@inline ϕ₄(ξ, η) =      ξ  *      η

# Assumption: if the time is represented by a number it is defined in seconds.
# TODO: extend these functions for `DateTime` times when these are supported in
# Oceananigans.
@inline simulation_day(time::Time{<:Number})      = time.time ÷ 86400
@inline seconds_in_day(time::Time{<:Number}, day) = time.time - day * 86400

@inline function stateindex(α::TabulatedAlbedo, i, j, k, grid, time, loc, Qs)
    FT = eltype(α)
    λ, φ, z = _node(i, j, k, grid, Center(), Center(), Center())

    φ = deg2rad(φ)
    λ = deg2rad(λ)

    day         = simulation_day(time)
    day2rad     = α.day_to_radians
    noon_in_sec = α.noon_in_seconds
    sec_of_day  = seconds_in_day(time, day)

    # Hour angle h
    h = (sec_of_day - noon_in_sec) * day2rad + λ

    # Declination angle δ
    march_first = 80
    δ = deg2rad((23 + 27/60) * sind(360 * (day - march_first) / 365.25))
    δ = convert(FT, δ)

    # Zenith angle of the sun (if smaller than 0 we are in the dark)
    cosθₛ = max(0, sin(φ) * sin(δ) + cos(h) * cos(δ) * cos(φ))

    # Maximum downwelling solar radiation for
    # a transparent atmosphere
    Qmax = α.S₀ * cosθₛ

    # Finding the transmissivity and capping it to 1
    𝓉 = ifelse(Qmax > 0, min(1, Qs / Qmax), 0)

    # finding the i-index in the table (depending on transmissivity)
    # we assume that the transmissivity is tabulated with a constant spacing
    𝓉₁ = @inbounds α.𝓉_values[1]
    Δ𝓉 = @inbounds α.𝓉_values[2] - 𝓉₁
    i⁻, i⁺, ξ = interpolator((𝓉 - 𝓉₁) / Δ𝓉)

    # finding the j-index in the table (depending on latitude)
    # we assume that the transmissivity is tabulated with a constant spacing
    φ₁ = @inbounds α.φ_values[1]
    Δφ = @inbounds α.φ_values[2] - φ₁
    j⁻, j⁺, η = interpolator((abs(φ) - φ₁) / Δφ)

    # Bilinear interpolation!
    α =  @inbounds ϕ₁(ξ, η) * getindex(α.α_table, i⁻, j⁻) +
                   ϕ₂(ξ, η) * getindex(α.α_table, i⁻, j⁺) +
                   ϕ₃(ξ, η) * getindex(α.α_table, i⁺, j⁻) +
                   ϕ₄(ξ, η) * getindex(α.α_table, i⁺, j⁺)

    return α
end

