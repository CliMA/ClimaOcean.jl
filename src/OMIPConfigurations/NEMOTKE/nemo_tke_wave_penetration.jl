# Surface-wave penetration source (Mellor & Blumberg 2004), nn_etau=1: exponential
# decay below mixed layer with latitude-dependent length scale (nn_htau=1).
# NEMO source: zdftke.F90 lines ~493-526, 777-784.

# Latitude-dependent decay length:
#   nn_htau = 0  →  h_τ = 10 m
#   nn_htau = 1  →  h_τ = max(0.5, min(30, 45·|sin(φ)|)) m
@inline function wave_decay_length(φ_degrees, p)
    FT = typeof(φ_degrees)
    if p.surface_length_scale_formulation == 0
        return FT(10)
    else
        s = abs(sind(φ_degrees))
        return max(FT(0.5), min(FT(30), 45 * s))
    end
end

# WP source: rn_efr · e_surf · exp(-z/h_τ) · (1 - ℵ)
@inline function wave_penetration_source(z_c, e_surf, h_τ, ℵ, p)
    FT = typeof(z_c)
    return p.Cᶠ * e_surf * exp(-z_c / max(h_τ, FT(1e-10))) * (one(FT) - ℵ)
end
