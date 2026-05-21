# Mixing-length scale algorithm. Default (nn_mxl=2): two-pass gradient limiter
# applied directly on the 3D ℓ field by the column kernel; both ℓ_m and ℓ_ε
# take the same final value.
# NEMO source: zdftke.F90 lines ~526-583.

# Natural buoyancy length: ℓ₀(k) = max(rn_mxl0, sqrt(2·e/N²)).
@inline function natural_length_scale(e_k, N²_k, p)
    FT    = typeof(e_k)
    Nsafe = max(N²_k, FT(1e-32))
    return max(p.minimum_mixing_length, sqrt(2 * e_k / Nsafe))
end
