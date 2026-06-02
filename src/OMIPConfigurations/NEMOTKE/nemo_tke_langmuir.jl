# Langmuir-cell source (Axell 2002), enabled by ln_lc.
# NEMO source: zdftke.F90 lines ~341-369.

# Stokes velocity proxy: u_s = Cˢ · √|τ|. NEMO uses Cˢ = 0.016.
@inline stokes_velocity(τ_magnitude, p) = p.Cˢ * sqrt(τ_magnitude)

# LC source at depth z_c (cell center): (Cᴸ · u_s · sin(π z / h_LC))³ / h_LC
# valid for z_c < h_LC; sin(0) = 0 makes the formula return zero outside that range.
@inline function langmuir_source(z_c, h_LC, u_s, p)
    FT  = typeof(z_c)
    h̃   = max(h_LC, FT(1e-10))
    arg = ifelse(z_c < h_LC, FT(π) * z_c / h̃, zero(FT))
    w   = p.Cᴸ * u_s * sin(arg)
    return w * w * w / h̃
end
