# KPP parameters with calibrated Large 1994 / MITgcm defaults.

"""
    KPPParameters{FT}

Parameters of the K-Profile Parameterization (Large, McWilliams, & Doney 1994)
as implemented in MITgcm/pkg/kpp.
"""
struct KPPParameters{FT}
    # Boundary-layer depth
    Riᶜ                           :: FT
    Cᵉᵏ                           :: FT
    Cᴹᴼ                           :: FT
    Cᶜᵛ                           :: FT
    boundary_layer_solar_fraction :: FT
    limit_hbl_stable              :: Bool

    # Velocity scales
    κᵥ                            :: FT
    ε                             :: FT
    Cˢᵗ                           :: FT

    # Monin-Obukhov universal functions — momentum
    Aᵐ                            :: FT
    Bᵐ                            :: FT
    Cᵐ                            :: FT
    ζᵐ                            :: FT

    # Monin-Obukhov universal functions — scalars
    Aˢ                            :: FT
    Bˢ                            :: FT
    Cˢ                            :: FT
    ζˢ                            :: FT

    # Interior shear instability
    Ri∞                           :: FT
    ν₀ˢʰ                          :: FT
    κ₀ˢʰ                          :: FT

    # Internal-wave background (MITgcm `viscArNr`, `diffusKzS`)
    νⁱʷ                           :: FT
    κⁱʷ                           :: FT

    # Interior convective instability
    N²ᶜᵒⁿ                         :: FT
    νᶜᵒⁿ                          :: FT
    κᶜᵒⁿ                          :: FT

    # Boundary-layer nonlocal transport
    C★                            :: FT

    # Numerical safeguards
    minimum_boundary_layer_depth  :: FT
    minimum_friction_velocity     :: FT
end

"""
    KPPParameters(FT = Float64; kwargs...)

Construct a `KPPParameters{FT}` with MITgcm defaults; override any field via keyword.
"""
function KPPParameters(FT::DataType = Float64;
                       Riᶜ                           = 0.3,
                       Cᵉᵏ                           = 0.7,
                       Cᴹᴼ                           = 1.0,
                       Cᶜᵛ                           = 1.8,
                       boundary_layer_solar_fraction = 1.0,
                       limit_hbl_stable              = true,

                       κᵥ                            = 0.4,
                       ε                             = 0.1,
                       Cˢᵗ                           = 5.0,

                       Aᵐ                            = 1.257,
                       Bᵐ                            = 8.380,
                       Cᵐ                            = 16.0,
                       ζᵐ                            = -0.2,

                       Aˢ                            = -28.86,
                       Bˢ                            = 98.96,
                       Cˢ                            = 16.0,
                       ζˢ                            = -1.0,

                       Ri∞                           = 0.7,
                       ν₀ˢʰ                          = 5e-3,
                       κ₀ˢʰ                          = 5e-3,

                       νⁱʷ                           = 5e-5,
                       κⁱʷ                           = 5e-6,

                       N²ᶜᵒⁿ                         = -0.2e-4,
                       νᶜᵒⁿ                          = 0.1,
                       κᶜᵒⁿ                          = 0.1,

                       C★                            = 10.0,

                       minimum_boundary_layer_depth  = 1.0,
                       minimum_friction_velocity     = 1e-6)

    return KPPParameters{FT}(convert(FT, Riᶜ),
                             convert(FT, Cᵉᵏ),
                             convert(FT, Cᴹᴼ),
                             convert(FT, Cᶜᵛ),
                             convert(FT, boundary_layer_solar_fraction),
                             limit_hbl_stable,

                             convert(FT, κᵥ),
                             convert(FT, ε),
                             convert(FT, Cˢᵗ),

                             convert(FT, Aᵐ),
                             convert(FT, Bᵐ),
                             convert(FT, Cᵐ),
                             convert(FT, ζᵐ),

                             convert(FT, Aˢ),
                             convert(FT, Bˢ),
                             convert(FT, Cˢ),
                             convert(FT, ζˢ),

                             convert(FT, Ri∞),
                             convert(FT, ν₀ˢʰ),
                             convert(FT, κ₀ˢʰ),

                             convert(FT, νⁱʷ),
                             convert(FT, κⁱʷ),

                             convert(FT, N²ᶜᵒⁿ),
                             convert(FT, νᶜᵒⁿ),
                             convert(FT, κᶜᵒⁿ),

                             convert(FT, C★),

                             convert(FT, minimum_boundary_layer_depth),
                             convert(FT, minimum_friction_velocity))
end
