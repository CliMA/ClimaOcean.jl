# NEMO 3.6 zdftke namelist parameters with OMIP-2 ORCA1 defaults.
# NEMO names appear as trailing comments for traceability.

"""
    NEMOTKEParameters{FT}

Parameters of the NEMO 3.6 TKE vertical-mixing closure
(Blanke & Delecluse 1993; Gaspar et al. 1990; Madec et al. 2017).
Defaults reproduce the OMIP-2 ORCA1 NEMO preset (Tsujino 2020, Iovino 2023).
"""
struct NEMOTKEParameters{FT}
    Cᴷ                                 :: FT     # rn_ediff  Kₘ = Cᴷ·ℓ·√e
    Cᴰ                                 :: FT     # rn_ediss  Kolmogorov dissipation
    Cᵇ                                 :: FT     # rn_ebb    surface-BC: e = Cᵇ·u★²
    Cᴸ                                 :: FT     # rn_lc     Langmuir coeff (Axell 2002)
    Cᶠ                                 :: FT     # rn_efr    wave-penetration fraction
    Cˢ                                 :: FT     # us = Cˢ · √|τ| (Stokes-drift proxy)
    κᶜⁿᵛ                               :: FT     # rn_avevd  EVD overwrite κ when N²<0
    νᵇ                                 :: FT     # avmb (rn_avm0) background ν
    κᵇ                                 :: FT     # avtb (rn_avt0) background κ

    minimum_TKE                        :: FT     # rn_emin
    minimum_surface_TKE                :: FT     # rn_emin0
    minimum_mixing_length              :: FT     # rn_mxl0

    mixing_length_formulation          :: Int    # nn_mxl   2 = gradient-limited
    wave_penetration_formulation       :: Int    # nn_etau  1 = exp below MLD
    surface_length_scale_formulation   :: Int    # nn_htau  1 = lat-dep

    apply_langmuir_circulation         :: Bool   # ln_lc
    apply_wave_penetration             :: Bool   # nn_etau > 0
    apply_enhanced_vertical_diffusion  :: Bool   # ln_zdfevd
    apply_evd_to_momentum              :: Bool   # nn_evdm == 1
    apply_prandtl_richardson           :: Bool   # nn_pdl == 1
end

"""
    NEMOTKEParameters(FT = Float64; kwargs...)

Construct `NEMOTKEParameters{FT}` with NEMO 3.6 ORCA1 OMIP-2 defaults; override any field via keyword.
"""
function NEMOTKEParameters(FT::DataType = Float64;
                           Cᴷ                                = 0.1,
                           Cᴰ                                = 0.7,
                           Cᵇ                                = 3.75,
                           Cᴸ                                = 0.15,
                           Cᶠ                                = 1.0,
                           Cˢ                                = 0.016,
                           κᶜⁿᵛ                              = 100.0,
                           νᵇ                                = 1.2e-4,
                           κᵇ                                = 1.2e-5,

                           minimum_TKE                       = sqrt(2) * 1e-6,
                           minimum_surface_TKE               = 1e-4,
                           minimum_mixing_length             = 0.04,

                           mixing_length_formulation         = 2,
                           wave_penetration_formulation      = 1,
                           surface_length_scale_formulation  = 1,

                           apply_langmuir_circulation        = true,
                           apply_wave_penetration            = true,
                           apply_enhanced_vertical_diffusion = true,
                           apply_evd_to_momentum             = true,
                           apply_prandtl_richardson          = false)

    return NEMOTKEParameters{FT}(convert(FT, Cᴷ),
                                 convert(FT, Cᴰ),
                                 convert(FT, Cᵇ),
                                 convert(FT, Cᴸ),
                                 convert(FT, Cᶠ),
                                 convert(FT, Cˢ),
                                 convert(FT, κᶜⁿᵛ),
                                 convert(FT, νᵇ),
                                 convert(FT, κᵇ),

                                 convert(FT, minimum_TKE),
                                 convert(FT, minimum_surface_TKE),
                                 convert(FT, minimum_mixing_length),

                                 mixing_length_formulation,
                                 wave_penetration_formulation,
                                 surface_length_scale_formulation,

                                 apply_langmuir_circulation,
                                 apply_wave_penetration,
                                 apply_enhanced_vertical_diffusion,
                                 apply_evd_to_momentum,
                                 apply_prandtl_richardson)
end
