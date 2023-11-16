struct Radiation{FT, SW, LW, OE, IE, OA, IA}
    downwelling_shortwave_radiation :: SW
    downwelling_longwave_radiation :: LW
    ocean_emissivity :: OE
    ice_emissivity :: IE
    ocean_albedo :: OA
    ice_albedo :: IA
    stefan_boltzmann_constant :: FT
    reference_temperature :: FT
end

function Radiation(FT=Float64;
                   downwelling_shortwave_radiation = nothing,
                   downwelling_longwave_radiation = nothing,
                   ocean_emissivity = 0.97,
                   ice_emissivity = 1.0,
                   ocean_albedo = 0.3,
                   ice_albedo = 0.7,
                   stefan_boltzmann_constant = 5.67e-8)

    if downwelling_shortwave_radiation isa AbstractArray
        # Replace FT with appropriate eltype
        FT = eltype(downwelling_shortwave_radiation)
    elseif downwelling_longwave_radiation isa AbstractArray
        FT = eltype(downwelling_longwave_radiation)
    end

    ocean_emissivity isa Number && (ocean_emissivity = convert(FT, ocean_emissivity))
    ice_emissivity   isa Number && (ice_emissivity   = convert(FT, ice_emissivity))
    ocean_albedo     isa Number && (ocean_albedo     = convert(FT, ocean_albedo))
    ice_albedo       isa Number && (ice_albedo       = convert(FT, ice_albedo))

    return Radiation(downwelling_shortwave_radiation,
                     downwelling_longwave_radiation,
                     ocean_emissivity,
                     ice_emissivity,
                     ocean_albedo,
                     ice_albedo,
                     stefan_boltzmann_constant)
end

