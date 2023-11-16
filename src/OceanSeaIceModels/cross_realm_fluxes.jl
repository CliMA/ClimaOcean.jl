struct RelativeAtmosphereOceanVelocity end
struct AtmosphereVelocity end

#####
##### Bulk formula
#####

struct BulkFormula{T, CD}
    transfer_velocity :: T
    transfer_coefficient :: CD
end

function BulkFormula(FT=Float64;
                     transfer_velocity = RelativeAtmosphereOceanVelocity(),
                     transfer_coefficient = convert(FT, 1e-3))

    return BulkFormula(transfer_velocity, transfer_coefficient)
end

#####
##### Abstraction for fluxes across the realms
#####

struct CrossRealmFlux{EQ, F}
    formula :: EQ
    flux :: F
end

"""
    CrossRealmFlux(flux_field; formula = nothing)

May the realms communicate.
"""
function CrossRealmFlux(flux_field; formula = nothing)

    if isnothing(formula) # constant coefficient then
        formula = BulkFormula(eltype(flux_field))
    end
                        
    return CrossRealmFlux(formula, flux_field)
end

