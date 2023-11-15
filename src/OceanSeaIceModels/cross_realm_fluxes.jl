struct CrossRealmFluxes{EQ, F}
    formula :: EQ
    fluxes :: F
end

struct RelativeSpeed end
struct AtmosphereOnlySpeed end

struct FluxFormula{T, CD}
    transfer_velocity_scale :: T
    transfer_coefficient :: CD
end

function AtmosphereOceanMomentumFlux(; transfer_velocity_scale = RelativeSpeed(),
                                     drag_coefficient = 1e-3)

    return AtmosphereOceanMomentumFlux(transfer_velocity_scale,
                                       drag_coefficient)
end



