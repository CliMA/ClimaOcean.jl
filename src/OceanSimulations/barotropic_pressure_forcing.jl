using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ
using Oceananigans.Forcings: MultipleForcings

struct XDirection end
struct YDirection end

struct BarotropicPressureForcing{D, P}
    direction :: D
    pressure :: P
end

const XDirectionBPF = BarotropicPressureForcing{<:XDirection}
const YDirectionBPF = BarotropicPressureForcing{<:YDirection}

@inline (bpf::XDirectionBPF)(i, j, k, grid, clock, fields) = - ∂xᶠᶜᶜ(i, j, k, grid, bpf.pressure)
@inline (bpf::YDirectionBPF)(i, j, k, grid, clock, fields) = - ∂yᶜᶠᶜ(i, j, k, grid, bpf.pressure)

forcing_barotropic_pressure(something) = nothing
forcing_barotropic_pressure(f::BarotropicPressureForcing) = f.pressure.data

function forcing_barotropic_pressure(mf::MultipleForcings)
    n = findfirst(f -> f isa BarotropicPressureForcing, mf.forcings)
    if isnothing(n)
        return nothing
    else
        return forcing_barotropic_pressure(mf.forcings[n])
    end
end

forcing_barotropic_pressure(sim::Simulation) =
    forcing_barotropic_pressure(sim.model)

function forcing_barotropic_pressure(model::HydrostaticFreeSurfaceModel)
    u_pressure = forcing_barotropic_pressure(model.forcings.v)
    v_pressure = forcing_barotropic_pressure(model.forcings.u)
    @assert u_pressure === v_pressure
    return u_pressure
end
    
