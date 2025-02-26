using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ
using Oceananigans.Forcings: MultipleForcings

struct XDirection end
struct YDirection end

struct BarotropicPressureForcing{D, P, FT}
    direction :: D
    pressure :: P
    reference_density :: FT
end

const XDirectionBPF = BarotropicPressureForcing{<:XDirection}
const YDirectionBPF = BarotropicPressureForcing{<:YDirection}

@inline (bpf::XDirectionBPF)(i, j, k, grid, clock, fields) = - ∂xᶠᶜᶜ(i, j, k, grid, bpf.pressure) / bpf.reference_density
@inline (bpf::YDirectionBPF)(i, j, k, grid, clock, fields) = - ∂yᶜᶠᶜ(i, j, k, grid, bpf.pressure) / bpf.reference_density

forcing_barotropic_pressure(something) = nothing
forcing_barotropic_pressure(f::BarotropicPressureForcing) = f.pressure.data

function forcing_barotropic_pressure(mf::MultipleForcings)
    n = findfirst(f -> f isa BarotropicPressureForcing, mf.forcing)
    if isnothing(n)
        return nothing
    else
        return forcing_barotropic_pressure(mf.forcing[n])
    end
end

forcing_barotropic_pressure(sim::Simulation) =
    forcing_barotropic_pressure(sim.model)

function forcing_barotropic_pressure(model::HydrostaticFreeSurfaceModel)
    u_pressure = forcing_barotropic_pressure(model.forcing.u)
    v_pressure = forcing_barotropic_pressure(model.forcing.v)
    @assert u_pressure === v_pressure
    return u_pressure
end
    
