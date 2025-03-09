using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ
using Oceananigans.Forcings: MultipleForcings

using Adapt

struct XDirection end
struct YDirection end

struct BarotropicPotentialForcing{D, P}
    direction :: D
    potential :: P
end

Adapt.adapt_structure(to, bpf::BarotropicPotentialForcing) =
    BarotropicPotentialForcing(adapt(to, bpf.direction),
                               adapt(to, bpf.potential))

const XDirectionBPF = BarotropicPotentialForcing{<:XDirection}
const YDirectionBPF = BarotropicPotentialForcing{<:YDirection}

@inline (bpf::XDirectionBPF)(i, j, k, grid, clock, fields) = - ∂xᶠᶜᶜ(i, j, k, grid, bpf.potential)
@inline (bpf::YDirectionBPF)(i, j, k, grid, clock, fields) = - ∂yᶜᶠᶜ(i, j, k, grid, bpf.potential)

forcing_barotropic_potential(something) = nothing
forcing_barotropic_potential(f::BarotropicPotentialForcing) = f.potential.data

function forcing_barotropic_potential(tf::Tuple)
    n = findfirst(f -> f isa BarotropicPotentialForcing, tf)
    if isnothing(n)
        return nothing
    else
        return forcing_barotropic_potential(tf[n])
    end
end

function forcing_barotropic_potential(mf::MultipleForcings)
    n = findfirst(f -> f isa BarotropicPotentialForcing, mf.forcing)
    if isnothing(n)
        return nothing
    else
        return forcing_barotropic_potential(mf.forcing[n])
    end
end

forcing_barotropic_potential(sim::Simulation) = forcing_barotropic_potential(sim.model)

function forcing_barotropic_potential(model::HydrostaticFreeSurfaceModel)
    u_potential = forcing_barotropic_potential(model.forcing.u)
    v_potential = forcing_barotropic_potential(model.forcing.v)
    @assert u_potential === v_potential
    return u_potential
end
    
