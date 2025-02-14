using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.ECCO: all_ECCO_dates
using Oceananigans
using Oceananigans.Fields: ZeroField, location
using Oceananigans.Grids: architecture
using Oceananigans.Models: AbstractModel, update_model_field_time_series!
using Oceananigans.Units
using Adapt

using Printf
using KernelAbstractions: @index, @kernel

#####
##### A Data structure that holds flux statistics
#####

struct FluxStatistics{F}
    avg :: F
    std :: F
    max :: F
    min :: F
end

function FluxStatistics(f::Field)
    a = similar(f)
    s = similar(f)
    p = similar(f)
    m = similar(f)

    fill!(a, 0)
    fill!(s, 0)
    fill!(p, 0)
    fill!(m, 0)

    return FluxStatistics(a, s, p, m)
end

Adapt.adapt_structure(to, f::FluxStatistics) = FluxStatistics(Adapt.adapt(to, f.avg),
                                                              Adapt.adapt(to, f.std),
                                                              Adapt.adapt(to, f.max),
                                                              Adapt.adapt(to, f.min))

function update_stats!(stats::FluxStatistics, flux, Δt, stop_time)
    grid = flux.grid
    arch = architecture(grid)
    launch!(arch, grid, :xy, _update_stats!, stats, flux, Δt, stop_time)
end

@kernel function _update_stats!(stats, flux, Δt, stop_time)
    i, j = @index(Global, NTuple)

    @inbounds begin
         f = flux[i, j, 1]
        stats.avg[i, j, 1] += f * Δt / stop_time
        stats.std[i, j, 1] += f^2 * Δt / stop_time
        stats.max[i, j, 1] = max(stats.max[i, j, 1], f)
        stats.min[i, j, 1] = min(stats.min[i, j, 1], f)
    end
end

finalize_std!(f::FluxStatistics) = @. f.std = sqrt(f.std - f.avg^2)

#####
##### A function to compute flux statistics
#####

function compute_flux_climatology(earth)
    net_fluxes = coupled_model.interfaces.net_fluxes.ocean_surface
    τx = FluxStatistics(net_fluxes.u)
    τy = FluxStatistics(net_fluxes.v)
    Jᵀ = FluxStatistics(net_fluxes.T)
    Jˢ = FluxStatistics(net_fluxes.S)

    stats = (; τx, τy, Jᵀ, Jˢ)

    function update_flux_stats!(earth)
        Δt = earth.Δt
        stop_time = earth.stop_time

        update_stats!(τx, net_fluxes.u, Δt, stop_time)
        update_stats!(τy, net_fluxes.v, Δt, stop_time)
        update_stats!(Jᵀ, net_fluxes.T, Δt, stop_time)
        update_stats!(Jˢ, net_fluxes.S, Δt, stop_time)

        return nothing
    end

    add_callback!(earth, update_flux_stats!, IterationInterval(1))

    run!(earth)

    finalize_std!(τx)
    finalize_std!(τy)
    finalize_std!(Jᵀ)
    finalize_std!(Jˢ)

    return stats
end

#####
##### A prescribed ocean...
#####

struct PrescribedOcean{A, G, C, U, T} <: AbstractModel{Nothing}
    architecture :: A       
    grid :: G        
    clock :: Clock{C}
    velocities :: U
    tracers :: T
end

PrescribedOcean(; grid, clock= Clock{Float64}(time = 0), velocities, tracers) = 
    PrescribedOcean(architecture(grid), grid, clock, velocities, tracers)

import Oceananigans.TimeSteppers: time_step!

function time_step!(model::PrescribedOcean, Δt)
    tick!(model.clock, Δt)
    update_model_field_time_series!(model, model.clock)
end

# ...with prescribed velocity and tracer fields
version = ECCO4Monthly()
dates   = all_ECCO_dates(version)[1:24]

u = ECCOFieldTimeSeries(:u_velocity,  version; dates)
v = ECCOFieldTimeSeries(:v_velocity,  version; dates)
T = ECCOFieldTimeSeries(:temperature, version; dates)
S = ECCOFieldTimeSeries(:salinity,    version; dates)

grid = u.grid

ocean_model = PrescribedOcean(; grid, velocities=(; u, v, w=ZeroField()), tracers=(; T, S))
ocean = Simulation(ocean_model, Δt=1hour, stop_time=365days)

#####
##### Need to extend a couple of methods
#####

import ClimaOcean.OceanSeaIceModels: reference_density, heat_capacity

reference_density(ocean::Simulation{<:PrescribedOcean}) = 1025.6
heat_capacity(ocean::Simulation{<:PrescribedOcean}) = 3995.6

import ClimaOcean.OceanSeaIceModels.InterfaceComputations: surface_flux

surface_flux(fts::FieldTimeSeries) = Field{location(fts)[1:2]..., Nothing}(fts.grid)

#####
##### A prescribed atmosphere...
#####

atmosphere = JRA55PrescribedAtmosphere(CPU(); backend=JRA55NetCDFBackend(10))

#####
##### A prescribed earth...
#####

earth_model = OceanSeaIceModel(ocean, nothing; atmosphere, radiation = Radiation())
earth = Simulation(earth_model, Δt=1hour, stop_time=365days)

# Just checking that the ocean evolves as expected
function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T

    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))

    umax = (maximum(abs, interior(u)),
            maximum(abs, interior(v)),
            maximum(abs, interior(w)))

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg = @sprintf("Iter: %d, time: %s, Δt: %s", iteration(sim), prettytime(sim), prettytime(sim.Δt))
    msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e) m s⁻¹, extrema(T): (%.2f, %.2f) ᵒC, wall time: %s",
                    umax..., Tmax, Tmin, prettytime(step_time))

    @info msg 

    wall_time[] = time_ns()
end

add_callback!(earth, progress, IterationInterval(100))

stats = compute_flux_climatology(earth)


