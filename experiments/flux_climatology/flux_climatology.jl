using ClimaOcean
using ClimaOcean.ECCO
using Dates
using Oceananigans
using Oceananigans.Utils
using Oceananigans.Fields: ZeroField, location, VelocityFields
using Oceananigans.Grids: architecture
using Oceananigans.Models: AbstractModel
using Oceananigans.OutputReaders: extract_field_time_series, update_field_time_series!
using Oceananigans.Units
using Adapt

using Printf
using KernelAbstractions: @index, @kernel

import Oceananigans.TimeSteppers: time_step!, update_state!, reset!, tick!
import Oceananigans.Models: timestepper, update_model_field_time_series!

import ClimaOcean.OceanSeaIceModels: reference_density, heat_capacity
import Oceananigans.Architectures: on_architecture

ClimaOcean.DataWrangling.dataset_defaults.FloatType = Float64

#####
##### A Data structure that holds flux statistics
#####

struct FluxStatistics{F}
    mean :: F
    meansq :: F
    std :: F
    max :: F
    min :: F
end

function FluxStatistics(f::Field)
    mean = similar(f)
    meansq = similar(f)
    std = similar(f)
    max = similar(f)
    min = similar(f)

    fill!(mean, 0)
    fill!(meansq, 0)
    fill!(std, 0)
    fill!(max, 0)
    fill!(min, 0)

    return FluxStatistics(mean, meansq, std, max, min)
end

Adapt.adapt_structure(to, f::FluxStatistics) = FluxStatistics(Adapt.adapt(to, f.mean),
                                                              Adapt.adapt(to, f.meansq),
                                                              Adapt.adapt(to, f.std),
                                                              Adapt.adapt(to, f.max),
                                                              Adapt.adapt(to, f.min))

on_architecture(arch, f::FluxStatistics) = FluxStatistics(on_architecture(arch, f.mean),
                                                          on_architecture(arch, f.meansq),
                                                          on_architecture(arch, f.std),
                                                          on_architecture(arch, f.max),
                                                          on_architecture(arch, f.min))

function update_stats!(stats::FluxStatistics, flux, iteration)
    grid = flux.grid
    arch = architecture(grid)
    launch!(arch, grid, :xy, _update_stats!, stats, flux, iteration)
    return nothing
end

@kernel function _update_stats!(stats, flux, iteration)
    i, j = @index(Global, NTuple)

    inverse_iteration = 1 / (iteration + 1)

    # use iterative averaging via
    # mean_n = (x1 + ... + xn) / n ->
    # -> mean_{n+1} = mean_n * (1 - 1/(n+1)) * x_{n+1} / (n+1)
    @inbounds begin
        f = flux[i, j, 1]
        stats.mean[i, j, 1] *= 1 - inverse_iteration
        stats.mean[i, j, 1] += f * inverse_iteration
        stats.meansq[i, j, 1] *= 1 - inverse_iteration
        stats.meansq[i, j, 1] += f^2 * inverse_iteration
        stats.std[i, j, 1] = sqrt(stats.meansq[i, j, 1] - stats.mean[i, j, 1]^2)
        stats.max[i, j, 1] = max(stats.max[i, j, 1], f)
        stats.min[i, j, 1] = min(stats.min[i, j, 1], f)
    end
end

#####
##### A function to compute flux statistics
#####

function compute_flux_climatology(earth)
    net_fluxes = earth.model.interfaces.net_fluxes.ocean_surface
    τx = FluxStatistics(net_fluxes.u)
    τy = FluxStatistics(net_fluxes.v)
    Jᵀ = FluxStatistics(net_fluxes.T)
    Jˢ = FluxStatistics(net_fluxes.S)

    atmos_ocean_fluxes = earth.model.interfaces.atmosphere_ocean_interface.fluxes
    Qc = FluxStatistics(atmos_ocean_fluxes.sensible_heat)
    Qv = FluxStatistics(atmos_ocean_fluxes.latent_heat)

    stats = (; τx, τy, Jᵀ, Jˢ, Qc, Qv)

    function update_flux_stats!(earth)
        iteration = earth.model.clock.iteration
        update_stats!(τx, net_fluxes.u, iteration)
        update_stats!(τy, net_fluxes.v, iteration)
        update_stats!(Jᵀ, net_fluxes.T, iteration)
        update_stats!(Jˢ, net_fluxes.S, iteration)
        update_stats!(Qc, atmos_ocean_fluxes.sensible_heat, iteration)
        update_stats!(Qv, atmos_ocean_fluxes.latent_heat, iteration)

        return nothing
    end

    add_callback!(earth, update_flux_stats!, IterationInterval(1))

    run!(earth)

    return stats
end

#####
##### A prescribed ocean...
#####

struct PrescribedOcean{A, G, C, U, T, F} <: AbstractModel{Nothing, A}
    architecture :: A
    grid :: G
    clock :: Clock{C}
    velocities :: U
    tracers :: T
    timeseries :: F
end

function PrescribedOcean(timeseries;
                         grid,
                         clock=Clock{Float64}(time = 0))

    τx = Field{Face, Center, Nothing}(grid)
    τy = Field{Center, Face, Nothing}(grid)
    Jᵀ = Field{Center, Center, Nothing}(grid)
    Jˢ = Field{Center, Center, Nothing}(grid)

    u = XFaceField(grid,  boundary_conditions=FieldBoundaryConditions(grid, (Face,   Center, Center), top = FluxBoundaryCondition(τx)))
    v = YFaceField(grid,  boundary_conditions=FieldBoundaryConditions(grid, (Center, Face,   Center), top = FluxBoundaryCondition(τy)))
    T = CenterField(grid, boundary_conditions=FieldBoundaryConditions(grid, (Center, Center, Center), top = FluxBoundaryCondition(Jᵀ)))
    S = CenterField(grid, boundary_conditions=FieldBoundaryConditions(grid, (Center, Center, Center), top = FluxBoundaryCondition(Jˢ)))

    PrescribedOcean(architecture(grid), grid, clock, (; u, v, w=ZeroField()), (; T, S), timeseries)
end

# ...with prescribed velocity and tracer fields
dataset = ECCO4Monthly()
arch    = CPU()

start_date = DateTime(1992, 1, 1)
end_date   = DateTime(1992, 12, 31)

stop_time = Day(end_date - start_date).value * Oceananigans.Units.days

time_indices_in_memory = 13

u_meta = Metadata(:u_velocity;  start_date, end_date, dataset)
v_meta = Metadata(:v_velocity;  start_date, end_date, dataset)
T_meta = Metadata(:temperature; start_date, end_date, dataset)
S_meta = Metadata(:salinity;    start_date, end_date, dataset)

u = FieldTimeSeries(u_meta, arch; time_indices_in_memory)
v = FieldTimeSeries(v_meta, arch; time_indices_in_memory)
T = FieldTimeSeries(T_meta, arch; time_indices_in_memory)
S = FieldTimeSeries(S_meta, arch; time_indices_in_memory)

grid = ECCO.ECCO_immersed_grid(arch)

ocean_model = PrescribedOcean((; u, v, T, S); grid)
ocean = Simulation(ocean_model; Δt=3hours, stop_time)

#####
##### Need to extend a couple of methods
#####

function time_step!(model::PrescribedOcean, Δt; callbacks=[], euler=true)
    tick!(model.clock, Δt)
    time = Oceananigans.Units.Time(model.clock.time)

    possible_fts = merge(model.velocities, model.tracers)
    time_series_tuple = extract_field_time_series(possible_fts)

    for fts in time_series_tuple
        update_field_time_series!(fts, time)
    end

    parent(model.velocities.u) .= parent(model.timeseries.u[time])
    parent(model.velocities.v) .= parent(model.timeseries.v[time])
    parent(model.tracers.T)    .= parent(model.timeseries.T[time])
    parent(model.tracers.S)    .= parent(model.timeseries.S[time])

    return nothing
end

update_state!(::PrescribedOcean) = nothing
timestepper(::PrescribedOcean) = nothing

reference_density(ocean::Simulation{<:PrescribedOcean}) = 1025.6
heat_capacity(ocean::Simulation{<:PrescribedOcean}) = 3995.6

#####
##### A prescribed atmosphere...
#####

atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(1000))

#####
##### A prescribed earth...
#####

earth_model = OceanSeaIceModel(ocean, nothing; atmosphere, radiation = Radiation(arch))
earth = Simulation(earth_model; Δt=3hours, stop_time)

wall_time = Ref(time_ns())

# Just checking that the ocean evolves as expected
function progress(sim)
    net_fluxes = sim.model.interfaces.net_fluxes.ocean_surface
    atmos_ocean_fluxes = sim.model.interfaces.atmosphere_ocean_interface.fluxes

    Q = net_fluxes.Q
    Qc = atmos_ocean_fluxes.sensible_heat
    Qv = atmos_ocean_fluxes.latent_heat

    τx = net_fluxes.u
    τy = net_fluxes.v

    Qmax = maximum(Q)
    Qmin = minimum(Q)
    Qcmax = maximum(Qc)
    Qcmin = minimum(Qc)
    Qvmax = maximum(Qv)
    Qvmin = minimum(Qv)
    τmax = (maximum(abs, τx), maximum(abs, τy))

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg = @sprintf("Iter: %d, time: %s", iteration(sim), prettytime(sim))
    msg *= @sprintf(", max|τ|: (%.2e, %.2e) m² s⁻², extrema(Q): (%.2f, %.2f) W m⁻², extrema(Qc): (%.2f, %.2f) W m⁻², extrema(Qv): (%.2f, %.2f) W m⁻², wall time: %s",
                    τmax..., Qmin, Qmax, Qcmin, Qcmax, Qvmin, Qvmax, prettytime(step_time))

    @info msg

    wall_time[] = time_ns()
end

add_callback!(earth, progress, IterationInterval(16))

stats = compute_flux_climatology(earth)
