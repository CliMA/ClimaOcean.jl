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

#####
##### A Data structure that holds flux statistics
#####

struct FluxStatistics{F}
    flux :: F
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

    return FluxStatistics(f, mean, meansq, std, max, min)
end

Adapt.adapt_structure(to, f::FluxStatistics) = FluxStatistics(Adapt.adapt(to, f.flux),
                                                              Adapt.adapt(to, f.mean),
                                                              Adapt.adapt(to, f.meansq),
                                                              Adapt.adapt(to, f.std),
                                                              Adapt.adapt(to, f.max),
                                                              Adapt.adapt(to, f.min))

on_architecture(arch, f::FluxStatistics) = FluxStatistics(on_architecture(arch, f.flux),
                                                          on_architecture(arch, f.mean),
                                                          on_architecture(arch, f.meansq),
                                                          on_architecture(arch, f.std),
                                                          on_architecture(arch, f.max),
                                                          on_architecture(arch, f.min))

function update_stats!(stats::FluxStatistics, iteration)
    grid = stats.flux.grid
    arch = architecture(grid)
    launch!(arch, grid, :xy, _update_stats!, stats, iteration)
    return nothing
end

function update_stats!(stats::NamedTuple, iteration)
    for stat in stats
        update_stats!(stat, iteration)
    end
    return nothing
end

@kernel function _update_stats!(stats, iteration)
    i, j = @index(Global, NTuple)

    inverse_iteration = 1 / (iteration + 1)

    # use iterative averaging via
    # mean_n = (x1 + ... + xn) / n ->
    # -> mean_{n+1} = mean_n * (1 - 1/(n+1)) * x_{n+1} / (n+1)
    @inbounds begin
        f = stats.flux[i, j, 1]
        stats.mean[i, j, 1]   *= 1 - inverse_iteration
        stats.mean[i, j, 1]   += f * inverse_iteration
        stats.meansq[i, j, 1] *= 1 - inverse_iteration
        stats.meansq[i, j, 1] += f^2 * inverse_iteration
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
    Qu = FluxStatistics(atmos_ocean_fluxes.upwelling_longwave)
    Qs = FluxStatistics(atmos_ocean_fluxes.downwelling_shortwave)
    Qℓ = FluxStatistics(atmos_ocean_fluxes.downwelling_longwave)

    stats = (; τx, τy, Jᵀ, Jˢ, Qc, Qv, Qu, Qs, Qℓ)

    function update_flux_stats!(earth)
        iteration = earth.model.clock.iteration
        update_stats!(stats, iteration)
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
##### A prescribed ocean...
#####

# ...with prescribed velocity and tracer fields
dataset = ECCO4Monthly()
arch    = CPU()

start_date = DateTime(1992, 1, 1)
end_date   = DateTime(1992, 1, 2)

stop_time = Day(end_date - start_date).value * Oceananigans.Units.days 

time_indices_in_memory = 2

u_meta = Metadata(:u_velocity;  start_date, end_date, dataset)
v_meta = Metadata(:v_velocity;  start_date, end_date, dataset)
T_meta = Metadata(:temperature; start_date, end_date, dataset)
S_meta = Metadata(:salinity;    start_date, end_date, dataset)

u = FieldTimeSeries(u_meta, arch; time_indices_in_memory)
v = FieldTimeSeries(v_meta, arch; time_indices_in_memory)
T = FieldTimeSeries(T_meta, arch; time_indices_in_memory)
S = FieldTimeSeries(S_meta, arch; time_indices_in_memory)

grid = u.grid

ocean_model = PrescribedOcean((; u, v, T, S); grid)
ocean = Simulation(ocean_model; Δt=3hours, stop_time)

set!(ocean_model.tracers.T,    first(T_meta))
set!(ocean_model.tracers.S,    first(S_meta))
set!(ocean_model.velocities.u, first(u_meta))
set!(ocean_model.velocities.v, first(v_meta))

#####
##### A prescribed atmosphere...
#####

atmosphere = ClimaOcean.ECCO.ECCOPrescribedAtmosphere(arch; start_date, end_date, time_indices_in_memory, dataset)

#####
##### A prescribed earth...
#####         

earth_model = OceanSeaIceModel(ocean, nothing; atmosphere, radiation = Radiation(arch))

Qtecco = Metadata(:net_heat_flux; start_date, end_date, dataset)
Qcecco = Metadata(:sensible_heat_flux; start_date, end_date, dataset)
Qvecco = Metadata(:latent_heat_flux; start_date, end_date, dataset)
Qℓecco = Metadata(:net_longwave; start_date, end_date, dataset)
Qsecco = Metadata(:downwelling_shortwave; start_date, end_date, dataset)

Qℓecco = FieldTimeSeries(Qℓecco, arch; time_indices_in_memory)
Qtecco = FieldTimeSeries(Qtecco, arch; time_indices_in_memory)
Qcecco = FieldTimeSeries(Qcecco, arch; time_indices_in_memory)
Qvecco = FieldTimeSeries(Qvecco, arch; time_indices_in_memory)
Qsecco = FieldTimeSeries(Qsecco, arch; time_indices_in_memory)

Qt = earth_model.interfaces.net_fluxes.ocean_surface.Q
Qc = earth_model.interfaces.atmosphere_ocean_interface.fluxes.sensible_heat
Qv = earth_model.interfaces.atmosphere_ocean_interface.fluxes.latent_heat
Qℓ = earth_model.interfaces.atmosphere_ocean_interface.fluxes.downwelling_longwave + 
     earth_model.interfaces.atmosphere_ocean_interface.fluxes.upwelling_longwave 
Qs = earth_model.interfaces.atmosphere_ocean_interface.fluxes.downwelling_shortwave

import Oceananigans.Fields: interior

interior(b::Oceananigans.AbstractOperations.BinaryOperation, idx...) = 
    interior(compute!(Field(b)), idx...)    

fig = Figure()
ax = Axis(fig[1, 1], title = "ECCO net heat flux")
heatmap!(ax, Qtecco[1], colorrange=(-150, 150), colormap=:balance)
ax = Axis(fig[1, 2], title = "CO net heat flux")
hm = heatmap!(ax, -Qt, colorrange=(-150, 150), colormap=:balance)
Colorbar(fig[1, 3], hm)
ax = Axis(fig[1, 4], title = "ECCO - CO")
hm = heatmap!(ax, interior(Qtecco[1], :, :, 1) .+ interior(Qt, :, :, 1), colorrange=(-20, 20), colormap=:balance)
Colorbar(fig[1, 5], hm)

ax = Axis(fig[2, 1], title = "ECCO sensible flux")
heatmap!(ax, Qcecco[1], colorrange=(-100, 100), colormap=:balance)
ax = Axis(fig[2, 2], title = "CO sensible flux")
hm = heatmap!(ax, -Qc, colorrange=(-100, 100), colormap=:balance)
Colorbar(fig[2, 3], hm)
ax = Axis(fig[2, 4], title = "ECCO - CO")
hm = heatmap!(ax, interior(Qcecco[1], :, :, 1) .+ interior(Qc, :, :, 1), colorrange=(-20, 20), colormap=:balance)
Colorbar(fig[2, 5], hm)

ax = Axis(fig[3, 1], title = "ECCO latent flux")
heatmap!(ax, Qvecco[1], colorrange=(-150, 150), colormap=:balance)
ax = Axis(fig[3, 2], title = "CO latent flux")
heatmap!(ax, -Qv, colorrange=(-150, 150), colormap=:balance)
hm = Colorbar(fig[3, 3], hm)
ax = Axis(fig[3, 4], title = "ECCO - CO")
hm = heatmap!(ax, interior(Qvecco[1], :, :, 1) .+ interior(Qv, :, :, 1), colorrange=(-20, 20), colormap=:balance)
Colorbar(fig[3, 5], hm)

ax = Axis(fig[4, 1], title = "ECCO longwave flux")
heatmap!(ax, Qℓecco[1], colorrange=(-100, 100), colormap=:balance)
ax = Axis(fig[4, 2], title = "CO longwave flux")
hm = heatmap!(ax, Qℓ, colorrange=(-100, 100), colormap=:balance)
Colorbar(fig[4, 3], hm)
ax = Axis(fig[4, 4], title = "ECCO - CO")
hm = heatmap!(ax, interior(Qℓecco[1], :, :, 1) .- interior(Qℓ, :, :, 1), colorrange=(-1, 1), colormap=:balance)
Colorbar(fig[4, 5], hm)
