struct PrescribedOcean{G, C, T, V} <: ClimaOcean.OceanSeaIceModels.AbstractPrescribedComponent
    grid :: G
    clock :: C
    tracers :: T
    velocities :: V
end

function default_ocean_velocities(grid, times)
    uo = FieldTimeSeries{Center, Center, Nothing}(grid, times)
    vo = FieldTimeSeries{Center, Center, Nothing}(grid, times)
    return (u=uo, v=vo)
end

function default_ocean_tracers(grid, times)
    To = FieldTimeSeries{Center, Center, Nothing}(grid, times)
    So = FieldTimeSeries{Center, Center, Nothing}(grid, times)
    parent(So) .= 35
    return (T=To, S=So)
end

function PrescribedOcean(grid, times=[zero(grid)];
                         clock       = Clock{Float64}(time = 0),
                         velocities  = default_ocean_velocities(grid, times),
                         tracers     = default_ocean_tracers(grid, times))

    return PrescribedOcean(grid, clock, velocities, tracers)
end

function ComponentExchanger(ocean::PrescribedOcean, grid) 

    
    state = (; u  = Field{Center, Center, Nothing}(grid),
               v  = Field{Center, Center, Nothing}(grid),
               T  = Field{Center, Center, Nothing}(grid),
               S  = Field{Center, Center, Nothing}(grid),
               q  = Field{Center, Center, Nothing}(grid),
               Qs = Field{Center, Center, Nothing}(grid),
               Qℓ = Field{Center, Center, Nothing}(grid),
               Mp = Field{Center, Center, Nothing}(grid))

    return ComponentExchanger(state, regridder)
end

#####
##### Extending ClimaOcean interface
#####

reference_density(ocean::PrescribedOcean) = convert(eltype(ocean.grid), 1025)
heat_capacity(ocean::PrescribedOcean) = convert(eltype(ocean.grid), 3995)

#####
##### Extend utility functions to grab the state of the ocean
#####

function ocean_surface_salinity(ocean::PrescribedOcean)
    kᴺ = size(ocean.model.grid, 3)
    time = Time(ocean.clock.time)
    return interior(ocean.tracers.S[time], :, :, kᴺ:kᴺ)
end

function ocean_surface_velocities(ocean::PrescribedOcean)
    kᴺ = size(ocean.model.grid, 3)
    return view(ocean.model.velocities.u, :, :, kᴺ), view(ocean.model.velocities.v, :, :, kᴺ)
end

# When using an Oceananigans simulation, we assume that the exchange grid is the ocean grid
# We need, however, to interpolate the surface pressure to the ocean grid
interpolate_state!(exchanger, grid, ::PrescribedOcean, coupled_model) = nothing

function ComponentExchanger(ocean::PrescribedOcean, grid) 
    ocean_grid = ocean.grid
    
    if ocean_grid == grid
        kᴺ = grid.Nz
        u = view(ocean.model.velocities.u, :, :, kᴺ:kᴺ)
        v = view(ocean.model.velocities.v, :, :, kᴺ:kᴺ)
        T = view(ocean.model.tracers.T,    :, :, kᴺ:kᴺ)
        S = view(ocean.model.tracers.S,    :, :, kᴺ:kᴺ)
    else
        u = Field{Center, Center, Nothing}(grid)
        v = Field{Center, Center, Nothing}(grid)
        T = Field{Center, Center, Nothing}(grid)
        S = Field{Center, Center, Nothing}(grid)
    end

    return ComponentExchanger((; u, v, T, S), nothing)
end
