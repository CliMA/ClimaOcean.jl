struct PrescribedOcean{G, C, T, V} <: ClimaOcean.OceanSeaIceModels.PrescribedComponent
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

function ComponentExchanger(atmosphere::PrescribedOcean, grid) 

    
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
