module ClimaOcean

export one_degree_tripolar_ocean, 
       half_degree_tripolar_ocean, 
       latitude_longitude_ocean, 
       sixth_degree_tripolar_ocean,
       Progress

using Reexport
using Oceananigans
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ

@reexport using NumericalEarth
@reexport using NumericalEarth.DataWrangling
@reexport using NumericalEarth.EarthSystemModels
@reexport using NumericalEarth.EarthSystemModels.InterfaceComputations

#####
##### Source code
#####

struct Progress{W} <: Function
    wall_time :: W
end

Progress() = Progress(Ref(time_ns()))

function (p::Progress)(sim)
    
    sea_ice = sim.model.sea_ice
    ocean   = sim.model.ocean

    msg1 = @sprintf("time: %s, iteration: %d, Δt: %s, ", prettytime(sim), iteration(sim), prettytime(sim.Δt))

    if sea_ice isa Simulation
        hmax = maximum(sea_ice.model.ice_thickness)
        ℵmax = maximum(sea_ice.model.ice_concentration)
        msg2 = @sprintf("max(h): %.2e m, max(ℵ): %.2e ", hmax, ℵmax)
    else
        msg2 = @sprintf("")
    end

    Tmax = maximum(ocean.model.tracers.T)
    Tmin = minimum(ocean.model.tracers.T)
    Smax = maximum(ocean.model.tracers.S)
    Smin = minimum(ocean.model.tracers.S)
    umax = maximum(ocean.model.velocities.u)
    vmax = maximum(ocean.model.velocities.v)
    wmax = maximum(ocean.model.velocities.w)

    step_time = 1e-9 * (time_ns() - p.wall_time[])

    msg4 = @sprintf("extrema(T, S): (%.2f, %.2f) ᵒC, (%.2f, %.2f) psu, ", Tmax, Tmin)
    msg5 = @sprintf("maximum(u): (%.2e, %.2e, %.2e) m/s, ", umax, vmax, wmax)
    msg6 = @sprintf("wall time: %s \n", prettytime(step_time))

    @info msg1 * msg2 * msg4 * msg5 * msg6

    p.wall_time[] = time_ns()

    return nothing
end

include("InitialConditions/InitialConditions.jl")
include("Diagnostics/Diagnostics.jl")
include("OceanConfigurations/OceanConfigurations.jl")
include("SeaIceConfigurations/SeaIceConfigurations.jl")

@reexport using NumericalEarth.DataWrangling: ETOPO, ECCO, GLORYS, EN4, JRA55
@reexport using NumericalEarth.Bathymetry
@reexport using NumericalEarth.EarthSystemModels
@reexport using NumericalEarth.Atmospheres
@reexport using NumericalEarth.Oceans
@reexport using NumericalEarth.SeaIces
using .InitialConditions
using .OceanConfigurations
using .SeaIceConfigurations

end # module
