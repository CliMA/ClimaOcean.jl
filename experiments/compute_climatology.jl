using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.ECCO: all_ECCO_dates
using Oceananigans
using Oceananigans.Models: AbstractModel, update_model_field_time_series!
using Oceananigans.Units

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
##### A prescribed atmosphere...
#####

atmosphere = JRA55PrescribedAtmosphere(CPU(); backend=JRA55NetCDFBackend(10))

#####
##### A prescribed earth...
#####

earth = OceanSeaIceModel(ocean, nothing; atmosphere, radiation = Radiation())