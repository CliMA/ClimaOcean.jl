using CondaPkg

using ClimaOcean.OceanSeaIceModels: reference_density, heat_capacity, SeaIceSimulation

import Oceananigans.Fields: set!
import Oceananigans.TimeSteppers: time_step!

import ClimaOcean.OceanSeaIceModels: OceanSeaIceModel

"""
    install_veros()

Install the Copernicus Marine CLI using CondaPkg.
Returns a NamedTuple containing package information if successful.
"""
function install_veros()
    CondaPkg.add("veros"; channel = "conda-forge")
    cli = CondaPkg.which("veros")
    @info "... the veros CLI has been installed at $(cli)."
    return cli
end

struct VerosOceanSimulation{S}
    setup :: S
end

time_step!(sim::VerosOceanSimulation, Δt) = sim.setup.step()

function remove_outputs(setup::Symbol)
    rm("$(setup).averages.nc", force=true)
    rm("$(setup).energy.nc", force=true)
    rm("$(setup).overturning.nc", force=true)
    rm("$(setup).snapshot.nc", force=true)
    return nothing
end

const Field2D = Field{<:Any, <:Any, <:Nothing}

function set!(field::Field2D, pyarray::Py, k=pyconvert(Int, pyarray.shape[2]))
    array = PyArray(pyarray)
    Nx, Ny, Nz = size(array)
    set!(field, view(array, 3:Nx-2, 3:Ny-2, k, 1))
    return field
end

function veros_ocean_simulation(setup, setup_name)
    setups = pyimport("veros.setups." * setup)
    setup  = @eval $setups.$setup_name()

    # instantiate the setup
    setup.setup()

    return VerosOceanSimulation(setup) 
end

function surface_grid(setup::VerosOceanSimulation)

    xf = Array(PyArray(setup.state.variables.xu))
    yf = Array(PyArray(setup.state.variables.yu))
    
    xc = Array(PyArray(setup.state.variables.xt))
    yc = Array(PyArray(setup.state.variables.yt))
    
    xf = xf[2:end-2]
    yf = yf[2:end-2]

    xc = xc[3:end-2]
    yc = yc[3:end-2]

    xf[1] = xf[2] - 2xc[1]
    yf[1] = sign(yf[2]) * (yf[2] - 2yc[1])

    TX = if xf[1] == 0 && xf[end] == 360
        Periodic
    else
        Bounded
    end

    Nx = length(xc) 
    Ny = length(yc) 

    return LatitudeLongitudeGrid(size=(Nx, Ny), longitude=xf, latitude=yf, topology=(TX, Bounded, Flat))
end

function veros_set!(ocean::VerosOceanSimulation, v, x)
    s = ocean.setup
    pyexec("""
       with setup.state.variables.unlock():
           setup.state.variables.__setattr__(y, t)
       """, Main, (y=v, t=x, setup=s))
end

function veros_settings_set!(ocean::VerosOceanSimulation, v, x)
    s = ocean.setup
    pyexec("""
       with setup.state.settings.unlock():
           setup.state.settings.__setattr__(y, t)
       """, Main, (y=v, t=x, setup=s))
end

function OceanSeaIceModel(ocean::VerosOceanSimulation, sea_ice=nothing;
                          atmosphere = nothing,
                          radiation = Radiation(),
                          clock = Clock(time=0),
                          ocean_reference_density = 1020.0,
                          ocean_heat_capacity = 3998.0, 
                          sea_ice_reference_density = reference_density(sea_ice),
                          sea_ice_heat_capacity = heat_capacity(sea_ice),
                          interfaces = nothing)

    if sea_ice isa SeaIceSimulation
        if !isnothing(sea_ice.callbacks)
            pop!(sea_ice.callbacks, :stop_time_exceeded, nothing)
            pop!(sea_ice.callbacks, :stop_iteration_exceeded, nothing)
            pop!(sea_ice.callbacks, :wall_time_limit_exceeded, nothing)
            pop!(sea_ice.callbacks, :nan_checker, nothing)
        end
    end

    # Contains information about flux contributions: bulk formula, prescribed fluxes, etc.
    if isnothing(interfaces) && !(isnothing(atmosphere) && isnothing(sea_ice))
        interfaces = ComponentInterfaces(atmosphere, ocean, sea_ice;
                                         ocean_reference_density,
                                         ocean_heat_capacity,
                                         sea_ice_reference_density,
                                         sea_ice_heat_capacity,
                                         radiation)
    end

    arch = CPU()

    ocean_sea_ice_model = OceanSeaIceModel(arch,
                                           clock,
                                           atmosphere,
                                           sea_ice,
                                           ocean,
                                           interfaces)

    # Make sure the initial temperature of the ocean
    # is not below freezing and above melting near the surface
    initialization_update_state!(ocean_sea_ice_model)

    return ocean_sea_ice_model
end