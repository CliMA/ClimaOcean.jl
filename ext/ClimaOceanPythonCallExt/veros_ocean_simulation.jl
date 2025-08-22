using CondaPkg

using Oceananigans.Grids: topology
using ClimaOcean.OceanSeaIceModels: reference_density, heat_capacity, SeaIceSimulation

import Oceananigans.Fields: set!
import Oceananigans.TimeSteppers: time_step!, initialize!

import ClimaOcean.OceanSeaIceModels: OceanSeaIceModel, default_nan_checker
import Oceananigans.Architectures: architecture

import Base: eltype

"""
    install_veros()

Install the Veros ocean model Marine CLI using CondaPkg.
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

default_nan_checker(model::OceanSeaIceModel{<:Any, <:Any, <:VerosOceanSimulation}) = nothing

initialize!(::ClimaOceanPythonCallExt.VerosOceanSimulation{Py}) = nothing
time_step!(ocean::VerosOceanSimulation, Δt) = ocean.setup.step(ocean.setup.state)
architecture(model::OceanSeaIceModel{<:Any, <:Any, <:VerosOceanSimulation}) = CPU()
eltype(model::OceanSeaIceModel{<:Any, <:Any, <:VerosOceanSimulation}) = Float64

function remove_outputs(setup::Symbol)
    rm("$(setup).averages.nc", force=true)
    rm("$(setup).energy.nc", force=true)
    rm("$(setup).overturning.nc", force=true)
    rm("$(setup).snapshot.nc", force=true)
    return nothing
end

const CCField2D = Field{<:Center, <:Center, <:Nothing}
const FCField2D = Field{<:Face,   <:Center, <:Nothing}
const CFField2D = Field{<:Center, <:Face,   <:Nothing}

function set!(field::CCField2D, pyarray::Py, k=pyconvert(Int, pyarray.shape[2]))
    array = PyArray(pyarray)
    Nx, Ny, Nz = size(array)
    set!(field, view(array, 3:Nx-2, 3:Ny-2, k, 1))
    return field
end

function set!(field::FCField2D, pyarray::Py, k=pyconvert(Int, pyarray.shape[2]))
    array = PyArray(pyarray)
    Nx, Ny, Nz = size(array)
    TX, TY, _  = topology(field.grid)
    i_indices  = TX == Periodic ? UnitRange(3, Nx-2) : UnitRange(2, Nx-2)
    set!(field, view(array, i_indices, 3:Ny-2, k, 1))
    return field
end

function set!(field::CFField2D, pyarray::Py, k=pyconvert(Int, pyarray.shape[2]))
    array = PyArray(pyarray)
    Nx, Ny, Nz = size(array)
    set!(field, view(array, 3:Nx-2, 2:Ny-2, k, 1))
    return field
end

"""
    VerosOceanSimulation(setup, setup_name::Symbol)

Creates and initializes a preconfigured Veros ocean simulation using the 
specified setup module and setup name.

Arguments
==========
- `setup::AbstractString`: The name of the Veros setup module to import (e.g., `"global_4deg"`).
- `setup_name::Symbol`: The name of the setup class or function within the module to instantiate (e.g., `:GlobalFourDegreeSetup`).
"""
function VerosOceanSimulation(setup, setup_name::Symbol)
    setups = pyimport("veros.setups." * setup)
    setup  = @eval $setups.$setup_name()

    # instantiate the setup
    setup.setup()

    return VerosOceanSimulation(setup) 
end

"""
    surface_grid(ocean::VerosOceanSimulation)

Constructs a `LatitudeLongitudeGrid` representing the surface grid of the given `VerosOceanSimulation` object.
Notes: Veros always uses a LatitudeLongitudeGrid with 2 halos in both the latitude and longitude directions.
Both latitude and longitude can be either stretched or uniform, depending on the setup, and while the meridional
direction (latitude) is always Bounded, the zonal direction (longitude) can be either Periodic or Bounded.

Arguments
==========
- `ocean::VerosOceanSimulation`: The ocean simulation object containing the grid state variables.
"""
function surface_grid(ocean::VerosOceanSimulation)

    xf = Array(PyArray(ocean.setup.state.variables.xu))
    yf = Array(PyArray(ocean.setup.state.variables.yu))
    
    xc = Array(PyArray(ocean.setup.state.variables.xt))
    yc = Array(PyArray(ocean.setup.state.variables.yt))
    
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

    return LatitudeLongitudeGrid(size=(Nx, Ny), longitude=xf, latitude=yf, topology=(TX, Bounded, Flat), halo=(2, 2))
end

"""
    set!(class::Py, v, x)

Set the `s` property of a Python `class` to the value of `x`.
"""
function set!(class::Py, s, x)
    pyexec("""
       with class.unlock():
           class.__setattr__(y, t)
       """, Main, (y=s, t=x, class=class))
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