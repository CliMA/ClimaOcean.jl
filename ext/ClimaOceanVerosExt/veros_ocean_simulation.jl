using CondaPkg

using Oceananigans.Grids: topology
using OffsetArrays: OffsetArray

import Oceananigans.Fields: set!
import Oceananigans.TimeSteppers: time_step!, initialize!

import Oceananigans.Architectures: architecture

import ClimaOcean.OceanSeaIceModels: default_nan_checker
import ClimaOcean.OceanSeaIceModels: reference_density, 
                                     heat_capacity, 
                                     ocean_temperature, 
                                     ocean_salinity, 
                                     ocean_surface_salinity,
                                     ocean_surface_velocities

import ClimaOcean.Oceans: ocean_simulation

import Base: eltype

"""
    install_veros()

Install the Veros ocean model Marine CLI using CondaPkg.
Returns a NamedTuple containing package information if successful.
"""
function install_veros()
    CondaPkg.add_pip("veros")
    cli = CondaPkg.which("veros")
    @info "... the veros CLI has been installed at $(cli)."
    return cli
end

struct VerosOceanSimulation{S}
    setup :: S
end

default_nan_checker(model::OceanSeaIceModel{<:Any, <:Any, <:VerosOceanSimulation}) = nothing

initialize!(::ClimaOceanVerosExt.VerosOceanSimulation{Py}) = nothing
time_step!(ocean::VerosOceanSimulation, Δt) = ocean.setup.step(ocean.setup.state)
architecture(model::OceanSeaIceModel{<:Any, <:Any, <:VerosOceanSimulation}) = CPU()
eltype(model::OceanSeaIceModel{<:Any, <:Any, <:VerosOceanSimulation}) = Float64

reference_density(ocean::VerosOceanSimulation) = pyconvert(eltype(ocean), ocean.setup.state.settings.rho_0)
heat_capacity(ocean::VerosOceanSimulation) = convert(eltype(ocean), 3995)

function ocean_surface_velocities(ocean::VerosOceanSimulation)
    u = PyArray(ocean.setup.state.variables.u)
    v = PyArray(ocean.setup.state.variables.v)
    Nxu, Nyu, Nzu = size(u)
    Nxv, Nyv, Nzv = size(v)
    grid = surface_grid(ocean)
    TX, TY, _  = topology(grid)
    i_indices  = TX == Periodic ? UnitRange(3, Nxu-2) : UnitRange(2, Nxu-2)

    u_view = view(u, i_indices, 3:Nyu-2, Nzu, 1)
    v_view = view(v, 3:Nxv-2,   2:Nyv-2, Nzv, 1)

    return u_view, v_view
end

function ocean_surface_salinity(ocean::VerosOceanSimulation) 
    S = ocean_salinity(ocean)
    Nx, Ny, Nz = size(S)
    return view(S, :, :, Nz)
end

# Veros hardcodes 2 halos in the x and y direction,
# and each prognostic variable is 4 dimensional, where the first three dimensions
# are x, y, z, and the last index differentiate between variable, tendency at n and tendency at n-1
function ocean_temperature(ocean::VerosOceanSimulation) 
    T = PyArray(ocean.setup.state.variables.temp)
    Nx, Ny, Nz = size(T)
    return view(T, 3:Nx-2, 3:Ny-2, 1:Nz, 1)
end

# Veros hardcodes 2 halos in the x and y direction,
# and each prognostic variable is 4 dimensional, where the first three dimensions
# are x, y, z, and the last index differentiate between variable, tendency at n and tendency at n-1
function ocean_salinity(ocean::VerosOceanSimulation) 
    S = PyArray(ocean.setup.state.variables.salt)
    Nx, Ny, Nz = size(S)
    return view(S, 3:Nx-2, 3:Ny-2, 1:Nz, 1)
end

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
function VerosOceanSimulation(setup::String, setup_name::Symbol)
    setups = pyimport("veros.setups." * setup)
    ocean  = @eval $setups.$setup_name()

    # instantiate the setup
    ocean.setup()

    return VerosOceanSimulation(ocean) 
end

# We assume that if we pass a python object, this is a veros simulation
function ocean_simulation(ocean::Py)

    # instantiate the setup
    ocean.setup()

    return VerosOceanSimulation(ocean)
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
    set!(ocean, v, x; path = :variable)

Set the `v` variable in the `ocean` model to the value of `x`.
the path corresponds to the path inside the class where to locate the 
variable `v` to set. It can be either `:variables` or `:settings`.
"""
function set!(ocean::VerosOceanSimulation, v, x; path = :variables)
    setup = ocean.setup
    if path == :variables
        pyexec("""
        with setup.state.variables.unlock():
            setup.state.variables.__setattr__(y, t)
        """, Main, (y=v, t=x, setup=setup))
    elseif path == :settings
        pyexec("""
        with setup.state.settings.unlock():
            setup.state.settings.__setattr__(y, t)
        """, Main, (y=v, t=x, setup=setup))
    else
        error("path must be either :variable or :settings.")
    end
end
