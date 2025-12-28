module Atmospheres

export atmosphere_simulation, PrescribedAtmosphere

using Oceananigans
using Oceananigans.Utils
using Oceananigans.Fields: Center
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: grid_name, architecture, topology, Flat, prettysummary
using Oceananigans.OutputReaders: FieldTimeSeries, update_field_time_series!, extract_field_time_series
using Oceananigans.TimeSteppers: Clock, tick!
using Oceananigans.Units: Time

using Adapt
using Thermodynamics.Parameters: AbstractThermodynamicsParameters
using KernelAbstractions: @kernel, @index
using ClimaOcean.OceanSeaIceModels.InterfaceComputations: interface_kernel_parameters

import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Fields: set!

import ClimaOcean.OceanSeaIceModels: interpolate_state!, 
                                     update_net_fluxes!, 
                                     thermodynamics_parameters, 
                                     surface_layer_height, 
                                     boundary_layer_height

import ClimaOcean.OceanSeaIceModels.InterfaceComputations: ComponentExchanger, initialize!, net_fluxes

# Can be extended by atmosphere models
function atmosphere_simulation end

include("thermodynamic_parameters.jl")
include("prescribed_atmosphere.jl")
include("prescribed_atmosphere_regridder.jl")
include("interpolate_atmospheric_state.jl")

net_fluxes(::PrescribedAtmosphere) = nothing

end # module Atmospheres
