module AtmosphereSimulations

export atmosphere_simulation

using Oceananigans
using Oceananigans.Utils
using Oceananigans.Fields: Center
using Oceananigans.Grids: grid_name, architecture, topology, Flat
using Oceananigans.OutputReaders: FieldTimeSeries, update_field_time_series!, extract_field_time_series
using Oceananigans.TimeSteppers: Clock, tick!
using Oceananigans.Utils: prettysummary, Time

using Adapt
using Thermodynamics.Parameters: AbstractThermodynamicsParameters
using KernelAbstractions: @kernel, @index
using ClimaOcean.OceanSeaIceModels.InterfaceComputations: interface_kernel_parameters

import Oceananigans.TimeSteppers: time_step!, update_state!

import ClimaOcean.OceanSeaIceModels: interpolate_state!, 
                                     compute_net_fluxes!, 
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

end # module AtmosphereSimulations
