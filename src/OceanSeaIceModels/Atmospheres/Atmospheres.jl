module Atmospheres

export PrescribedAtmosphere, PrognosticAtmosphere

using KernelAbstractions: @kernel, @index

include("atmospheric_parameters.jl")
include("prescribed_atmospheres.jl")

end