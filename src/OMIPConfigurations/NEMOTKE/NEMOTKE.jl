# NEMO 3.6 TKE vertical-mixing closure (zdftke + zdfevd) — Blanke & Delecluse 1993,
# Gaspar et al. 1990, Madec et al. 2017. OMIP-2 ORCA1 preset; vendored as a submodule.

module NEMOTKE

export NEMOTKEVerticalDiffusivity, NEMOTKEParameters

using Adapt
using KernelAbstractions: @index, @kernel

using Oceananigans
using Oceananigans.BuoyancyFormulations: ∂z_b
using Oceananigans.Fields: Field
using Oceananigans.Grids: Center, Face, znode, peripheral_node, φnode
using Oceananigans.Operators
using Oceananigans.Operators: ℑxᶜᵃᵃ, ℑyᵃᶜᵃ, ∂zᶠᶜᶠ, ∂zᶜᶠᶠ, Δzᶜᶜᶜ
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: AbstractScalarDiffusivity, VerticalFormulation,
                                       VerticallyImplicitTimeDiscretization,
                                       getclosure
using Oceananigans.Utils: launch!, prettysummary

import Oceananigans.TurbulenceClosures: viscosity, diffusivity,
                                        viscosity_location, diffusivity_location,
                                        with_tracers,
                                        compute_closure_fields!, build_closure_fields
import Oceananigans.TimeSteppers: step_closure_prognostics!, time_discretization

include("nemo_tke_parameters.jl")
include("nemo_tke_vertical_diffusivity.jl")
include("nemo_tke_surface_forcing.jl")
include("nemo_tke_mixing_length.jl")
include("nemo_tke_diffusivities.jl")
include("nemo_tke_langmuir.jl")
include("nemo_tke_wave_penetration.jl")
include("nemo_tke_evd.jl")
include("nemo_tke_compute_closure_fields.jl")

end # module NEMOTKE
