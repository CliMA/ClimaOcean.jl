# # OMIP evaluation figures
#
# The figures compare the output of one or more OMIP simulations with observations. Each run is read with
# `OMIPOutput`, whose `averaging_window` selects the records that time means and climatologies average; here the
# last ten years of the first 60-year forcing cycle. The observations are downloaded on first use.

using ClimaOcean
using Oceananigans.Units
using CairoMakie
using GeoMakie

averaging_window = (50, 60) .* 365days

one_degree = OMIPOutput("one_degree_omip"; prefix = "one_degree", label = "1ᵒ", averaging_window)
quarter_degree = OMIPOutput("quarter_degree_omip"; prefix = "quarter_degree", label = "1/4ᵒ", averaging_window)

# One figure, comparing the two runs side by side

sea_surface_temperature_bias(one_degree, quarter_degree)

# and every figure, saved to `omip_figures/`

omip_figures(one_degree, quarter_degree; directory = "omip_figures")
