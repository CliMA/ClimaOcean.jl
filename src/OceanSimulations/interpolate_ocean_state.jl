# When using an Oceananigans simulation, we assume that the exchange grid is the ocean grid
# We need, however, to interpolate the surface pressure to the ocean grid
interpolate_ocean_state!(interfaces, ::Simulation{<:HydrostaticFreeSurfaceModel}, coupled_model) = nothing