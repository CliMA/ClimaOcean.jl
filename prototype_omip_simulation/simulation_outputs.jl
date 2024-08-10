
function set_outputs!(coupled_simulation)

    ocean = coupled_simulation.model.ocean
    sea_ice = coupled_simulation.model.sea_ice

    fluxes = (u = ocean.model.velocities.u.boundary_conditions.top.condition,
              v = ocean.model.velocities.v.boundary_conditions.top.condition,
              T = ocean.model.tracers.T.boundary_conditions.top.condition,
              S = ocean.model.tracers.S.boundary_conditions.top.condition)
    
    sea_ice_outputs = merge(sea_ice.model.velocities, 
                           (; h = sea_ice.model.ice_thickness, a = sea_ice.model.ice_concentration),
                           sea_ice.model.ice_dynamics.auxiliary_fields)

    ocean.output_writers[:fluxes] = JLD2OutputWriter(ocean.model, fluxes,
                                                     schedule = TimeInterval(0.1days),
                                                     overwrite_existing = false,
                                                     array_type = Array{Float32},
                                                     filename = "surface_fluxes")
    
    ocean.output_writers[:surface] = JLD2OutputWriter(ocean.model, merge(model.tracers, model.velocities),
                                                      schedule = TimeInterval(0.1days),
                                                      overwrite_existing = false,
                                                      array_type = Array{Float32},
                                                      filename = "surface",
                                                      indices = (:, :, grid.Nz))
    
    ocean.output_writers[:snapshots] = JLD2OutputWriter(ocean.model, merge(model.tracers, model.velocities),
                                                        schedule = TimeInterval(10days),
                                                        overwrite_existing = false,
                                                        array_type = Array{Float32},
                                                        filename = "snapshots")
    
    ocean.output_writers[:checkpoint] = Checkpointer(ocean.model, 
                                                     schedule = TimeInterval(60days),
                                                     overwrite_existing = true,
                                                     prefix = "checkpoint_ocean")

    sea_ice.output_writers[:snapshots] = JLD2OutputWriter(sea_ice.model, sea_ice_outputs,
                                                          schedule = TimeInterval(0.1days),
                                                          overwrite_existing = false,
                                                          array_type = Array{Float32},
                                                          filename = "snapshots_sea_ice")
    
#     sea_ice.output_writers[:checkpoint] = Checkpointer(model, 
#                                                        schedule = TimeInterval(60days),
#                                                        overwrite_existing = true,
#                                                        prefix = "checkpoint_sea_ice")

    return nothing
end
