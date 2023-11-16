ice_grid = LatitudeLongitudeGrid(arch,
                                 size = (Nx, Ny),
                                 longitude = (0, 360),
                                 halo = (7, 7),
                                 latitude = (southern_limit, northern_limit),
                                 topology = (Periodic, Bounded, Flat))

ice_grid = ImmersedBoundaryGrid(ice_grid, GridFittedBottom(bottom_height))

Nz = size(grid, 3)
So = ocean_model.tracers.S
ocean_surface_salinity = view(So, :, :, Nz)
bottom_bc = IceWaterThermalEquilibrium(ocean_surface_salinity)

u, v, w = ocean_model.velocities
ocean_surface_velocities = (u = view(u, :, :, Nz), #interior(u, :, :, Nz),
                            v = view(v, :, :, Nz), #interior(v, :, :, Nz),
                            w = ZeroField())

ice_model = SlabSeaIceModel(ice_grid;
                            velocities = ocean_surface_velocities,
                            advection = nothing,
                            ice_consolidation_thickness = 0.05,
                            ice_salinity = 4,
                            internal_heat_flux = ConductiveFlux(conductivity=2),
                            top_heat_flux = ConstantField(0), # W m⁻²
                            top_heat_boundary_condition = PrescribedTemperature(0),
                            bottom_heat_boundary_condition = bottom_bc,
                            bottom_heat_flux = ice_ocean_heat_flux)

set!(ice_model, h=ℋᵢ) 

ice = Simulation(ice_model, Δt=5minutes, verbose=false)

