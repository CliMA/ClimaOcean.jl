using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Printf
using Dates
using CUDA
using JLD2

arch = GPU()

Nx = 720 # longitudinal direction
Ny = 360 # meridional direction
Nz = 100

z_faces = ExponentialDiscretization(Nz, -6000, 0; scale=1800)

const z_surf = z_faces(Nz)

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             latitude = (-80, 80),
                             longitude = (0, 360),
                             z = z_faces,
                             halo = (7, 7, 7))

bottom_height = regrid_bathymetry(grid; minimum_depth=20, major_basins=1, interpolation_passes=25)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

#####
##### A Propgnostic Ocean model
#####

tracer_advection   = WENO(order=7)
momentum_advection = WENOVectorInvariant(order=5)
free_surface       = SplitExplicitFreeSurface(grid; cfl=0.8, fixed_Δt=40minutes)

@inline Δ²ᵃᵃᵃ(i, j, k, grid, lx, ly, lz) =  2 * (1 / (1 / Δx(i, j, k, grid, lx, ly, lz)^2 + 1 / Δy(i, j, k, grid, lx, ly, lz)^2))
@inline geometric_νhb(i, j, k, grid, lx, ly, lz, clock, fields, λ) = Δ²ᵃᵃᵃ(i, j, k, grid, lx, ly, lz)^2 / λ

horizontal_viscosity = HorizontalScalarBiharmonicDiffusivity(ν=geometric_νhb, discrete_form=true, parameters=25days)
catke_closure = ClimaOcean.OceanSimulations.default_ocean_closure() 
closure = (catke_closure, horizontal_viscosity)

dataset = ECCO4Monthly()
start_date = DateTime(1992, 1, 1)
@inline mask(x, y, z, t) = z ≥ z_surf - 1
Smetadata = Metadata(:salinity; dataset, start_date=date)
FS = DatasetRestoring(Smetadata, arch; rate = 1/30days, mask, time_indices_in_memory=10)

ocean = ocean_simulation(grid; Δt=1minutes,
                         momentum_advection,
                         tracer_advection,
                         timestepper = :SplitRungeKutta3,
                         free_surface,
                         forcing = (; S = FS),
                         closure)

set!(ocean.model, T=Metadatum(:temperature; dataset, date=start_date),
                  S=Metadatum(:salinity;    dataset, date=start_date))

#####
##### A Prognostic Sea-ice model
#####

# Default sea-ice dynamics and salinity coupling are included in the defaults
sea_ice = sea_ice_simulation(grid, ocean; advection=WENO(order=7))

set!(sea_ice.model, h=Metadatum(:sea_ice_thickness;     dataset, date=start_date),
                    ℵ=Metadatum(:sea_ice_concentration; dataset, date=start_date))

#####
##### A Prescribed Atmosphere model
#####

dir = "./"
dataset = MultiYearJRA55()
backend = JRA55NetCDFBackend(30)

atmosphere = JRA55PrescribedAtmosphere(arch; dir, dataset, backend, include_rivers_and_icebergs=true, start_date)
radiation  = Radiation()

#####
##### An ocean-sea ice coupled model
#####

omip = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
omip = Simulation(omip, Δt=30minutes, stop_time=20*365years)


ocean.output_writers[:surface] = JLD2Writer(ocean.model, ocean_outputs;
                                            schedule = TimeInterval(180days),
                                            filename = "halfdegree_ocean_surface_fields",
                                            indices = (:, :, grid.Nz),
                                            overwrite_existing = true)

sea_ice.output_writers[:surface] = JLD2Writer(ocean.model, sea_ice_outputs;
                                              schedule = TimeInterval(180days),
                                              filename = "halfdegree_sea_ice_surface_fields",
                                              overwrite_existing = true)

ocean.output_writers[:sample_average] = JLD2Writer(ocean.model, ocean_outputs;
                                        schedule = TimeInterval(365days),
                                        filename = "halfdegree_ocean_complete_fields",
                                        overwrite_existing = true)

run!(omip)