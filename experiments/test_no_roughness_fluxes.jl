using ClimaOcean
using Oceananigans

using SeawaterPolynomials: TEOS10EquationOfState

using ClimaOcean
using ClimaOcean.OceanSeaIceModels: Radiation
using ClimaOcean.DataWrangling.JRA55: JRA55_prescribed_atmosphere
using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: SimilarityScales, compute_atmosphere_ocean_fluxes!

#####
##### Single column grid with one vertical level (2.5 meters deep)
#####

Nz = 1
Lz = 2.5

# In the middle of the Atlantic (one degree spacing)
φ₁ = - 1
φ₂ = + 0

λ₁ = 160
λ₂ = 161

grid = LatitudeLongitudeGrid(size = (1, 1, 1), longitude = (λ₁, λ₂), latitude = (φ₁, φ₂), z = (-Lz, 0))

#####
##### Define the simulation
#####

Qᵀ = Field{Center, Center, Nothing}(grid)
Fˢ = Field{Center, Center, Nothing}(grid)
τˣ = Field{Face, Center, Nothing}(grid)
τʸ = Field{Center, Face, Nothing}(grid)

ocean_boundary_conditions = (u = FieldBoundaryConditions(top=FluxBoundaryCondition(τˣ)),
                             v = FieldBoundaryConditions(top=FluxBoundaryCondition(τʸ)),
                             T = FieldBoundaryConditions(top=FluxBoundaryCondition(Qᵀ)),
                             S = FieldBoundaryConditions(top=FluxBoundaryCondition(Fˢ)))

atmosphere = JRA55_prescribed_atmosphere(CPU(), 1:2)

buoyancy  = SeawaterBuoyancy(equation_of_state = TEOS10EquationOfState())

# The two ocean models
ocean_1 = HydrostaticFreeSurfaceModel(; grid, buoyancy, boundary_conditions = ocean_boundary_conditions)
ocean_2 = HydrostaticFreeSurfaceModel(; grid, buoyancy, boundary_conditions = deepcopy(ocean_boundary_conditions))

set!(ocean_1, u = 1, v = 2, T = 1, S = 30)
set!(ocean_2, u = 1, v = 2, T = 1, S = 30)

ocean_simulation_1 = Simulation(ocean_1; Δt = 1)
ocean_simulation_2 = Simulation(ocean_2; Δt = 1)

radiation = Radiation()

# Setting up the comparison: the two coupled models
coupled_model_1 = OceanSeaIceModel(ocean_simulation_1, nothing; 
                                   atmosphere, radiation, 
                                   roughness_lengths = SimilarityScales(1e-3, 1e-3))
                             
coupled_model_2 = OceanSeaIceModel(ocean_simulation_2, nothing; 
                                   atmosphere, radiation, 
                                   roughness_lengths = SimilarityScales(1e-3, 1e-3, 1e-3))
                             
# Computing the fluxes:
compute_atmosphere_ocean_fluxes!(coupled_model_1)
compute_atmosphere_ocean_fluxes!(coupled_model_2)
