using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids
using ClimaOcean.Oceans
using ClimaOcean.ECCO
using ClimaOcean.JRA55
using ClimaOcean.DataWrangling
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium
using ClimaOceanBiogeochemistry
using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients #, transfer_surface_atmospheric_state_for_bgc!
using ClimaOceanBiogeochemistry.CarbonSystemSolvers.UniversalRobustCarbonSolver: UniversalRobustCarbonSystem
using ClimaOceanBiogeochemistry.CarbonSystemSolvers: CarbonSystemParameters, CarbonSolverOptions, CarbonCoefficientParameters
using KernelAbstractions: @kernel, @index
using Oceananigans.Architectures: architecture
using Oceananigans: TendencyCallsite, TimeStepCallsite, UpdateStateCallsite
using Oceananigans.Utils: launch!
using Oceananigans.Grids: inactive_cell
using Printf
using Dates
using CUDA
using JLD2

function synch!(clock1::Clock, clock2)
    # Synchronize the clocks
    clock1.time = clock2.time
    clock1.iteration = clock2.iteration
    clock1.last_Δt = clock2.last_Δt
end

synch!(model1, model2) = synch!(model1.clock, model2.clock)
restart_iteration = "460000" 

arch = GPU()

Nx = 720 # longitudinal direction 
Ny = 360 # meridional direction 
Nz = 100

z_faces = ExponentialDiscretization(Nz, -6000, 0; scale=1800, mutable=true)

const z_surf = z_faces.cᵃᵃᶠ(Nz)

grid = TripolarGrid(arch;
                    size = (Nx, Ny, Nz),
                    z = z_faces,
                    halo = (7, 7, 7))

if isfile("bottom_height.jld2")
    bottom_height = on_architecture(arch, jldopen("bottom_height.jld2")["bottom"])
else
    bottom_height = regrid_bathymetry(grid; minimum_depth=20, major_basins=1, interpolation_passes=25)
    jldsave("bottom_height.jld2", bottom=bottom_height)
end

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

#####
##### A Propgnostic Ocean model
#####

using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization, FivePointHorizontalFilter, DiffusiveFormulation, AdvectiveFormulation, IsopycnalSkewSymmetricDiffusivity
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity, CATKEMixingLength, CATKEEquation

buffer_scheme      = Centered(order=2)
buffer_scheme      = WENO(order=3; buffer_scheme)
buffer_scheme      = WENO(order=5; buffer_scheme)
tracer_advection   = WENO(order=7; buffer_scheme)
momentum_advection = WENOVectorInvariant(order=5)
free_surface       = SplitExplicitFreeSurface(grid; substeps=150) 

@inline νhb(i, j, k, grid, ℓx, ℓy, ℓz, clock, fields, λ) = Oceananigans.Operators.Az(i, j, k, grid, ℓx, ℓy, ℓz)^2 / λ
@inline zerofunc(i, j, k, grid, ℓx, ℓy, ℓz, clock, fields, λ) = zero(grid)

horizontal_viscosity = HorizontalScalarBiharmonicDiffusivity(ν=νhb, discrete_form=true, parameters=40days) 
catke_closure = ClimaOcean.Oceans.default_ocean_closure() 
eddy_closure = IsopycnalSkewSymmetricDiffusivity(κ_skew=500, κ_symmetric=100) 
closure = (catke_closure, eddy_closure, horizontal_viscosity, VerticalScalarDiffusivity(ν=1e-5, κ=2e-6))

dataset = EN4Monthly()
start_date = DateTime(1993, 1, 1)
end_date   = DateTime(1994, 1, 1)
@inline mask(x, y, z, t) = z ≥ z_surf - 1
Smetadata = Metadata(:salinity; dataset, start_date, end_date)
FS = DatasetRestoring(Smetadata, grid; rate = 1/30days, mask, time_indices_in_memory=2) 

Fmetadata = Metadata(:dissolved_iron; dataset = ECCO4DarwinMonthly(), start_date, end_date)
FF = DatasetRestoring(Fmetadata, grid; rate = 1/30days, mask, time_indices_in_memory = 2)

Dustmetadata = Metadata(:aeolian_iron_deposition; dataset = ECCO4DarwinMonthly(), start_date, end_date)
AeolianIron_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(20.0))
surface_forcing = (; S=FS, Fe=FF)

FC  = Field{Center, Center, Nothing}(grid)
FA  = Field{Center, Center, Nothing}(grid)
FF2 = Field{Center, Center, Nothing}(grid)

# Initialize carbon surface boundary condition (CO2 fluxes)
BC  = FieldBoundaryConditions(
    top = FluxBoundaryCondition(FC),
    )
BA = FieldBoundaryConditions(
    top = FluxBoundaryCondition(FA),
    )
BF = FieldBoundaryConditions(
    top = FluxBoundaryCondition(FF2),
    )
surface_boundary_conditions = (; DIC=BC, ALK=BA, Fe=BF)

# Demonstrate how to use carbon solver options to set the maximum number 
#  of iterations for the pH and pCO₂ solvers.
# Could also set some of the carbon system parameters here if we wanted to, 
#  e.g. for ionic strength dependent carbon chemistry, we could set the 
#  coefficient parameters:
#     parms = CarbonSystemParameters(
#                 Pᵘˢ = CarbonCoefficientParameters(
#                     a₀ = 0.02,
#     ))
chem_params = CarbonSystemParameters(
    Sᵒᵖᵗˢ=CarbonSolverOptions(
        Iᴴ⁺ₘₐₓ=250,
))

bgc = CarbonAlkalinityNutrients(;
            grid,
            maximum_net_community_production_rate = 0.5/365.25days,
            atmospheric_pCO₂                      = 380e-6,
            carbon_system_params                  = chem_params,
)

ocean = ocean_simulation(grid; Δt=1minutes,
                         momentum_advection,
                         tracer_advection,
                         timestepper = :SplitRungeKutta3,
                         free_surface,
                         radiative_forcing = nothing,
                         closure,
                         tracers = (:T, :S, :DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                         biogeochemistry=bgc,
                         forcing = surface_forcing,
                         boundary_conditions=surface_boundary_conditions,
)

set!(ocean.model, T=EN4Metadatum(:temperature; date=start_date),
                  S=EN4Metadatum(:salinity;    date=start_date),
                  DIC=Metadatum(:dissolved_inorganic_carbon;     dataset = ECCO4DarwinMonthly(), date=start_date),
                  ALK=Metadatum(:alkalinity;                     dataset = ECCO4DarwinMonthly(), date=start_date),
                  PO₄=Metadatum(:phosphate;                      dataset = ECCO4DarwinMonthly(), date=start_date),
                  NO₃=Metadatum(:nitrate;                        dataset = ECCO4DarwinMonthly(), date=start_date),
                  DOP=Metadatum(:dissolved_organic_phosphorus;   dataset = ECCO4DarwinMonthly(), date=start_date),
                  POP=Metadatum(:particulate_organic_phosphorus; dataset = ECCO4DarwinMonthly(), date=start_date),
                  Fe =Metadatum(:dissolved_iron;                 dataset = ECCO4DarwinMonthly(), date=start_date),
)

#####
##### A Prognostic Sea-ice model
#####

# Default sea-ice dynamics and salinity coupling are included in the defaults
sea_ice = sea_ice_simulation(grid, ocean; advection=WENO(order=7, minimum_buffer_upwind_order=1)) 

set!(sea_ice.model, h=Metadatum(:sea_ice_thickness;     dataset=ECCO4Monthly()),
                    ℵ=Metadatum(:sea_ice_concentration; dataset=ECCO4Monthly()))

#####
##### A Prescribed Atmosphere model
#####

dir = "./forcing_data"
dataset = MultiYearJRA55()
backend = JRA55NetCDFBackend(30)

atmosphere = JRA55PrescribedAtmosphere(arch; dir, dataset, backend, include_rivers_and_icebergs=true, start_date)
radiation  = Radiation()

#####
##### An ocean-sea ice coupled model
#####

sim = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
sim = Simulation(sim, Δt=20minutes, stop_time=1*365days) 

#@show sim.model.interfaces.sea_ice_ocean_interface.flux_formulation

# Figure out the outputs....

synch!(sim.model, ocean.model)

wall_time = Ref(time_ns())

function check_salinity(sim)
   if minimum(sim.model.ocean.model.tracers.S) < 10
       @info "minimum of salinity < 10 at iteration $(sim.model.clock.iteration)"
   end
end

function save_datafile(sim; checkpoint=false)
    ocean  = sim.model.ocean
    seaice = sim.model.sea_ice
    suffix = "halfdegree_iteration$(ocean.model.clock.iteration)" 

    MyFloat32 = checkpoint ? Float64 : Float32

    uo = MyFloat32.(Array(interior(ocean.model.velocities.u)))  
    vo = MyFloat32.(Array(interior(ocean.model.velocities.v)))
    wo = MyFloat32.(Array(interior(ocean.model.velocities.w)))
    To = MyFloat32.(Array(interior(ocean.model.tracers.T)))
    So = MyFloat32.(Array(interior(ocean.model.tracers.S)))
    
    eo = if haskey(ocean.model.tracers, :e)
        MyFloat32.(Array(interior(ocean.model.tracers.e)))
    else
        0
    end
    
    ηo = MyFloat32.(Array(interior(ocean.model.free_surface.displacement))) 
    ui = MyFloat32.(Array(interior(seaice.model.velocities.u)))
    vi = MyFloat32.(Array(interior(seaice.model.velocities.v)))
    hi = MyFloat32.(Array(interior(seaice.model.ice_thickness)))
    ℵi = MyFloat32.(Array(interior(seaice.model.ice_concentration)))

    if checkpoint
        suffix *= "_checkpoint"
    end

    jldsave(suffix * ".jld2"; uo, vo, wo, To, clock = ocean.model.clock, So, eo, ηo, ui, vi, hi, ℵi)
end

using Statistics

function progress(sim)
    sea_ice = sim.model.sea_ice
    ocean   = sim.model.ocean
    bgc     = ocean.model.biogeochemistry
    hmax = maximum(sea_ice.model.ice_thickness)
    ℵmax = maximum(sea_ice.model.ice_concentration)
    Tmax = maximum(sim.model.interfaces.atmosphere_sea_ice_interface.temperature)
    Tmin = minimum(sim.model.interfaces.atmosphere_sea_ice_interface.temperature)
    umax = maximum(ocean.model.velocities.u)
    vmax = maximum(ocean.model.velocities.v)
    wmax = maximum(ocean.model.velocities.w)
    cmin = minimum(ocean.model.tracers.DIC)
    cmax = maximum(ocean.model.tracers.DIC)
    amin = minimum(ocean.model.tracers.ALK)
    amax = maximum(ocean.model.tracers.ALK)
    pmin = minimum(ocean.model.tracers.PO₄)
    pmax = maximum(ocean.model.tracers.PO₄)
    pco2ave = mean(bgc.ocean_pCO₂)*1e6
    phave = mean(bgc.pH)
    iave = mean(bgc.PAR)
    step_time = 1e-9 * (time_ns() - wall_time[])

    msg1 = @sprintf("time: %s, iteration: %d, Δt: %s, ", prettytime(sim), iteration(sim), prettytime(sim.Δt))
    msg2 = @sprintf("extrema(DIC): (%.2f, %.2f), mol C/m³, ", cmin, cmax)
    msg3 = @sprintf("extrema(ALK): (%.2f, %.2f) mol eq/m³, ", amin, amax)
    msg4 = @sprintf("extrema(PO₄): (%.2e, %.2e) mol P/m³, ", pmin, pmax)
    msg5 = @sprintf("average pH: %.2f, ", phave)
    msg6 = @sprintf("average pCO₂: %.2f ppm, ", pco2ave)
    msg7 = @sprintf("average I: %.2f, ", iave)
    msg8 = @sprintf("extrema(T): (%.2f, %.2f) ᵒC, ", Tmax, Tmin)
    msg9 = @sprintf("wall time: %s \n", prettytime(step_time))

    @info msg1 * msg2 * msg3 * msg4 * msg5 * msg6 * msg7 * msg8 * msg9

     wall_time[] = time_ns()

     return nothing
end

# And add it as a callback to the simulation.
add_callback!(sim, progress, IterationInterval(10))
add_callback!(sim, sim -> save_datafile(sim), IterationInterval(500))
add_callback!(sim, sim -> save_datafile(sim; checkpoint=true), IterationInterval(20000))

#function update_boundary_conditions!(sim, params)
#    transfer_surface_atmospheric_state_for_bgc!(sim, params)
#    return nothing
#end
#
#incident_PAR = Field{Center, Center, Nothing}(grid)
#
#"""
#Set the PAR to 40% of the downwelling shortwave radiation
#    at the ocean surface.
#"""
#@kernel function calculate_par_from_surface_Qs!(grid, parfac, incident_par, Qs)
#    FT = eltype(grid)
#    i, j = @index(Global, NTuple)
#    ks    = size(grid, 3)
#    inactive = inactive_cell(i, j, ks, grid)
#
#    @inbounds incident_par[i, j, 1] = ifelse(
#        inactive,
#        zero(FT),
#        parfac * Qs[i, j, 1]
#    )
#end
#"""
#Use salt forcing to inform DIC and ALK forcing using surface average values
#    to convert salinity flux to carbon and alkalinity 'virtual' fluxes.
#"""
#@kernel function calculate_bgc_fw_forcing!(grid, C₀, A₀, S₀, FS, FC, FA)
#    FT = eltype(grid)
#    i, j = @index(Global, NTuple)
#    ks    = size(grid, 3)
#    inactive = inactive_cell(i, j, ks, grid)
#
#    # Add carbon flux to the CO2 flux
#    @inbounds begin
#        FC[i, j, 1] = ifelse(
#            inactive,
#            zero(FT),
#            (C₀[i, j, ks]/S₀[i, j, ks]) * FS[i, j, 1],
#        )
#        FA[i, j, 1] = ifelse(
#            inactive,
#            zero(FT),
#            (A₀[i, j, ks]/S₀[i, j, ks]) * FS[i, j, 1],
#        )
#    end
#end
#
#"""
#    Transfer the surface shortwave radiation and surface salinity, to the incident PAR field
#        (par is ~40% of Qs) and surface DIC and ALK forcing fields respectively.
#"""
#@inline function transfer_surface_atmospheric_state_for_bgc!(simulation::Simulation, params)
#    grid = simulation.model.ocean.model.grid
#
#    kernel_args = (
#        grid,
#        params.par_fraction,
#        params.incident_PAR,
#        simulation.model.interfaces.atmosphere_ocean_interface.fluxes.downwelling_shortwave,
#        )
##
#    launch!(
#	    architecture(grid),
#        grid,
#        :xy,
#        calculate_par_from_surface_Qs!,
#        kernel_args...
#    )
#
#    kernel_args =(
#        grid,
#        simulation.model.ocean.model.tracers.DIC,
#        simulation.model.ocean.model.tracers.ALK,
#        simulation.model.ocean.model.tracers.S,
#        simulation.model.interfaces.net_fluxes.ocean.S, # The FW forcing
#        simulation.model.ocean.model.tracers.DIC.boundary_conditions.top.condition,
#        simulation.model.ocean.model.tracers.ALK.boundary_conditions.top.condition,
#    )
#    launch!(
#        architecture(grid),
#        grid,
#        :xy,
#        calculate_bgc_fw_forcing!,
#        kernel_args...,
#    )
#    return nothing
#end
add_callback!(sim, 
    ClimaOceanBiogeochemistry.transfer_surface_atmospheric_state_for_bgc!, 
    callsite = TimeStepCallsite(), 
    IterationInterval(1),
)
add_callback!(sim, 
    ClimaOceanBiogeochemistry.calculate_air_sea_carbon_exchange!, 
    callsite = TimeStepCallsite(), 
    IterationInterval(1),
)
run!(sim)
