#####
##### OMIP configurations
#####

"""
    OMIPConfiguration

Every physical and numerical choice of an OMIP simulation, grouped into the `grid`, the `ocean`, the `sea_ice`, the
`interfaces` between them and the atmosphere, and the `forcing`, together with the `name` that prefixes the output
files and the `time_step`. Built by [`OneDegree`](@ref), [`QuarterDegree`](@ref) and [`TwelfthDegree`](@ref) with
the defaults of the production runs, or directly through the keyword constructor below.
"""
struct OMIPConfiguration{G, O, S, I, F, T}
    name :: String
    grid :: G
    ocean :: O
    sea_ice :: S
    interfaces :: I
    forcing :: F
    time_step :: T
end

const vertically_implicit_advection = AdaptiveVerticallyImplicitDiscretization(cfl = 1)

"""
    OMIPConfiguration(; kwargs...)

An OMIP configuration with the settings common to every resolution as defaults. The keyword arguments are

Grid
====

- `name`: prefix of the output files. Default: `"omip"`.
- `dataset`: the NEMO eORCA mesh and bathymetry. Default: `ORCAOne()`.
- `vertical_coordinate`: the vertical discretization. Default: 100 levels stretching from 1.5 m to 100 m.
- `Nz`: number of vertical levels. Default: 100.
- `minimum_depth`: minimum depth [m] of the wet columns. Default: 0.
- `polar_cap_latitude`: latitude north of which the initial temperature and salinity are horizontally uniform,
  or `nothing`. Default: `nothing`.

Ocean
=====

- `vertical_mixing`: the vertical mixing closure. Default: CATKE with `Cᵇ = 0.01` and diffusivities capped at 3 m² s⁻¹.
- `background_diffusivity`: vertical background tracer diffusivity [m² s⁻¹], a number or a function of
  `(λ, φ, z, t)`. Default: 10⁻⁵.
- `background_viscosity`: vertical background viscosity [m² s⁻¹]. Default: 3 × 10⁻⁵.
- `eddy_closure`: the Gent-McWilliams and Redi closure, or `nothing`. Default: `nothing`.
- `biharmonic_timescale`: timescale on which the biharmonic viscosity damps the grid scale, or `nothing`.
  Default: `nothing`.
- `momentum_advection`: Default: `WENOVectorInvariant()`, vertically implicit where the Courant number exceeds one.
- `tracer_advection`: advection of temperature and salinity. Default: seventh-order `WENO`, vertically implicit
  where the Courant number exceeds one.
- `barotropic_substeps`: substeps of the split-explicit free surface per time step. Default: 300.
- `reference_density`: [kg m⁻³]. Default: 1026.
- `chlorophyll_dataset`: the chlorophyll climatology attenuating the penetrating shortwave. Default: `SeaWiFSMonthly()`.
- `initial_condition_dataset`: the temperature and salinity climatology the ocean starts from. Default: `WOAAnnual()`.
- `salinity_restoring_dataset`: the surface salinity climatology the ocean is restored to. Default: `WOAMonthly()`.
- `salinity_restoring_piston_velocity`: [m s⁻¹]. Default: 0.254 m day⁻¹.
- `total_water_normalization`: whether to hold the water of the ocean, sea ice and snow at its initial value.
  Default: `true`.

Sea ice
=======

- `sea_ice_advection`: Default: `IncrementalRemapping()`.
- `sea_ice_ocean_drag_reference_depth`: depth [m] over which the ocean velocity of the ice--ocean drag is averaged.
  Default: 6.
- `sea_ice_conductivity`: [W m⁻¹ K⁻¹]. Default: 2.
- `sea_ice_thickness_categories`: number of uniform thickness categories through which snow and ice conduct heat.
  Default: 4.
- `sea_ice_initial_condition`: the sea ice thickness and concentration the sea ice starts from. Default: the ECCO
  state estimate of January 1993.

Interfaces
==========

- `atmosphere_ocean_fluxes`: Default: the bulk formulae of [Large and Yeager (2009)](@cite large2009global),
  iterated five times.
- `atmosphere_sea_ice_fluxes`: Default: Monin--Obukhov similarity with the SHEBA stability functions, roughness
  lengths of 5 × 10⁻⁴ m for momentum and 5 × 10⁻⁵ m for heat and moisture, and a gustiness of at least 0.2 m s⁻¹.
- `atmosphere_ocean_interface_specific_humidity`, `atmosphere_sea_ice_interface_specific_humidity`: Default: the
  saturation specific humidity enhanced by the moist-air factor of Gill (1982), 98% of it over the ocean.
- `sea_ice_ocean_heat_transfer_coefficient`: Stanton number of the three-equation ice--ocean heat flux.
  Default: 0.0057.
- `ocean_albedo`: Default: 0.06.

Forcing
=======

- `forcing_dataset`: the JRA55-do reanalysis. Default: `MultiYearJRA55()`.
- `river_spread_radius`, `maximum_river_spread_cells`: the footprint over which each river mouth discharges.
  Default: 1.2ᵒ and no cap.

Time stepping
=============

- `time_step`: a [`TimeStepSchedule`](@ref). Default: 30 minutes.
"""
function OMIPConfiguration(; name = "omip",
                             dataset = ORCAOne(),
                             vertical_coordinate = ReferenceToStretchedDiscretization(extent = 5500,
                                                                                      constant_spacing = 1.5,
                                                                                      constant_spacing_extent = 1.5,
                                                                                      maximum_spacing = 100,
                                                                                      stretching = LinearStretching(0.0749),
                                                                                      rounding_digits = 4,
                                                                                      mutable = true),
                             Nz = 100,
                             minimum_depth = 0,
                             polar_cap_latitude = nothing,
                             vertical_mixing = CATKEVerticalDiffusivity(VerticallyImplicitTimeDiscretization();
                                                                        mixing_length = CATKEMixingLength(Cᵇ = 0.01),
                                                                        maximum_viscosity = 3,
                                                                        maximum_tracer_diffusivity = 3,
                                                                        maximum_tke_diffusivity = 3,
                                                                        negative_tke_damping_time_scale = 10),
                             background_diffusivity = 1e-5,
                             background_viscosity = 3e-5,
                             eddy_closure = nothing,
                             biharmonic_timescale = nothing,
                             momentum_advection = WENOVectorInvariant(time_discretization = vertically_implicit_advection),
                             tracer_advection = WENO(order = 7, time_discretization = vertically_implicit_advection),
                             barotropic_substeps = 300,
                             reference_density = 1026,
                             chlorophyll_dataset = SeaWiFSMonthly(),
                             initial_condition_dataset = WOAAnnual(),
                             salinity_restoring_dataset = WOAMonthly(),
                             salinity_restoring_piston_velocity = 0.254 / days,
                             total_water_normalization = true,
                             sea_ice_advection = IncrementalRemapping(),
                             sea_ice_ocean_drag_reference_depth = 6,
                             sea_ice_conductivity = 2,
                             sea_ice_thickness_categories = 4,
                             sea_ice_initial_condition = MetadataSet((:sea_ice_thickness, :sea_ice_concentration);
                                                                     dataset = ECCO4Monthly(),
                                                                     date = DateTime(1993, 1, 1)),
                             atmosphere_ocean_fluxes = CoefficientBasedFluxes(transfer_coefficients = LargeYeagerTransferCoefficients(),
                                                                              solver_stop_criteria = FixedIterations(5)),
                             atmosphere_sea_ice_fluxes = SimilarityTheoryFluxes(stability_functions = atmosphere_sea_ice_stability_functions(),
                                                                                similarity_form = COARELogarithmicSimilarityProfile(),
                                                                                subgrid_velocities = ConvectiveGustiness(minimum_gustiness = 0.2),
                                                                                momentum_roughness_length = 5e-4,
                                                                                temperature_roughness_length = 5e-5,
                                                                                water_vapor_roughness_length = 5e-5),
                             atmosphere_ocean_interface_specific_humidity =
                                 ImpureSaturationSpecificHumidity(AtmosphericThermodynamics.Liquid(), 0.98;
                                                                  saturation_enhancement = GillSaturationEnhancement()),
                             atmosphere_sea_ice_interface_specific_humidity =
                                 ImpureSaturationSpecificHumidity(AtmosphericThermodynamics.Ice();
                                                                  saturation_enhancement = GillSaturationEnhancement()),
                             sea_ice_ocean_heat_transfer_coefficient = 0.0057,
                             ocean_albedo = 0.06,
                             forcing_dataset = MultiYearJRA55(),
                             river_spread_radius = 1.2,
                             maximum_river_spread_cells = nothing,
                             time_step = TimeStepSchedule(30minutes))

    grid = (; dataset, vertical_coordinate, Nz, minimum_depth, polar_cap_latitude)

    ocean = (; vertical_mixing, background_diffusivity, background_viscosity, eddy_closure, biharmonic_timescale,
               momentum_advection, tracer_advection, barotropic_substeps, reference_density, chlorophyll_dataset,
               initial_condition_dataset, salinity_restoring_dataset, salinity_restoring_piston_velocity,
               total_water_normalization)

    sea_ice = (; advection = sea_ice_advection, ocean_drag_reference_depth = sea_ice_ocean_drag_reference_depth,
                 conductivity = sea_ice_conductivity, thickness_categories = sea_ice_thickness_categories,
                 initial_condition = sea_ice_initial_condition)

    interfaces = (; atmosphere_ocean_fluxes, atmosphere_sea_ice_fluxes, atmosphere_ocean_interface_specific_humidity,
                    atmosphere_sea_ice_interface_specific_humidity, sea_ice_ocean_heat_transfer_coefficient,
                    ocean_albedo)

    forcing = (; dataset = forcing_dataset, river_spread_radius, maximum_river_spread_cells)

    return OMIPConfiguration(name, grid, ocean, sea_ice, interfaces, forcing, time_step)
end

"""
    OneDegree(; kwargs...)

The OMIP configuration on the 1ᵒ NEMO eORCA1 mesh: Gent-McWilliams and Redi diffusivities of 500 and 800 m² s⁻¹, a
biharmonic viscosity damping the grid scale in 50 days, river discharge spread over the 8 closest wet cells, and a
90-minute time step. The keyword arguments of [`OMIPConfiguration`](@ref) change any setting.
"""
OneDegree(; kwargs...) = OMIPConfiguration(; merge((; name = "one_degree",
                                                      dataset = ORCAOne(),
                                                      eddy_closure = IsopycnalSkewSymmetricDiffusivity(κ_skew = 500, κ_symmetric = 800),
                                                      biharmonic_timescale = 50days,
                                                      momentum_advection = WENOVectorInvariant(order = 5, time_discretization = vertically_implicit_advection),
                                                      river_spread_radius = nothing,
                                                      maximum_river_spread_cells = 8,
                                                      time_step = TimeStepSchedule(90minutes)), kwargs)...)

"""
    henyey_diffusivity(λ, φ, z, t)

The latitude-dependent internal-wave diffusivity of [Henyey et al. (1986)](@cite henyey1986energy),
``κ = \\max(2 × 10⁻⁶, 10⁻⁵ |\\sin φ|)`` m² s⁻¹.
"""
@inline henyey_diffusivity(λ, φ, z, t) = max(2e-6, 1e-5 * abs(sind(φ)))

# Exponential vertical spacing over `extent` whose top cell is `surface_spacing` thick
function surface_refined_exponential_discretization(Nz, extent, surface_spacing)
    top_spacing(scale) = extent * expm1(extent / Nz / scale) / expm1(extent / scale)
    lower, upper = extent / Nz / 1000, 1000 * extent

    scale = (lower + upper) / 2
    for _ in 1:200
        scale = (lower + upper) / 2
        top_spacing(scale) < surface_spacing ? (lower = scale) : (upper = scale)
    end

    return ExponentialDiscretization(Nz, -extent, 0; scale, mutable = true)
end

"""
    QuarterDegree(; kwargs...)

The OMIP configuration on the 1/4ᵒ NEMO eORCA025 mesh: 100 levels whose spacing grows exponentially from 1.5 m at
the surface, columns at least 20 m deep, resolved eddies with no Gent-McWilliams or Redi diffusivity and no
explicit horizontal viscosity, the background diffusivity of [`henyey_diffusivity`](@ref), a horizontally uniform
initial state north of 88ᵒN, and a time step of 4 minutes for the first 60 days, 30 minutes until the end of the
first year, and 45 minutes afterwards. The keyword arguments of [`OMIPConfiguration`](@ref) change any setting.
"""
QuarterDegree(; kwargs...) = OMIPConfiguration(; merge((; name = "quarter_degree",
                                                          dataset = ORCAQuarter(),
                                                          vertical_coordinate = surface_refined_exponential_discretization(100, 5500, 1.5),
                                                          minimum_depth = 20,
                                                          polar_cap_latitude = 88,
                                                          background_diffusivity = henyey_diffusivity,
                                                          time_step = TimeStepSchedule(4minutes, 60days => 30minutes, 365days => 45minutes)), kwargs)...)

"""
    TwelfthDegree(; kwargs...)

The OMIP configuration on the 1/12ᵒ NEMO eORCA12 mesh, with the settings of [`QuarterDegree`](@ref), 200 barotropic
substeps, and a time step of 6 minutes for the first 181 days and 15 minutes afterwards. The keyword arguments of
[`OMIPConfiguration`](@ref) change any setting.
"""
TwelfthDegree(; kwargs...) = OMIPConfiguration(; merge((; name = "twelfth_degree",
                                                          dataset = ORCATwelfth(),
                                                          vertical_coordinate = surface_refined_exponential_discretization(100, 5500, 1.5),
                                                          minimum_depth = 20,
                                                          polar_cap_latitude = 88,
                                                          background_diffusivity = henyey_diffusivity,
                                                          barotropic_substeps = 200,
                                                          time_step = TimeStepSchedule(6minutes, 181days => 15minutes)), kwargs)...)

timescale_summary(::Nothing) = "nothing"
timescale_summary(timescale) = prettytime(timescale)

function Base.show(io::IO, configuration::OMIPConfiguration)
    grid, ocean, forcing = configuration.grid, configuration.ocean, configuration.forcing
    print(io, "OMIPConfiguration ", configuration.name, " on ", summary(grid.dataset), " with ", grid.Nz, " levels", '\n',
          "├── eddy closure: ", summary(ocean.eddy_closure), '\n',
          "├── biharmonic timescale: ", timescale_summary(ocean.biharmonic_timescale), '\n',
          "├── background diffusivity: ", prettysummary(ocean.background_diffusivity), '\n',
          "├── minimum depth: ", grid.minimum_depth, " m, polar cap latitude: ", prettysummary(grid.polar_cap_latitude), '\n',
          "├── forcing: ", summary(forcing.dataset), '\n',
          "└── time step: ", configuration.time_step)
end
