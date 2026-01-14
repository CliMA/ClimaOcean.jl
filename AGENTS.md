# ClimaOcean.jl rules for agent-assisted development

## Project Overview

ClimaOcean.jl is a Julia package for realistic ocean-only and coupled ocean + sea-ice simulations driven by prescribed atmospheres. It is built on top of [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) and [ClimaSeaIce.jl](https://github.com/CliMA/ClimaSeaIce.jl).

Key features include:
- **OceanSeaIceModel**: Core abstraction for coupling ocean, sea ice, and prescribed atmosphere
- **PrescribedAtmosphere**: Interface for driving ocean simulations with reanalysis data (JRA55)
- **Data Wrangling**: Utilities for downloading and processing bathymetry (ETOPO), ocean state (ECCO, GLORYS, EN4), and atmospheric forcing (JRA55)
- **Surface Flux Computation**: Similarity theory and coefficient-based turbulent flux calculations
- **ocean_simulation**: Convenience constructor for realistic ocean simulations with sensible defaults

### Relationship to Oceananigans

ClimaOcean provides the "realistic modeling" layer on top of Oceananigans:
- Oceananigans handles the fluid dynamics (solvers, grids, fields, operators)
- ClimaOcean handles coupling, forcing, and data wrangling for realistic simulations
- Users should be proficient with Oceananigans to use ClimaOcean effectively

## Language & Environment

- **Language**: Julia 1.10+
- **Architectures**: CPU and GPU (via Oceananigans/KernelAbstractions.jl)
- **Key Dependencies**:
  - `Oceananigans.jl` - Core fluid dynamics
  - `ClimaSeaIce.jl` - Sea ice model
  - `SeawaterPolynomials.jl` - TEOS-10 equation of state
  - `SurfaceFluxes.jl` - Surface flux parameterizations (via Thermodynamics.jl)
- **Testing**: Test suite covering data downloading, surface fluxes, and coupled models

## Code Style & Conventions

### Julia Best Practices

1. **Explicit Imports**: Use explicit imports in source code
   - Import from Oceananigans modules explicitly
   - Example: `using Oceananigans.Utils: launch!`
   
2. **Type Stability**: Prioritize type-stable code for GPU performance
   - All structs must be concretely typed
   
3. **Kernel Functions**: For GPU compatibility:
   - Use KernelAbstractions.jl syntax (`@kernel`, `@index`)
   - Keep kernels type-stable and allocation-free
   - Use `ifelse` instead of short-circuiting if-statements when possible
   - No error messages inside kernels
   - Models _never_ go inside kernels
   
4. **Documentation**:
   - Use DocStringExtensions.jl for consistent docstrings
   - Include `$(SIGNATURES)` for automatic signature documentation
   - Add examples in docstrings when helpful

5. **Memory Efficiency**:
   - Favor inline computations over temporary allocations
   - Design solutions that work within the existing Oceananigans framework

### Naming Conventions

- **Files**: snake_case (e.g., `ocean_sea_ice_model.jl`, `similarity_theory_turbulent_fluxes.jl`)
- **Types**: PascalCase (e.g., `OceanSeaIceModel`, `PrescribedAtmosphere`, `SimilarityTheoryFluxes`)
- **Functions**: snake_case (e.g., `ocean_simulation`, `compute_atmosphere_ocean_fluxes!`)
- **Kernels**: May be prefixed with underscore (e.g., `_compute_flux_kernel`)
- **Variables**: Use _either_ an English long name or mathematical notation with readable unicode

### Module Structure

```
src/
├── ClimaOcean.jl                 # Main module, exports
├── Bathymetry.jl                 # Bathymetry regridding utilities
├── SeaIceSimulations.jl          # Sea ice simulation setup
├── OceanSimulations/             # Ocean simulation constructors
│   ├── OceanSimulations.jl
│   ├── ocean_simulation.jl       # Hydrostatic ocean simulation setup
│   └── nonhydrostatic_ocean_simulation.jl
├── OceanSeaIceModels/            # Coupled model implementation
│   ├── OceanSeaIceModels.jl      # Module definition
│   ├── ocean_sea_ice_model.jl    # Core model type
│   ├── PrescribedAtmospheres.jl  # Atmosphere data structures
│   ├── time_step_ocean_sea_ice_model.jl
│   └── InterfaceComputations/    # Flux calculations
│       ├── interface_states.jl
│       ├── similarity_theory_turbulent_fluxes.jl
│       ├── coefficient_based_turbulent_fluxes.jl
│       ├── atmosphere_ocean_fluxes.jl
│       └── ...
├── DataWrangling/                # Data downloading and processing
│   ├── DataWrangling.jl          # Module definition
│   ├── metadata.jl               # Metadata types
│   ├── inpainting.jl             # Gap-filling algorithms
│   ├── restoring.jl              # Dataset restoring forcing
│   ├── ECCO/                     # ECCO data utilities
│   ├── JRA55/                    # JRA55 atmosphere data
│   ├── ETOPO/                    # Bathymetry data
│   ├── Copernicus/               # GLORYS data (via extension)
│   └── EN4/                      # EN4 data
├── InitialConditions/            # Initial condition utilities
│   └── diffuse_tracers.jl
└── Diagnostics/                  # Diagnostic utilities
    └── mixed_layer_depth.jl
```

## Testing Guidelines

### Running Tests

```julia
# All tests
Pkg.test("ClimaOcean")

# Run specific test groups by setting TEST_GROUP environment variable
ENV["TEST_GROUP"] = "fluxes"  # or "JRA55", "ecco2_monthly", "bathymetry", etc.
Pkg.test("ClimaOcean")

# Available test groups:
# - init: Download test data
# - JRA55: JRA55 atmosphere tests
# - ecco2_monthly, ecco2_daily, ecco4_en4: Dataset tests
# - downloading, copernicus_downloading: Download tests
# - fluxes: Surface flux tests
# - bathymetry: Bathymetry regridding tests
# - ocean_sea_ice_model: Coupled model tests
# - distributed: MPI/distributed tests
# - reactant: Reactant.jl integration tests
```

### Writing Tests

- Place tests in `test/` directory
- Follow the existing test group structure in `runtests.jl`
- Use `include("runtests_setup.jl")` for common setup
- Test on both CPU and GPU when possible
- Name test files `test_<feature>.jl`

### Data Dependencies

- Tests download data from external sources (ECCO, JRA55, ETOPO)
- Use `@get_scratch!` for caching downloaded files
- Be mindful of data download times in CI

## Common Development Tasks

### Adding New Data Sources

1. Create a new subdirectory in `src/DataWrangling/` (e.g., `NewDataset/`)
2. Define dataset types and metadata structures
3. Implement required interface functions:
   - `download_dataset(metadata)`
   - `z_interfaces(dataset)`
   - `native_grid(dataset)`
4. Add to `DataWrangling.jl` includes
5. Export from `ClimaOcean.jl` if needed
6. Add tests for downloading and usage

### Modifying Surface Flux Calculations

- Flux calculations are in `src/OceanSeaIceModels/InterfaceComputations/`
- `SimilarityTheoryFluxes` uses Monin-Obukhov similarity theory
- `CoefficientBasedFluxes` uses simple transfer coefficients
- Roughness lengths: `roughness_lengths.jl`
- State interpolation: `interpolate_atmospheric_state.jl`

### Adding New Coupled Model Features

1. Core model: `src/OceanSeaIceModels/ocean_sea_ice_model.jl`
2. Time stepping: `time_step_ocean_sea_ice_model.jl`
3. Interface fluxes: `InterfaceComputations/`
4. Update exports in `OceanSeaIceModels.jl` and `ClimaOcean.jl`

## Documentation

### Building Docs Locally

```sh
julia --project=docs/ docs/make.jl
```

### Viewing Docs

```julia
using LiveServer
serve(dir="docs/build")
```

### Documentation Style

- Use Documenter.jl syntax for cross-references
- Include code examples in documentation pages
- Add references to papers in `climaocean.bib`
- In examples, use `using ClimaOcean` and avoid explicit imports of exported names

### Writing Examples

- Explain at the top of the file what a simulation is doing
- Use Literate.jl style - let code speak for itself
- Common patterns:
  - Single-column simulations (e.g., `single_column_os_papa_simulation.jl`)
  - Regional simulations (e.g., `mediterranean_simulation_with_ecco_restoring.jl`)
  - Near-global simulations (e.g., `near_global_ocean_simulation.jl`)
- Demonstrate data wrangling, model setup, and visualization

## Important Files to Know

### Core Implementation

- `src/ClimaOcean.jl` - Main module, all exports
- `src/OceanSeaIceModels/ocean_sea_ice_model.jl` - Core coupled model type
- `src/OceanSeaIceModels/PrescribedAtmospheres.jl` - Atmosphere data structures
- `src/OceanSimulations/ocean_simulation.jl` - Ocean simulation constructor
- `src/OceanSeaIceModels/InterfaceComputations/` - Surface flux calculations

### Data Wrangling

- `src/DataWrangling/metadata.jl` - Metadata types and interface
- `src/DataWrangling/JRA55/JRA55_prescribed_atmosphere.jl` - JRA55 atmosphere
- `src/DataWrangling/ECCO/ECCO.jl` - ECCO data utilities
- `src/Bathymetry.jl` - Bathymetry regridding

### Configuration

- `Project.toml` - Package dependencies and compat bounds
- `test/runtests.jl` - Test configuration

### Examples

- `examples/single_column_os_papa_simulation.jl` - Single column at Ocean Station Papa
- `examples/mediterranean_simulation_with_ecco_restoring.jl` - Regional with restoring
- `examples/near_global_ocean_simulation.jl` - Near-global hydrostatic
- `examples/diurnal_large_eddy_simulation.jl` - LES with diurnal forcing
- `examples/idealized_single_column_simulation.jl` - Idealized single column

## Physics Domain Knowledge

### Air-Sea Fluxes

- Momentum flux (wind stress): τ = ρₐ Cᴰ |Δu| Δu
- Sensible heat flux: Qₛ = ρₐ cₚ Cᴴ |Δu| ΔT
- Latent heat flux: Qₗ = ρₐ Lᵥ Cᴱ |Δu| Δq
- Monin-Obukhov similarity theory for stability-dependent transfer coefficients
- Roughness length parameterizations (Charnock, COARE, etc.)

### Ocean-Sea Ice Coupling

- Freezing-limited ocean temperature
- Sea ice thermodynamics via ClimaSeaIce.jl
- Ice-ocean heat and salt fluxes
- Ice concentration and thickness tracking

### Radiation

- Downwelling shortwave and longwave radiation
- Albedo parameterizations (latitude-dependent, tabulated)
- Surface temperature feedback

### Datasets

- **JRA55**: 3-hourly atmospheric reanalysis (1958-present)
- **ECCO**: Monthly/daily ocean state estimates
- **GLORYS**: Copernicus ocean reanalysis (via extension)
- **ETOPO**: Global bathymetry
- **EN4**: Subsurface ocean temperature/salinity

## Common Pitfalls

1. **Type Instability**: Especially in kernel functions - check with `@code_warntype`
2. **Missing Data**: ECCO/GLORYS data needs inpainting for land points
3. **Coordinate Systems**: Be careful with longitude conventions (0-360 vs -180-180)
4. **Time Zones**: JRA55 uses UTC; be consistent with dates
5. **Grid Staggering**: Remember Oceananigans C-grid locations for fluxes
6. **Memory**: Large FieldTimeSeries can consume significant memory; use `InMemory(chunk_size)` or `OnDisk()`

## Git Workflow

- Follow ColPrac (Collaborative Practices for Community Packages)
- Create feature branches for new work
- Write descriptive commit messages
- Update tests and documentation with code changes
- Check CI passes before merging

## Helpful Resources

- ClimaOcean docs: https://clima.github.io/ClimaOceanDocumentation/stable/
- Oceananigans docs: https://clima.github.io/OceananigansDocumentation/stable/
- ClimaSeaIce: https://github.com/CliMA/ClimaSeaIce.jl
- Discussions: https://github.com/CliMA/Oceananigans.jl/discussions
- JRA55 documentation: https://jra.kishou.go.jp/JRA-55/
- ECCO data portal: https://ecco-group.org/
- ETOPO: https://www.ncei.noaa.gov/products/etopo-global-relief-model

## When Unsure

1. Check existing examples in `examples/` directory
2. Look at similar implementations in the codebase
3. Review tests for usage patterns
4. Ask in GitHub discussions
5. Check Oceananigans documentation for underlying functionality
6. Review the ClimaOcean preprint for high-level overview

## AI Assistant Behavior

- Prioritize type stability and GPU compatibility
- Follow established patterns in existing ClimaOcean code
- Add tests for new functionality
- Update exports in `ClimaOcean.jl` when adding public API
- Consider data downloading and caching implications
- Reference physics equations in comments when implementing flux calculations
- Maintain consistency with both ClimaOcean and Oceananigans coding styles

## Current Development Focus

Active areas of development to be aware of:

- **Reactant.jl integration**: XLA compilation for performance (see `ext/ClimaOceanReactantExt.jl`)
- **Copernicus/GLORYS support**: Alternative to ECCO for ocean state (see `ext/ClimaOceanCopernicusMarineExt.jl`)
- **Nonhydrostatic ocean simulations**: LES with realistic forcing
- **Improved sea ice coupling**: Better ice-ocean interactions
- **Additional datasets**: EN4, multi-year JRA55, etc.
- **Surface flux improvements**: Better roughness length parameterizations
- **Distributed computing**: MPI support for large simulations
