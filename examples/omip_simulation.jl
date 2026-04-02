# # OMIP simulation
#
# Global ocean–sea ice simulation following the OMIP protocol (Griffies et al. 2016).
# Uses `omip_simulation` for turnkey setup: 0.5° tripolar grid with 100 vertical
# levels, JRA55 atmospheric forcing, EN4 salinity restoring, and ECCO4 sea ice
# initial conditions. OMIP-protocol diagnostics are attached automatically.

using ClimaOcean
using ClimaOcean.OMIPConfigurations
using Oceananigans.Units
using CUDA

# ### Build the simulation
#
# `omip_simulation(:half_degree)` constructs the full coupled model:
# - Ocean: 720×360×100 TripolarGrid with CATKE, GM, biharmonic viscosity
# - Sea ice: prognostic ClimaSeaIce model
# - Atmosphere: JRA55 prescribed reanalysis with rivers and icebergs
# - Diagnostics: surface fields, 3-D fields, zonal means, checkpoints

arch = GPU()

simulation = omip_simulation(:half_degree;
    arch,
    forcing_dir = "forcing_data",
    restoring_dir = "climatology",
    Δt = 30minutes,
    stop_time = 2 * 365days,
    output_dir = "omip_output",
    filename_prefix = "halfdegree",
)

# The simulation comes with OMIP diagnostics pre-attached:
# - `halfdegree_surface_partN.jld2`   — SST, SSS, SSH, MLD, fluxes, sea ice (2D, ~15-day averages)
# - `halfdegree_fields_partN.jld2`    — T, S, u, v, w, TKE (3D, ~15-day averages)
# - `halfdegree_scalars_partN.jld2`   — zonal means, global means
# - `halfdegree_checkpoint_*`         — full restart files every 90 days
#
# Files split yearly (24 time instances per file).

# ### Run!

run!(simulation)

# ### Diagnostic report

simulation_report(simulation.model.ocean, filename = "omip_report.png")

# ![](omip_report.png)

# ### Custom diagnostics example
#
# You can also build the simulation without diagnostics and attach them manually:
#
# ```julia
# sim = omip_simulation(:half_degree; arch=GPU(), diagnostics=false)
#
# # Attach with custom intervals
# add_omip_diagnostics!(sim;
#     surface_averaging_interval = 30days,
#     field_averaging_interval = 30days,
#     file_splitting_interval = 365days,
# )
#
# run!(sim)
# ```

# ### ORCA grid variant
#
# For the ORCA (NEMO-compatible) grid:
#
# ```julia
# sim = omip_simulation(:orca; arch=GPU(), stop_time=300*365days)
# run!(sim)
# ```
