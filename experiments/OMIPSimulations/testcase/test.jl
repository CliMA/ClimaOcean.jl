using ClimaOcean
using Oceananigans
using Oceananigans.Units
using CUDA
using Dates

run_name = "orca_corrected_snow_kskew800_ksymm800_bih50days_2"

sim = omip_simulation(:orca;
                      arch = CPU(),
                      Nz = 70,
                      depth = 5500,
                      κ_skew = 800,
                      κ_symmetric = 800,
                      biharmonic_timescale = 50days,
                      Δt = 30minutes,
                      Δz_top = 1.5,
                      flux_configuration = :corrected,
                      with_snow = true,
                      end_date = DateTime(1958, 12, 1),
                      diagnostics = false,
                      forcing_dir = "./forcing_data",
                      with_ice_dynamics = true,
                      output_dir = "$(run_name)_run",
                      filename_prefix = run_name)

sim.stop_iteration = 10

run!(sim)

