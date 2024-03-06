using Preferences
const iscray = parse(Bool, load_preference(Base.UUID("3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"), "iscray", "false"))
@debug "Preloading GTL library" iscray
if iscray
    import Libdl
    Libdl.dlopen_e("libmpi_gtl_cuda", Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)
end

using MPI
MPI.Init()

using Oceananigans
using Oceananigans: architecture
using ClimaOcean
using ClimaOcean.ECCO2
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity
using Oceananigans.Coriolis: ActiveCellEnstrophyConserving
using Oceananigans.Units
using ClimaOcean.OceanSimulations
using ClimaOcean.OceanSeaIceModels
using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: Radiation
using ClimaOcean.VerticalGrids: exponential_z_faces
using ClimaOcean.JRA55
using ClimaOcean.JRA55: JRA55NetCDFBackend, JRA55_prescribed_atmosphere
using Printf
using JLD2

#####
##### Global Ocean at 1/12th of a degree
#####

include("correct_oceananigans.jl") 
bathymetry_file = "bathymetry12.jld2"

# 100 vertical levels
z_faces = exponential_z_faces(Nz=100, depth=6000)

Nx = 4320
Ny = 1800
Nz = length(z_faces) - 1

bottom = zeros(Nx, Ny, 1)

arch = Distributed(GPU(), partition = Partition(8))

grid = load_balanced_regional_grid(arch; 
                                   size = (Nx, Ny, Nz), 
                                   z = z_faces, 
                                   latitude  = (-75, 75),
                                   longitude = (0, 360),
                                   halo = (7, 7, 7),
                                   interpolation_passes = 20,
                                   maximum_size = 650,
                                   minimum_depth = 10,
                                   connected_regions_allowed = 3, # We allow the oceans, the med, the bering sea
                                   bathymetry_file)
 
@show grid                                   
                              
#####
##### The Ocean component
#####                             

free_surface = SplitExplicitFreeSurface(; grid, cfl=0.7, fixed_Δt=270)

ocean = ocean_simulation(grid; Δt = 10, free_surface, closure = RiBasedVerticalDiffusivity())
model = ocean.model

# Initializing the model
set!(model, 
     T = ECCO2Metadata(:temperature), 
     S = ECCO2Metadata(:salinity))

#####
##### The atmosphere
#####

backend    = JRA55NetCDFBackend(5) 
atmosphere = JRA55_prescribed_atmosphere(arch; backend)
radiation  = Radiation()

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)

function profiled_time_steps!(model, Δt; gc_steps = 100, profiled_steps = 10)
    # initial time steps
    for step in 1:10
        time_step!(model, Δt)
    end

    nranks = MPI.Comm_size(MPI.COMM_WORLD)
    rank   = MPI.Comm_rank(MPI.COMM_WORLD)

    if rank == 0
      @info "start profiling"
    end
    elapsed_time = Float64[0]
    # Perform profiling
    for step in 1:gc_steps
        for nogc in 1:profiled_steps
            elapsed_time[1] += @elapsed begin
                    NVTX.@range "one time step" begin
                    time_step!(model, Δt)
                end
            end
        end
        GC.gc()
    end

    elapsed_time[1] = elapsed_time[1] / nranks / gc_steps / profiled_steps
    MPI.Allreduce!(elapsed_time, +, MPI.COMM_WORLD)

    if rank == 0
        file = "time_twelth.jld2"
        while isfile(file)
                file = "new_" * file
        end
        jldsave(file, elapsed_time=elapsed_time[1])
    end

    MPI.Barrier(MPI.COMM_WORLD)

    return nothing
end

profiled_time_steps!(coupled_model, 0.1)