include("runtests_setup.jl")

using MPI

MPI.Init()
atexit(MPI.Finalize)  

using Oceananigans.Units
using Oceananigans.DistributedComputations
using Oceananigans.Architectures: on_architecture
using Dates
using ClimaSeaIce
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium

archs = [Distributed(CPU(); partition = Partition(y = DistributedComputations.Equal()), synchronized_communication=true),
         Distributed(GPU(); partition = Partition(y = DistributedComputations.Equal()), synchronized_communication=true)]

function analytical_immersed_tripolar_grid(underlying_grid::TripolarGrid;
                                           radius = 5, # degrees
                                           active_cells_map = false)

    λp = underlying_grid.conformal_mapping.first_pole_longitude
    φp = underlying_grid.conformal_mapping.north_poles_latitude
    φm = underlying_grid.conformal_mapping.southernmost_latitude

    Lz = underlying_grid.Lz

    # We need a bottom height field that ``masks'' the singularities
    bottom_height(λ, φ) = ((abs(λ - λp) < radius)       & (abs(φp - φ) < radius)) |
                          ((abs(λ - λp - 180) < radius) & (abs(φp - φ) < radius)) | (φ < φm) ? 0 : - Lz

    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height); active_cells_map)

    return grid
end

@testset "Distributed Models" begin
    for arch in archs
        @info "Testing on architecture: $arch"

        Nx, Ny, Nz = 100, 100, 30
        underlying_grid = TripolarGrid(arch; size = (Nx, Ny, Nz), z = (-6000, 0), halo = (7, 7, 4))
        grid = analytical_immersed_tripolar_grid(underlying_grid; active_cells_map=true)
        free_surface = SplitExplicitFreeSurface(grid; cfl=0.7, fixed_Δt=10minutes)

        Δt = 10

        ocean = ocean_simulation(grid; Δt, free_surface, timestepper = :SplitRungeKutta3)
        sea_ice = sea_ice_simulation(grid, ocean; advection=WENO(order=7))

        set!(sea_ice.model, h=Metadatum(:sea_ice_thickness;     dataset=ECCO4Monthly()),
                            ℵ=Metadatum(:sea_ice_concentration; dataset=ECCO4Monthly()))

        radiation  = Radiation(arch)
        atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(5))

        coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

        stop_iteration = 5
        simulation = Simulation(coupled_model; Δt, verbose=false, stop_time=stop_iteration * Δt)

        run!(simulation)

        @test coupled_model.clock.iteration == stop_iteration
    end
end
