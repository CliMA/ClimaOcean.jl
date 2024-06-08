using Oceananigans.Grids: architecture, location, node, with_halo
using ClimaOcean.DataWrangling.ECCO: ecco_mask

function ecco_immersed_grid(metadata)
    mask = ecco_mask(metadata)
    grid = with_halo((3, 3, 3), mask.grid)
    
    Nx, Ny, Nz = size(grid)
    bathymetry = zeros(Nx, Ny)
    
    for i in 1:Nx, j in 1:Ny
        for k in Nz:-1:1
            if mask[i, j, k]
                bathymetry[i, j] = 0.5 .* (grid.zᵃᵃᶜ[k] + grid.zᵃᵃᶠ[k+1])
                break;
            elseif k == 1
                bathymetry[i, j] = - grid.Lz - 1
            end
        end
    end

    return ImmersedBoundaryGrid(grid, GridFittedBottom(bathymetry))
end