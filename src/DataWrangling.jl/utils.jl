using GlobalShenanigans: set!
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Dierckx 

struct SplineInterpolation end
struct LinearInterpolation end
struct SpectralInterpolation end 

@inline function interpolate_one_level(old_array, old_grid, new_grid, loc; interpolation_method = LinearInterpolation())
    old_field = Field{Center, Center, Center}(old_grid)
    set!(old_field, old_array)
    fill_halo_regions!(old_field)
    Nx_new, Ny_new = size(new_grid)[[1, 2]]
    new_array = zeros(Nx_new, Ny_new)

    if interpolation_method isa SplineInterpolation
        data = old_field.data.parent[:, :, 1]
        spline = Spline2D(parent(old_grid.xᶜᵃᵃ), parent(old_grid.yᵃᶜᵃ), data)
        
        for i in 1:Nx_new, j in 1:Ny_new
            if loc == Face
                new_array[i, j] = evaluate(spline, new_grid.xᶜᵃᵃ[i], new_grid.yᵃᶠᵃ[j])
            else
                new_array[i, j] = evaluate(spline, new_grid.xᶜᵃᵃ[i], new_grid.yᵃᶠᵃ[j])
            end
        end
    else
        for i in 1:Nx_new, j in 1:Ny_new
            new_array[i, j] = interpolate(old_field, new_grid.xᶜᵃᵃ[i], new_grid.yᵃᶠᵃ[j], old_grid.zᵃᵃᶜ[1])
        end
    end

    return new_array
end

@inline function propagate_step(vec, Nz_mine)
    Nx, Ny, Nz = size(vec)
    vec2 = deepcopy(vec)
    @inbounds begin
        for k in 1:Nz_mine
            for i in 1:Nx, j in 1:Ny
                neigh_west  = i == 1  ? vec[Nx, j, k] : vec[i - 1, j, k]
                neigh_north = j == Ny ? 0.0 : vec[i, j + 1, k]
                neigh_east  = i == Nx ? vec[1, j, k]  : vec[i + 1, j, k]
                neigh_south = j == 1  ? 0.0 : vec[i, j - 1, k]
                neigbours = [neigh_west, neigh_east, neigh_north, neigh_south]
                non_null  = Int.(neigbours .!= 0)
                if (vec[i, j, k] == 0) & (sum(non_null) > 0)
                    vec2[i, j, k] = dot(non_null, neigbours) / sum(non_null) 
                end
            end
        end
    end
    return vec2
end

@inline function horizonthal_filter(vec, iᵢ, iₑ, jᵢ, jₑ, Nz_mine)
    vec2 = deepcopy(vec)
    for k in 1:Nz_mine
        for i in iᵢ:iₑ, j in jᵢ:jₑ
            neigbours = [vec[i, j, k], vec[i + 1, j, k], vec[i - 1, j, k], vec[i, j + 1, k], vec[i, j - 1, k]]
            non_null  = Int.(neigbours .!= 0)
            if sum(non_null) > 0
                vec2[i, j, k] = sum(neigbours) / sum(non_null)
            end
        end
    end
    return vec2
end

@inline function filter_boundary_values(old_vector, Nz_mine, passes)

    Nx, Ny, Nz = size(old_vector)
    new_vector = deepcopy(old_vector)
    grid = RectilinearGrid(size = (Nx, Ny), 
                              y = (-1, 1),
                              x = (-1, 1),
                       topology = (Periodic, Bounded, Flat), halo = (20, 20))

    field = Field{Center, Center, Nothing}(grid)
    for k = 1:Nz_mine
        set!(field, old_vector[:, :, k])
        fill_halo_regions!(field, architecture(grid))
        @info "level $k in filter boundary"
        for pass = 1:passes
            field = horizonthal_filter(field,   -15,    30, 2, Ny-1, 1)
            field = horizonthal_filter(field, Nx-30, Nx+15, 2, Ny-1, 1)

            new_vector[:, :, k] = interior(field)[:, :, 1]
        end
    end
    return new_vector
end

@inline function copy_boundary_values(old_vector, Nz_mine)

    Nx, Ny, Nz = size(old_vector)
    new_vector = deepcopy(old_vector)
    grid = RectilinearGrid(size = (Nx, Ny), 
                              y = (-1, 1),
                              x = (-1, 1),
                       topology = (Periodic, Bounded, Flat), halo = (20, 20))

    field = Field{Center, Center, Nothing}(grid)
    for k = 1:Nz_mine
        set!(field, old_vector[:, :, k])
        fill_halo_regions!(field, architecture(grid))
        @info "level $k in filter boundary"
        for pass = 1:100
            field = horizonthal_filter(field,   -15,    30, 2, Ny-1, 1)
            field = horizonthal_filter(field, Nx-30, Nx+15, 2, Ny-1, 1)

            new_vector[:, :, k] = interior(field)[:, :, 1]
        end
    end
    return new_vector
end

function interpolate_one_level_in_passes(array_old, Nxₒ, Nyₒ, Nxₙ, Nyₙ, passes; interpolation_method = LinearInterpolation())

    ΔNx = floor((Nxₒ - Nxₙ) / passes)
    ΔNy = floor((Nyₒ - Nyₙ) / passes)

    Nx = deepcopy(Nxₒ)
    Ny = deepcopy(Nyₒ)

    @assert Nxₒ == Nxₙ + passes * ΔNx
    @assert Nyₒ == Nyₙ + passes * ΔNy
    
    array = deepcopy(array_old)
    latitude = 75.0

    for pass = 1:passes
        array_full = deepcopy(array)
        Nxₒ = Nx
        Nyₒ = Ny
        Nx -= Int(ΔNx) 
        Ny -= Int(ΔNy)
        if pass == 1
            oldlat = 89.9999999999999999
        else
            oldlat = latitude
        end
        old_grid = RectilinearGrid(size = (Nxₒ, Nyₒ), y = (-oldlat,   oldlat),   x = (-180, 180), topology = (Periodic, Bounded, Flat))
        new_grid = RectilinearGrid(size = (Nx,  Ny ), y = (-latitude, latitude), x = (-180, 180), topology = (Periodic, Bounded, Flat))
    
        @show Nxₒ, Nyₒ, Nx, Ny, pass
        array = interpolate_one_level(array_full, old_grid, new_grid, Center; interpolation_method)
    end

    return array
end
