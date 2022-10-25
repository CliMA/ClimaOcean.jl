using GlobalShenanigans: set!
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Dierckx 

"""
struct specifying the interpolation methods 
"""
struct SplineInterpolation end
struct LinearInterpolation end
struct SpectralInterpolation end 


@inline function z_face_from_center(z_center)
    z_faces = zeros(length(z_center) + 1)
    for k = 2:length(z_center) + 1
        z_faces[k] = z_faces[k-1] + 2*(z_center[k-1] - z_faces[k-1])
    end
    
    return reverse(z_faces)
end

@inline z_center_from_face(z_faces) = 0.5 .* (z_faces[2:end] .+ z_faces[1:end-1])

"""
    interpolate_one_level(old_array, old_grid, new_grid, loc; interpolation_method = LinearInterpolation())

interpolates an `old_array` living on `old_grid` to a `new_grid`

Arguments
=========

* `old_array`: array to be interpolated

* `old_grid`: grid of `size(old_array)`

* `new_grid`: grid to interpolate to

KeywordsArguments
=================

* `interpolation_method`: either `SplineInterpolation` or `LinearInterpolation`

"""
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

"""
    propagate_step(vec, Nz_mine)

propagate the field in immersed location
"""
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

"""
    horizonthal_filter(vec, iᵢ, iₑ, jᵢ, jₑ, Nz_mine)

horizontal box filter of vector `vec` from `iᵢ`, `iₑ` and `jᵢ`, `jₑ`.

"""
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

"""
    interpolate_one_level_in_passes(array_old, Nxₒ, Nyₒ, Nxₙ, Nyₙ, passes; interpolation_method = LinearInterpolation())

applies passes to one_level interpolation
""" 
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