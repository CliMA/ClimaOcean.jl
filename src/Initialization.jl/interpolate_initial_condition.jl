using Dierckx 
using Oceananigans.Fields: interpolate

function interp_initial_condition(old_vector, Nx_old, Ny_old, Nx_new, Ny_new, Nz; interpolation_method = LinearInterpolation())
    new_vector = zeros(Nx_new, Ny_new, Nz)
    passes = 1
    for k in 1:Nz
        new_vector[:, :, k] = interpolate_one_level_in_passes(old_vector[:, :, k], Nx_old, Ny_old, Nx_new, Ny_new, passes; interpolation_method)
    end
    return new_vector
end

function check_zeros(bathymetry, z_faces, old_array; ensure_positivity = false)
    Nx, Ny, Nz = size(old_array)

    zc = z_center_from_face(z_faces)

    if ensure_positivity 
        old_array[old_array .< 0] .= 0.0
    end
    
    count = 0
    for k in 1:Nz
        @info "we are at k = $k"
        for i in 1:Nx, j in 1:Ny
            if (zc[k] > bathymetry[i, j]) && (old_array[i, j, k] == 0)
                count += 1
                @info "There is a zero! at $i, $j with z = $(zc[k]) and bat = $(bathymetry[i, j])"
            end
        end
    end

    return count
end

function substitute_min_value(old_array, min_val) 
    old_array[old_array.< min_val] .= min_val 
    return old_array
end

function read_and_interpolate_var(folder, month, new_size; interpolation_method = SplineInterpolation())

    path_to_file = "/Volumes/files/Version4/Release4/interp_monthly/$folder/1992/$(folder)_1992_$(month).nc"

    @info path_to_file
    interp_old = ncread(path_to_file, folder)

    interp_old = reverse(interp_old[:, :, :, 1], dims = 3)
    
    passes     = 200

    Nx, Ny, Nz = size(interp_old)

    for pass = 1:passes
        @info "pass $pass in propagate step"
        interp_old = propagate_step(interp_old, Nz)
    end

    interp_new  = interp_initial_condition(interp_old, Nx, Ny, new_size...; interpolation_method)

    return interp_new
end

function read_and_interpolate_var_from_array(old_array, new_size)

    passes = 200

    interp_old = deepcopy(old_array)

    Nx, Ny, Nz = size(interp_old)

    for pass = 1:passes
        @info "pass $pass in propagate step"
        interp_old = propagate_step(interp_old, Nz)
    end

    interp_new = interp_initial_condition(interp_old, Nx, Ny, new_size...)

    return interp_new
end

months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]

function calculate_initial_conditions(months, interpolation_method = SplineInterpolation())

    bathymetry = jldopen("/Users/simonesilvestri/global-solution/bathymetries/bathymetry-one-degree.jld2")["bathymetry"]
    z_faces    = jldopen("/Users/simonesilvestri/development/GlobalShenanigans.jl/data/zgrid.jld2")["z"][5:end-4]

    arrays = [zeros(360, 150, 48), zeros(360, 150, 48)]
    for month in months
        @info "month number $month"
        for (var, folder) in zip(arrays, ["THETA", "SALT"])
            @info "interpolating $folder"
            var .= read_and_interpolate_var(folder, month, (360, 150, 48); interpolation_method)
            ensure_positivity = folder == "SALT"
            check_zeros(bathymetry, z_faces, var; ensure_positivity)
        end

        jldsave("initial_condition_month_$(month).jld2", T = arrays[1], S = arrays[2])
    end

end

    
