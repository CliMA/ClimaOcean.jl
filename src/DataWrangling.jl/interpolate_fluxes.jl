using Dierckx
using NetCDF

# grab interpolate_fluxes
function read_and_interpolate_flux(var, folder, yr, dims, degree; passes = 1, filter_passes = 5)

    lat = 75

    Nxₙ  = Int(360 / degree)
    Nyₙ  = Int(2lat/ degree)

    flux = zeros(Nxₙ, Nyₙ, 12)
    months  = ("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

    for (idx, month) in enumerate(months)
        println("$var, $yr, $month")
        
        path_to_file = "/Volumes/files/Version4/Release4/interp_monthly/$folder/$yr/$(folder)_$(yr)_$(month).nc"

        @show path_to_file
        interp_old = ncread(path_to_file, folder)
        if dims == 3
            data = interp_old[:, :, 1:1, 1] 
        else
            data = interp_old
        end

        for i in 1:100
            @info "propagating step $i"
            data = propagate_step(data, 1)
        end

        Nxₒ  = 720
        Nyₒ  = 360
        ΔNx = floor((Nxₒ - Nxₙ) / passes)
        ΔNy = floor((Nyₒ - Nyₙ) / passes)

        @assert Nxₒ == Nxₙ + passes * ΔNx
        @assert Nyₒ == Nyₙ + passes * ΔNy
        
        Nx = deepcopy(Nxₒ)
        Ny = deepcopy(Nyₒ)

        for pass = 1:passes
            data_full = deepcopy(data)

            Nxₒ = Nx
            Nyₒ = Ny
            Nx -= Int(ΔNx) 
            Ny -= Int(ΔNy)
            if pass == 1
                oldlat = 89.9999999999999999
            else
                oldlat = lat
            end
            old_grid = RectilinearGrid(size = (Nxₒ, Nyₒ), y = (-oldlat,   oldlat), x = (-180, 180), topology = (Periodic, Bounded, Flat))
            new_grid = RectilinearGrid(size = (Nx,  Ny ), y = (-lat, lat),         x = (-180, 180), topology = (Periodic, Bounded, Flat))
            
            @show Nxₒ, Nyₒ, Nx, Ny, pass
            data = interpolate_one_level(data_full, old_grid, new_grid, Center)

            for i in 1:filter_passes
                @info "horizontal filter pass $i"
                data = horizonthal_filter(data, 2, Nx-1, 2, Ny-1, 1)
            end
    
        end

        flux[:, :, idx] .= data
    end

    @show size(flux)
    return flux
end

function save_fluxes(year, degree)
    vars    = ("tau_x", "tau_y", "Q_flux", "S_flux", "temp", "salt")
    folders = ("oceTAUE", "oceTAUN", "oceQnet", "oceFWflx", "THETA", "SALT")

    τˣ = read_and_interpolate_flux(vars[1], folders[1], year, 2, degree, passes = 4, filter_passes = 8)
    τʸ = read_and_interpolate_flux(vars[2], folders[2], year, 2, degree, passes = 4, filter_passes = 8)
    Qᶠ = read_and_interpolate_flux(vars[3], folders[3], year, 2, degree, passes = 4, filter_passes = 8)
    Sᶠ = read_and_interpolate_flux(vars[4], folders[4], year, 2, degree, passes = 4, filter_passes = 8)
    Tₛ = read_and_interpolate_flux(vars[5], folders[5], year, 3, degree, passes = 4, filter_passes = 8)
    Sₛ = read_and_interpolate_flux(vars[6], folders[6], year, 3, degree, passes = 4, filter_passes = 8)
    
    jldsave("boundary_conditions_$(degree).jld2"; τˣ, τʸ, Qᶠ, Sᶠ, Tₛ, Sₛ)
end
