using ImageMagick
using PyCall
using FFTW
using FastSphericalHarmonics
using Dierckx

"""
    interpolate_bathymetry_from_file(resolution, latitude; 
                                     filename = "data/bathymetry-ice-21600x10800.jld2", 
                                     interpolation_method = LinearInterpolation(), 
                                     minimum_depth = 6)

Generate a latitude longitude bathymetry array a

"""
abstract type AbstractInterpolation end 

struct SplineInterpolation <: AbstractInterpolation end

struct LinearInterpolation <: AbstractInterpolation 
    passes :: Int
end

struct SpectralInterpolation{F, C} <: AbstractInterpolation 
    filter_function :: F
    spectral_coeff  :: C
end 

LinearInterpolation(; passes = 5) = LinearInterpolation(passes)
SpectralInterpolation(; filter_func = (l) -> exp(-l * (l+1)/ 180 / 240), spectral_coeff = nothing) = 
            SpectralInterpolation(filter_func, spectral_coeff)


function interpolate_bathymetry_from_file(resolution, latitude; 
                                          filename = "data/bathymetry-ice-21600x10800.jld2", 
                                          interpolation_method = LinearInterpolation()
                                          minimum_depth = 6)

    file = jldopen(filename)
    bathy_old = Float64.(file["bathymetry"])

    Nx  = Int(360 / resolution)
    Ny  = Int(2latitude / resolution)

    Nn = (Nx, Ny)  
    bathy = twodimensional_interpolation(bathy_old, interpolation_method, Nn)

    if interpolation_method isa SpectralInterpolation
        if calc_coeff isa nothing
            spectral_coeff = etopo1_to_spherical_harmonics(bathy_old, Nyₒ)
        else 
            spectral_coeff = calc_coeff
        end

        bathy = bathymetry_from_etopo1(Nxₙ, Nyₙ, spher_harm_coeff, filter_func)
    else 
        bathy = interpolate_one_level_in_passes(bathy_old, Nxₒ, Nyₒ, Nxₙ, Nyₙ, passes; interpolation_method)
    end

    # apparently bathymetry is reversed in the longitude direction, therefore we have to swap it
    bathy = reverse(bathy, dims = 2)
    
    bathy[bathy .> 0] .= ABOVE_SEA_LEVEL

    fixed_bathymetry = remove_connected_regions(bathy)

    fixed_bathymetry[bathy .> - minimum_depth] .= ABOVE_SEA_LEVEL
    
    return fixed_bathymetry
end

function etopo1_to_spherical_harmonics(etopo1, Nmax)

    lmax = Nmax - 2

    # latitude interpolate and transpose
    etopo1_center = 0.5 *(etopo1[:, 2:end] .+ etopo1[:, 1:end-1])
    etopo1_center = etopo1_center'

    # Drop the 360 point
    etopo1_center = etopo1_center[:, 1:end-1]

    # longitude interpolate
    fft_etopo1_center = rfft(etopo1_center, 2)
    fft_etopo1_center = fft_etopo1_center[:, 1:end-1]
    etopo1_interp     = irfft(fft_etopo1_center, 2lmax +1, 2)

    # Sperical harmonic filtering
    return sph_transform(etopo1_interp)
end

function bathymetry_from_etopo1(Nφ, Nλ, spher_harm_coeff, filter)

    lmax_interp = Nφ - 1

    spher_harm_coeff_filter = zeros(lmax_interp + 1, 2lmax_interp +1)
    for l = 0:lmax_interp, m = -l:l
        spher_harm_coeff_filter[sph_mode(l,m)] = spher_harm_coeff[sph_mode(l, m)] * filter(l)
    end

    etopo1_filter = sph_evaluate(spher_harm_coeff_filter) # dimensions are Nφ, 2Nφ - 1 

    # longitude interpolate 2
    fft_etopo1_filter = rfft(etopo1_filter, 2) # dimensions are Nφ, Nφ 

    mmax_interp = Nλ ÷ 2 + 1

    fft_etopo1_interp = zeros(Complex{Float64}, Nφ, mmax_interp)

    if mmax_interp <= size(fft_etopo1_filter, 2)
        fft_etopo1_interp .= fft_etopo1_filter[:, 1:mmax_interp]
    else
        fft_etopo1_interp[:, 1:size(fft_etopo1_filter, 2)] .= fft_etopo1_filter
    end

    etopo1_final = irfft(fft_etopo1_interp, Nλ, 2)

    return etopo1_final
end

function remove_connected_regions(bat)

    bathymetry = deepcopy(bat)
    batneg     = deepcopy(bathymetry)

    batneg[batneg.>0] .= 0
    batneg[batneg.<0] .= 1

    labels = sckikitimage.label(batneg)
    try
        total_elements = zeros(maximum(labels))

        for i in 1:lastindex(total_elements)
            total_elements[i] = sum(labels[labels.==i])
        end

        ocean_idx = findfirst(x -> x == maximum(x), total_elements)
        second_maximum = maximum(filter((x) -> x != total_elements[ocean_idx], total_elements))

        bering_sea_idx = findfirst(x -> x == second_maximum, total_elements)

        labels = Float64.(labels)
        labels[labels.==0] .= NaN

        for i in 1:length(total_elements)
            if (i != ocean_idx) && (i != bering_sea_idx)
                labels[labels.==i] .= NaN
            end
        end

        bathymetry .+= labels
        bathymetry[isnan.(bathymetry)] .= ABOVE_SEA_LEVEL
    catch err
        println("this is the error $err")
    end
    return bathymetry
end

function write_bathymetry_to_file(prefix, bathy, lat)
    Nxₙ, Nyₙ = size(bathy)
    output_prefix = prefix * "-$(Int(Nxₙ))x$(Int(Nyₙ))-latitude-$(lat)"
    jldsave(output_prefix * ".jld2", bathymetry = bathy)
end
