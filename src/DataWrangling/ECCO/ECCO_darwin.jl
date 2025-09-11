using MeshArrays
using JLD2
using Glob

struct ECCO2DarwinMonthly <:SomeECCODataset end
struct ECCO4DarwinMonthly <:SomeECCODataset end

# URLs for the ECCO datasets specific to each version
const ECCO4Darwin_url = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC90/ECCO-Darwin/"
const ECCO2Darwin_url = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/ECCO-Darwin_extension/"

Base.size(data::Metadata{<:ECCO4DarwinMonthly}) = (720,  360, 50, length(data.dates))
Base.size(::Metadatum{<:ECCO4DarwinMonthly})    = (720,  360, 50, 1)
Base.size(data::Metadata{<:ECCO2DarwinMonthly}) = (1440, 720, 50, length(data.dates))
Base.size(::Metadatum{<:ECCO2DarwinMonthly})    = (1440, 720, 50, 1)

# The whole range of dates in the different dataset datasets
all_dates(::ECCO4DarwinMonthly, name) = DateTime(1992, 1, 1) : Month(1) : DateTime(2023, 3, 1)
all_dates(::ECCO2DarwinMonthly, name) = DateTime(1992, 1, 1) : Month(1) : DateTime(2025, 5, 1)

ECCO_Darwin_timestep(::Metadatum{<:ECCO4DarwinMonthly}) = 3600
ECCO_Darwin_timeref(::Metadatum{<:ECCO4DarwinMonthly}) = DateTimeProlepticGregorian(1992, 1, 1, 12, 0, 0)

ECCO_Darwin_timestep(::Metadatum{<:ECCO2DarwinMonthly}) = 1200
ECCO_Darwin_timeref(::Metadatum{<:ECCO2DarwinMonthly}) = DateTimeProlepticGregorian(1992, 1, 1, 0, 0, 0)

# File name generation specific to each Dataset dataset
function metadata_filename(metadata::Metadatum{<:Union{ECCO2DarwinMonthly, ECCO4DarwinMonthly}})
    shortname = dataset_variable_name(metadata)
    
    reference_date = ECCO_Darwin_timeref(metadata)
    timestep_size  = ECCO_Darwin_timestep(metadata)

    # Explicitly convert to Int to avoid return of a float
    iternum = Int(Dates.value((metadata.dates - reference_date) / (timestep_size * 1e3)))
    iterstr = string(iternum, pad=10)

    return shortname * "." * iterstr * ".data"
end

# Convenience functions
default_mask_value(::ECCO4DarwinMonthly) = 0
default_mask_value(::ECCO2DarwinMonthly) = 0

dataset_variable_name(data::Metadata{<:Union{ECCO2DarwinMonthly,ECCO4DarwinMonthly}}) = ECCO_darwin_dataset_variable_names[data.name]

location(::Metadata{<:Union{ECCO2DarwinMonthly, ECCO4DarwinMonthly}}) = (Center, Center, Center)

variable_is_three_dimensional(::Metadata{<:Union{ECCO2DarwinMonthly, ECCO4DarwinMonthly}}) = true

ECCO_darwin_dataset_variable_names = Dict(
    :temperature => "THETA",
    :salinity    => "SALTanom",
    :DIC => "DIC",
    :ALK => "ALK",
    :PO₄ => "PO4",
    :NO₃ => "NO3",
    :DOP => "DOP",
    :POP => "POP",
    :Fe  => "FeT",
    :Siᵀ => "SiO2",
)

ECCO_darwin_scale_factor = Dict(
    :temperature => 1,
    :salinity    => 1,
    :DIC => 1e-3,
    :ALK => 1e-3,
    :PO₄ => 1e-3,
    :NO₃ => 1e-3,
    :DOP => 1e-3,
    :POP => 1e-3,
    :Fe  => 1e-3,
    :Siᵀ => 1e-3,
)

ECCO_darwin_offset_factor = Dict(
    :temperature => 0,
    :salinity    => 35,
    :DIC => 0,
    :ALK => 0,
    :PO₄ => 0,
    :NO₃ => 0,
    :DOP => 0,
    :POP => 0,
    :Fe  => 0,
    :Siᵀ => 0,
)

function default_download_directory(::ECCO4DarwinMonthly)
    path = joinpath(download_ECCO_cache, "v4_darwin", "monthly")
    return mkpath(path)
end

function default_download_directory(::ECCO2DarwinMonthly)
    path = joinpath(download_ECCO_cache, "v2_darwin", "monthly")
    return mkpath(path)
end

metadata_url(m::Metadata{<:ECCO4DarwinMonthly}) = ECCO4Darwin_url * "monthly/" * dataset_variable_name(m) * "/" * metadata_filename(m)
metadata_url(m::Metadata{<:ECCO2DarwinMonthly}) = ECCO2Darwin_url * "monthly/" * dataset_variable_name(m) * "/" * metadata_filename(m)

ECCO_darwin_native_grid(::ECCO4DarwinMonthly) = GridSpec(ID=:LLC90)
ECCO_darwin_native_size(::ECCO4DarwinMonthly) = (90, 1170, 50)
ECCO_darwin_native_grid(::ECCO2DarwinMonthly) = GridSpec(ID=:LLC270)
ECCO_darwin_native_size(::ECCO2DarwinMonthly) = (270, 3510, 50)

longitude_interfaces(::ECCO4DarwinMonthly) = (-180, 180)

"""
    retrieve_data(metadata::Metadatum{<:ECCO4DarwinMonthly})

Read a ECCO4DarwinMonthly data file and regrid using MeshArrays on to regular lat-lon grid
"""
function retrieve_data(metadata::Metadatum{<:Union{ECCO4DarwinMonthly, ECCO2DarwinMonthly}})
    native_size = ECCO_darwin_native_size(metadata.dataset)
    native_grid = ECCO_darwin_native_grid(metadata.dataset)
    native_data = zeros(Float32, prod(native_size)) # Native LLC90 grid at precision of the input binary file

    read!(metadata_path(metadata), native_data)
    native_data = bswap.(native_data)

    meshed_data   = read(reshape(native_data, native_size...), native_grid)
    Nx, Ny, Nz, _ = size(metadata)
    data          = zeros(Float32, Nx, Ny, Nz) # Native LLC90 grid at precision of the input binary file
    
    # Download the native grid data from MeshArrays repo (only if not in already in datadeps)
    native_grid_coords = GridLoad(native_grid; option="full")

    # Check if the interpolation coefficients are already calculated
    interp_file = joinpath(dirname(metadata_path(metadata)),"native_interp_coeffs.jld2")
    if !isfile(interp_file)
        # Calculate coefficients to interpolate from native grid to regular lat-lon grid (as in the ECCO netcdf files)
        resolution_X = 360/Nx
        resolution_Y = 180/Ny

        # Regular lat-lon grid
        longitudes = longitude_interfaces(metadata.dataset)
        latitudes  = latitude_interfaces(metadata.dataset)
        lon = [i for i = longitudes[1]+resolution_X/2:resolution_X:longitudes[2]-resolution_X/2, 
                     j = latitudes[1]+resolution_Y/2:resolution_Y:latitudes[2]-resolution_Y/2]
        lat = [j for i = longitudes[1]+resolution_X/2:resolution_X:longitudes[2]-resolution_X/2, 
                     j = latitudes[1]+resolution_Y/2:resolution_Y:latitudes[2]-resolution_Y/2]
        
        # Interpolation factors for the native grid (writes out to a file "interp_file" for later use)
        coeffs = interpolation_setup(; Γ=native_grid_coords, lat, lon, filename=interp_file)
    else
        # Read the coefficients from the file that was previously calculated
        coeffs = interpolation_setup(interp_file)
    end

    # Read continental mask on the native model grid 
    native_grid_fac_center = GridLoadVar("hFacC", native_grid)

    # Interpolate each masked layer on to the native lat lon grid
    for k in 1:Nz
        i, j, c = MeshArrays.Interpolate(
            meshed_data[:, k] * land_mask(native_grid_fac_center[:, k]), 
            coeffs,
        )
        data[:, :, k] = c
    end
    
    # Reverse the z-axis
    data = reverse(data, dims=3)

    # Fill NaNs in Antarctica with zeros
    data[isnan.(data)] .= 0.f0

    # Scale data according to metadata.scale_factor
    return data .* ECCO_darwin_scale_factor[metadata.name] .+ ECCO_darwin_offset_factor[metadata.name]
end
