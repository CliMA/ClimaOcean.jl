using MeshArrays
using JLD2
using Glob

struct ECCO2DarwinMonthly <:ECCODataset end
struct ECCO4DarwinMonthly <:ECCODataset end

# URLs for the ECCO datasets specific to each version
const ECCO4Darwin_url = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC90/ECCO-Darwin/"
const ECCO2Darwin_url = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/ECCO-Darwin_extension/"

Base.size(data::Metadata{<:ECCO4DarwinMonthly}) = (720,  360, 50, length(data.dates))
Base.size(::Metadatum{<:ECCO4DarwinMonthly})    = (720,  360, 50, 1)
Base.size(data::Metadata{<:ECCO2DarwinMonthly}) = (1440, 720, 50, length(data.dates))
Base.size(::Metadatum{<:ECCO2DarwinMonthly})    = (1440, 720, 50, 1)

metadata_time_step(::Metadatum{<:ECCO4DarwinMonthly}) = 3600
metadata_epoch(::Metadatum{<:ECCO4DarwinMonthly}) = DateTime(1992, 1, 1, 12, 0, 0)

metadata_time_step(::Metadatum{<:ECCO2DarwinMonthly}) = 1200
metadata_epoch(::Metadatum{<:ECCO2DarwinMonthly}) = DateTime(1992, 1, 1, 0, 0, 0)

# The whole range of dates in the different dataset datasets
all_dates(dataset::ECCO4DarwinMonthly, name) = metadata_epoch(dataset) : Month(1) : DateTime(2023, 3, 1)
all_dates(dataset::ECCO2DarwinMonthly, name) = metadata_epoch(dataset) : Month(1) : DateTime(2025, 5, 1)

# File name generation specific to each Dataset dataset
"""
    metadata_filename(metadata::Metadatum{<:Union{ECCO2DarwinMonthly, ECCO4DarwinMonthly}})

Generate the filename for a given ECCO Darwin dataset and date.

The filename is constructed using the dataset variable name, and the iteration number is calculated
from the date and epoch.
"""
function metadata_filename(metadata::Metadatum{<:Union{ECCO2DarwinMonthly, ECCO4DarwinMonthly}})
    shortname = dataset_variable_name(metadata)

    reference_date = metadata_epoch(metadata)
    timestep_size  = metadata_time_step(metadata)

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
    :temperature                    => "THETA",
    :salinity                       => "SALTanom",
    :dissolved_inorganic_carbon     => "DIC",
    :alkalinity                     => "ALK",
    :phosphate                      => "PO4",
    :nitrate                        => "NO3",
    :dissolved_organic_phosphorus   => "DOP",
    :particulate_organic_phosphorus => "POP",
    :dissolved_iron                 => "FeT",
    :dissolved_silicate             => "SiO2",
    :dissolved_oxygen               => "O2",
)

"""
    concentration_units(metadatum::Metadatum{<:Union{ECCO2DarwinMonthly, ECCO4DarwinMonthly}}) 

Set up conversion from the ECCODarwin output data to standard units
  -  salinity = SALTanom + 35
  -  biogeochemical tracer concentrations are in uL => umol/L in the output files from Darwin
"""
function concentration_units(metadatum::Metadatum{<:Union{ECCO2DarwinMonthly, ECCO4DarwinMonthly}}) 
    if dataset_variable_name(metadatum) == "SALTanom"
        return GramPerKilogramMinus35()
    elseif dataset_variable_name(metadatum) != "THETA"
        return MicromolePerLiter()
    else
        return nothing
    end
end

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

# Functions for reading the ECCO binary files using MeshArrays
binary_data_grid(::ECCO4DarwinMonthly) = GridSpec(ID=:LLC90)
binary_data_size(::ECCO4DarwinMonthly) = (90, 1170, 50)
binary_data_grid(::ECCO2DarwinMonthly) = GridSpec(ID=:LLC270)
binary_data_size(::ECCO2DarwinMonthly) = (270, 3510, 50)

longitude_interfaces(::ECCO4DarwinMonthly) = (-180, 180)

"""
    retrieve_data(metadata::Metadatum{<:ECCO4DarwinMonthly})

Read a ECCO4DarwinMonthly data file and regrid using MeshArrays on to regular lat-lon grid
"""
function retrieve_data(metadata::Metadatum{<:Union{ECCO4DarwinMonthly, ECCO2DarwinMonthly}})
    native_size = binary_data_size(metadata.dataset)
    native_grid = binary_data_grid(metadata.dataset)
    native_data = zeros(Float32, prod(native_size)) # Native LLC grid at precision of the input binary file

    read!(metadata_path(metadata), native_data)
    native_data = bswap.(native_data)

    meshed_data   = read(reshape(native_data, native_size...), native_grid)
    Nx, Ny, Nz, _ = size(metadata)
    data          = zeros(Float32, Nx, Ny, Nz) # Native LLC grid at precision of the input binary file
    mask          = zeros(Float32, Nx, Ny, Nz)

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
            meshed_data[:, k], 
            coeffs,
        )
        data[:, :, k] = c
        i, j, c = MeshArrays.Interpolate(
            land_mask(native_grid_fac_center[:, k]), 
            coeffs,
        )
        mask[:, :, k] = c
    end
    
    # Reverse the z-axis
    data = reverse(data, dims=3)
    mask = reverse(mask, dims=3)

    # Fill NaNs in Antarctica with zeros
    data[isnan.(data)] .= default_mask_value(metadata.dataset)

    return data .* mask
end
