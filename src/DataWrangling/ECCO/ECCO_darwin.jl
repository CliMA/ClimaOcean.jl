using MeshArrays

# URLs for the ECCO datasets specific to each version
const ECCO4Darwin_url = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC90/ECCO-Darwin/"

Base.size(::Metadata{<:AnyDateTime, <:ECCO4DarwinMonthly}) = (720,  360, 50, 1)
Base.size(data::Metadata{<:Any, <:ECCO4DarwinMonthly}) = (720,  360, 50, length(data.dates))

# The whole range of dates in the different dataset datasets
all_dates(::ECCO4DarwinMonthly, name) = DateTime(1992, 1, 1) : Month(1) : DateTime(2023, 12, 1)

# File name generation specific to each Dataset dataset
function metadata_filename(metadata::ECCOMetadata{<:AnyDateTime, <:ECCO4DarwinMonthly})
    shortname = short_name(metadata)
    
    reference_date = DateTimeProlepticGregorian(1992, 1, 1, 12, 0, 0)
    timestep_size  = 3600

    iternum = Dates.value((metadata.dates - reference_date) / (timestep_size * 1e3))
    iterstr = string(iternum, pad=10)

    return shortname * "." * iterstr * ".data"
end

# Convenience functions
short_name(data::Metadata{<:Any, <:ECCO4DarwinMonthly}) = ECCO_darwin_short_names[data.name]

location(::ECCOMetadata{<:Any, <:ECCO4DarwinMonthly}) = (Center, Center, Center)

variable_is_three_dimensional(::ECCOMetadata{<:Any, <:ECCO4DarwinMonthly}) = true

ECCO_darwin_short_names = Dict(
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
    :DIC => 1e-3,
    :ALK => 1e-3,
    :PO₄ => 1e-3,
    :NO₃ => 1e-3,
    :DOP => 1e-3,
    :POP => 1e-3,
    :Fe  => 1e-3,
    :Siᵀ => 1e-3,
)

metadata_url(m::Metadata{<:Any, <:ECCO4DarwinMonthly}) = ECCO4Darwin_url * "monthly/" * short_name(m) * "/" * metadata_filename(m)

ECCO_darwin_native_grid(::ECCO4DarwinMonthly) = GridSpec("LatLonCap", MeshArrays.GRID_LLC90)
ECCO_darwin_native_size(::ECCO4DarwinMonthly) = (90, 1170, 50)

"""
    ECCO_darwin_model_data(metafile)

Read a ECCO4DarwinMonthly data file and regrid using MeshArrays on to regular lat-lon grid
"""
function ECCO_darwin_model_data(metadata, path)
    native_size = ECCO_darwin_native_size(metadata.dataset)
    native_grid = ECCO_darwin_native_grid(metadata.dataset)
    native_data = zeros(Float32, prod(native_size)) # Native LLC90 grid at precision of the input binary file

    read!(path, native_data)
    native_data = bswap.(native_data)

    meshed_data   = read(reshape(native_data, native_size...), native_grid)
    Nx, Ny, Nz, _ = size(metadata)
    data          = zeros(Float32, Nx, Ny, Nz) # Native LLC90 grid at precision of the input binary file

    # Download the native grid data from MeshArrays repo (only if not in already in datadeps)
    native_grid_coords = GridLoad(native_grid; option="full")

    # Calculate coefficients to interpolate from native grid to regular lat-lon grid (as in the ECCO netcdf files)
    resolution_X = 360/Nx
    resolution_Y = 180/Ny

    lon = [i for i = -180+resolution_X/2:resolution_X:180-resolution_X/2, 
                 j = -90+resolution_Y/2:resolution_Y:90-resolution_Y/2]
    lat = [j for i = -180+resolution_X/2:resolution_X:180-resolution_X/2, 
                 j = -90+resolution_Y/2:resolution_Y:90-resolution_Y/2]
    
    (f, i, j, w, _, _, _) = InterpolationFactors(native_grid_coords, vec(lon), vec(lat))
    coeffs = (lon=lon, lat=lat, f=f, i=i, j=j, w=w)

    # Read continental mask on the native model grid 
    native_grid_fac_center = GridLoadVar("hFacC", native_grid)

    # Interpolate each masked layer
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
    return data .* ECCO_darwin_scale_factor[metadata.name]
end

retrieve_data(metadata::ECCOMetadata{<:Any, <:ECCO4DarwinMonthly}, path) = ECCO_darwin_model_data(metadata, path)
