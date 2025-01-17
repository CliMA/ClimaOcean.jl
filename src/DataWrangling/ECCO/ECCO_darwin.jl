using MeshArrays

Base.size(::ECCOMetadata{<:AbstractCFDateTime, <:ECCO4DarwinMonthly})   = (720, 360, 50, 1)
Base.size(data::ECCOMetadata{<:Any, <:ECCO4DarwinMonthly}) = (720, 360, 50, length(data.dates))

# File name generation specific to each Dataset version
function metadata_filename(metadata::ECCOMetadata{<:AbstractCFDateTime, <:ECCO4DarwinMonthly})
    shortname = short_name(metadata)
    
    reference_date = DateTimeProlepticGregorian(1992, 1, 1, 12, 0, 0)
    timestep_size  = 3600

    iternum = Dates.value((metadata.dates - reference_date) / (timestep_size * 1e3))
    iterstr = string(iternum, pad=10)

    return shortname * "." * iterstr * ".data"
end

# Convenience functions
short_name(data::ECCOMetadata{<:Any, <:ECCO4DarwinMonthly}) = ECCO_darwin_short_names[data.name]

metadata_url(prefix, m::ECCOMetadata{<:Any, <:ECCO4DarwinMonthly}) = prefix * "/" * short_name(m) * "/" * metadata_filename(m)

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

# URLs for the ECCO datasets specific to each version
urls(::ECCOMetadata{<:Any, <:ECCO4DarwinMonthly}) = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC90/ECCO-Darwin/monthly/"

ECCO_darwin_native_grid(::ECCO4DarwinMonthly) = GridSpec("LatLonCap", MeshArrays.GRID_LLC90)
ECCO_darwin_native_size(::ECCO4DarwinMonthly) = (90, 1170, 50)

"""
    ECCO_darwin_model_data(metafile)

Read a ECCO4DarwinMonthly data file and regrid using MeshArrays on to regular lat-lon grid
"""
function ECCO_darwin_model_data(metadata, path)
    native_size = ECCO_darwin_native_size(metadata.version)
    native_grid = ECCO_darwin_native_grid(metadata.version)
    native_data = zeros(Float32, prod(native_size)) # Native LLC90 grid at precision of the input binary file

    read!(path, native_data)
    native_data = bswap.(native_data)

    meshed_data = read(reshape(native_data, native_size...), native_grid)
    Nx, Ny, Nz, _ = size(metadata)
    coeffs = interpolation_setup()
    data = zeros(Float32, Nx, Ny, Nz) # Native LLC90 grid at precision of the input binary file

    # Interpolate each layer
    for k in 1:Nz
        i, j, c = Interpolate(meshed_data[:, k], coeffs)
        data[:, :, k] = c
    end
    
    # Reverse the z-axis
    data = reverse(data, dims=3)

    # Fill NaNs in Antarctica with zeros
    data[isnan.(data)] .= 0.f0

    scale_factor = ECCO_darwin_scale_factor[metadata.name]
    # Scale data according to metadata.scale_factor
    return data .* scale_factor
end

retrieve_data(metadata::ECCOMetadata{<:Any, <:ECCO4DarwinMonthly}, path) = ECCO_darwin_model_data(metadata, path)
