using CFTime
using Dates
using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: netrc_downloader
using ClimaOcean.InitialConditions: three_dimensional_regrid!, interpolate!

import Dates: year, month, day

using Base: @propagate_inbounds
using Downloads

import Oceananigans.Fields: set!, location
import Base

using MeshArrays

struct ECCO270DarwinMonthly end
struct ECCO4DarwinMonthly end

"""
    struct ECCODarwinMetadata{D, V}

Metadata information for an ECCO dataset:
- `name`: The name of the dataset.
- `dates`: The dates of the dataset, in an `AbstractCFDateTime` format.
- `version`: The version of the dataset, could be `ECCO270DarwinMonthly`, or `ECCO4DarwinMonthly`.
- `dir`: The directory where the dataset is stored.
"""
struct ECCODarwinMetadata{D, R, V}
    name    :: Symbol
    dates   :: D
    date_ref:: R
    dt      :: Int
    version :: V
    dir     :: String
end

Base.show(io::IO, metadata::ECCODarwinMetadata) = 
    print(io, "ECCODarwinMetadata:", '\n',
    "├── name:    $(metadata.name)", '\n',
    "├── dates:   $(metadata.dates)", '\n',
    "├── date_ref:$(metadata.date_ref)", '\n',
    "├── dt:      $(metadata.dt)", '\n',
    "├── version: $(metadata.version)", '\n',
    "└── dir:     $(metadata.dir)")

Base.summary(md::ECCODarwinMetadata{<:Any, <:ECCO4DarwinMonthly})   = "ECCO4DarwinMonthly $(md.name) metadata ($(first(md.dates))--$(last(md.dates)))"
Base.summary(md::ECCODarwinMetadata{<:Any, <:ECCO270DarwinMonthly}) = "ECCO270DarwinMonthly $(md.name) metadata ($(first(md.dates))--$(last(md.dates)))"

Base.summary(md::ECCODarwinMetadata{<:AbstractCFDateTime, <:ECCO4DarwinMonthly})   = "ECCO4DarwinMonthly $(md.name) metadata at $(md.dates)"
Base.summary(md::ECCODarwinMetadata{<:AbstractCFDateTime, <:ECCO270DarwinMonthly}) = "ECCO270DarwinMonthly $(md.name) metadata at $(md.dates)"

"""
    ECCODarwinMetadata(name::Symbol; 
                dates    = DateTimeProlepticGregorian(1993, 1, 1),
                date_ref = DateTimeProlepticGregorian(1993, 1, 1, 12, 0, 0),
                dt       = 3600,
                version  = ECCO4DarwinMonthly(),
                dir      = download_ECCO_cache)

Construct an `ECCODarwinMetadata` object with the specified parameters.

Arguments
=========
- `name::Symbol`: The name of the metadata.

Keyword Arguments
=================
- `dates`: The date(s) of the metadata. Note this can either be a single date,
           representing a snapshot, or a range of dates, representing a time-series.
           Default: `DateTimeProlepticGregorian(1993, 1, 1)`.


- `date_ref`: The reference date. Default: `DateTimeProlepticGregorian(1993, 1, 1, 12, 0, 0)`.

- `dt`: The time step. Default: `3600`.
           
- `version`: The data version. Supported versions are `ECCO270DarwinMonthly()` or `ECCO4DarwinMonthly()`.

- `dir`: The directory of the data file. Default: `download_ECCO_cache`.
"""
function ECCODarwinMetadata(name::Symbol; 
    dates    = DateTimeProlepticGregorian(1993, 1, 1),
    date_ref = DateTimeProlepticGregorian(1992, 1, 1, 12, 0, 0),
    dt       = 3600,
    version  = ECCO4DarwinMonthly(),
    dir      = download_ECCO_cache,
)

    return ECCODarwinMetadata(name, dates, date_ref, dt, version, dir)
end

function ECCODarwinMetadata(name::Symbol, 
    dates,
    version::ECCO4DarwinMonthly;
    dir      = download_ECCO_cache,
)
    date_ref = DateTimeProlepticGregorian(1992, 1, 1, 12, 0, 0)
    dt       = 3600
    return ECCODarwinMetadata(name, dates, date_ref, dt, version, dir)
end

function ECCODarwinMetadata(name::Symbol, 
    dates, 
    version::ECCO270DarwinMonthly; 
    dir=download_ECCO_cache)

    date_ref = DateTimeProlepticGregorian(1992, 1, 1, 0, 0, 0)
    dt = 1200
    ECCODarwinMetadata(name, dates, date_ref, dt, version, dir)
end

# Treat ECCODarwinMetadata as an array to allow iteration over the dates.
Base.eltype(metadata::ECCODarwinMetadata) = Base.eltype(metadata.dates)

@propagate_inbounds Base.getindex(m::ECCODarwinMetadata, i::Int) = ECCODarwinMetadata(m.name, m.dates[i],   m.date_ref, m.dt, m.version, m.dir)
@propagate_inbounds Base.first(m::ECCODarwinMetadata)            = ECCODarwinMetadata(m.name, m.dates[1],   m.date_ref, m.dt, m.version, m.dir)
@propagate_inbounds Base.last(m::ECCODarwinMetadata)             = ECCODarwinMetadata(m.name, m.dates[end], m.date_ref, m.dt, m.version, m.dir)

@inline function Base.iterate(m::ECCODarwinMetadata, i=1)
    if (i % UInt) - 1 < length(m)
        return ECCODarwinMetadata(m.name, m.dates[i], m.date_ref, m.dt, m.version, m.dir), i + 1
    else
        return nothing
    end
end


Base.axes(metadata::ECCODarwinMetadata{<:AbstractCFDateTime})    = 1
Base.first(metadata::ECCODarwinMetadata{<:AbstractCFDateTime})   = metadata
Base.last(metadata::ECCODarwinMetadata{<:AbstractCFDateTime})    = metadata
Base.iterate(metadata::ECCODarwinMetadata{<:AbstractCFDateTime}) = (metadata, nothing)
Base.iterate(::ECCODarwinMetadata{<:AbstractCFDateTime}, ::Any)  = nothing

Base.length(metadata::ECCODarwinMetadata) = length(metadata.dates)
Base.size(data::ECCODarwinMetadata{<:Any, <:Any, <:ECCO270DarwinMonthly}) = (1080, 540, 50, length(data.dates))
Base.size(data::ECCODarwinMetadata{<:Any, <:Any, <:ECCO4DarwinMonthly})   = (360,  180, 50, length(data.dates))

Base.length(metadata::ECCODarwinMetadata{<:AbstractCFDateTime}) = 1
Base.size(::ECCODarwinMetadata{<:AbstractCFDateTime, <:Any, <:ECCO270DarwinMonthly}) = (1080, 540, 50, 1)
Base.size(::ECCODarwinMetadata{<:AbstractCFDateTime, <:Any, <:ECCO4DarwinMonthly})   = (360,  180, 50, 1)

# The whole range of dates in the different dataset versions
all_ECCO_dates(::ECCO4DarwinMonthly)   = DateTimeProlepticGregorian(1992, 1, 1) : Month(1) : DateTimeProlepticGregorian(2023, 12, 31)
all_ECCO_dates(::ECCO270DarwinMonthly) = DateTimeProlepticGregorian(1992, 1, 1) : Month(1) : DateTimeProlepticGregorian(2023, 12, 31)

# File names of metadata containing multiple dates
metadata_filename(metadata) = [metadata_filename(metadatum) for metadatum in metadata]

# File name generation specific to each Dataset version
function metadata_filename(metadata::ECCODarwinMetadata{<:AbstractCFDateTime})
    shortname = short_name(metadata)
    
    iternum = Dates.value((metadata.dates-metadata.date_ref)/(metadata.dt*1e3))
    iterstr = string(iternum, pad=10)

    return shortname * "." * iterstr * ".data"
end

# Convenience functions
metadata_path(metadata) = joinpath(metadata.dir, metadata_filename(metadata))
short_name(data::ECCODarwinMetadata{<:Any, <:Any, <:ECCO270DarwinMonthly}) = ECCO270Darwin_short_names[data.name]
short_name(data::ECCODarwinMetadata{<:Any, <:Any, <:ECCO4DarwinMonthly})   = ECCO4Darwin_short_names[data.name]

metadata_url(prefix, m::ECCODarwinMetadata{<:Any, <:Any, <:ECCO4DarwinMonthly}) = prefix * "/" * short_name(m) * "/" * metadata_filename(m)
metadata_url(prefix, m::ECCODarwinMetadata{<:Any, <:Any, <:ECCO270DarwinMonthly}) = prefix * "/" * short_name(m) * "/" * metadata_filename(m)

location(data::ECCODarwinMetadata) = ECCODarwin_location[data.name]

variable_is_three_dimensional(data::ECCODarwinMetadata) =
    data.name == :DIC ||
    data.name == :ALK ||
    data.name == :PO₄ ||
    data.name == :NO₃ ||
    data.name == :DOP ||
    data.name == :POP ||
    data.name == :Fe  ||
    data.name == :Siᵀ

ECCO4Darwin_short_names = Dict(
    :DIC => "DIC",
    :ALK => "ALK",
    :PO₄ => "PO4",
    :NO₃ => "NO3",
    :DOP => "DOP",
    :POP => "POP",
    :Fe  => "FeT",
    :Siᵀ => "SiO2",
    :PAR => "PAR",
)

ECCO270Darwin_short_names = Dict(
    :DIC => "DIC",
    :ALK => "ALK",
    :PO₄ => "PO4",
    :NO₃ => "NO3",
    :DOP => "DOP",
    :POP => "POP",
    :Fe  => "FeT",
    :Siᵀ => "SiO2",
    :PAR => "PAR",
)

ECCODarwin_location = Dict(
    :DIC => (Center, Center, Center),
    :ALK => (Center, Center, Center),
    :PO₄ => (Center, Center, Center),
    :NO₃ => (Center, Center, Center),
    :DOP => (Center, Center, Center),
    :POP => (Center, Center, Center),
    :Fe  => (Center, Center, Center),
    :Siᵀ => (Center, Center, Center),
    :PAR => (Center, Center, Nothing),
)

# URLs for the ECCO datasets specific to each version
urls(::ECCODarwinMetadata{<:Any, <:Any, <:ECCO270DarwinMonthly}) = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/ECCO-Darwin_extension/monthly"
urls(::ECCODarwinMetadata{<:Any, <:Any, <:ECCO4DarwinMonthly})   = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC90/ECCO-Darwin/monthly/"

"""
    download_dataset(metadata::ECCODarwinMetadata; url = urls(metadata))

Download the dataset specified by `ECCODarwinMetadata`. If `ECCODarwinMetadata.dates` is a single date, 
the dataset is downloaded directly. If `ECCODarwinMetadata.dates` is a vector of dates, each date
is downloaded individually.
The data download requires a username and password to be provided in the `ECCO_USERNAME` and
`ECCO_PASSWORD` environment variables. This can be done by exporting the environment variables
in the shell before running the script, or by launching julia with 

```
ECCO_USERNAME=myusername ECCO_PASSWORD=mypassword julia 
```

Arguments
=========
- `metadata::ECCODarwinMetadata`: The metadata specifying the dataset to be downloaded.
"""
function download_dataset(metadata::ECCODarwinMetadata; url = urls(metadata))
    username = get(ENV, "ECCO_USERNAME", nothing)
    password = get(ENV, "ECCO_PASSWORD", nothing)
    dir = metadata.dir
    
    # Create a temporary directory to store the .netrc file
    # The directory will be deleted after the download is complete
    @root mktempdir(dir) do tmp

        # Write down the username and password in a .netrc file
        downloader = netrc_downloader(username, password, "ecco.jpl.nasa.gov", tmp)

        asyncmap(metadata, ntasks=10) do metadatum # Distribute the download among tasks

            fileurl  = metadata_url(url, metadatum) 
            filepath = metadata_path(metadatum)

            if !isfile(filepath)
                instructions_msg = "\n See ClimaOcean.jl/src/ECCO/README.md for instructions."
                if isnothing(username)
                    msg = "Could not find the ECCO_PASSWORD environment variable. \
                           See ClimaOcean.jl/src/ECCO/README.md for instructions on obtaining \
                           and setting your ECCO_USERNAME and ECCO_PASSWORD." * instructions_msg
                    throw(ArgumentError(msg))
                elseif isnothing(password)
                    msg = "Could not find the ECCO_PASSWORD environment variable. \
                           See ClimaOcean.jl/src/ECCO/README.md for instructions on obtaining \
                           and setting your ECCO_USERNAME and ECCO_PASSWORD." * instructions_msg
                    throw(ArgumentError(msg))
                end
                # Get the data file
                Downloads.download(fileurl, filepath; downloader, progress=download_progress)
                # Get the meta file
                Downloads.download(replace(fileurl, "data" => "meta"), 
                                   replace(filepath, "data" => "meta"); 
                                   downloader, 
                                   progress=download_progress)
            end
        end
    end
    
    return nothing
end

ECCODarwinMeshGrid(::ECCO4DarwinMonthly)   = GridSpec("LatLonCap", MeshArrays.GRID_LLC90)
ECCODarwinMeshGrid(::ECCO270DarwinMonthly) = GridSpec("LatLonCap", MeshArrays.GRID_LLC90)

"""
    ECCODarwinModelMeta(metafile)

Read ECCODarwin metadata file and return as a NamedTuple

"""
function ECCODarwinModelMeta(fileName)
    meta = read(fileName,String)
    meta = split(meta,";\n")
    meta = meta[isempty.(meta).==false]
    meta = replace.(meta,Ref(",\n"=>";"))
    meta = replace.(meta,Ref("\n"=>""))
    meta = replace.(meta,Ref("}"=>"]"))
    meta = replace.(meta,Ref("{"=>"["))
    meta = replace.(meta,Ref("'"=>"\""))
    meta = replace.(meta,Ref(";]"=>"]"))
    meta = replace.(meta,Ref(","=>" "))

    meta = split.(meta,"=")
    meta = [[replace(x[1]," "=>"") x[2]] for x in meta]

    metaDict = Dict{String,Any}(m[1] => m[2] for m in meta)

    for k in keys(metaDict)
        val = eval(Meta.parse(metaDict[k]))
        if isa(val[1],String)
            val = replace.(val,Ref(" "=>""))
        end
        if length(val) == 1
            val = val[1]
        end
        metaDict[k] = val
    end
    metaDict["dataprec"] = titlecase(metaDict["dataprec"])

    meta = (; zip(Symbol.(keys(metaDict)),values(metaDict))...)
    return meta
end

"""
    ECCODarwinModelData(metafile)

Read a ECCODarwin data file and regrid using MeshArrays on to regular lat-lon grid

"""
function ECCODarwinModelData(file_name::String, meta_data, mesh_grid)
    T     = eval(Symbol(meta_data.dataprec))
    fid   = open(file_name)
    ndims = meta_data.nDims
    sizes = meta_data.dimList[:,1]
    bin   = Array{T,1}(undef,*(meta_data.dimList[:,1]...));

    read!(fid,bin)
    bin   = hton.(bin)
    close(fid)

    ma    = read(
                reshape(bin, sizes...),
                mesh_grid,
    )

    #grid  = GridLoad(mesh_grid; option="full")
    #lat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
    #lon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
    #(f,i,j,w)=InterpolationFactors(grid, vec(lon), vec(lat))
    #coeffs=(lon=lon, lat=lat, f=f, i=i, j=j, w=w)
    coeffs = interpolation_setup()
    data   = Array{T,ndims}(
                undef,
                (size(coeffs.lon)...,size(ma,2)),
    )

    # Interpolate each layer
    for z in 1:size(ma,2)
        i, j, c = Interpolate(ma[:,z], coeffs)
        data[:,:,z] = c
    end

    return data
end

"""
    ECCO_field(metadata::ECCODarwinMetadata;
               architecture = CPU(),
               inpainting = nothing,
               mask = nothing,
               horizontal_halo = (7, 7),
               cache_inpainted_data = false)

Return a `Field` on `architecture` described by `ECCOMetadata` with
`horizontal_halo` size.
If not `nothing`, the `inpainting` method is used to fill the cells
within the specified `mask`. `mask` is set to `ECCO_mask` for non-nothing
`inpainting`.
"""
function ECCO_field(metadata::ECCODarwinMetadata;
                    architecture = CPU(),
                    inpainting = NearestNeighborInpainting(Inf),
                    mask = nothing,
                    horizontal_halo = (7, 7),
                    cache_inpainted_data = true)

    z = ClimaOcean.DataWrangling.ECCO.ECCO_z
    input_data_grid = LatitudeLongitudeGrid(architecture; 
                z,
                size      = (720, 360, 50),
                latitude  = (-90, 90),
                longitude = (0, 360),
    )
    field_input_data = Field{Center, Center, Center}(input_data_grid)
    inpainted_path = inpainted_metadata_path(metadata)

    if !isnothing(inpainting) && isfile(inpainted_path)
        file = jldopen(inpainted_path, "r")
        maxiter = file["inpainting_maxiter"]

        # read data if generated with the same inpainting
        if maxiter == inpainting.maxiter
            data = file["data"]
            close(file)
            copyto!(parent(field_input_data), data)
            return field_input_data
        end

        close(file)
    end

    download_dataset(metadata)
    path = metadata_path(metadata)

    # Open ECCODarwin binary data file
    data = ECCODarwinModelData(
                path, 
                ECCODarwinModelMeta(
                    replace(path,"data"=>"meta"),
                ),
                ECCODarwinMeshGrid(
                    metadata.version,
                ),
    )

    if variable_is_three_dimensional(metadata)
        data = reverse(data, dims=3)
    end        
    
    # Convert data from Union(FT, missing} to FT
    FT = eltype(field_input_data)
    data[ismissing.(data)] .= 1e10 # Artificially large number!
    data = if location(field_input_data)[2] == Face # ?
        new_data = zeros(FT, size(field_input_data))
        new_data[:, 1:end-1, :] .= data
        new_data    
    else
        data = Array{FT}(data)
    end
    
    # ECCO4 data is on a -180, 180 longitude grid as opposed to ECCO2 data that
    # is on a 0, 360 longitude grid. To make the data consistent, we shift ECCO4
    # data by 180 degrees in longitude
    if metadata.version isa ECCO4DarwinMonthly 
        Nx = size(data, 1)
        if variable_is_three_dimensional(metadata)
            shift = (Nx ÷ 2, 0, 0)
        else
            shift = (Nx ÷ 2, 0)
        end
        data = circshift(data, shift)
    end
    set!(field_input_data, data)
    fill_halo_regions!(field_input_data)

    if !isnothing(inpainting)
        # Respect user-supplied mask, but otherwise build default ECCO mask.
        if isnothing(mask)
            mask = ECCO_mask(metadata, architecture; data_field=field_input_data)
        end

        # Make sure all values are extended properly
        name = string(metadata.name)
        date = string(metadata.dates)
        version = summary(metadata.version)
        @info string("Inpainting ", version, " ", name, " data from ", date, "...")
        start_time = time_ns()
        
        inpaint_mask!(field_input_data, mask; inpainting)
        fill_halo_regions!(field_input_data)

        elapsed = 1e-9 * (time_ns() - start_time)
        @info string(" ... (", prettytime(elapsed), ")")
    
        # We cache the inpainted data to avoid recomputing it
        @root if cache_inpainted_data
            file = jldopen(inpainted_path, "w+")
            file["data"] = on_architecture(CPU(), parent(field_input_data))
            file["inpainting_maxiter"] = inpainting.maxiter
            close(file)
        end
    end

    return field_input_data
end