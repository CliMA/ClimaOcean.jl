using CFTime
using Dates
using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: AnyDateTime
using Oceananigans.DistributedComputations

import Dates: year, month, day
using Downloads

import Oceananigans.Fields: set!, location
import Base
import ClimaOcean.DataWrangling: all_dates, metadata_filename, download_dataset, default_download_directory
using ZipFile

struct EN4Monthly end

const EN4Metadata{D} = Metadata{<:EN4Monthly, D}
const EN4Metadatum   = Metadatum{<:EN4Monthly}

const EN4_url = "http://www.metoffice.gov.uk/hadobs/en4/data/en4-2-1/EN.4.2.2/EN.4.2.2.analyses.g10."

"""
    EN4Metadatum(name; 
                  date = first_date(EN4Monthly()), 
                  dir = download_EN4_cache)

an alias to construct a [`Metadatum`](@ref) of [`EN4Montly`](@ref)
"""
function EN4Metadatum(name; 
                       date = first_date(EN4Monthly()), 
                       dir = download_EN4_cache)
  
    return Metadatum(name; date, dir, dataset=EN4Monthly())
end

default_download_directory(::EN4Monthly) = download_EN4_cache

datasetstr(md::EN4Metadata) = string(md.dataset)

datestr(md::EN4Metadata) = string(first(md.dates), "--", last(md.dates))
datestr(md::EN4Metadatum) = string(md.dates)

Base.summary(md::EN4Metadata) = string("EN4Metadata{", datasetstr(md), "} of ",
                                        md.name, " for ", datestr(md))

Base.size(data::Metadata{<:EN4Monthly}) = (360, 173, 42, length(data.dates))

Base.size(::Metadatum{<:EN4Monthly}) = (360, 173, 42, 1)

# The whole range of dates in the different dataset datasets
all_dates(::EN4Monthly, name) = DateTime(1900, 1, 1) : Month(1) : DateTime(2024, 12, 1)

# Fallback, actually, we do not really need the name for EN4 since all
# variables have the same frequency and the same time-range, differently from JRA55
# all_dates(dataset::Union{<:EN4Monthly, <:EN42Monthly, <:EN42Daily}) = all_dates(dataset, :temperature)

# File name generation specific to each Dataset dataset
function metadata_filename(metadata::Metadatum{<:EN4Monthly})
    shortname = short_name(metadata)
    yearstr  = string(Dates.year(metadata.dates))
    monthstr = string(Dates.month(metadata.dates), pad=2)
    return shortname * "_" * yearstr * "_" * monthstr * ".nc"
end

# Convenience functions
short_name(data::Metadata{<:EN4Monthly}) = EN4_short_names[data.name]

location(data::EN4Metadata) = EN4_location[data.name]

variable_is_three_dimensional(data::EN4Metadata) =
    data.name == :temperature ||
    data.name == :salinity

EN4_short_names = Dict(
    :temperature           => "temperature",
    :salinity              => "salinity"
)

EN4_location = Dict(
    :temperature           => (Center, Center, Center),
    :salinity              => (Center, Center, Center)
)

function metadata_url(m::Metadata{<:EN4Monthly})
    year = string(Dates.year(m.dates))
    return EN4_url * year * ".zip"
end

function metadata_path(m::Metadata{<:EN4Monthly})
    year = string(Dates.year(m.dates))
    month = string(Dates.month(m.dates))
    zipfile = m.dir * "EN4_" * year * ".zip"
    extracted_file = m.dir * "EN4_" * year * "_" * month * ".nc"
    return zipfile, extracted_file
end

"""
    download_dataset(metadata::EN4Metadata; url = urls(metadata))

Download the dataset specified by the `metadata::EN4Metadata`. If `metadata.dates` is a single date,
the dataset is downloaded directly. If `metadata.dates` is a vector of dates, each date
is downloaded individually.

The data download requires a username and password to be provided in the `EN4_USERNAME` and
`EN4_PASSWORD` environment variables respectively. This can be done by exporting the
environment variables in the shell before running the script, or by launching julia with

```
EN4_USERNAME=myusername EN4_PASSWORD=mypassword julia
```

or by invoking

```julia
julia> ENV["EN4_USERNAME"] = "myusername"

julia> ENV["EN4_PASSWORD"] = "mypassword"
```

within julia.


Arguments
=========
- `metadata::EN4Metadata`: The metadata specifying the dataset to be downloaded.
"""
function download_dataset(metadata::Metadata{<:EN4Monthly})
    dir = metadata.dir

    # Create a temporary directory to store the .netrc file
    # The directory will be deleted after the download is complete
    @root mktempdir(dir) do tmp
        ntasks = Threads.nthreads()

        asyncmap(metadata; ntasks) do metadatum # Distribute the download among tasks
            fileurl  = metadata_url(metadatum)
            zippath = metadata_path(metadatum)[1]
            extracted_file = metadata_path(metadatum)[2]
            if !isfile(zippath) & !isfile(extracted_file)
                instructions_msg = "\n See ClimaOcean.jl/src/DataWrangling/EN4/README.md for instructions."
                @info "Downloading EN4 data: $(metadatum.name) in $(metadatum.dir)..."
                Downloads.download(fileurl, zippath; progress=download_progress)
            end
        ZipFile.extract(metadata.dir*"*.zip", metadata.dir)
        end
    end

    return nothing
end
