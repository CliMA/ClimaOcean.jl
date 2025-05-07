using CFTime
using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: AnyDateTime, Celsius, Kelvin
using Dates
using Downloads
using Oceananigans.DistributedComputations
using ZipFile

import Dates: year, month, day

import Base
import Oceananigans.Fields: location
import ClimaOcean.DataWrangling: all_dates, metadata_filename, download_dataset,
                                 default_download_directory, metadata_path,
                                 short_name, dataset_latitude_extent

struct EN4Monthly end

const EN4Metadata{D} = Metadata{<:EN4Monthly, D}
const EN4Metadatum   = Metadatum{<:EN4Monthly}

const EN4_url_pre2021  = "http://www.metoffice.gov.uk/hadobs/en4/data/en4-2-1/EN.4.2.2/EN.4.2.2.analyses.g10."
const EN4_url_post2021 = "http://www.metoffice.gov.uk/hadobs/en4/data/en4-2-1/EN.4.2.2.analyses.g10."

"""
    EN4Metadatum(name;
                 date = first_date(EN4Monthly()),
                 dir = download_EN4_cache)

An alias to construct a [`Metadatum`](@ref) of `EN4Monthly`.
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

# File name generation specific to each Dataset dataset
function metadata_filename(metadata::Metadatum{<:EN4Monthly})
    shortname = short_name(metadata)
    yearstr  = string(Dates.year(metadata.dates))
    monthstr = string(Dates.month(metadata.dates), pad=2)
    return "EN.4.2.2.f.analysis.g10." * yearstr * lpad(string(monthstr), 2, '0') * ".nc"
end

# Convenience functions
short_name(data::Metadata{<:EN4Monthly}) = EN4_short_names[data.name]

location(data::EN4Metadata) = EN4_location[data.name]

dataset_latitude_extent(data::Metadata{<:EN4Monthly, <:Any}) = (-83.5, 89.5)

variable_is_three_dimensional(data::EN4Metadata) =
    data.name == :temperature ||
    data.name == :salinity

EN4_short_names = Dict(
    :temperature => "temperature",
    :salinity    => "salinity"
)

EN4_location = Dict(
    :temperature => (Center, Center, Center),
    :salinity    => (Center, Center, Center)
)

function metadata_url(m::EN4Metadata)
    year = string(Dates.year(m.dates))
    if Dates.year(m.dates) < 2021
        return EN4_url_pre2021 * year * ".zip"
    else
        return EN4_url_post2021 * year * ".zip"
    end
end

## This function is explicitly for the downloader to check if the zip file/extracted file exists,
## then to download the relevant URL (from above)

function metadata_path(m::Metadata{V, <:Union{AbstractCFDateTime, Dates.AbstractDateTime}}) where V<:EN4Monthly
    year = string(Dates.year(m.dates))
    month = string(Dates.month(m.dates))
    filename = "EN.4.2.2.f.analysis.g10." * year * lpad(string(month), 2, '0') * ".nc"
    unzipped_filepath = joinpath(m.dir, filename)
    return unzipped_filepath
end

function metadata_zippath(m::Metadata{V, <:Union{AbstractCFDateTime, Dates.AbstractDateTime}}) where V<:EN4Monthly
    year = string(Dates.year(m.dates))
    month = string(Dates.month(m.dates))
    zippath = joinpath(m.dir, "EN4_" * year * ".zip")
    return zippath
end

function unzip(file, exdir="")
    filepath = isabspath(file) ? file : joinpath(pwd(), file)
    basepath = dirname(filepath)
    outpath = (exdir == "" ? basepath : (isabspath(exdir) ? exdir : joinpath(pwd(), exdir)))
    isdir(outpath) ? "" : mkdir(outpath)
    zarchive = ZipFile.Reader(filepath)
    for f in zarchive.files
        filepath = joinpath(outpath, f.name)
        if endswith(f.name, "/") || endswith(f.name, "\\")
            mkdir(filepath)
        else
            write(filepath, read(f))
        end
    end
    close(zarchive)
end

function download_dataset(metadata::Metadata{<:EN4Monthly})
    dir = metadata.dir
    missingzips = []

    @root begin
        ntasks = Threads.nthreads()

        asyncmap(metadata; ntasks) do metadatum # Distribute the download among tasks
            fileurl = metadata_url(metadatum)
            zippath = metadata_zippath(metadatum)
            extracted_file = metadata_path(metadatum)
            if !isfile(extracted_file) & !isfile(zippath)
                push!(missingzips, zippath)
                @info "Downloading EN4 data: $(metadatum.name) in $(metadatum.dir)..."
                Downloads.download(fileurl, zippath; progress=download_progress)
            elseif !isfile(extracted_file) & isfile(zippath)
                push!(missingzips, zippath)
            end
        end

        for zips in unique(missingzips)
            unzip(zips, dir)
        end
    end

    return nothing
end
