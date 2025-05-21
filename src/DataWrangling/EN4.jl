module EN4

export EN4Metadatum, EN4_immersed_grid, adjusted_EN4_tracers, initialize!
export EN4Monthly

using ClimaOcean
using Oceananigans
using NCDatasets
using JLD2
using Downloads: download
using Adapt
using Scratch

using ClimaOcean.DataWrangling:
    Metadata,
    Metadatum,
    BoundingBox,
    inpaint_mask!,
    NearestNeighborInpainting,
    download_progress,
    compute_native_date_range,
    Kelvin,
    Celsius

using KernelAbstractions: @kernel, @index

using Oceananigans.Architectures: architecture

using Dates: year, month, day
using Oceananigans.DistributedComputations: @root

using Dates
import Downloads
import ZipFile

import ClimaOcean.DataWrangling:
    all_dates,
    metadata_filename,
    download_dataset,
    default_download_directory,
    metadata_path,
    temperature_units,
    dataset_variable_name,
    metaprefix,
    z_interfaces,
    longitude_interfaces,
    latitude_interfaces,
    is_three_dimensional,
    reversed_vertical_axis,
    inpainted_metadata_path,
    available_variables

import Oceananigans.Fields: location

download_EN4_cache::String = ""
function __init__()
    global download_EN4_cache = @get_scratch!("EN4")
end

EN4_dataset_variable_names = Dict(
    :temperature => "temperature",
    :salinity    => "salinity"
)

struct EN4Monthly end

default_download_directory(::EN4Monthly) = download_EN4_cache
Base.size(::EN4Monthly, variable) = (360, 173, 42)
all_dates(::EN4Monthly, variable) = DateTime(1900, 1, 1) : Month(1) : DateTime(2024, 12, 1)
temperature_units(::EN4Monthly) = Kelvin()
reversed_vertical_axis(::EN4Monthly) = true

longitude_interfaces(::EN4Monthly) = (0.5, 360.5)
latitude_interfaces(::EN4Monthly) = (-83.5, 89.5)
available_variables(::EN4Monthly) = EN4_dataset_variable_names

z_interfaces(::EN4Monthly) = [
    -5500.0,
    -5200.5986,
    -4901.459,
    -4602.695,
    -4304.469,
    -4007.0154,
    -3710.6636,
    -3415.8828,
    -3123.3313,
    -2833.926,
    -2548.9243,
    -2270.0168,
    -1999.4135,
    -1739.8875,
    -1494.7288,
    -1267.542,
    -1061.846,
    -880.5102,
    -725.1709,
    -595.8501,
    -490.9494,
    -407.6244,
    -342.3651,
    -291.5638,
    -251.9186,
    -220.6418,
    -195.5097,
    -174.8189,
    -157.3019,
    -142.0345,
    -128.3522,
    -115.7825,
    -103.9913,
    -92.7437,
    -81.8753,
    -71.2711,
    -60.8509,
    -50.5586,
    -40.3553,
    -30.214,
    -20.1158,
    -10.0475,
      0.0,
]

const EN4Metadata{D} = Metadata{<:EN4Monthly, D}
const EN4Metadatum   = Metadatum{<:EN4Monthly}

const EN4_url_pre2021  = "http://www.metoffice.gov.uk/hadobs/en4/data/en4-2-1/EN.4.2.2/EN.4.2.2.analyses.g10."
const EN4_url_post2021 = "http://www.metoffice.gov.uk/hadobs/en4/data/en4-2-1/EN.4.2.2.analyses.g10."

function inpainted_metadata_filename(metadata::EN4Metadata)
    original_filename = metadata_filename(metadata)
    without_extension = original_filename[1:end-3]
    var = string(metadata.name)
    return without_extension * "_" * var *"_inpainted.jld2"
end

inpainted_metadata_path(metadata::EN4Metadata) = joinpath(metadata.dir, inpainted_metadata_filename(metadata))

"""
    EN4Metadatum(name;
                 date = first_date(EN4Monthly()),
                 dir = download_EN4_cache)

An alias to construct a [`Metadatum`](@ref) of `EN4Monthly`.
"""
function EN4Metadatum(name;
                      date = first_date(EN4Monthly(), name),
                      dir = download_EN4_cache)

    return Metadatum(name; date, dir, dataset=EN4Monthly())
end

metaprefix(::EN4Metadata) = "EN4Metadata"
metaprefix(::EN4Metadatum) = "EN4Metadatum"

# Note, EN4 files contain all variables, so the filenames do not
# depend on metadata.name.
function metadata_filename(metadata::Metadatum{<:EN4Monthly})
    yearstr  = string(Dates.year(metadata.dates))
    monthstr = string(Dates.month(metadata.dates), pad=2)
    return "EN.4.2.2.f.analysis.g10." * yearstr * lpad(string(monthstr), 2, '0') * ".nc"
end

# Convenience functions
dataset_variable_name(data::EN4Metadata) = EN4_dataset_variable_names[data.name]
location(::EN4Metadata) = (Center, Center, Center)
is_three_dimensional(::EN4Metadata) = true

## This function is explicitly for the downloader to check if the zip file/extracted file exists,
## then to download the relevant URL (from above)
function metadata_zippath(m::EN4Metadata)
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

function metadata_url(m::EN4Metadata)
    year = string(Dates.year(m.dates))
    if Dates.year(m.dates) < 2021
        return EN4_url_pre2021 * year * ".zip"
    else
        return EN4_url_post2021 * year * ".zip"
    end
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

end # Module
