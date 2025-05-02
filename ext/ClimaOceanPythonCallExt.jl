module ClimaOceanPythonCallExt

using ClimaOcean
using CondaPkg
using PythonCall

using Dates: DateTime
using ClimaOcean.DataWrangling.Copernicus: CopernicusMetadata

import ClimaOcean.DataWrangling: download_dataset

"""
    install_copernicusmarine()

Install the Copernicus Marine CLI using CondaPkg.
Returns a NamedTuple containing package information if successful.
"""
function install_copernicusmarine()
    @info "Installing the copernicusmarine CLI..."
    CondaPkg.add("copernicusmarine"; channel = "conda-forge")
    cli = CondaPkg.which("copernicusmarine")
    @info "... the copernicusmarine CLI has been installed at $(cli)."
    return cli
end

function download_dataset(meta::CopernicusMetadata, grid=nothing; skip_existing = true, additional_kw...)
    output_directory = meta.dir
    output_filename = ClimaOcean.DataWrangling.metadata_filename(meta)
    output_path = joinpath(output_directory, output_filename)
    isfile(output_path) && return output_path
    # rm(output_path, force=true)

    toolbox = try 
        pyimport("copernicusmarine")
    catch
        install_copernicusmarine()
        pyimport("copernicusmarine")
    end

    variables = [
        ClimaOcean.DataWrangling.Copernicus.copernicus_dataset_variable_names[meta.name]
    ]

    dataset_id = ClimaOcean.DataWrangling.Copernicus.copernicusmarine_dataset_id(meta.dataset)

    kw = (; coordinates_selection_method = "outside",
          skip_existing,
          dataset_id,
          variables,
          output_filename,
          output_directory)

    if !(meta.dataset isa ClimaOcean.DataWrangling.Copernicus.GLORYSStatic)
        start_datetime = ClimaOcean.DataWrangling.Copernicus.start_date_str(meta.dates)
        end_datetime = ClimaOcean.DataWrangling.Copernicus.end_date_str(meta.dates)
        kw = merge(kw, (; start_datetime, end_datetime))
    end

    kw = with_longitude_bounds(kw, meta.bounding_box)
    kw = with_latitude_bounds(kw, meta.bounding_box)
    kw = with_depth_bounds(kw, meta.bounding_box)

    @show kw

    toolbox.subset(; kw..., additional_kw...)

    return output_path
end

with_longitude_bounds(kw, ::Nothing) = kw
with_latitude_bounds(kw, ::Nothing) = kw
with_depth_bounds(kw, ::Nothing) = kw

const BBOX = ClimaOcean.DataWrangling.BoundingBox

with_longitude_bounds(kw, bounding_box::BBOX) = with_longitude_bounds(kw, bounding_box.longitude)
with_latitude_bounds(kw, bounding_box::BBOX) = with_latitude_bounds(kw, bounding_box.latitude)
with_depth_bounds(kw, bounding_box::BBOX) = with_depth_bounds(kw, bounding_box.z)

function with_longitude_bounds(kw, longitude)
    minimum_longitude = longitude[1]
    maximum_longitude = longitude[2]
    return merge(kw, (; minimum_longitude, maximum_longitude))
end

function with_latitude_bounds(kw, latitude)
    minimum_latitude = latitude[1]
    maximum_latitude = latitude[2]
    return merge(kw, (; minimum_latitude, maximum_latitude))
end

function with_depth_bounds(kw, z)
    minimum_depth = - z[2]
    maximum_depth = - z[1]
    return merge(kw, (; minimum_depth, maximum_depth))
end

end # module ClimaOceanPythonCallExt 
