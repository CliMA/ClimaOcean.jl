module ClimaOceanPythonCallExt

using ClimaOcean
using CondaPkg
using PythonCall
using Oceananigans.DistributedComputations: @root

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

    toolbox = try 
        pyimport("copernicusmarine")
    catch
        install_copernicusmarine()
        pyimport("copernicusmarine")
    end

    variable_name = ClimaOcean.DataWrangling.Copernicus.copernicus_dataset_variable_names[meta.name]
    variables = PythonCall.pylist([variable_name])

    dataset_id = ClimaOcean.DataWrangling.Copernicus.copernicusmarine_dataset_id(meta.dataset)
    datetime_kw = if meta.dataset isa ClimaOcean.DataWrangling.Copernicus.GLORYSStatic
        NamedTuple()
    else
        start_datetime = ClimaOcean.DataWrangling.Copernicus.start_date_str(meta.dates)
        end_datetime = ClimaOcean.DataWrangling.Copernicus.end_date_str(meta.dates)
        (; start_datetime, end_datetime)
    end

    lon_kw = longitude_bounds_kw(meta.bounding_box)
    lat_kw = latitude_bounds_kw(meta.bounding_box)
    z_kw = depth_bounds_kw(meta.bounding_box)

    kw = (; coordinates_selection_method = "outside",
          skip_existing,
          dataset_id,
          variables,
          output_filename,
          output_directory)

    additional_kw = NamedTuple(name => value for (name, value) in additional_kw)
    kw = merge(kw, datetime_kw, lon_kw, lat_kw, z_kw, additional_kw)

    @root toolbox.subset(; kw...)

    return output_path
end

longitude_bounds_kw(::Nothing) = NamedTuple()
latitude_bounds_kw(::Nothing) = NamedTuple()
depth_bounds_kw(::Nothing) = NamedTuple()

const BBOX = ClimaOcean.DataWrangling.BoundingBox

longitude_bounds_kw(bounding_box::BBOX) = longitude_bounds_kw(bounding_box.longitude)
latitude_bounds_kw(bounding_box::BBOX) = latitude_bounds_kw(bounding_box.latitude)
depth_bounds_kw(bounding_box::BBOX) = depth_bounds_kw(bounding_box.z)

function longitude_bounds_kw(longitude)
    minimum_longitude = longitude[1]
    maximum_longitude = longitude[2]
    return (; minimum_longitude, maximum_longitude)
end

function latitude_bounds_kw(latitude)
    minimum_latitude = latitude[1]
    maximum_latitude = latitude[2]
    return (; minimum_latitude, maximum_latitude)
end

function depth_bounds_kw(z)
    minimum_depth = - z[2]
    maximum_depth = - z[1]
    return (; minimum_depth, maximum_depth)
end

end # module ClimaOceanPythonCallExt 
