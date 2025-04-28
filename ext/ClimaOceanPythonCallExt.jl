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
    CondaPkg.add("copernicusmarine"; channel = "conda-forge")
    cli = CondaPkg.which("copernicusmarine")
    @info "The copernicusmarine CLI has been installed at $(cli)."
    return cli
end


function download_dataset(meta::CopernicusMetadata, grid=nothing; skip_existing = true, kw...)
    start_datetime = ClimaOcean.DataWrangling.Copernicus.start_date_str(meta.dates)
    end_datetime = ClimaOcean.DataWrangling.Copernicus.end_date_str(meta.dates)
    dataset_id = ClimaOcean.DataWrangling.Copernicus.copernicusmarine_dataset_id(meta.dataset)

    toolbox = pyimport("copernicusmarine")
    output_filename = ClimaOcean.DataWrangling.metadata_filename(meta)
    output_directory = meta.dir
    variables = [
        ClimaOcean.DataWrangling.Copernicus.copernicus_short_names[meta.name]
    ]

    toolbox.subset(; skip_existing, kw...,
                   dataset_id,
                   variables,
                   start_datetime,
                   end_datetime,
                   output_filename,
                   output_directory)
end

end # module ClimaOceanPythonCallExt 
