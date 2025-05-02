include("runtests_setup.jl")

using PythonCall
using CondaPkg

@testset "Downloading Copernicus data" begin
    variables = (:temperature, :salinity, :u_velocity, :v_velocity)
    bounding_box = ClimaOcean.DataWrangling.BoundingBox(longitude=(200, 202), latitude=(35, 37))
    dataset = ClimaOcean.DataWrangling.Copernicus.GLORYSDaily()
    for variable in variables
        metadata = Metadata(variable; dataset, bounding_box)
        filepath = ClimaOcean.DataWrangling.metadata_path(metadata)
        isfile(filepath) && rm(filepath; force=true)
        ClimaOcean.DataWrangling.download_dataset(metadata)
    end
end
