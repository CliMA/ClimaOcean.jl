include("runtests_setup.jl")

using CopernicusMarine

@testset "Downloading Copernicus data" begin
    bounding_box = ClimaOcean.DataWrangling.BoundingBox(longitude=(200, 202), latitude=(35, 37))

    # Physics datasets
    variables = (:temperature, :salinity, :u_velocity, :v_velocity)
    dataset = ClimaOcean.DataWrangling.Copernicus.GLORYSDaily()
    for variable in variables
        metadatum = Metadatum(variable; dataset, bounding_box)
        filepath = ClimaOcean.DataWrangling.metadata_path(metadatum)
        isfile(filepath) && rm(filepath; force=true)
        ClimaOcean.DataWrangling.download_dataset(metadatum)
    end

    # Biogeochemistry datasets
    variables = (:nitrate, :phosphate, :dissolved_silicate)
    dataset = ClimaOcean.DataWrangling.Copernicus.GLORYSBGCDaily()
    for variable in variables
        metadatum = Metadatum(variable; dataset, bounding_box)
        filepath = ClimaOcean.DataWrangling.metadata_path(metadatum)
        isfile(filepath) && rm(filepath; force=true)
        ClimaOcean.DataWrangling.download_dataset(metadatum)
    end
end
