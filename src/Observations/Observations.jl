"""
Observational climatologies registered with NumericalEarth's `DataWrangling`, so that `Metadata` and `Metadatum` of
them download, read and interpolate onto any grid like every other dataset.
"""
module Observations

export BoyerMontegutMonthly, HadISSTMonthly, NCEPMonthly

using Dates: Dates, DateTime, Month, month, year
using Downloads: Downloads
using NCDatasets: Dataset
using Oceananigans.DistributedComputations: @root
using Oceananigans.Fields: Center
using NumericalEarth.DataWrangling: DataWrangling, Metadata, Metadatum, metadata_path, download_with_retries

observations_cache::String = ""

function __init__()
    global observations_cache = DataWrangling.download_cache("ClimaOceanObservations")
    return nothing
end

# Climatologies are dated in this nominal year
const climatology_year = 2000
const climatology_dates = [DateTime(climatology_year, m, 1) for m in 1:12]

abstract type MonthlyObservation end

DataWrangling.default_download_directory(::MonthlyObservation) = observations_cache
DataWrangling.dataset_location(::MonthlyObservation, name) = (Center, Center, Nothing)
DataWrangling.averaging_window(metadatum::Metadatum{<:MonthlyObservation}) = DataWrangling.calendar_month_window(metadatum)
DataWrangling.is_three_dimensional(::Metadata{<:MonthlyObservation}) = false
DataWrangling.longitude_interfaces(::Metadata{<:MonthlyObservation}) = (-180, 180)
DataWrangling.latitude_interfaces(::Metadata{<:MonthlyObservation}) = (-90, 90)
DataWrangling.dataset_variable_name(metadata::Metadata{<:MonthlyObservation}) =
    DataWrangling.available_variables(metadata.dataset)[metadata.name]

# Every record of a dataset lives in one file, so the inpainted fields are told apart by their date
DataWrangling.inpainted_metadata_path(metadatum::Metadatum{<:MonthlyObservation}) =
    joinpath(metadatum.dir, string(metadatum.dataset, "_", metadatum.name, "_",
                                   Dates.format(metadatum.dates, "yyyy-mm"), "_inpainted.jld2"))

function Downloads.download(metadata::Metadata{<:MonthlyObservation}; kwargs...)
    for metadatum in metadata
        path = metadata_path(metadatum)
        @root isfile(path) || download_with_retries(observation_url(metadatum), path; kwargs...)
    end
    return nothing
end

# The record of `name` at the time index `n` of the file of `metadatum`, from -180ᵒ eastward and from the south pole
function regular_record(metadatum, n)
    Dataset(metadata_path(metadatum)) do dataset
        longitudes = mod.(Float64.(dataset["lon"][:]) .+ 180, 360) .- 180
        latitudes = Float64.(dataset["lat"][:])
        record = Float32.(coalesce.(dataset[DataWrangling.dataset_variable_name(metadatum)][:, :, n], NaN32))
        record[sortperm(longitudes), sortperm(latitudes)]
    end
end

#####
##### Mixed-layer depth climatology of de Boyer Montégut et al. (2004)
#####

"""
    BoyerMontegutMonthly()

The monthly climatology of the mixed-layer depth of [de Boyer Montégut et al. (2004)](@cite deboyermontegut2004mixed),
with a density threshold of 0.03 kg m⁻³ from 10 m (DR003), on a 2ᵒ grid. Provides `:mixed_layer_depth` in m.
"""
struct BoyerMontegutMonthly <: MonthlyObservation end

DataWrangling.all_dates(::BoyerMontegutMonthly, args...) = climatology_dates
DataWrangling.available_variables(::BoyerMontegutMonthly) = Dict(:mixed_layer_depth => "mld")
Base.size(::BoyerMontegutMonthly, variable) = (180, 90, 1)
DataWrangling.metadata_filename(::BoyerMontegutMonthly, name, date, region) = "mld_DR003_c1m_reg2.0.nc"
observation_url(::Metadatum{<:BoyerMontegutMonthly}) = "https://mld.ifremer.fr/data/mld_DR003_c1m_reg2.0.nc"

function DataWrangling.retrieve_data(metadatum::Metadatum{<:BoyerMontegutMonthly})
    record = regular_record(metadatum, month(metadatum.dates))
    record[record .> 1f8] .= NaN32
    return record
end

#####
##### HadISST sea ice concentration
#####

"""
    HadISSTMonthly()

Monthly sea ice concentration of the Hadley Centre Sea Ice and Sea Surface Temperature data set
[Rayner et al. (2003)](@cite rayner2003global) on a 1ᵒ grid, from 1870. Provides `:sea_ice_concentration`.
"""
struct HadISSTMonthly <: MonthlyObservation end

DataWrangling.all_dates(::HadISSTMonthly, args...) = DateTime(1870, 1, 1):Month(1):DateTime(2024, 12, 1)
DataWrangling.available_variables(::HadISSTMonthly) = Dict(:sea_ice_concentration => "sic")
Base.size(::HadISSTMonthly, variable) = (360, 180, 1)
DataWrangling.metadata_filename(::HadISSTMonthly, name, date, region) = "HadISST_ice.nc"
observation_url(::Metadatum{<:HadISSTMonthly}) = "https://www.metoffice.gov.uk/hadobs/hadisst/data/HadISST_ice.nc.gz"

function Downloads.download(metadata::Metadata{<:HadISSTMonthly}; kwargs...)
    path = metadata_path(first(metadata))
    @root if !isfile(path)
        download_with_retries(observation_url(first(metadata)), path * ".gz"; kwargs...)
        run(pipeline(`gunzip -c $(path * ".gz")`, stdout = path))
    end
    return nothing
end

function DataWrangling.retrieve_data(metadatum::Metadatum{<:HadISSTMonthly})
    Dataset(metadata_path(metadatum)) do dataset
        dates = dataset["time"][:]
        n = findfirst(d -> year(d) == year(metadatum.dates) && month(d) == month(metadatum.dates), dates)
        longitudes = mod.(Float64.(dataset["longitude"][:]) .+ 180, 360) .- 180
        latitudes = Float64.(dataset["latitude"][:])
        concentration = Float32.(coalesce.(dataset["sic"][:, :, n], NaN32))
        concentration[(concentration .< 0) .| (concentration .> 1)] .= NaN32
        concentration[sortperm(longitudes), sortperm(latitudes)]
    end
end

#####
##### NCEP reanalysis wind stress
#####

"""
    NCEPMonthly()

The monthly long-term mean wind stress on the surface of the NCEP/NCAR reanalysis [Kalnay et al. (1996)](@cite kalnay1996ncep),
positive eastward and northward. Provides `:zonal_wind_stress` and `:meridional_wind_stress` in N m⁻². The Gaussian
latitudes of the reanalysis are taken as equally spaced, which displaces them by less than 0.1ᵒ.
"""
struct NCEPMonthly <: MonthlyObservation end

const ncep_server = "https://psl.noaa.gov/thredds/fileServer/Datasets/ncep.reanalysis.derived/surface_gauss"

DataWrangling.all_dates(::NCEPMonthly, args...) = climatology_dates
DataWrangling.available_variables(::NCEPMonthly) = Dict(:zonal_wind_stress => "uflx", :meridional_wind_stress => "vflx")
Base.size(::NCEPMonthly, variable) = (192, 94, 1)
DataWrangling.metadata_filename(::NCEPMonthly, name, date, region) =
    string(DataWrangling.available_variables(NCEPMonthly())[name], ".sfc.mon.ltm.nc")
observation_url(metadatum::Metadatum{<:NCEPMonthly}) = string(ncep_server, "/", DataWrangling.metadata_filename(metadatum))

# The reanalysis fluxes are positive upward, out of the ocean
DataWrangling.retrieve_data(metadatum::Metadatum{<:NCEPMonthly}) = - regular_record(metadatum, month(metadatum.dates))

end # module
