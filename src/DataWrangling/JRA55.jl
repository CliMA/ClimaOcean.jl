module JRA55

using Oceananigans
using Oceananigans.Units
 
using Oceananigans.Architectures: arch_array
using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: child_architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: λnodes, φnodes, on_architecture
using Oceananigans.Fields: interpolate!
using Oceananigans.OutputReaders: Cyclical, TotallyInMemory, AbstractInMemoryBackend, FlavorOfFTS, time_indices

using ClimaOcean.OceanSeaIceModels:
    PrescribedAtmosphere,
    TwoStreamDownwellingRadiation

using CUDA: @allowscalar

using NCDatasets
using JLD2 
using Dates

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend, update_field_time_series!
using Downloads: download

# A list of all variables provided in the JRA55 dataset:
JRA55_variable_names = (:river_freshwater_flux,
                        :rain_freshwater_flux,
                        :snow_freshwater_flux,
                        :iceberg_freshwater_flux,
                        :specific_humidity,
                        :sea_level_pressure,
                        :relative_humidity,
                        :downwelling_longwave_radiation,
                        :downwelling_shortwave_radiation,
                        :temperature,
                        :eastward_velocity,
                        :northward_velocity)

filenames = Dict(
    :river_freshwater_flux           => "RYF.friver.1990_1991.nc",   # Freshwater fluxes from rivers
    :rain_freshwater_flux            => "RYF.prra.1990_1991.nc",     # Freshwater flux from rainfall
    :snow_freshwater_flux            => "RYF.prsn.1990_1991.nc",     # Freshwater flux from snowfall
    :iceberg_freshwater_flux         => "RYF.licalvf.1990_1991.nc",  # Freshwater flux from calving icebergs
    :specific_humidity               => "RYF.huss.1990_1991.nc",     # Surface specific humidity
    :sea_level_pressure              => "RYF.psl.1990_1991.nc",      # Sea level pressure
    :relative_humidity               => "RYF.rhuss.1990_1991.nc",    # Surface relative humidity
    :downwelling_longwave_radiation  => "RYF.rlds.1990_1991.nc",     # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "RYF.rsds.1990_1991.nc",     # Downwelling shortwave radiation
    :temperature                     => "RYF.tas.1990_1991.nc",      # Near-surface air temperature
    :eastward_velocity               => "RYF.uas.1990_1991.nc",      # Eastward near-surface wind
    :northward_velocity              => "RYF.vas.1990_1991.nc",      # Northward near-surface wind
)

jra55_short_names = Dict(
    :river_freshwater_flux           => "friver",   # Freshwater fluxes from rivers
    :rain_freshwater_flux            => "prra",     # Freshwater flux from rainfall
    :snow_freshwater_flux            => "prsn",     # Freshwater flux from snowfall
    :iceberg_freshwater_flux         => "licalvf",  # Freshwater flux from calving icebergs
    :specific_humidity               => "huss",     # Surface specific humidity
    :sea_level_pressure              => "psl",      # Sea level pressure
    :relative_humidity               => "rhuss",    # Surface relative humidity
    :downwelling_longwave_radiation  => "rlds",     # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "rsds",     # Downwelling shortwave radiation
    :temperature                     => "tas",      # Near-surface air temperature
    :eastward_velocity               => "uas",      # Eastward near-surface wind
    :northward_velocity              => "vas",      # Northward near-surface wind
)

field_time_series_short_names = Dict(
    :river_freshwater_flux           => "Fri", # Freshwater fluxes from rivers
    :rain_freshwater_flux            => "Fra", # Freshwater flux from rainfall
    :snow_freshwater_flux            => "Fsn", # Freshwater flux from snowfall
    :iceberg_freshwater_flux         => "Fic", # Freshwater flux from calving icebergs
    :specific_humidity               => "qa",  # Surface specific humidity
    :sea_level_pressure              => "pa",  # Sea level pressure
    :relative_humidity               => "rh",  # Surface relative humidity
    :downwelling_longwave_radiation  => "Ql",  # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "Qs",  # Downwelling shortwave radiation
    :temperature                     => "Ta",  # Near-surface air temperature
    :eastward_velocity               => "ua",  # Eastward near-surface wind
    :northward_velocity              => "va",  # Northward near-surface wind
)

urls = Dict(
    :shortwave_radiation => "https://www.dropbox.com/scl/fi/z6fkvmd9oe3ycmaxta131/" *
                            "RYF.rsds.1990_1991.nc?rlkey=r7q6zcbj6a4fxsq0f8th7c4tc&dl=0",

    :river_freshwater_flux => "https://www.dropbox.com/scl/fi/21ggl4p74k4zvbf04nb67/" * 
                              "RYF.friver.1990_1991.nc?rlkey=ny2qcjkk1cfijmwyqxsfm68fz&dl=0",

    :rain_freshwater_flux => "https://www.dropbox.com/scl/fi/5icl1gbd7f5hvyn656kjq/" *
                             "RYF.prra.1990_1991.nc?rlkey=iifyjm4ppwyd8ztcek4dtx0k8&dl=0",

    :snow_freshwater_flux => "https://www.dropbox.com/scl/fi/1r4ajjzb3643z93ads4x4/" *
                             "RYF.prsn.1990_1991.nc?rlkey=auyqpwn060cvy4w01a2yskfah&dl=0",

    :iceberg_freshwater_flux => "https://www.dropbox.com/scl/fi/44nc5y27ohvif7lkvpyv0/" *
                                "RYF.licalvf.1990_1991.nc?rlkey=w7rqu48y2baw1efmgrnmym0jk&dl=0",

    :specific_humidity => "https://www.dropbox.com/scl/fi/66z6ymfr4ghkynizydc29/" *
                          "RYF.huss.1990_1991.nc?rlkey=107yq04aew8lrmfyorj68v4td&dl=0",

    :sea_level_pressure => "https://www.dropbox.com/scl/fi/0fk332027oru1iiseykgp/" *
                           "RYF.psl.1990_1991.nc?rlkey=4xpr9uah741483aukok6d7ctt&dl=0",

    :relative_humidity => "https://www.dropbox.com/scl/fi/1agwsp0lzvntuyf8bm9la/" *
                          "RYF.rhuss.1990_1991.nc?rlkey=8cd0vs7iy1rw58b9pc9t68gtz&dl=0",

    :downwelling_longwave_radiation  => "https://www.dropbox.com/scl/fi/y6r62szkirrivua5nqq61/" *
                                        "RYF.rlds.1990_1991.nc?rlkey=wt9yq3cyrvs2rbowoirf4nkum&dl=0",

    :downwelling_shortwave_radiation => "https://www.dropbox.com/scl/fi/z6fkvmd9oe3ycmaxta131/" *
                                        "RYF.rsds.1990_1991.nc?rlkey=r7q6zcbj6a4fxsq0f8th7c4tc&dl=0",

    :temperature => "https://www.dropbox.com/scl/fi/fpl0npwi476w635g6lke9/" *
                    "RYF.tas.1990_1991.nc?rlkey=0skb9pe6lgbfbiaoybe7m945s&dl=0",

    :eastward_velocity => "https://www.dropbox.com/scl/fi/86wetpqla2x97isp8092g/" *
                          "RYF.uas.1990_1991.nc?rlkey=rcaf18sh1yz0v9g4hjm1249j0&dl=0",

    :northward_velocity => "https://www.dropbox.com/scl/fi/d38sflo9ddljstd5jwgml/" *
                           "RYF.vas.1990_1991.nc?rlkey=f9y3e57kx8xrb40gbstarf0x6&dl=0",
)

const AA = Oceananigans.Architectures.AbstractArchitecture

JRA55_prescribed_atmosphere(time_indices=Colon(); kw...) =
    JRA55_prescribed_atmosphere(CPU(), time_indices; kw...)

JRA55_prescribed_atmosphere(arch::Distributed, time_indices=Colon(); kw...) =
    JRA55_prescribed_atmosphere(child_architecture(arch), time_indices; kw...)

# TODO: allow the user to pass dates
function JRA55_prescribed_atmosphere(architecture::AA, time_indices=Colon();
                                     backend = nothing,
                                     time_indexing = Cyclical(),
                                     measurement_height = 10,  # meters
                                     other_kw...)

    if isnothing(backend) # apply a default
        Ni = try
            length(time_indices)
        catch
            Inf
        end

        # Manufacture a default for the number of fields to keep InMemory
        Nf = min(24, Ni)
        backend = NetCDFBackend(Nf)
    end

    kw = (; time_indices, time_indexing, backend, architecture)
    kw = merge(kw, other_kw) 

    ua  = reanalysis_field_time_series(:eastward_velocity;               kw...)
    va  = reanalysis_field_time_series(:northward_velocity;              kw...)
    Ta  = reanalysis_field_time_series(:temperature;                     kw...)
    qa  = reanalysis_field_time_series(:specific_humidity;               kw...)
    ra  = reanalysis_field_time_series(:relative_humidity;               kw...)
    pa  = reanalysis_field_time_series(:sea_level_pressure;              kw...)
    Fra = reanalysis_field_time_series(:rain_freshwater_flux;            kw...)
    Fsn = reanalysis_field_time_series(:snow_freshwater_flux;            kw...)
    Fri = reanalysis_field_time_series(:river_freshwater_flux;           kw...)
    Fic = reanalysis_field_time_series(:iceberg_freshwater_flux;         kw...)
    Ql  = reanalysis_field_time_series(:downwelling_longwave_radiation;  kw...)
    Qs  = reanalysis_field_time_series(:downwelling_shortwave_radiation; kw...)

    times = ua.times

    velocities = (u = ua,
                  v = va)

    tracers = (T = Ta,
               q = qa,
               r = ra)

    freshwater_flux = (rain = Fra,
                       snow = Fsn,
                       rivers = Fri,
                       icebergs = Fic)
                       
    pressure = pa

    downwelling_radiation = TwoStreamDownwellingRadiation(shortwave=Qs, longwave=Ql)

    FT = eltype(ua)
    measurement_height = convert(FT, measurement_height)

    atmosphere = PrescribedAtmosphere(times, FT;
                                      velocities,
                                      freshwater_flux,
                                      tracers,
                                      downwelling_radiation,
                                      measurement_height,
                                      pressure)

    return atmosphere
end

end # module

