# Supported datasets

ClimaOcean currently ships connectors for the following data products:

| Dataset            | Supported Variables                                      | Documentation Link                                                                                 |
|--------------------|-----------------------------------------------------------|-----------------------------------------------------------------------------------------------------|
| `ETOPO2022`        | [Supported variables](@ref dataset-etopo2022-vars)        | [NOAA ETOPO 2022 overview](https://www.ncei.noaa.gov/products/etopo-global-relief-model)           |
| `ECCO2Monthly`     | [Supported variables](@ref dataset-ecco2monthly-vars)     | [ECCO2 documentation](https://ecco.jpl.nasa.gov/products/all/)                                     |
| `ECCO2Daily`       | [Supported variables](@ref dataset-ecco2daily-vars)       | [ECCO2 documentation](https://ecco.jpl.nasa.gov/products/all/)                                     |
| `ECCO4Monthly`     | [Supported variables](@ref dataset-ecco4monthly-vars)     | [ECCO V4r4 product guide](https://ecco-group.org/products-ECCO-V4r4.htm)                           |
| `EN4Monthly`       | [Supported variables](@ref dataset-en4monthly-vars)       | [Met Office EN4 overview](https://www.metoffice.gov.uk/hadobs/en4/)                                |
| `GLORYSDaily`      | [Supported variables](@ref dataset-glorysdaily-vars)      | [Copernicus GLORYS product page](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description) |
| `GLORYSMonthly`    | [Supported variables](@ref dataset-glorysmonthly-vars)    | [Copernicus GLORYS product page](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description) |
| `RepeatYearJRA55`  | [Supported variables](@ref dataset-repeatyearjra55-vars)  | [JRA-55 Reanalysis](https://jra.kishou.go.jp/JRA-55/index_en.html)                                 |
| `MultiYearJRA55`   | [Supported variables](@ref dataset-multiyearjra55-vars)   | [JRA-55 Reanalysis](https://jra.kishou.go.jp/JRA-55/index_en.html)                                 |

## [Supported variables for ETOPO2022](@id dataset-etopo2022-vars)

- `:bottom_height` - Global bathymetry/topography on a 1 arc-minute grid.

## [Supported variables for ECCO2Monthly](@id dataset-ecco2monthly-vars)

- `:temperature` - Potential temperature on the ECCO2 three-dimensional grid (degC).
- `:salinity` - Practical salinity (psu).
- `:u_velocity` - Zonal velocity component (m s^-1).
- `:v_velocity` - Meridional velocity component (m s^-1).
- `:free_surface` - Sea surface height anomaly (m).
- `:sea_ice_thickness` - Effective sea-ice thickness (m).
- `:sea_ice_concentration` - Sea-ice area fraction (dimensionless).
- `:net_heat_flux` - Net surface heat flux into the ocean (W m^-2).

## [Supported variables for ECCO2Daily](@id dataset-ecco2daily-vars)

- `:temperature` - Potential temperature on the ECCO2 three-dimensional grid (degC).
- `:salinity` - Practical salinity (psu).
- `:u_velocity` - Zonal velocity component (m s^-1).
- `:v_velocity` - Meridional velocity component (m s^-1).
- `:free_surface` - Sea surface height anomaly (m).
- `:sea_ice_thickness` - Effective sea-ice thickness (m).
- `:sea_ice_concentration` - Sea-ice area fraction (dimensionless).
- `:net_heat_flux` - Net surface heat flux into the ocean (W m^-2).

## [Supported variables for ECCO4Monthly](@id dataset-ecco4monthly-vars)

- `:temperature` - Potential temperature on the ECCO V4 grid (degC).
- `:salinity` - Practical salinity (psu).
- `:u_velocity` - Zonal velocity component (m s^-1).
- `:v_velocity` - Meridional velocity component (m s^-1).
- `:free_surface` - Sea surface height anomaly (m).
- `:sea_ice_thickness` - Effective sea-ice thickness (m).
- `:sea_ice_concentration` - Sea-ice area fraction (dimensionless).
- `:net_heat_flux` - Net surface heat flux into the ocean (W m^-2).
- `:sensible_heat_flux` - Surface sensible heat flux (W m^-2).
- `:latent_heat_flux` - Surface latent heat flux (W m^-2).
- `:net_longwave` - Net longwave radiation at the surface (W m^-2).
- `:downwelling_shortwave` - Downward shortwave radiation at the surface (W m^-2).
- `:downwelling_longwave` - Downward longwave radiation at the surface (W m^-2).

## [Supported variables for EN4Monthly](@id dataset-en4monthly-vars)

- `:temperature` - Objective analyses of ocean temperature (degC).
- `:salinity` - Objective analyses of ocean salinity (psu).

## [Supported variables for GLORYSDaily](@id dataset-glorysdaily-vars)

- `:temperature` - Potential temperature (degC).
- `:salinity` - Practical salinity (psu).
- `:u_velocity` - Zonal velocity component (m s^-1).
- `:v_velocity` - Meridional velocity component (m s^-1).
- `:sea_ice_concentration` - Sea-ice area fraction (dimensionless).
- `:sea_ice_thickness` - Sea-ice thickness (m).
- `:sea_ice_u_velocity` - Zonal sea-ice drift (m s^-1).
- `:sea_ice_v_velocity` - Meridional sea-ice drift (m s^-1).
- `:free_surface` - Sea surface height (m).
- `:depth` - Static bathymetry/depth (m).

## [Supported variables for GLORYSMonthly](@id dataset-glorysmonthly-vars)

- `:temperature` - Potential temperature (degC).
- `:salinity` - Practical salinity (psu).
- `:u_velocity` - Zonal velocity component (m s^-1).
- `:v_velocity` - Meridional velocity component (m s^-1).
- `:sea_ice_concentration` - Sea-ice area fraction (dimensionless).
- `:sea_ice_thickness` - Sea-ice thickness (m).
- `:sea_ice_u_velocity` - Zonal sea-ice drift (m s^-1).
- `:sea_ice_v_velocity` - Meridional sea-ice drift (m s^-1).
- `:free_surface` - Sea surface height (m).
- `:depth` - Static bathymetry/depth (m).

## [Supported variables for RepeatYearJRA55](@id dataset-repeatyearjra55-vars)

- `:temperature` - 2 m air temperature (K).
- `:specific_humidity` - 2 m specific humidity (kg kg^-1).
- `:eastward_velocity` - 10 m eastward wind (m s^-1).
- `:northward_velocity` - 10 m northward wind (m s^-1).
- `:sea_level_pressure` - Sea-level pressure (Pa).
- `:downwelling_shortwave_radiation` - Downward shortwave radiation (W m^-2).
- `:downwelling_longwave_radiation` - Downward longwave radiation (W m^-2).
- `:rain_freshwater_flux` - Liquid precipitation flux (kg m^-2 s^-1).
- `:snow_freshwater_flux` - Solid precipitation flux (kg m^-2 s^-1).
- `:river_freshwater_flux` - River discharge flux (kg m^-2 s^-1).
- `:iceberg_freshwater_flux` - Iceberg calving flux (kg m^-2 s^-1).

## [Supported variables for MultiYearJRA55](@id dataset-multiyearjra55-vars)

- `:temperature` - 2 m air temperature (K).
- `:specific_humidity` - 2 m specific humidity (kg kg^-1).
- `:eastward_velocity` - 10 m eastward wind (m s^-1).
- `:northward_velocity` - 10 m northward wind (m s^-1).
- `:sea_level_pressure` - Sea-level pressure (Pa).
- `:downwelling_shortwave_radiation` - Downwelling shortwave radiation (W m^-2).
- `:downwelling_longwave_radiation` - Downwelling longwave radiation (W m^-2).
- `:rain_freshwater_flux` - Precipitation flux from liquid (kg m^-2 s^-1).
- `:snow_freshwater_flux` - Precipitation flux from snow/ice (kg m^-2 s^-1).
- `:river_freshwater_flux` - River discharge flux (kg m^-2 s^-1).
- `:iceberg_freshwater_flux` - Iceberg calving flux (kg m^-2 s^-1).
