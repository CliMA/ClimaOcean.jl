#####
##### OMIP evaluation figures, implemented in the extension loaded with Makie and GeoMakie
#####
##### Every figure takes one or more `OMIPOutput`s, drawing each as a column of panels or a line so that runs
##### compare side by side, and returns a `Figure`.
#####

"""
    omip_figures(outputs...; directory = ".")

Draw every OMIP figure for `outputs` and save them as PNG files in `directory`. Requires Makie and GeoMakie.
"""
function omip_figures end

"""
    sea_surface_temperature_bias(outputs...)

Time-mean sea surface temperature minus the World Ocean Atlas climatology.
"""
function sea_surface_temperature_bias end

"""
    sea_surface_salinity_bias(outputs...)

Time-mean sea surface salinity minus the World Ocean Atlas climatology.
"""
function sea_surface_salinity_bias end

"""
    sea_surface_height(outputs...)

Time-mean sea surface height, and its difference from the ECCO state estimate of 1992--2012 with both global means
removed.
"""
function sea_surface_height end

"""
    mixed_layer_depth(outputs...)

Summer minimum and winter maximum of the monthly mixed-layer depth, with the climatology of de Boyer Montégut et al.
(2004).
"""
function mixed_layer_depth end

"""
    sea_ice_concentration(outputs...)

March and September sea ice concentration in both hemispheres.
"""
function sea_ice_concentration end

"""
    sea_ice_concentration_bias(outputs...)

March and September sea ice concentration minus the HadISST climatology of 1979--2007, with the observed 15% edge.
"""
function sea_ice_concentration_bias end

"""
    surface_fluxes(outputs...)

Time-mean net surface heat and freshwater fluxes into the ocean.
"""
function surface_fluxes end

"""
    wind_stress(outputs...)

Time-mean zonal and meridional wind stress on the ocean, and their difference from the NCEP reanalysis climatology.
"""
function wind_stress end

"""
    sea_surface_height_variability(outputs...)

Root-mean-square variability of the sea surface height about its time mean.
"""
function sea_surface_height_variability end

"""
    sea_ice_extent(outputs...)

Monthly climatology of the Arctic and Antarctic sea ice extent, with NSIDC over the same forcing years.
"""
function sea_ice_extent end

"""
    sea_ice_area(outputs...)

Monthly climatology of the Arctic and Antarctic sea ice area, with NSIDC over the same forcing years.
"""
function sea_ice_area end

"""
    arctic_sea_ice_volume(outputs...)

Monthly climatology of the Arctic sea ice volume, with PIOMAS over the same forcing years.
"""
function arctic_sea_ice_volume end

"""
    sea_ice_area_time_series(outputs...)

Arctic and Antarctic sea ice area over the whole simulation.
"""
function sea_ice_area_time_series end

"""
    arctic_sea_ice_volume_time_series(outputs...)

Arctic sea ice volume over the whole simulation.
"""
function arctic_sea_ice_volume_time_series end

"""
    kinetic_energy(outputs...)

Global-mean kinetic energy over the whole simulation.
"""
function kinetic_energy end

"""
    global_drift(outputs...)

Drift of the global-mean temperature and salinity over the whole simulation.
"""
function global_drift end

"""
    mean_profiles(outputs...)

Time-mean horizontal-mean profiles of temperature and salinity.
"""
function mean_profiles end

"""
    zonal_mean_sections(outputs...)

Zonal-mean temperature and salinity with the World Ocean Atlas contours, and buoyancy with the initial contours.
"""
function zonal_mean_sections end

"""
    zonal_mean_drift(outputs...)

Zonal-mean temperature and salinity minus the World Ocean Atlas, and buoyancy minus its initial value.
"""
function zonal_mean_drift end

"""
    zonal_mean_mixed_layer_depth(outputs...)

Zonal mean of the summer minimum and winter maximum mixed-layer depth, with the climatology of de Boyer Montégut et
al. (2004).
"""
function zonal_mean_mixed_layer_depth end

"""
    drift_profiles(outputs...)

Drift of the horizontal-mean temperature and salinity as a function of depth and time.
"""
function drift_profiles end

"""
    strait_transports(outputs...)

Annual-mean volume transport through Bering Strait, Drake Passage and the Indonesian Throughflow, with observations.
"""
function strait_transports end

"""
    atlantic_overturning(outputs...)

Time-mean Atlantic meridional overturning streamfunction, including the Gent-McWilliams transport where present.
"""
function atlantic_overturning end

"""
    atlantic_overturning_at_rapid(outputs...)

Atlantic overturning streamfunction at 26.5ᵒN as a profile and as the six-monthly maximum below 500 m, with RAPID.
"""
function atlantic_overturning_at_rapid end

"""
    near_surface_currents(outputs...)

Time-mean speed, eastward and northward velocity of the top model cell.
"""
function near_surface_currents end

"""
    equatorial_undercurrent(outputs...)

Time-mean eastward velocity along 140ᵒW and along the equator in the Pacific.
"""
function equatorial_undercurrent end

"""
    sea_level_drift(outputs...)

Global-mean sea surface height over the whole simulation.
"""
function sea_level_drift end

"""
    content_conservation(outputs...)

Drift of the salt in the ocean, in the sea ice and in both, of the heat and the volume of the ocean, and of the water
in the ocean, sea ice and snow.
"""
function content_conservation end

"""
    barotropic_streamfunction(outputs...)

Time-mean barotropic streamfunction of the global ocean and of the North Atlantic.
"""
function barotropic_streamfunction end

"""
    arctic_freshwater_transports(outputs...)

Annual-mean liquid and solid freshwater transport through Fram and Davis Straits, with observations.
"""
function arctic_freshwater_transports end

"""
    zonal_rms_bias(outputs...)

Zonal root-mean-square of the time-mean temperature and salinity minus the World Ocean Atlas, and the part of it that
cancels in the zonal mean.
"""
function zonal_rms_bias end

"""
    rms_bias_profiles(outputs...)

Horizontal root-mean-square of the temperature and salinity minus the World Ocean Atlas, as a function of depth and
time.
"""
function rms_bias_profiles end
