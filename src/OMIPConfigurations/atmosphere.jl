# Shared atmosphere setup for OMIP simulations

"""
    omip_atmosphere(arch; forcing_dir, start_date, end_date, backend_size=30)

Set up JRA55 prescribed atmosphere with river and iceberg forcing.
"""
function omip_atmosphere(arch;
                         forcing_dir,
                         start_date,
                         end_date,
                         backend_size = 30)

    dataset = MultiYearJRA55()
    backend = JRA55NetCDFBackend(backend_size)

    atmosphere = JRA55PrescribedAtmosphere(arch;
                                           dir = forcing_dir,
                                           dataset,
                                           backend,
                                           include_rivers_and_icebergs = true,
                                           start_date,
                                           end_date)

    radiation = Radiation()

    return atmosphere, radiation
end
