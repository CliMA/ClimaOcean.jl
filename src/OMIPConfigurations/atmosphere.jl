# Shared atmosphere / land / radiation setup for OMIP simulations

"""
    omip_atmosphere(arch; forcing_dir, start_date, end_date, backend_size=30)

Set up JRA55 prescribed atmosphere, land (river runoff + iceberg calving),
and radiation. Returns `(atmosphere, land, radiation)`.
"""
function omip_atmosphere(arch;
                         forcing_dir,
                         start_date,
                         end_date,
                         backend_size = 30)

    dataset = MultiYearJRA55()
    kw = (; dir = forcing_dir, dataset, start_date, end_date,
            time_indices_in_memory = backend_size)

    atmosphere = JRA55PrescribedAtmosphere(arch; kw...)
    land       = JRA55PrescribedLand(arch; kw...)
    radiation  = JRA55PrescribedRadiation(arch; kw...)

    return atmosphere, land, radiation
end
