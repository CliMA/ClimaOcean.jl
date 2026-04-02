"""
    add_omip_transport_diagnostics!(simulation; output_dir, filename_prefix, field_averaging_interval)

Add opt-in transport diagnostics: MOC streamfunction, meridional heat/salt
transports, and key strait mass transports.

These are Priority 1 in the OMIP protocol (Griffies et al. 2016) but require
basin masks, section definitions, and integration infrastructure.
"""
function add_omip_transport_diagnostics!(simulation;
                                         output_dir = ".",
                                         filename_prefix = "omip",
                                         field_averaging_interval = 15days)

    @warn "Transport diagnostics (transports=true) not yet implemented. " *
          "MOC, heat transports, and strait transports will be added in a future release."

    return nothing
end
