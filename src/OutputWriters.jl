module OutputWriters

using ClimaOcean.OceanSeaIceModels: OceanSeaIceModel

using Oceananigans.OutputWriters: default_checkpointed_properties

import Oceananigans: Checkpointer

function Checkpointer(coupled_model::OceanSeaIceModel; schedule,
                      dir = ".",
                      prefix = "checkpoint",
                      overwrite_existing = false,
                      verbose = false,
                      cleanup = false,
                      properties = default_checkpointed_properties(coupled_model.ocean.model))

    @info "I went in your new method"

    #=
    # Certain properties are required for `set!` to pickup from a checkpoint.
    required_properties = [:grid, :particles, :clock]

    if has_ab2_timestepper(model)
        push!(required_properties, :timestepper)
    end

    for rp in required_properties
        if rp ∉ properties
            @warn "$rp is required for checkpointing. It will be added to checkpointed properties"
            push!(properties, rp)
        end
    end

    for p in properties
        p isa Symbol || error("Property $p to be checkpointed must be a Symbol.")
        p ∉ propertynames(model) && error("Cannot checkpoint $p, it is not a model property!")

        if (p ∉ required_properties) && has_reference(Function, getproperty(model, p))
            @warn "model.$p contains a function somewhere in its hierarchy and will not be checkpointed."
            filter!(e -> e != p, properties)
        end
    end

    mkpath(dir)

    return Checkpointer(schedule, dir, prefix, properties, overwrite_existing, verbose, cleanup)
    =#
    return nothing
end

end # module