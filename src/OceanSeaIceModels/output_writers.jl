using JLD2

using Oceananigans.Utils: pretty_filesize
using Oceananigans.OutputWriters: checkpoint_path,
                                  cleanup_checkpoints,
                                  validate_checkpointed_properties

import Oceananigans.Fields: set!
import Oceananigans.OutputWriters: write_output!

wall_time = Ref(time_ns())

function set!(model::OSIMPA, checkpoint_file_path::AbstractString)
    addr = checkpointer_address(model)

    jldopen(checkpoint_file_path, "r") do file
        checkpointed_clock = file["$addr/clock"]

        # Update model clock
        set_clock!(model, checkpointed_clock)
    end

    # deal with model components

    for component in (model.ocean.model, model.atmosphere)
        set!(component, checkpoint_file_path)
    end

    return nothing
end


function write_output!(c::Checkpointer, model::OSIM)

    filepath = checkpoint_path(model.clock.iteration, c)
    c.verbose && @info "Checkpointing to file $filepath..."

    t1 = time_ns()

    # write OceanSeaIceModel
    mode = "w"
    write_output!(c, model, filepath, mode,
                  properties = default_checkpointed_properties(model))

    # write the OceanSeaIceModel components
    components = (model.atmosphere, model.ocean.model)

    mode = "a" # to append in the file already written above
    for component in components
        write_output!(c, component, filepath, mode,
                      properties = default_checkpointed_properties(component))
    end

    t2, sz = time_ns(), filesize(filepath)

    c.verbose && @info "Checkpointing done: time=$(prettytime((t2 - t1) * 1e-9)), size=$(pretty_filesize(sz))"

    c.cleanup && cleanup_checkpoints(c)

    return nothing
end
