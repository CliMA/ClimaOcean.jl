module Simulations

import Oceananigans: run!

function run!(coupled_model::OceanSeaIceModel; schedule,
    dir = ".",
    prefix = "checkpoint",
    overwrite_existing = false,
    verbose = false,
    cleanup = false,
    properties = default_checkpointed_properties(coupled_model.ocean.model))

    function run!(sim; pickup=false)

        start_run = time_ns()
    
        if we_want_to_pickup(pickup)
            set!(sim, pickup)
        end
    
        sim.initialized = false
        sim.running = true
        sim.run_wall_time = 0.0
    
        while sim.running
            time_step!(sim)
        end
    
        for callback in values(sim.callbacks) 
            finalize!(callback, sim)
        end
    
        # Increment the wall clock
        end_run = time_ns()
        sim.run_wall_time += 1e-9 * (end_run - start_run)
    
        return nothing
    end
    
