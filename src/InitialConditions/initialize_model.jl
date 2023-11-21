using ClimaOcean.ECCO2: ecco2_tracer_fields, adjusted_ecco_field, ecco2_location

"""
    initialize!(model;
                filename = "./data/adjusted_ecco_tracers.nc", 
                kwargs...)

Initialize `model`. The keyword arguments `kwargs...` take the form `name=data`,
where `name` refers to one of the fields of `model.velocities` or `model.tracers`, 
and the `data` may be
 (1) an array
 (2) a function with arguments `(x, y, z)` for 3D fields, `(x, y)` for 2D fields and `(x)` for 1D fields
 (3) a symbol corresponding to a fldname of `ecco2_tracer_fields`
 (4) a checkpoint path containing the checkpointed field
The keyword argument `filename` is the path to a netcdf file containing the adjusted ecco fields.
If the file does not exist it will be created
"""
function initialize!(model;
                     filename = "./data/adjusted_ecco_tracers.nc", 
                     kwargs...)
    
    grid = model.grid
    arch = architecture(grid)
    
    custom_fields     = Dict()
    ecco2_fields      = Dict()
    checkpoint_fields = Dict()
    
    # Differentiate between custom initialization and ECCO2 initialization
    for (fldname, value) in kwargs
        if value ∈ keys(ecco2_tracer_fields)
            ecco2_fields[fldname] = ecco2_tracer_fields[value]
        elseif value isa String
            checkpoint_fields[fldname] = value
        else
            custom_fields[fldname] = value
        end
    end

    # Custom fields not present in the ECCO dataset
    if !isempty(custom_fields)
        set!(model; custom_fields...)
    end

    # Fields initialized from checkpointer
    if isempty!(checkpoint_fields)
        for (fldname, path) in checkpoint_fields
            data = jldopen(path)[string(fldname) * "/data"]
            set!(model, fldname => data)
        end
    end

    # Fields initialized from ECCO2
    if !isempty(ecco2_fields)
        mask = ecco2_center_mask(architecture(grid))
    
        for (fldname, variable_name) in ecco2_fields
            f = adjusted_ecco_field(variable_name; 
                                    architecture = arch,
                                    filename,
                                    mask)

            f_grid = Field(ecco2_location[variable_name], grid)   

            three_dimensional_regrid!(f_grid, f)

            set!(model; fldname => f_grid)
        end
    end

    return nothing
end
