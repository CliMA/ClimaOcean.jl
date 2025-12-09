function ComponentExchanger(atmosphere::PrescribedAtmosphere, grid) 

    regridder = atmosphere_regridder(atmosphere, grid)

    state = (; u  = Field{Center, Center, Nothing}(grid),
               v  = Field{Center, Center, Nothing}(grid),
               T  = Field{Center, Center, Nothing}(grid),
               p  = Field{Center, Center, Nothing}(grid),
               q  = Field{Center, Center, Nothing}(grid),
               Qs = Field{Center, Center, Nothing}(grid),
               Qℓ = Field{Center, Center, Nothing}(grid),
               Mp = Field{Center, Center, Nothing}(grid))

    return ComponentExchanger(state, regridder)
end

# Note that Field location can also affect fractional index type.
# Here we assume that we know the location of Fields that will be interpolated.
fractional_index_type(FT, Topo) = FT
fractional_index_type(FT, ::Flat) = Nothing

function atmosphere_regridder(atmosphere::PrescribedAtmosphere, exchange_grid)
    atmos_grid = atmosphere.grid
    arch = architecture(exchange_grid)
    Nx, Ny, Nz = size(exchange_grid)

    # Make a NamedTuple of fractional indices
    # Note: we could use an array of FractionalIndices. Instead, for compatbility
    # with Reactant we construct FractionalIndices on the fly in `interpolate_atmospheric_state`.
    FT = eltype(atmos_grid)
    TX, TY, TZ = topology(exchange_grid)
    fi = TX() isa Flat ? nothing : Field{Center, Center, Nothing}(exchange_grid, FT)
    fj = TY() isa Flat ? nothing : Field{Center, Center, Nothing}(exchange_grid, FT)
    frac_indices = (i=fi, j=fj) # no k needed, only horizontal interpolation

    return frac_indices
end

function initialize!(exchanger::ComponentExchanger, grid, atmosphere::PrescribedAtmosphere)

    frac_indices = exchanger.regridder
    atmos_grid = atmosphere.grid
    kernel_parameters = interface_kernel_parameters(grid)
    launch!(arch, grid, kernel_parameters, _compute_fractional_indices!, frac_indices, grid, atmos_grid)

    return nothing
end

@kernel function _compute_fractional_indices!(indices_tuple, exchange_grid, atmos_grid)
    i, j = @index(Global, NTuple)
    kᴺ = size(exchange_grid, 3) # index of the top ocean cell
    X = _node(i, j, kᴺ + 1, exchange_grid, Center(), Center(), Face())
    if topology(atmos_grid) == (Flat, Flat, Flat)
        fractional_indices_ij = FractionalIndices(nothing, nothing, nothing)
    else
        fractional_indices_ij = FractionalIndices(X, atmos_grid, Center(), Center(), Center())
    end
    fi = indices_tuple.i
    fj = indices_tuple.j
    @inbounds begin
        if !isnothing(fi)
            fi[i, j, 1] = fractional_indices_ij.i
        end

        if !isnothing(fj)
            fj[i, j, 1] = fractional_indices_ij.j
        end
    end
end
