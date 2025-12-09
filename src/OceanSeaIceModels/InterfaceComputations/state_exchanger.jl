
mutable struct StateExchanger{G, AS, OS, SS, AEX, OEX, SIE}
    exchange_grid :: G
    exchange_atmosphere_state :: AST
    exchange_ocean_state :: OS
    exchange_sea_ice_state :: SS
    atmosphere_exchanger :: AEX
    ocean_exchanger :: OEX
    sea_ice_exchanger :: SIE
end

mutable struct ExchangeAtmosphereState{FC, CF, CC}
    u  :: FC
    v  :: CF
    T  :: CC
    q  :: CC
    p  :: CC
    Qs :: CC
    Qℓ :: CC
    Mp :: CC
end

mutable struct ExchangeOceanState{FC, CF, CC}
    u :: FC
    v :: CF
    T :: CC
    S :: CC
end

mutable struct ExchangeSeaIceState{FC, CF, CC}
    u :: FC
    v :: CF
    h :: CC
    ℵ :: CC
end

ExchangeAtmosphereState(grid) = ExchangeAtmosphereState(Field{Center, Center, Nothing}(grid),
                                                        Field{Center, Center, Nothing}(grid),
                                                        Field{Center, Center, Nothing}(grid),
                                                        Field{Center, Center, Nothing}(grid),
                                                        Field{Center, Center, Nothing}(grid),
                                                        Field{Center, Center, Nothing}(grid),
                                                        Field{Center, Center, Nothing}(grid),
                                                        Field{Center, Center, Nothing}(grid))

# Note that Field location can also affect fractional index type.
# Here we assume that we know the location of Fields that will be interpolated.
fractional_index_type(FT, Topo) = FT
fractional_index_type(FT, ::Flat) = Nothing

StateExchanger(ocean::Simulation, ::Nothing) = nothing

function StateExchanger(ocean::Simulation, atmosphere)
    # TODO: generalize this
    exchange_grid = ocean.model.grid
    exchange_atmosphere_state = ExchangeAtmosphereState(exchange_grid)
    exchanger = atmosphere_exchanger(atmosphere, exchange_grid, exchange_atmosphere_state)

    return StateExchanger(ocean.model.grid, exchange_atmosphere_state, exchanger)
end

function atmosphere_exchanger(atmosphere::PrescribedAtmosphere, exchange_grid, exchange_atmosphere_state)
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

initialize!(exchanger::StateExchanger, ::Nothing) = nothing

function initialize!(exchanger::StateExchanger, atmosphere::PrescribedAtmosphere)
    atmos_grid = atmosphere.grid
    exchange_grid = exchanger.exchange_grid
    arch = architecture(exchange_grid)
    frac_indices = exchanger.atmosphere_exchanger
    kernel_parameters = interface_kernel_parameters(exchange_grid)
    launch!(arch, exchange_grid, kernel_parameters,
            _compute_fractional_indices!, frac_indices, exchange_grid, atmos_grid)
    return nothing
end

@kernel function _compute_fractional_indices!(indices_tuple, exchange_grid, atmos_grid)
    i, j = @index(Global, NTuple)
    kᴺ = size(exchange_grid, 3) # index of the top ocean cell
    X = _node(i, j, kᴺ + 1, exchange_grid, c, c, f)
    if topology(atmos_grid) == (Flat, Flat, Flat)
        fractional_indices_ij = FractionalIndices(nothing, nothing, nothing)
    else
        fractional_indices_ij = FractionalIndices(X, atmos_grid, c, c, c)
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
