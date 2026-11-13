#####
##### World Ocean Atlas fields on the TEOS-10 variables of the ocean model
#####

host_grid(grid) = on_architecture(cpu_architecture(architecture(grid)), grid)

# Sea pressure [dbar] of a cell centered at depth -z, to the accuracy of the hydrostatic 1 dbar per meter
@inline sea_pressure(z) = max(0, -z)

"""
    convert_to_teos10!(T, S, atlas = SAARAtlas())

Convert the in-situ temperature `T` [ᵒC] and Practical Salinity `S` of the World Ocean Atlas, set on an ocean grid,
to Conservative Temperature [ᵒC] and Absolute Salinity [g kg⁻¹] in place.
"""
function convert_to_teos10!(T, S, atlas = SAARAtlas())
    grid = host_grid(T.grid)
    Nx, Ny, Nz = size(grid)
    Tʰ = Array(interior(T))
    Sʰ = Array(interior(S))

    Threads.@threads for j in 1:Ny
        for k in 1:Nz, i in 1:Nx
            Tⁱ, Sᴾ = Tʰ[i, j, k], Sʰ[i, j, k]
            (isnan(Tⁱ) || isnan(Sᴾ)) && continue
            λ = λnode(i, j, k, grid, Center(), Center(), Center())
            φ = φnode(i, j, k, grid, Center(), Center(), Center())
            p = sea_pressure(znode(i, j, k, grid, Center(), Center(), Center()))
            Sᴬ = absolute_salinity(atlas, Sᴾ, p, λ, φ)
            Tʰ[i, j, k] = Θ_from_T(Sᴬ, Tⁱ, p)
            Sʰ[i, j, k] = Sᴬ
        end
    end

    set!(T, Tʰ)
    set!(S, Sʰ)

    return T, S
end

"""
    convert_to_absolute_salinity!(S, atlas = SAARAtlas())

Convert every record of the Practical Salinity time series `S` to Absolute Salinity [g kg⁻¹] in place.
"""
function convert_to_absolute_salinity!(S::FieldTimeSeries, atlas = SAARAtlas())
    grid = host_grid(S.grid)
    Nx, Ny, Nz = size(grid)

    for n in 1:length(S.times)
        Sʰ = Array(interior(S[n]))
        Threads.@threads for j in 1:Ny
            for k in 1:Nz, i in 1:Nx
                Sᴾ = Sʰ[i, j, k]
                isnan(Sᴾ) && continue
                λ = λnode(i, j, k, grid, Center(), Center(), Center())
                φ = φnode(i, j, k, grid, Center(), Center(), Center())
                p = sea_pressure(znode(i, j, k, grid, Center(), Center(), Center()))
                Sʰ[i, j, k] = absolute_salinity(atlas, Sᴾ, p, λ, φ)
            end
        end
        set!(S[n], Sʰ)
    end

    return S
end

"""
    flatten_polar_cap!(T, S, latitude)

Replace `T` and `S` north of `latitude` with their mean over the wet cells of each level there. Does nothing
when `latitude` is `nothing`.
"""
flatten_polar_cap!(T, S, ::Nothing) = nothing

function flatten_polar_cap!(T, S, latitude)
    grid = host_grid(T.grid)
    Nx, Ny, Nz = size(grid)
    Tʰ = Array(interior(T))
    Sʰ = Array(interior(S))

    cap = [φnode(i, j, 1, grid, Center(), Center(), Center()) > latitude for i in 1:Nx, j in 1:Ny]

    for k in 1:Nz
        wet = [cap[i, j] && !inactive_node(i, j, k, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
        any(wet) || continue
        view(Tʰ, :, :, k)[wet] .= sum(view(Tʰ, :, :, k)[wet]) / count(wet)
        view(Sʰ, :, :, k)[wet] .= sum(view(Sʰ, :, :, k)[wet]) / count(wet)
    end

    set!(T, Tʰ)
    set!(S, Sʰ)

    return nothing
end
