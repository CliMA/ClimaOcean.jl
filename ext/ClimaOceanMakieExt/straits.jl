#####
##### Transports through straits
#####

"""
A section of constant `j` (`axis = :v`, crossed by the meridional transport) or constant `i` (`axis = :u`, crossed by
the zonal transport) over the index ranges `i` and `j`.
"""
struct StraitSection
    i :: UnitRange{Int}
    j :: UnitRange{Int}
    axis :: Symbol
end

# Sections on eORCA1 whose wet cells form one run from coast to coast
const eorca1_sections = (bering = StraitSection(111:117, 251:251, :v),
                         drake = StraitSection(220:220, 53:71, :u),
                         indonesian_throughflow = StraitSection(38:57, 130:130, :v),
                         fram = StraitSection(267:277, 275:275, :v),
                         davis = StraitSection(233:240, 254:254, :v))

# Elsewhere, sections are the model rows and columns closest to these boxes, in degrees East and North
const strait_boxes = (bering = (longitudes = (187, 196), latitude = 66.0, axis = :v),
                      drake = (longitude = 293, latitudes = (-68, -54), axis = :u),
                      indonesian_throughflow = (longitudes = (110, 130), latitude = -8.5, axis = :v),
                      fram = (longitudes = (340, 375), latitude = 79.0, axis = :v),
                      davis = (longitudes = (298, 307), latitude = 66.5, axis = :v))

in_longitude_range(λ, (λ₁, λ₂)) = mod(λ - λ₁, 360) ≤ λ₂ - λ₁

function geographic_section(output, box)
    λ, φ = center_longitudes(output), center_latitudes(output)
    Nx, Ny = size(λ)

    if box.axis == :v
        columns(j) = filter(i -> in_longitude_range(λ[i, j], box.longitudes), 1:Nx)
        j = argmin(j -> isempty(columns(j)) ? Inf : mean(abs.(φ[columns(j), j] .- box.latitude)), 1:Ny)
        i = columns(j)
        return StraitSection(minimum(i):maximum(i), j:j, :v)
    else
        i = argmin(i -> longitude_distance(λ[i, Ny ÷ 2], box.longitude), 1:Nx)
        j = filter(j -> box.latitudes[1] ≤ φ[i, j] ≤ box.latitudes[2], 1:Ny)
        return StraitSection(i:i, minimum(j):maximum(j), :u)
    end
end

"""
The strait sections of the grid of `output`: those of `eorca1_sections` on eORCA1, and the rows and columns
closest to `strait_boxes` on other grids, where they are not checked to close from coast to coast.
"""
strait_sections(output) = cached(output, :strait_sections) do
    size(model_grid(output))[1:2] == (360, 296) ? eorca1_sections :
        map(box -> geographic_section(output, box), strait_boxes)
end

section_transport(transport, section) = sum(transport[i, j, k] for i in section.i, j in section.j, k in axes(transport, 3))

"""
The volume transport [Sv] through every strait section at every record of the three-dimensional output, positive
northward or eastward.
"""
strait_transport_series(output) = cached(output, :strait_transport_series) do
    u = timeseries(output, "fields", "uvol")
    v = timeseries(output, "fields", "vvol")
    sections = strait_sections(output)
    Nt = length(v.times)
    transports = map(_ -> zeros(Nt), sections)

    for n in 1:Nt
        uⁿ, vⁿ = record_array(u, n), record_array(v, n)
        for (name, section) in pairs(sections)
            transports[name][n] = section_transport(section.axis == :v ? vⁿ : uⁿ, section) / 1e6
        end
    end

    merge(transports, (; years = in_years(v.times)))
end

const reference_salinity = 34.8
const cubic_kilometers_per_year = seconds_per_year / 1e9

"""
The liquid freshwater transport ``∫ v (S★ - S) / S★ dA`` and the solid one ``∫ ℵ h vⁱ ρⁱ/ρʷ (S★ - Sⁱ) / S★ dx``
[km³ yr⁻¹] through Fram and Davis Straits, relative to the Arctic reference salinity ``S★ = 34.8``.
"""
arctic_freshwater_series(output) = cached(output, :arctic_freshwater_series) do
    sections = (; fram = strait_sections(output).fram, davis = strait_sections(output).davis)
    grid = model_grid(output)
    Nz = size(grid, 3)
    mask = ocean_mask(output)

    v = timeseries(output, "fields", "vvol")
    S = timeseries(output, "fields", "so")
    liquid = map(_ -> zeros(length(v.times)), sections)

    for n in eachindex(v.times)
        vⁿ, Sⁿ = record_array(v, n), record_array(S, n)
        for (name, section) in pairs(sections)
            liquid[name][n] = sum(vⁿ[i, j, k] * (reference_salinity - (Sⁿ[i, j-1, k] + Sⁿ[i, j, k]) / 2) / reference_salinity
                                  for i in section.i, j in section.j, k in 1:Nz
                                  if mask[i, j-1, k] && mask[i, j, k]) * cubic_kilometers_per_year
        end
    end

    vⁱ = timeseries(output, "surface", "siv")
    ℵ = timeseries(output, "surface", "siconc")
    h = timeseries(output, "surface", "sithick")
    freshwater_fraction = sea_ice_density / 1000 * (reference_salinity - sea_ice_salinity) / reference_salinity
    solid = map(_ -> zeros(length(vⁱ.times)), sections)

    for n in eachindex(vⁱ.times)
        vⁱⁿ, ℵⁿ, hⁿ = record_array(vⁱ, n), record_array(ℵ, n), record_array(h, n)
        for (name, section) in pairs(sections)
            solid[name][n] = sum(vⁱⁿ[i, j, 1] * (ℵⁿ[i, j-1, 1] * hⁿ[i, j-1, 1] + ℵⁿ[i, j, 1] * hⁿ[i, j, 1]) / 2 *
                                 Δxᶜᶠᶜ(i, j, Nz, grid) * freshwater_fraction
                                 for i in section.i, j in section.j
                                 if mask[i, j-1, Nz] && mask[i, j, Nz]) * cubic_kilometers_per_year
        end
    end

    (; liquid, solid, liquid_years = in_years(v.times), solid_years = in_years(vⁱ.times))
end

function annual_means(years, values)
    year_index = floor.(Int, years)
    centers = sort(unique(year_index))
    return centers .+ 0.5, [mean(values[year_index .== y]) for y in centers]
end
