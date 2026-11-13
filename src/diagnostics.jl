#####
##### OMIP output
#####

using NumericalEarth.Diagnostics: frazil_temperature_flux, net_ocean_temperature_flux, sea_ice_ocean_temperature_flux,
                                  atmosphere_ocean_temperature_flux, net_ocean_salinity_flux, sea_ice_ocean_salinity_flux

# Volume and tracer transports through the cell faces, for streamfunctions and section transports offline
@inline zonal_volume_flux(i, j, k, grid, u) = @inbounds u[i, j, k] * Axᶠᶜᶜ(i, j, k, grid)
@inline meridional_volume_flux(i, j, k, grid, v) = @inbounds v[i, j, k] * Ayᶜᶠᶜ(i, j, k, grid)
@inline zonal_tracer_transport(i, j, k, grid, u, c) = @inbounds u[i, j, k] * Axᶠᶜᶜ(i, j, k, grid) * ℑxᶠᵃᵃ(i, j, k, grid, c)
@inline meridional_tracer_transport(i, j, k, grid, v, c) = @inbounds v[i, j, k] * Ayᶜᶠᶜ(i, j, k, grid) * ℑyᵃᶠᵃ(i, j, k, grid, c)

@inline squared(i, j, k, grid, ψ) = @inbounds ψ[i, j, k]^2
@inline u²ᶜᶜᶜ(i, j, k, grid, u) = ℑxᶜᵃᵃ(i, j, k, grid, squared, u)
@inline v²ᶜᶜᶜ(i, j, k, grid, v) = ℑyᵃᶜᵃ(i, j, k, grid, squared, v)
@inline unit_volume(i, j, k, grid) = one(grid)

# The meridional velocity of the Gent-McWilliams skew flux, v★ = - ∂z(κ ϵ Sʸ)
@inline function bolus_meridional_velocity(i, j, k, grid, clock, closure, buoyancy, fields)
    return - δzᵃᵃᶜ(i, j, k, grid, κ_ϵSyᶜᶠᶠ, clock, closure, closure.κ_skew, buoyancy, fields) * Δz⁻¹ᶜᶠᶜ(i, j, k, grid)
end

bolus_meridional_volume_flux(model, ::Nothing) = nothing

function bolus_meridional_volume_flux(model, closure::IsopycnalSkewSymmetricDiffusivity)
    grid = model.grid
    v★ = KernelFunctionOperation{Center, Face, Center}(bolus_meridional_velocity, grid, model.clock, closure,
                                                        model.buoyancy, fields(model))
    return KernelFunctionOperation{Center, Face, Center}(meridional_volume_flux, grid, v★)
end

function isopycnal_closure(model)
    n = findfirst(c -> c isa IsopycnalSkewSymmetricDiffusivity, model.closure)
    return isnothing(n) ? nothing : model.closure[n]
end

"""
    add_omip_diagnostics!(simulation;
                          output_directory = ".",
                          prefix = "omip",
                          surface_averaging_interval = 5days,
                          field_averaging_interval = 15days,
                          averages_interval = 5days,
                          checkpoint_interval = 720days,
                          file_splitting_interval = 360days)

Write the OMIP output of the coupled `simulation` to `output_directory`, in files named after `prefix`:

- `<prefix>_surface`: surface fields averaged over `surface_averaging_interval`: temperature `tos`, salinity `sos`,
  free surface `zos`, velocities `uos` and `vos`, the squares `tossq`, `sossq` and `zossq`, mixed-layer depth `mlotst`,
  wind stress `tauuo` and `tauvo`, the net heat and freshwater fluxes `hfds` and `wfo`, sensible and latent heat
  `hfss` and `hfls`, the components of the surface temperature and salinity fluxes, and the sea ice concentration
  `siconc`, thickness `sithick`, velocities `siu` and `siv`, surface temperature `sitemptop` and snow thickness
  `sisnthick`;
- `<prefix>_fields`: three-dimensional fields averaged over `field_averaging_interval`: temperature `to`, salinity
  `so`, velocities `uo`, `vo` and `wo`, buoyancy `bo`, turbulent kinetic energy `tke`, the squares `uosq` and `vosq`,
  the volume transports `uvol` and `vvol`, the heat and salt transports `uvolto`, `vvolto`, `uvolso` and `vvolso`,
  and, with a Gent-McWilliams closure, its meridional volume transport `vvolgm`;
- `<prefix>_averages`: global quantities averaged over `averages_interval`: mean sea level `zosga`, horizontal-mean
  profiles `to_h`, `so_h` and `bo_h`, and the integrals of volume `voco`, salinity `soco` and temperature `hoco` of
  the ocean and of the volumes of sea ice `sivol` and snow `snvol`;
- `<prefix>_checkpoint`: checkpoints every `checkpoint_interval`, from which `run!(simulation; pickup = true)` restarts.
"""
function add_omip_diagnostics!(simulation;
                               output_directory = ".",
                               prefix = "omip",
                               surface_averaging_interval = 5days,
                               field_averaging_interval = 15days,
                               averages_interval = 5days,
                               checkpoint_interval = 720days,
                               file_splitting_interval = 360days)

    coupled_model = simulation.model
    ocean = coupled_model.ocean.model
    sea_ice = coupled_model.sea_ice.model
    grid = ocean.grid
    Nz = size(grid, 3)

    T, S = ocean.tracers.T, ocean.tracers.S
    u, v, w = ocean.velocities
    η = ocean.free_surface.displacement
    hⁱ = sea_ice.ice_thickness
    hˢ = sea_ice.snow_thickness
    ℵ = sea_ice.ice_concentration
    net_fluxes = coupled_model.interfaces.net_fluxes.ocean
    atmosphere_ocean_fluxes = coupled_model.interfaces.atmosphere_ocean_interface.fluxes

    tos = view(T, :, :, Nz)
    sos = view(S, :, :, Nz)

    surface_outputs = Dict{Symbol, Any}(:tos => tos,
                                        :sos => sos,
                                        :zos => η,
                                        :uos => view(u, :, :, Nz),
                                        :vos => view(v, :, :, Nz),
                                        :tossq => tos * tos,
                                        :sossq => sos * sos,
                                        :zossq => Field(η * η),
                                        :mlotst => MixedLayerDepthField(ocean.buoyancy, grid, ocean.tracers),
                                        :tauuo => net_fluxes.u,
                                        :tauvo => net_fluxes.v,
                                        :hfds => net_fluxes.T,
                                        :wfo => net_fluxes.S,
                                        :hfss => atmosphere_ocean_fluxes.sensible_heat,
                                        :hfls => atmosphere_ocean_fluxes.latent_heat,
                                        :JTf => frazil_temperature_flux(coupled_model),
                                        :JTn => net_ocean_temperature_flux(coupled_model),
                                        :JTio => sea_ice_ocean_temperature_flux(coupled_model),
                                        :JTao => atmosphere_ocean_temperature_flux(coupled_model),
                                        :JSn => net_ocean_salinity_flux(coupled_model),
                                        :JSio => sea_ice_ocean_salinity_flux(coupled_model),
                                        :siconc => ℵ,
                                        :sithick => hⁱ,
                                        :siu => sea_ice.velocities.u,
                                        :siv => sea_ice.velocities.v,
                                        :sitemptop => sea_ice.snow_thermodynamics.top_surface_temperature,
                                        :sisnthick => hˢ)

    field_outputs = Dict{Symbol, Any}(:to => T,
                                      :so => S,
                                      :uo => u,
                                      :vo => v,
                                      :wo => w,
                                      :bo => buoyancy_operation(ocean),
                                      :tke => ocean.tracers.e,
                                      :uosq => KernelFunctionOperation{Center, Center, Center}(u²ᶜᶜᶜ, grid, u),
                                      :vosq => KernelFunctionOperation{Center, Center, Center}(v²ᶜᶜᶜ, grid, v),
                                      :uvol => KernelFunctionOperation{Face, Center, Center}(zonal_volume_flux, grid, u),
                                      :vvol => KernelFunctionOperation{Center, Face, Center}(meridional_volume_flux, grid, v),
                                      :uvolto => KernelFunctionOperation{Face, Center, Center}(zonal_tracer_transport, grid, u, T),
                                      :vvolto => KernelFunctionOperation{Center, Face, Center}(meridional_tracer_transport, grid, v, T),
                                      :uvolso => KernelFunctionOperation{Face, Center, Center}(zonal_tracer_transport, grid, u, S),
                                      :vvolso => KernelFunctionOperation{Center, Face, Center}(meridional_tracer_transport, grid, v, S))

    vvolgm = bolus_meridional_volume_flux(ocean, isopycnal_closure(ocean))
    isnothing(vvolgm) || (field_outputs[:vvolgm] = vvolgm)

    average_outputs = Dict{Symbol, Any}(:zosga => Average(η),
                                        :to_h => Average(T, dims = (1, 2)),
                                        :so_h => Average(S, dims = (1, 2)),
                                        :bo_h => Average(buoyancy_operation(ocean), dims = (1, 2)),
                                        :voco => Integral(KernelFunctionOperation{Center, Center, Center}(unit_volume, grid)),
                                        :soco => Integral(S),
                                        :hoco => Integral(T),
                                        :sivol => Integral(hⁱ * ℵ),
                                        :snvol => Integral(hˢ * ℵ))

    # The writers serialize against `simulation.model`, a coupled model that carries none of the ocean properties
    writer_kw = (; including = Symbol[], dir = output_directory, overwrite_files = true)

    compressed = (; jld2_kw = Dict(:compress => ZstdFilter()))

    simulation.output_writers[:surface] = JLD2Writer(ocean, surface_outputs; writer_kw..., compressed...,
                                                     schedule = AveragedTimeInterval(surface_averaging_interval),
                                                     file_splitting = TimeInterval(file_splitting_interval),
                                                     filename = prefix * "_surface")

    simulation.output_writers[:fields] = JLD2Writer(ocean, field_outputs; writer_kw..., compressed...,
                                                    schedule = AveragedTimeInterval(field_averaging_interval),
                                                    file_splitting = TimeInterval(file_splitting_interval),
                                                    filename = prefix * "_fields")

    simulation.output_writers[:averages] = JLD2Writer(ocean, average_outputs; writer_kw...,
                                                      schedule = AveragedTimeInterval(averages_interval),
                                                      file_splitting = TimeInterval(file_splitting_interval),
                                                      filename = prefix * "_averages")

    simulation.output_writers[:checkpointer] = Checkpointer(coupled_model;
                                                            schedule = TimeInterval(checkpoint_interval),
                                                            prefix = joinpath(output_directory, prefix * "_checkpoint"),
                                                            cleanup = false)

    return simulation
end
