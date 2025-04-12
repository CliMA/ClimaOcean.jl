using ClimaOcean
using Oceananigans
import Thermodynamics as AtmosphericThermodynamics

# Ocean state parameters
T₀ = 20   # Surface temperature, ᵒC
S₀ = 35   # Surface salinity
Tᵢ = T₀

# Atmospheric state parameters
Tₐ = 273.15 + T₀
u₁₀ = 10 # wind at 10 m, m/s
qₐ = 0.01 # specific humidity
Qℓ = 400 # shortwave radiation (W m⁻², positive means heating right now)
pₐ = 1e5

# Build the atmosphere
radiation = Radiation(ocean_albedo=0.0)
atmosphere_grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))
atmosphere_times = range(0, 86400.0, length=3)
atmosphere = PrescribedAtmosphere(atmosphere_grid, atmosphere_times, boundary_layer_height=600)

# Build ocean model at rest with initial temperature stratification
grid = RectilinearGrid(size=3, z=(-1, 0), topology=(Flat, Flat, Bounded))
ocean = ocean_simulation(grid, coriolis=FPlane(f=0))

parent(atmosphere.tracers.T) .= Tₐ     # K
parent(atmosphere.velocities.u) .= u₁₀ # m/s
parent(atmosphere.tracers.q) .= qₐ     # mass ratio
parent(atmosphere.downwelling_radiation.shortwave) .= Qℓ # W
parent(atmosphere.pressure) .= pₐ

model = OceanSeaIceModel(ocean; atmosphere, radiation)

Qv  = model.interfaces.atmosphere_ocean_interface.fluxes.latent_heat
Qc  = model.interfaces.atmosphere_ocean_interface.fluxes.sensible_heat
ρτx  = model.interfaces.atmosphere_ocean_interface.fluxes.x_momentum
ρτy  = model.interfaces.atmosphere_ocean_interface.fluxes.y_momentum
Fv  = model.interfaces.atmosphere_ocean_interface.fluxes.water_vapor

# Bu1ld thermodynamic and dynamic states in the atmosphere and interface.
# Notation:
#   ⋅ 𝒬 ≡ thermodynamic state vector
#   ⋅ 𝒰 ≡ "dynamic" state vector (thermodynamics + reference height + velocity)
ℂₐ = atmosphere.thermodynamics_parameters
𝒬ₐ = thermodynamic_atmospheric_state = AtmosphericThermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₐ, qₐ)
zₐ = 10

@show ℂₐ
@show 𝒬ₐ.ρ Tᵢ Sᵢ
Tᵢ += 273.15 # Convert to Kelvin
q_formulation = model.interfaces.atmosphere_ocean_interface.properties.specific_humidity_formulation
qₛ = ClimaOcean.OceanSeaIceModels.InterfaceComputations.saturation_specific_humidity(q_formulation, ℂₐ, 𝒬ₐ.ρ, Float64(Tᵢ), Float64(Sᵢ))

@show ρₐ = AtmosphericThermodynamics.air_density(ℂₐ, 𝒬ₐ)
@show qₛ
ΔU = u₁₀
τx = ρτx[1, 1, 1] / ρₐ
τy = ρτy[1, 1, 1] / ρₐ
u★ = sqrt(sqrt(τx^2 + τy^2))
Cd = (u★ / ΔU)^2
