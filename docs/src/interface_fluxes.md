# Turbulent fluxes at component interfaces

To motivate this tutorial, we first note that `ClimaOcean`'s [`OceanSeaIceModel`](@ref) has essentially two goals:

1. Manage time-stepping multiple component models forward simultaneously,
2. Compute and communicate fluxes between the component models.

This tutorial therefore touches on one of the two main purposes of `OceanSeaIceModel`:
computing turbulent fluxes between model components.

## Component interfaces we consider

The `OceanSeaIceModel` has atmosphere, ocean, and sea ice components (and will someday also have land and radiation components).
This tutorial will hopefully eventually cover all turbulent flux computations, but for the time being we
focus on atmosphere-ocean fluxes.
Future expansions of this tutorial should cover atmosphere-sea-ice fluxes, ocean-sea-ice fluxes, ocean-land fluxes,
and surface optical computations for radiation.

## Turbulent exchanges between the atmosphere and underlying media

Exchanges of properties like momentum, heat, water vapor, and trace gases between the fluid atmosphere and its underlying surfaces ---
ocean, sea ice, snow, land --- mediate the evolution of the Earth system.
Microscopic property exchange is mediated by a complex panolpy of processes including heat conduction, viscous and pressure form drag over rough surface elements, plunging breakers, and more.
To represent atmosphere-surface exchanges, we construct a model of the near-surface atmosphere that connects a turbulent "similarity layer", usually a few meters thick, with a "constant flux layer"
that buffers free atmospheric turbulence from microscopic surface exchange processes beneath.
The problem of modeling property exchange then turns to the task of modeling turbulent atmospheric fluxes just above the constant flux layer.

## Bulk formula and similarity theory

Within each grid cell, we represent atmosphere-surface turbulent fluxes of some quantity $\psi$ as

```math
J^ψ = \overline{w' \psi'}
```

where $w$ is the atmospheric vertical velocity, the overline $\overline{( \, )}$ denotes a horizontal average over a grid cell,
and primes denote deviations from the horizontal average.

!!! note
    Arguably, the averaging operator $\overline{( \, )}$ should also represent an average in time,
    which is implicit in the context of typical global Earth system modeling.
    Including time-averaging in the averaging operator is explicit required for analyzing obvservations,
    however, and may also be needed for very high resolution coupled modeling, and should be the subject
    of future research.

The essential turbulent fluxes that couple the ocean and atmosphere are

1. Momentum fluxes $\rho_a \overline{\bm{u}'w'}$,
   where $\rho_a$ is the atmosphere density at the air-sea interface, $\bm{u}$ is horizontal velocity, and $w$ is vertical velocity.

1. Sensible heat fluxes $\rho_a c_{a} \overline{w'\theta'}$ due to fluid dynamical heat transport, 
   where $\rho_a$ is the atmosphere density at the air-sea interface,
   $c_a$ is the atmosphere specific heat at constant pressure,
   $\theta$ is the atmosphere potential temperature, and the atmosphere specific heat at constant pressure.

1. Water vapor fluxes $\overline{w' q'}$ due to evaporation and condensation,
   where $q$ is the atmosphere specific humidity at the air-sea interface (the ratio between the mass of water and the total mass of an air parcel).

1. Latent heat fluxes $\rho_a \mathscr{L}_v \overline{w' q'}$ due to the conversion of liquid ocean water into
   water vapor during evaporation, and vice versa during condensation, where
   $\mathscr{L}_v$ is the latent heat of vaporization at the air-sea interface.

There are two ways by which turbulent fluxes may be computed: by specifying "transfer coefficients",
or by using Monin-Obukhov similarity theory.

### Coefficient-based fluxes

The simplest method for computing fluxes prescribes "transfer coefficients" that relate differences
between the near-surface atmosphere and the ocean surface to fluxes via scaling arguments,

```math
\overline{\bm{u}' w'} ≈ C_D      \, Δ \bm{u} \, U \\
\overline{w' \theta'} ≈ C_\theta \, Δ \theta \, U \\
\overline{w' q'}      ≈ C_q      \, Δ q \, U
```
where $Δ \bm{u}$ is the velocity difference between the atmosphere and the ocean surface,
$Δ \theta$ is the difference between the atmosphere and ocean temperature,
and $Δ q$ is the
and $q$ is the difference between the atmosphere specific humidity and the
saturation specific humidity computed using the ocean sea surface temperature and
accounting for ocean sea surface salinity via [Raoult's law](https://en.wikipedia.org/wiki/Raoult%27s_law).

The variable $U$ is a characteristic velocity scale, which is most simply just $U = | Δ \bm{u}|$.
However, an important class of parameterizations introduce an additional model for $U$ that
produces non-vanishing heat and moisture fluxes in zero-mean-wind conditions.
Usually these parameterizations are formulated as models for "gustiness" associated with atmospheric convection;
but more generally a common-thread is that $U$ may include contributions from turbulent motions
in addition to the relative mean velocity, $Δ \bm{u}$.

The variable $C_D$ is often called the drag coefficient, while $C_\theta$ and $C_q$ are the heat transfer
coefficient and vapor flux coefficient.
The simplest method for computing fluxes is merely to prescribe $C_D$, $C_\theta$, and $C_q$
as constants --- typically with a magnitude around $5 × 10^{-4}$--$2 × 10^{-3}$.
A comprehensive example will be given below, but we note briefly here that
`ClimaOcean` supports the computation of turbulent fluxes with constant coefficients via

```@example interface_fluxes
using ClimaOcean

coefficient_fluxes = CoefficientBasedFluxes(drag_coefficient=2e-3,
                                            heat_transfer_coefficient=2e-3,
                                            vapor_flux_coefficient=1e-3)
```

### Similarity theory for neutral boundary layers

Another way to compute fluxes is to use Monin-Obukhov similarity theory.
The formulation of similarity theory, which is based on dimensional arguments,
begins with the definition of "characteristic scales" which are related to momentum,
heat, and vapor fluxes through

```math
u_\star^2 ≡ | \overline{\bm{u}' w'} |  \\
u_\star \theta_\star ≡ \overline{w' \theta'} \\
u_\star q_\star ≡ \overline{w' q'}
```

where $u_\star$, often called the "friction velocity", is the characteristic scale for velocity,
$\theta_\star$ is the characteristic scale for temperature, and $q_\star$ is the characteristic scale
for water vapor.

To introduce similarity theory, we first consider momentum fluxes in "neutral" conditions,
or with zero buoyancy flux.
We further simplify the situation by considering unidirectional flow with $\bm{u} = u \, \bm{\hat x}$.
(To generalize our results to any flow direction, we simpliy rotate fluxes into the direction of the
relative velocity $Δ \bm{u}$.)
The fundamental supposition of similarity theory is that the vertical shear depends only on
height above the boundary, such that by dimensional analysis,

```math
\partial_z u \sim \frac{u_\star}{z} \, ,
\qquad \text{and thus} \qquad
\partial_z u = \frac{u_\star}{ϰ z} \, ,
```

where the second expression forms an equality by introducing the "Von Karman constant" $\varkappa$,
which is placed in the denominator by convention.
We can then integrate this expression from an inner scale $z=\ell$ up to $z=h$ to obtain

```math
u_a(h) - u_a(\ell_u) = \frac{u_⋆}{\varkappa} \log \left ( \frac{h}{\ell_u} \right )  
```

The inner length scale $\ell_u$, which is called the "momentum roughness length",
can be interpted as producing a characteristic upper value for the boundary layer shear, $u_⋆ / \ell_u$
in the region where similarity theory must be matched with the inner boundary layer (such as a viscous sublayer)
below.
Note that we take the inner velocity scale $u_a(\ell_u)$ as being equal to the velocity of the interface,
so $u_a(\ell_u) = u_o(0)$.

!!! note
    We currently assume that the input to the surface flux computation is the 
    atmospheric velocity at $z=h$. However, in coupled modeling context we are typically
    instead given the atmosphere velocity _averaged_ over the height of the first layer,
    or $⟨u_a⟩_h = \frac{1}{h} \int_0^h \, u_a \, \mathrm{d} z$.
    Formulating the flux computation in terms of $⟨u_a⟩_h$ rather than $u_a(h)$ 
    (e.g. [Nishizawa and Kitamura 2018](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018MS001534))
    is a minor modification to the algorithm and an important avenue for future work.

The roughness length in general depends on the physical nature of the surface.
For smooth, no-slip walls, experiments (cite) found agreement with a viscous sublayer model

```math
\ell_ν = \mathbb{C}_\nu \frac{\nu}{u_\star} \, ,
```

where $\nu$ is the kinematic viscosity of the fluid (the air in our case) and $\mathbb{C}_\nu$ is a free
parameter which was found to be $≈ 0.11$.
For air-water interfaces that develop a wind-forced spectrum of surface gravity waves, the alternative scaling

```math
\ell_g = \mathbb{C}_g \frac{g}{u_\star^2} \, ,
```

where $g$ is gravitational acceleration, has been repeatedly (and perhaps shockingly due to its simplicity) confirmed by field campaigns.
The free parameter $\mathbb{C}_g$ is often called the "Charnock parameter" and takes typical values
between $0$ and $0.03$ (cite Edson et al 2013).

```@example
using ClimaOcean
using CairoMakie

charnock_length = MomentumRoughnessLength(gravity_wave_parameter = 0.02,
                                          smooth_wall_parameter = 0,
                                          maximum_roughness_length = Inf)

smooth_wall_length = MomentumRoughnessLength(gravity_wave_parameter = 0,
                                             smooth_wall_parameter = 0.11)

default_roughness_length = MomentumRoughnessLength()
modified_default_length = MomentumRoughnessLength(gravity_wave_parameter = 0.011)

u★ = 1e-2:1e-2:3e-1
ℓg = charnock_length.(u★)
ℓν = smooth_wall_length.(u★)
ℓd = default_roughness_length.(u★)
ℓ2 = modified_default_length.(u★)

fig = Figure(size=(800, 400))
ax1 = Axis(fig[1, 1], xlabel="Friction velocity, u★ (m s⁻¹)", ylabel="Momentum roughness length ℓᵤ (m)")
lines!(ax1, u★, ℓd, label="ClimaOcean default")
lines!(ax1, u★, ℓg, label="Charnock")
lines!(ax1, u★, ℓν, label="Smooth wall")
lines!(ax1, u★, ℓ2, color=:black, label="ClimaOcean default with ℂg = 0.011")

ax2 = Axis(fig[1, 2], xlabel="Friction velocity, u★ (m s⁻¹)", ylabel="Momentum roughness length, ℓᵤ (m)")
u★ = 0.1:0.1:10
ℓd = default_roughness_length.(u★)
ℓ2 = modified_default_length.(u★)
lines!(ax2, u★, ℓd)
lines!(ax2, u★, ℓ2, color=:black)

Legend(fig[0, 1:2], ax1, orientation=:horizontal)

fig
```

!!! note
    The roughness length $\ell$ should not be interpreted as a physical length scale,
    a fact made clear by its submillimeter magnitude under (albeit calm) air-sea flux conditions.

## Computing fluxes and the effective similarity drag coefficient

ClimaOcean's default roughness length for air-sea fluxes is a function of the
friction velocity $u_\star$.
This formulation produces a nonlinear equation for $u_\star$, in terms of $Δ u = u_a(h) - u_o$,
which we emphasize by rearranging the similarity profile 

```math
u_\star = \frac{ϰ \, Δ u}{\log \left [ h / \ell_u(u_\star) \right ]} \, . 
```

This equation is solved for $u_\star$ using fixed-point iteration initialized with a reasonable
guess for $u_\star$.
Once $u_\star$ is obtained, the similarity drag coeffient may then be computed via

```math
C_D(h) ≡ \frac{u_\star^2}{|Δ u(h)|^2} = \frac{\varkappa^2}{\left ( \log \left [ h / \ell_u \right ] \right )^2} \,
```

where we have used the simple bulk velocity scale $U = Δ u$.
We have also indicated that, the effective similarity drag "coefficient" depends on the height $z=h$
at which the atmosphere velocity is computed to form the relative velocity $Δ u = u_a(h) - u_o$.
Most observational campaigns use $h = 10 \, \mathrm{m}$ and most drag coefficients reported in the
literature pertain to $h=10 \, \mathrm{m}$.

To compute fluxes with ClimaOcean, we build an `OceanSeaIceModel` with an atmosphere an ocean state
concocted such that we can evaluate fluxes over a range of relative atmosphere and oceanic conditions.
For this we use a $200 × 200$ horizontal grid and start with atmospheric winds that vary from
the relatively calm $u_a(10 \, \mathrm{m}) = 0.5 \, \mathrm{m \, s^{-1}}$ to a
blustery $u_a(10 \, \mathrm{m}) = 40 \, \mathrm{m \, s^{-1}}$.
We also initialize the ocean at rest with surface temperature $T_o = 20 \, \mathrm{{}^∘ C}$ and 
surface salinity $S_o = 35 $g \, kg^{-1}$ --- but the surface temperature and salinity won't matter until later.

```@example interface_fluxes
using Oceananigans
using ClimaOcean

# Atmosphere velocities
Nx = Ny = 200
uₐ = range(0.5, stop=40, length=Nx) # winds at 10 m, m/s

# Ocean state parameters
T₀ = 20   # Surface temperature, ᵒC
S₀ = 35   # Surface salinity

x = y = (0, 1)
z = (-1, 0)
atmos_grid = RectilinearGrid(size=(Nx, Ny); x, y, topology=(Periodic, Periodic, Flat))
ocean_grid = RectilinearGrid(size=(Nx, Ny, 1); x, y, z, topology=(Periodic, Periodic, Bounded))

# Build the atmosphere
atmosphere = PrescribedAtmosphere(atmos_grid, surface_layer_height=10)
interior(atmosphere.tracers.T) .= 273.15 + T₀ # K
interior(atmosphere.velocities.u, :, :, 1, 1) .= uₐ # m/s

kw = (momentum_advection=nothing, tracer_advection=nothing, closure=nothing)
ocean = ocean_simulation(ocean_grid; kw...)
set!(ocean.model, T=T₀, S=S₀)
```

Next we build two models with different flux formulations --- the default  "similarity model"
that uses similarity theory with "Charnock" gravity wave parameter $\mathbb{C}_g = 0.02$,
and a "coefficient model" with a constant drag coefficient $C_D = 2 × 10^{-3}$:

```@example interface_fluxes
neutral_similarity_fluxes = SimilarityTheoryFluxes(stability_functions=nothing)
interfaces = ComponentInterfaces(atmosphere, ocean; atmosphere_ocean_fluxes=neutral_similarity_fluxes)
default_model = OceanSeaIceModel(ocean; atmosphere, interfaces)

momentum_roughness_length = MomentumRoughnessLength(gravity_wave_parameter=0.04)
neutral_similarity_fluxes = SimilarityTheoryFluxes(stability_functions=nothing; momentum_roughness_length)
interfaces = ComponentInterfaces(atmosphere, ocean; atmosphere_ocean_fluxes=neutral_similarity_fluxes)
increased_roughness_model = OceanSeaIceModel(ocean; atmosphere, interfaces)

coefficient_fluxes = CoefficientBasedFluxes(drag_coefficient=2e-3)
interfaces = ComponentInterfaces(atmosphere, ocean; atmosphere_ocean_fluxes=coefficient_fluxes)
coefficient_model = OceanSeaIceModel(ocean; atmosphere, interfaces)
```

Note that `OceanSeaIceModel` computes fluxes upon instantiation, so after constructing
the two models we are ready to analyze the results.
We first verify that the similarity model friction velocity has been computed successfully,

```@example interface_fluxes
u★ = default_model.interfaces.atmosphere_ocean_interface.fluxes.friction_velocity
u★ = interior(u★, :, 1, 1)
extrema(u★)
```

and it seems that we've obtained a range of friction velocities, which is expected
given that our atmospheric winds varied from $0.5$ to $40 \, \mathrm{m \, s^{-1}}$.
Computing the drag coefficient for the similarity model is as easy as

```@example interface_fluxes
Cᴰ_default = @. (u★ / uₐ)^2
extrema(Cᴰ_default)
```

We'll also re-compute the drag coefficient for the coefficent model
(which we specified as constant), which verifies that the coefficient was correctly 
specified:

```@example interface_fluxes
u★_coeff = coefficient_model.interfaces.atmosphere_ocean_interface.fluxes.friction_velocity
u★_coeff = interior(u★_coeff, :, 1, 1)
Cᴰ_coeff = @. (u★_coeff / uₐ)^2
extrema(Cᴰ_coeff)
```

We'll compare the computed fluxes and drag coefficients from our two models with
a polynomial expression due to Large and Yeager (2009), and
an expression reported in Edson et al. (2013) that was developed at ECMWF,

```@example interface_fluxes
# From Large and Yeager (2009), equation XX
c₁ = 0.0027
c₂ = 0.000142
c₃ = 0.0000764
u★_LY = @. sqrt(c₁ * uₐ + c₂ * uₐ^2 + c₃ * uₐ^3)
Cᴰ_LY = @. (u★_LY / uₐ)^2

# From Edson et al. (2013), equation 20
c₁ = 1.03e-3
c₂ = 4e-5
p₁ = 1.48
p₂ = 0.21
Cᴰ_EC = @. (c₁ + c₂ * uₐ^p₁) / uₐ^p₂
u★_EC = @. sqrt(Cᴰ_EC) * uₐ
extrema(u★_EC)
```

Finally, we plot the results to compare the estimated friction velocity and effective
drag coefficient from the polynomials expressions with the two `OceanSeaIceModel`s:

```@example interface_fluxes
using CairoMakie
set_theme!(Theme(fontsize=14, linewidth=2))

# Extract u★ and compute Cᴰ for increased roughness model
u★_rough = increased_roughness_model.interfaces.atmosphere_ocean_interface.fluxes.friction_velocity
u★_rough = interior(u★_rough, :, 1, 1)
Cᴰ_rough = @. (u★_rough / uₐ)^2

fig = Figure(size=(800, 400))
axu = Axis(fig[1:2, 1], xlabel="uₐ (m s⁻¹) at 10 m", ylabel="u★ (m s⁻¹)")
lines!(axu, uₐ, u★, label="ClimaOcean default")
lines!(axu, uₐ, u★_rough, label="Increased roughness model")
lines!(axu, uₐ, u★_LY, label="Large and Yeager (2009) polynomial fit")
lines!(axu, uₐ, u★_EC, label="ECMWF polynomial fit (Edson et al 2013)")

axd = Axis(fig[1:2, 2], xlabel="uₐ (m s⁻¹) at 10 m", ylabel="1000 × Cᴰ")
lines!(axd, uₐ, 1000 .* Cᴰ_default, label="ClimaOcean default")
lines!(axd, uₐ, 1000 .* Cᴰ_rough, label="ClimaOcean default")
lines!(axd, uₐ, 1000 .* Cᴰ_LY, label="Large and Yeager (2009) polynomial fit")
lines!(axd, uₐ, 1000 .* Cᴰ_EC, label="ECMWF polynomial fit (Edson et al 2013)")

Legend(fig[3, 1:2], axd, nbanks = 2)

fig
```

## Non-neutral boundary layers and stability functions

The relationship between the relative air-sea state and turbulent fluxes
is modified by the presence of buoyancy fluxes --- "destabilizing" fluxes, which stimulate convection,
tend to increase turbulent exchange, while stabilizing fluxes suppress turbulence and turbulent exchange.
Monin-Obhukhov stability theory provides a scaling-argument-based framework
for modeling the effect of buoyancy fluxes on turbulent exchange.

### Buoyancy flux and stability of the near-surface atmosphere

Our next objective is to characterize the atmospheric statbility in terms of the buoyancy flux, $\overline{w' b'}$,
which requires a bit of thermodynamics background to define the buoyancy perturbation, $b'$.

#### The ideal gas law

To define the buoyancy perturbation $b'$, we start by writing down the ideal gas law for a mixture of 
dry air (subscript $d$, itself a composite of gases) and water vapor (subscript $v$) that coexist in equilibrium at the same temperature $T$,

```math
p = \rho_d R_d T + \rho_v R_v T \, ,
```

where $R_d$ is the gas constant for dry air,
$R_v$ is the gas constant for water vapor, $\rho_d$ is the density of the dry air component,
$\rho_v$ is the density of the water vapor.
Note that atmospheric thermodynamics constants are stored in `atmosphere_properties` of `OceanSeaIceModel.interfaces`,

```@example interface_fluxes
default_model.interfaces.atmosphere_properties
```

The total specific humidity, $q$, is defined as the mass ratio of total water,
including the vapor density $\rho_v$ and density of liquid and ice condensate $\rho_c$, 
to the total mixture density $\rho$,

```math
q = \frac{\rho_v + \rho_c}{\rho} ≈ \frac{\rho_v}{\rho}
```

The second expression above introduces the approximation that $q_c ≡ \rho_c / ρ ≪ q$,
which we make hereafter for the purpose of evaluating stability during surface flux computations.

!!! note
    It'd be nice to comment further on the validity of this approximation, and its practicality when working with atmospheric data products.)


Putting the pieces together, we write the ideal gas law as

```math
p = ρ R_m T
\qquad \text{where} \qquad
R_m(q) ≈ R_d \left (1 - q \right ) + R_v q = R_d \left ( 1 - δ q \right ) \, 
```

with $δ ≡ \frac{R_v}{R_d} - 1$. For water vapor and a standard dry air composite, $\delta ≈ 0.61$.

#### The buoyancy perturbation

Our next objective is to compute the buoyancy perturbation $b'$ in terms of
surface temperature $T$, specific humidity $q$, and thermodynamic constants.
The subgrid buoyancy perturbation $b'$ is defined 

```math
b' ≡ - g \left ( \rho - \bar{\rho} \right ) = g \frac{\alpha - \bar{\alpha}}{\bar{\alpha}}
\qquad \text{where} \qquad
α = \frac{1}{\rho} = \frac{R_m T}{p}
```

We then neglect pressure perturbations so that $p = \bar{p}$, which simplifies the specific volume perturbation to,

```math
α - \bar{\alpha} = \frac{1}{p} \left [ R_d (1 - δ q ) T' - δ q' T \right ]
```


in terms of the dry air gas constant $R_d$ (assuming dry air is a composite of nitrogen, oxygen, argon, and other trace gases),
and $R_v$, the gas constant for water vapor, and $q_c$, the mass fraction of liquid and ice water condensates.
This formula assumes also that the condenstate volume is negligible, so that condensates do not contribute to the mixture total pressure.
Conditions in Earth's atmosphere are such that usually $q_c ≪ q$, which motivates the approximation

```math
R_m(q) ≈ R_d \left (1 - q \right ) + R_v q = R_d \left ( 1 - δ q \right ) \, 
\qquad \text{where} \qquad δ ≡ \frac{R_v}{R_d} - 1 ≈ 0.61 \, .
```

We've thus shown why the number "0.61" (or 1.61) appears frequently in texts on turbulent fluxes.
Another reason to neglect $q_c$ is that atmospheric data products often do not
provide $q_c$ separately from $q$.

If we neglect pressure perturbations so that $p = \bar{p}$ and use $\bar{\alpha} = \bar{p} / {R_m(\bar{q}) \bar{T}}$, then
the buoyancy perturbation may be computed via

```math
b' ≡ - g \frac{\tilde{T}' - \overline{\tilde{T}}}{\overline{\tilde{T}}}
```

The characteristic buoyancy scale is then

```math
b_⋆ ≡ \frac{g}{T̃₀} \left [ \theta_⋆ \left ( 1 + δ q₀ \right ) + δ θ₀ q_⋆ \right ]
```

where

```math
T̃ = \frac{R_m(q)}{R_d} T ≈ T + δ q T \, ,
```

Virtual temperature plays the role of temperature in the computation of buoyancy
for moist air. More specifically, virtual temperature accounts for the fact that moist air
with $q > 0$ is heavier than perfectly dry air (because water vapor is heavier than dry air).

We begin by defining the characteristic buoyancy scale $b_\star$ in termms of the air-sea buoyancy flux,

```math
u_⋆ b_⋆ = \overline{w' b'}
```



To measure the stability of the atmosphere at $z=h$, we consider the ratio between shear production
and buoyancy production, which oceanographers usually call the "flux Richardson number",

```math
Ri_f(z=h) ≡ - \frac{\overline{w' b'}}{\partial_z \bar{\bm{u}} \, ⋅ \, \overline{\bm{u}' w'}} = \frac{ϰ \, h \, b_⋆}{u_\star^2} = \frac{h}{L}
\qquad \text{where} \qquad
L ≡ - \frac{u_\star^2}{ϰ b_\star} \, ,
```

is called the "Monin-Obhukhov length scale".
In the air-sea flux literature, the symbol $ζ = h / L$ is used in place of $Ri_f$, and referred to as the "stability parameter".

### The Monin-Obhukhov "stability functions"

The fundamental premise of Monin-Obhkhov stability theory is that the relationship
between air-surface differences and turbulent fluxes can be expressed by

```math
\frac{ϰ z}{u_\star} \partial_z \bar{u} = \phi_u(\zeta) \, ,
```

where $\phi_u(\zeta)$ is called the "stability function", and we recall that $\zeta$ is the flux Richardson number.
$\phi(\zeta)$ must be $> 1$ when turbulent fluxes are destabilizing and $ζ > 0$, representing enhanced fluxes,
and $< 1$ in the presence of stabilizing buoyancy fluxes with $ζ > 0$, representing a suppression of fluxes.
For convenience we introduce the transformation

```math
\phi_u(\zeta) = 1 - ζ \psi_u'(\zeta) \, ,
```

where $\psi'$ denote the derivative of $\psi$.
Then integrating from $z=\ell_u$ to $z=h$ as for the neutral case, we find that

```math
u_a(h) - u_a(\ell_u) = Δ u = \frac{u_\star}{\varkappa} \left [ \log \left (\frac{h}{\ell_u} \right ) + ψ_u \left ( \frac{h}{L} \right ) - ψ_u \left (\frac{\ell_u}{L} \right ) \right ]
```

Note that the last term, $\psi_u(\ell_u / L)$, is often neglected because $\ell_u / L$ is miniscule and $\psi_u(0) = 0$.
Similar formula hold for temperature and water vapor,

```math
Δ \theta = \frac{\theta_\star}{\varkappa} \left [ \log \left (\frac{h}{\ell_\theta} \right ) + ψ_\theta \left ( \frac{h}{L} \right ) - ψ_\theta \left (\frac{\ell_\theta}{L} \right ) \right ] \\
Δ q = \frac{q_\star}{\varkappa} \left [ \log \left (\frac{h}{\ell_q} \right ) + ψ_q \left ( \frac{h}{L} \right ) - ψ_q \left (\frac{\ell_q}{L} \right ) \right ]
```

We rearrange these into expresssions for $u_\star$, $\theta_\star$, and $q_\star$ in terms of $Δ u$, $Δ \theta$, and $Δ q$, and solve
the resulting system of three equations using fixed-point iteration.
At each iterate, both the roughness lengths and the Monin-Obhukhov length $L = - u_\star^2 / ϰ b_\star$ are estimated using the values of $u_\star$, $\theta_\star$, and $q_\star$
from the previous iterate.

```@example interface_fluxes
using ClimaOcean.OceanSeaIceModels.InterfaceComputations: saturation_specific_humidity
ρₐ = 1.2 # guess
Tₒ = 273.15 + 20 # in Kelvin
Sₒ = 35
interfaces = default_model.interfaces
ℂₐ = interfaces.atmosphere_properties
q_formulation = interfaces.atmosphere_ocean_interface.properties.specific_humidity_formulation
qₛ = saturation_specific_humidity(q_formulation, ℂₐ, ρₐ, Tₒ, Sₒ)
```

We then set the atmos state:

```@example interface_fluxes
interior(atmosphere.pressure) .= 101352
interior(atmosphere.tracers.q) .= qₛ

Tₐ = 273.15 .+ range(-40, stop=40, length=Ny)
Tₐ = reshape(Tₐ, 1, Ny)
interior(atmosphere.tracers.T) .= Tₐ

Oceananigans.TimeSteppers.update_state!(default_model)
u★ = default_model.interfaces.atmosphere_ocean_interface.fluxes.friction_velocity
θ★ = default_model.interfaces.atmosphere_ocean_interface.fluxes.temperature_scale

fig = Figure(size=(800, 400))
axu = Axis(fig[1, 1], xlabel="Wind speed uₐ (m s⁻¹)", ylabel="Air-sea temperature difference (K)")
axθ = Axis(fig[1, 2], xlabel="Wind speed uₐ (m s⁻¹)", ylabel="Air-sea temperature difference (K)")
axC = Axis(fig[2, 1:2], xlabel="Wind speed uₐ (m s⁻¹)", ylabel="Cᴰ / neutral Cᴰ")

ΔT = Tₐ .- 293.15
ΔT = reshape(ΔT, Ny, 1)
ΔT = dropdims(ΔT, dims=2)

heatmap!(axu, uₐ, ΔT, interior(u★, :, :))
heatmap!(axθ, uₐ, ΔT, interior(θ★, :, :))

u★ = interior(u★, :, :)
uₐ = reshape(uₐ, Nx, 1)
Cᴰ = @. (u★ / uₐ)^2

uₐ = uₐ[:]
#lines!(ax★, uₐ, Cᴰ[:, 1:10:end])
lines!(axC, uₐ, Cᴰ[:, 1]   ./ Cᴰ_default)
lines!(axC, uₐ, Cᴰ[:, 20]  ./ Cᴰ_default)
lines!(axC, uₐ, Cᴰ[:, 50]  ./ Cᴰ_default)
lines!(axC, uₐ, Cᴰ[:, 100] ./ Cᴰ_default)
lines!(axC, uₐ, Cᴰ[:, 150] ./ Cᴰ_default)
lines!(axC, uₐ, Cᴰ[:, 200] ./ Cᴰ_default)

xlims!(axC, 0, 10)
ylims!(axC, 0, 4)

fig
```

The coefficient-based formula then take the form
 
```math
u_\star = \sqrt{C_D | Δ \bm{u} | \, U} \\
\theta_\star = \frac{C_θ}{\sqrt{C_D}} \, Δ θ \, \sqrt{\frac{U}{|Δ \bm{u} |}} \\ 
q_\star = \frac{C_q}{\sqrt{C_D}} \, Δ q \, \sqrt{\frac{U}{| Δ \bm{u} |}} \\ 
```
