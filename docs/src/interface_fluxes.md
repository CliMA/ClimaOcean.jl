# Turbulent fluxes at component interfaces

To motivate this tutorial, we first note that `ClimaOcean`'s [`OceanSeaIceModel`](@ref) has essentially two goals:

1. Manage time-stepping multiple component models forward simultaneously,
2. Compute and communicate fluxes between the component models.

This tutorial therefore touches on the latter of the two main purposes of `OceanSeaIceModel`:
computing turbulent fluxes between model components.

## Component interfaces we consider

The `OceanSeaIceModel` has atmosphere, ocean, and sea ice components (and will someday also have land and radiation components).
We envision that the tutorial will eventually cover all turbulent flux computations; for the time being we
focus on atmosphere-ocean fluxes.
Future expansions of this tutorial should cover atmosphere-sea ice fluxes, ocean-sea ice fluxes, ocean-land fluxes,
and surface optical computations for radiation.

## Turbulent exchanges between the atmosphere and underlying media

Exchanges of properties like momentum, heat, water vapor, and trace gases between the fluid atmosphere and its underlying surfaces --
ocean, sea ice, snow, land -- mediate the evolution of the Earth system.
Microscopic property exchange is mediated by a complex panoply of processes including heat conduction, viscous and pressure form drag over rough surface elements, plunging breakers, and more.
To represent atmosphere-surface exchanges, we construct a model of the near-surface atmosphere that connects a turbulent "similarity layer",
which is usually a few meters thick, with a "constant flux layer" that buffers free atmospheric turbulence from microscopic surface exchange processes beneath.
The problem of modeling property exchange then turns to the task of modeling turbulent atmospheric fluxes just above the constant flux layer.

## Bulk formula and similarity theory

Within in each grid cell at horizontal position ``x, y, t``, the atmosphere-surface
turbulent fluxes of some quantity ``\psi`` -- at the bottom of the similarity layer, and thus throughout
the constant flux layer and across the surface -- is defined as

```math
J_\psi(x, y, t) = \overline{w' \psi'}
```

where ``w`` is the atmospheric vertical velocity, the overline ``\overline{( \; )}`` denotes a horizontal average over a grid cell,
and primes denote deviations from the horizontal average.

!!! note
    Arguably, the averaging operator ``\overline{( \; )}`` should also represent an average in time,
    which is implicit in the context of typical global Earth system modeling.
    Explicit time-averaging is required to evaluate flux observations, however,
    and may also be warranted for high-resolution coupled modeling.
    Flux computations in ClimaOcean currently compute fluxes in terms of the instantaneous states
    of its components, but spatial coarse-graining and time-averaging for computing fluxes at high
    resolution should be the subject of future research.

The essential turbulent fluxes that couple the ocean and atmosphere are

1. Momentum fluxes ``\rho_a \overline{\bm{u}'w'}``,
   where ``\rho_a`` is the atmosphere density at the air-sea interface and ``\bm{u}`` is horizontal velocity.

2. Sensible heat fluxes ``\rho_a c_{a} \overline{w'\theta'}`` due to fluid dynamical heat transport,
   where ``\rho_a`` is the atmosphere density at the air-sea interface,
   ``c_a`` is the atmosphere specific heat at constant pressure, and
   ``\theta`` is the atmosphere potential temperature.

3. Water vapor fluxes ``\overline{w' q'}`` due to evaporation and condensation,
   where ``q`` is the atmosphere specific humidity at the air-sea interface (the ratio between the mass of water and the total mass of an air parcel).

4. Latent heat fluxes ``\rho_a \mathscr{L}_v \overline{w' q'}`` due to the conversion of liquid ocean water into
   water vapor during evaporation, and vice versa during condensation, where
   ``\mathscr{L}_v`` is the latent heat of vaporization at the air-sea interface.

There are two ways by which turbulent fluxes may be computed: by specifying "transfer coefficients",
or by using Monin--Obukhov similarity theory.
In both cases, computing turbulent fluxes requires:

1. Atmosphere-surface differences in horizontal velocity, ``\Delta \bm{u}``,
2. Atmosphere-surface differences in temperature, ``\Delta \theta``,
3. The skin surface temperature ``T_s``, which is used to compute the surface specific humidity ``q_s`` and the
   atmosphere-surface specific humidity difference ``\Delta q``,
4. Additional atmosphere-surface trace gas differences for computing trace gas fluxes,
5. Possibly, additional "bulk" properties of the surface media and radiation fluxes
   in order to compute an equilibrium "skin" surface temperature that differs from the
   bulk temperature below the surface.
