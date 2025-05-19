```math
b' ≡ - g \left ( \rho - \bar{\rho} \right ) = g \frac{\alpha - \bar{\alpha}}{\bar{\alpha}}
\qquad \text{where} \qquad
α = \frac{1}{\rho} = \frac{R_m T}{p}
```

To define the buoyancy perturbation $b'$, we start by writing down the ideal gas law for a mixture of 
dry air (subscript $d$, itself a composite of gases) and water vapor (subscript $v$) that coexist in equilibrium at the same temperature $T$,

```math
p = \rho_d R_d T + \rho_v R_v T \, ,
```

where $R_d$ is the gas constant for dry air,
$R_v$ is the gas constant for water vapor, $\rho_d$ is the density of the dry air component,
$\rho_v$ is the density of the water vapor.
Note that atmospheric thermodynamics constants are stored in `atmosphere_properties` of `OceanSeaIceModel.interfaces`,


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
α - \bar{\alpha} = \frac{R_d}{p} \left [ T' - \mathscr{M} (q T - \overline{q T} ) \right ]
```

The characteristic buoyancy scale is then

```math
b_⋆ ≡ \frac{g}{T̃₀} \left [ \theta_⋆ \left ( 1 + δ q₀ \right ) + δ θ₀ q_⋆ \right ]
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
with $q > 0$ is lighter than perfectly dry air (because water vapor is lighter than dry air).
