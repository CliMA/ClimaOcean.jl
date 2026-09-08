#!/bin/bash
# Submit an OMIP simulation to SLURM.
#
# Usage:
#   ./launch.sh orca                           # ORCA with default fluxes
#   NCAR=true ./launch.sh orca                 # ORCA with NCAR bulk formulae
#   NCAR=true SNOW=true ./launch.sh orca       # ORCA + NCAR + snow
#   CB=0.1 NCAR=true ./launch.sh orca          # ORCA + NCAR + Cᵇ=0.1
#   KSKEW=1000 KSYMM=500 ./launch.sh orca      # ORCA with custom eddy diffusivities
#   PROFILE=true ./launch.sh orca              # nsys-profile run
#
# Credentials (e.g. ECCO_USERNAME, ECCO_WEBDAV_PASSWORD) are NOT set
# here. Export them in your shell or source a private file before
# launching, e.g.:
#
#   source ~/.ecco_credentials && ./launch.sh orca

set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: ./launch.sh <config> [extra sbatch args...]

Configurations:
  halfdegree      Half-degree TripolarGrid
  quarterdegree   1/4-degree TripolarGrid (2 GPUs)
  orca            ORCA grid
  twelfthdegree     1/10-degree TripolarGrid (4 GPUs)

Environment variables (physics):
  NCAR          Set to "true" for OMIP-2/NCAR bulk formulae
  CORRECTED     Set to "true" for corrected COARE 3.6 fluxes
  SNOW          Set to "true" to enable snow thermodynamics
  SNOW_CATEGORIES
                Sub-grid categories for the snow conductivity, independently of ICE_CATEGORIES.
                Snow and ice are both scaled by the same Fichefet factor by default, which assumes
                snow depth and ice thickness co-vary in lockstep in the sub-grid. Snow dominates the
                series resistance hs/ks + hi/ki because ks is ~6x smaller, so this is the larger of
                the two knobs at typical snow depths. Default: whatever ICE_CATEGORIES is. Adds
                "_snowcat<value>" to the run name.
  ICE_DYNAMICS  Set to "false" to disable sea-ice dynamics (thermo-only ice).
  RUN_NAME_OVERRIDE
                Use this exact run name instead of the one built from the options, to resume a run
                whose directory predates the 2026-09-03 rename. The physics is still whatever the
                options say, so pass the old defaults too (ICE_DRAGREF=none ICE_LIQUIDUS=linear).
  ICE_LATERAL   Sea-ice lateral boundary condition: "no_slip" (default) applies the viscous
                wall stress -2 eta u / Delta on coastlines; "free_slip" leaves them stress-free.
                The old quadratic "side drag" was inert (unit mismatch) and has been removed.
  ICE_BASAL     Set to "false" to disable the Lemieux et al. (2015) landfast basal stress on
                grounded keels. Default: "true".
  ICE_DRAG      Ice-ocean drag coefficient. Default: 5.5e-3 (Hibler/McPhee, also ClimaSeaIce
                own default). The previous value 3.24e-3 is GFDL SIS2 CDW.
  ICE_PSTAR     Ice compressive strength P* in N/m^2, the EVP rheology's strength parameter.
                Default: 27500 (CICE/ClimaSeaIce). Raising it stiffens the pack and slows the
                export; the literature range is roughly 20000-45000.
  IC_CONDITIONS Named preset for the initial state.
                  default    WOA Annual T and S everywhere and ECCO4Monthly sea ice for January
                             1993. Reproduces the previous model exactly.
                  summerice  the summer pack in BOTH hemispheres -- ECCO4Monthly September 1993
                             north of the equator, January 1993 south of it. A 1 January start
                             otherwise takes the Arctic at its seasonal maximum and the first
                             spring melts that excess into the Labrador, the delivery C15-18 found
                             to trigger the convection cap. Adds "_summerice".
                  blended    "summerice" plus the January ocean: WOA Monthly for January blended
                             into WOA Annual above 500 m, i.e. IC_BLEND=500, which IC_BLEND still
                             overrides. Adds "_icblend<val>_summerice".
                The two arms are separable so that "summerice" is a single-parameter twin of the
                control and "blended" of "summerice".
  IC_BLEND      Depth in metres above which the initial T and S come from WOA Monthly for the
                month of the start date, instead of WOA Annual, tapering to the annual field at
                the monthly climatology's reach (~1500 m). WOA Annual averages a seasonal cycle
                whose winter half has km-deep mixed layers, so a January start from it begins
                with a seasonal thermocline the season does not have. Unset uses WOA Annual
                throughout and reproduces the previous model exactly. Adds "_icblend<val>".
                Requires woa_t_monthly_<MM>.nc and woa_s_monthly_<MM>.nc in the climatology dir.
  ICE_DRAGREF   Depth in metres over which the ocean velocity is averaged to give the reference
                of the ice-ocean drag. McPhee's Cio = 5.5e-3 is defined against the under-ice
                boundary layer; the topmost cell is 1.5 m and is dragged by the ice itself, so
                referencing the drag there under-brakes the pack -- and, because the same relative
                velocity sets u* in the ice-ocean heat flux, under-melts it from below as well.
                *** DEFAULT CHANGED 2026-09-03: now 6, which takes Fram ice export from 2911 to
                2331 km3/yr against an observed 2300 (C11-24). Set ICE_DRAGREF=none to restore the
                topmost cell and reproduce runs launched before this date. *** Tagged
                "_dragref<val>" only when it differs from 6, so "_dragrefnone" marks the old
                behaviour.
  ICE_ITD_SHAPE Gamma-distribution sub-grid thickness closure for the conductivity, as
                "<minimum_shape>,<maximum_shape>,<transition_thickness>". The factor is s/(s-1)
                with s running from maximum_shape for thin level ice to minimum_shape for thick
                deformed pack over transition_thickness metres. Unlike the Fichefet sum this is
                finite and needs no truncation, so there is no arbitrary category count. Use with
                ICE_CATEGORIES=1 so the stored conductivity is the bare material value. Adds
                "_itd<min>-<max>-<h>" to the run name.
  ICE_Z0        Aerodynamic momentum roughness of the ice surface, in metres, used by the corrected
                flux configuration. Sea ice carries no gravity waves, so this is a geometric constant
                set by ridges, floe edges and sastrugi rather than a Charnock relation. The default
                5e-4 is the SHEBA multiyear-pack value (Andreas et al. 2010); smooth first-year ice
                sits nearer 1e-4, which cuts the neutral drag coefficient from 1.63e-3 to 1.21e-3 and
                the free-drift speed by 14%. Adds "_icez0<value>" to the run name.
  TRACER_ORDER  Order of the WENO tracer advection. Default: 7. Adds "_tracer<val>".
  BUFFER_ORDER  Order the WENO reconstruction bottoms out at when its stencil touches a boundary --
                the DOMAIN buffer, and on an ImmersedBoundaryGrid any stencil containing an inactive
                node (Advection/immersed_advective_fluxes.jl reduces the order there; it does NOT
                simply zero the flux). Default: 3, which is still non-monotone. Set BUFFER_ORDER=1 to
                let the recursion reach FIRST-ORDER UPWIND, which is monotone, in exactly those
                boundary-touching cells while the interior keeps the full order -- targeted at
                topography-induced extrema without paying accuracy everywhere. Adds "_buford<val>".
                Only 1 and 3 are accepted: with boundary_scheme replacing minimum_buffer_upwind_order
                in Oceananigans, BUFFER_ORDER=1 is now the old spelling of BOUNDARY_SCHEME=upwind and
                keeps the same run name, while a chain truncated above third order is no longer
                expressible.
  BOUNDARY_SCHEME
                The reconstruction the WENO buffer chain terminates in, i.e. what is used in the one
                boundary-touching cell described under BUFFER_ORDER.
                  default  Centered(order=2) for tracers, UpwindBiased(order=1) for momentum.
                  upwind   first-order upwind everywhere the stencil does not fit; = BUFFER_ORDER=1.
                  cwenoz   the third-order central-WENO reconstruction of Semplice, Travaglia and
                           Puppo (2022), whose stencil extends only inwards. It blends an inward
                           parabola, a linear polynomial and a constant with Z-weights, so it keeps
                           third-order accuracy on smooth data and falls to the constant only where
                           the data is genuinely rough -- unlike "upwind", which pays first order in
                           every boundary cell. Adds "_cwenoz".
                Sets both tracers and momentum. TRACER_BOUNDARY_SCHEME and MOMENTUM_BOUNDARY_SCHEME
                take the same three values and override it one component at a time, which is how to
                tell whether the boundary treatment acts through the tracers or through the momentum.
                A scheme shared by both keeps the plain tag; split, they add "_tr<scheme>" and
                "_mom<scheme>", so the two halves never collide in one run directory.
  TRACER_BOUNDARY_SCHEME, MOMENTUM_BOUNDARY_SCHEME
                BOUNDARY_SCHEME for the tracer and for the momentum reconstructions separately.
                Default: whatever BOUNDARY_SCHEME is.
  TEMPERATURE_VARIATION, SALINITY_VARIATION, MOMENTUM_VARIATION
                CWENOZ reference variations, used only where the boundary scheme is cwenoz. Each sets
                eps = variation^2: the cell-to-cell variation of the field below which the reconstruction
                reads the stencil as noise and keeps third order. The constant candidate takes over above
                a variation of about 8.6 * variation between adjacent averages.

                It carries the units of the field and no grid spacing, so one value serves every direction
                and every cell thickness. All default to 0, which reads eps off the stencil: a nonzero
                value acts only on the horizontal and is inert on the vertical. Each adds its own tag to
                the run name.

                MOMENTUM_VARIATION is a speed in m/s and applies to the vertical momentum reconstruction
                alone: the horizontal terms reconstruct a vorticity, a divergence flux and a squared
                velocity, so one constant cannot carry their units.
  ICE_LIQUIDUS  Freezing-point relation. "teos10" (default) is the linear fit to the TEOS-10
                freezing point expressed in CONSERVATIVE temperature, which is what the ocean
                carries: Tm = -0.054523 S, accurate to 0.013 K over S = 28-35.5, against 0.032 K
                for ClimaSeaIce's own 0.054 -- which is too warm at EVERY salinity and so biases the
                ice-ocean heat flux one way everywhere. The intercept stays 0 because fresh water
                freezes at 0 C and masked (bathymetry) cells carry S = 0. The pressure dependence,
                -7.53e-4 K/dbar and worth 0.5 K at 660 m, is applied separately where depth is known.
                "linear" restores ClimaSeaIce's own Tm = -0.054 S, which is up to 0.032 K TOO WARM
                and biases the ice-ocean heat flux rho c alpha_h u* (Theta - Tm) directly.
                *** DEFAULT CHANGED 2026-09-03. *** Adds "_liq<val>" when not "teos10".
  ICE_TILT      Set to "true" to add the ocean surface tilt term -g grad(eta) to the sea ice
                momentum equation. The term is f x u_geostrophic, so without it the ice cannot ride
                the ocean's dynamic topography and the uncompensated Coriolis force is absorbed by
                the ice-ocean drag. Default: false, which reproduces the previous model exactly.
                Adds "_icetilt" to the run name. Requires a ClimaSeaIce carrying the free-surface
                term.
  ICE_SALINITY  Bulk sea-ice salinity in psu (ClimaSeaIce ConstantField). Sets the salt returned
                per unit melt, so the freshwater a melting cell delivers goes as (S_ocean - S_ice)
                / S_ocean: 0.885 at 4 psu against 0.828 at 6. Multi-year Arctic ice is 2-4 psu,
                first-year ice 5-8. Default: 4. Adds "_sice<val>" to the run name.
  ICE_HEAT_TRANSFER
                Ice-ocean turbulent heat transfer coefficient alpha_h. Default: 0.0057 (McPhee
                Stanton number, the value consistent with a COMPUTED friction velocity). The
                previous 0.0095 is calibrated by Shi et al. (2021) against a FIXED u* = 0.002 m/s,
                so pairing it with the computed u* ~ 0.006 m/s inflated the exchange velocity ~3x.
                The salt transfer coefficient follows at R = alpha_h / alpha_s = 35.
                Default: true.
  ICE_CATEGORIES
                Number of equal-area sub-grid ice thickness categories used for the effective
                conductivity (Fichefet & Morales Maqueda 1997). Conduction through the cell-mean
                thickness underestimates growth because the flux is ~1/h; N categories placed at
                (2i-1)h/N raise the effective conductivity by sum(1/(2i-1)), i.e. 1.53 for N=3 and
                1.79 for N=5. Default: 1 (conduct through the mean, no enhancement). Adds
                "_ncat<value>" to the run name.
  KSKEW         Isopycnal skew diffusivity κ_skew (default: per-config; 0 = off;
                "nemo" = NEMO's Treguier et al. 1997 coefficient, Ro² × baroclinic
                growth rate, recomputed each step and depth-uniform;
                "cesm" = CESM's Danabasoglu & Marshall 2007 coefficient,
                κ_ref · clamp(N²/N²_ref, 0.1, 1), recomputed each step with
                full vertical structure;
                "hybrid" = Treguier horizontal × Danabasoglu-Marshall vertical
                shape, with the Rossby-radius cap tightened to 20 km)
  KSYMM         Isopycnal symmetric diffusivity κ_symmetric (default: per-config; 0 = off;
                "nemo" as above, but rising to its reference value in the tropics
                with a floor of one fifth of it, per NEMO nn_aht_ijk_t = 21;
                "cesm" as above)
  BIHARMONIC    Biharmonic viscosity timescale (default: per-config; "nothing" = off)
  BIHVISC       Constant biharmonic viscosity ν in m^4/s (default: unset).
                When set, overrides BIHARMONIC and uses ν directly instead of
                the grid-area-scaled νhb = Az^2 / λ form.
  CB            CATKE bottom-distance coefficient for the shear length scale Cᵇ (default: 0.28).
                It enters as min(Cˢ*depth, Cᵇ*height_above_bottom, ℓᴺ), so it caps the stable
                mixing length through the whole column, not only near the bottom.
  CUNB          CATKE weight on |N²| where the stratification is unstable, Cᵘⁿᵇ (default: 0.0).
                Continues the stratification length w★/√N² onto the unstable branch instead of
                letting it go unbounded and fall back to the geometric min(Cˢ*depth, Cᵇ*hab).
  CF            CATKE weight on the geometric floor where the column is actively convecting and
                weakly sheared, Cᶠ (default: 1.0 = unmodified). The floor max(ℓ★, ℓʰ) becomes
                max(ℓʰ, Cᶠ*ℓ★) only where ℓʰ > 0, so CF acts on convecting cells while CUNB acts
                on parked ones (ℓʰ = 0) — the two populations are disjoint.
  CF0           S²/|N²| below which the CF weight applies in full, Cᶠ⁰ (default: 1e9, i.e. at every
                shear). Lower it to restore the floor in sheared cells; the measured separation between
                convecting and shear-entraining cells is only ~1.4x, so there is no supported setting.
  CFD           Width in S²/|N²| over which the floor ramps back to full strength, Cᶠᵟ (default: 0.75).
                Inert while CF0 is large. Must stay finite.
                THRESHOLD, NOT A KNOB: below Cᵘⁿc/CᵘⁿD = 0.6197 it is exactly a no-op; above it the
                parked TKE equilibrium disappears and e decays to minimum_tke. Use 1.0, not 0.5.
  CP            CATKE convective penetration length coefficient for tracers Cᵉc (default: 0.112).
                Sets how far a convective plume entrains below the unstable layer.
  ICE_FW        Fraction of the sea ice-ocean mass exchange delivered to the ocean, volume and salt
                alike (default: 1). The withheld water leaves the ocean+ice+snow total and
                NORMALIZE_FRESHWATER returns it globally through the free surface, so the global
                budget still closes while the local delivery is scaled.
                Adds "_icefw<value>" to the run name.
  ICE_MELT_MIX  Set to "true" for extra vertical tracer diffusivity where the sea ice is melting into
                the ocean, the same device river runoff already gets. The ice-ocean exchange lands in
                the surface cell, so a melt event leaves a one-cell lid the closure must erode; in
                reality it is stirred through the ice-ocean boundary layer under keels of 1-3 m.
                Adds "_icemix" to the run name.
  ICE_MELT_K    That diffusivity in m^2/s (default: 5e-4). NOT the river value of 0.1, which mixes
                sqrt(2*k*t) ~ 131 m in a day and is convective adjustment rather than boundary-layer
                stirring; 5e-4 gives ~9 m/day.
  ICE_MELT_DEPTH      Depth over which it is applied, in m (default: 10).
  ICE_MELT_THRESHOLD  Melt rate below which it stays off, in m/s (default: 1e-9, about 0.03 m/yr),
                so a column near no net exchange does not flicker it on and off between steps.
  UNDER_ICE_NU  Floor on the vertical viscosity beneath sea ice, in m^2/s, over the top
                UNDER_ICE_NU_DEPTH metres and weighted by the ice concentration, the closure otherwise
                untouched. Bounds how deep the ice stress reaches: under ice CATKE's Ekman layer is
                3-16 m against 11-196 m in open water. McPhee's u*·l with u* ~ 1 cm/s and l ~ 1 m
                gives 1e-2. Off by default; adds "_icenu<nu>" to the run name.
  UNDER_ICE_NU_DEPTH  Depth in m over which it is applied (default: 20). Adds "d<depth>".
  ICE_ARCH      A seasonal ice arch: a basal stress arresting the ice inside a longitude-latitude box.
                "nares" is Nares Strait (77.5-82.5N, 78-58W) from December through July, the observed
                arch season (Kwok 2005, 2010) - the model exports 450 km3/yr there against ~130
                observed. "davis" is Baffin Bay and Davis Strait (65-80N, 80-48W) year round, which
                arrests the Davis delivery outright; psi is linear in the removed Davis export and
                blind to Fram (C21-3, C21-5), and no run has yet varied Davis while holding the
                global ice state fixed. Any other value is read as a raw box "l1,l2,p1,p2".
                Adds "_naresarch", "_davisarch" or "_arch<box>". Off by default.
                ! Davis Strait is far too deep for a basal stress to be physical there: read "davis"
                as a mechanism probe, like ICE_DYNAMICS=false, not as a tuning.
  ICE_ARCH_STRESS  The arresting stress in N/m^2 (default: 100). Adds "<stress>" to the tag.
  ICE_ARCH_MONTHS  First and last month it acts, inclusive, wrapping the year, as "m1,m2".
                Defaults to "12,7" for "nares" and "1,12" for "davis". When set, adds "m<m1>-<m2>" to
                the tag, so a seasonal arch never collides with the year-round run's directory.
  ICE_MELTWATER_TB
                Deliver ice meltwater at the interface temperature Tb instead of at Conservative
                Temperature 0. *** DEFAULT CHANGED 2026-09-03: now "true". Set ICE_MELTWATER_TB=false
                to reproduce runs launched before this date; that is tagged "_nomeltTb". ***
                WARNING: applied to the WHOLE ice mass flux, and top melt is produced near 0, so it
                over-corrects that part -- read it as an upper bound on the correction. ClimaSeaIce produces basal meltwater at
                Tb, so the default hands the ocean ~0.23 W/m2 per (m/yr) of basal melt it never paid
                for. NOTE: applied to the whole ice mass flux, so it over-corrects top melt, which is
                produced near 0. An upper bound on the correction, not the exact treatment.
  ICE_VSF       Set to "true" to deliver that exchange as a virtual salt flux at fixed ocean volume
                instead of as a real volume flux, isolating the volume pathway from the freshwater
                amount. Does not conserve total salt. Overrides ICE_FW.
                Adds "_icevsf" to the run name.
  CLOSURE       Ocean vertical closure: "catke" (default), "simple", "nori", "rbvd",
                "kpp", or "nemo_tke"
                ("simple" = ConvectiveAdjustment + depth-stepped background κ/ν;
                 "nori"   = NORi Richardson-number closure
                            (xkykai/NORiOceanParameterization.jl, vendored);
                 "rbvd"   = Oceananigans' built-in RiBasedVerticalDiffusivity
                            (Richardson-number-based, with built-in κ-clip and
                             time-averaging smoothing);
                 "kpp"    = K-Profile Parameterization (Large 1994 / MITgcm,
                            vendored in `KPP/`);
                 "nemo_tke" = NEMO 3.6 TKE scheme (Blanke & Delecluse 1993,
                            Madec 2017), OMIP-2 ORCAOne preset: prognostic e,
                            gradient-limited length scale, Langmuir +
                            Mellor-Blumberg wave penetration + EVD.
                            Vendored in `NEMOTKE/`;
                 all ignore CB)
  PARTIAL_CELLS Set to "true" to use partial bottom cells (PartialCellBottom) instead of
                full-cell GridFittedBottom bathymetry. Resolves sill depths and slopes
                continuously; targets the too-shallow NADW from staircased overflows.
                Adds "_pcells" to the run name.
  BBL_KAPPA     Diffusive bottom boundary layer coefficient in m² s⁻¹ (NEMO rn_ahtbbl;
                its ORCA reference value is 1000). Dense water upslope of a deeper
                neighbour is diffused along the bottom, mimicking the gravity current a
                z-coordinate model otherwise mixes away within a grid cell or two.
                Off by default; adds "_bbl<κ>" to the run name.
  BBL_GAMMA     Advective bottom boundary layer coefficient in seconds (NEMO rn_gambbl,
                reference value 10; Campin & Goosse 1999, NEMO nn_bbl_adv = 2). Opens an
                overturning circuit down the step: dense shelf water enters the bottom of
                the deep column, the water it displaces rises through that column and
                returns to the shelf. The deep cell is therefore flushed towards the shelf
                density rather than mixed towards a mean, which is the ceiling BBL_KAPPA
                cannot pass. Downslope speed is u = γ g Δρ/ρ₀, so the transport shuts off
                as the contrast is consumed. Combines with BBL_KAPPA; the two act on
                different failures. Adds "_cg<γ>" to the run name.
  OVERFLOW_RESTORE
                Diagnostic only, in DAYS. Pins T and S on the East Greenland slope
                (36-26 W, 62-66.5 N, below 1500 m) to observed Denmark Strait Overflow
                Water, Θ = 2.0 and Sᴬ = 35.065, σθ = 27.892. This is not a
                parameterisation: it forces the delivered density so the question "does
                the overflow deficit set the AMOC?" can be answered without first solving
                how to get dense water down a staircase. Both BBL schemes left the
                delivered density unchanged, so neither tested that. Adds "_dsow<days>"
                to the run name.
  LAB_RESTORE   Diagnostic only, in DAYS. Restores SALINITY ONLY in the deep Labrador
                interior (65-40 W, 52-66 N, columns whose bottom is below 2000 m) above
                200 m toward WOA Annual Absolute Salinity. Campaign 27 priced the
                convection: removing d metres of freshwater from that column takes the
                autumn resistance D(1000 m) from orca's 1.170 to noicedyn's 0.750 at
                d = 1.54-1.65 m, and the model's measured freshwater-content excess over
                WOA there is 1.638 m - the same number. This pins the box to observation
                so "does the Labrador freshwater bias set the AMOC?" can be answered
                without first solving what causes the bias. It must be a restoring, not a
                one-off correction: davisarch stopped adding 177 g/kg m of freshwater to
                that box over sixteen years and kept 2, so a fixed flux is consumed by
                the lateral compensation within a year. Salinity only - the year-1
                resistance change is 100% haline (haline 0.785 against a full 0.783,
                thermal 0.998). Adds "_labrest<days>" to the run name.
  ML_TAPER      Set to "true" to ramp the isopycnal-closure slopes linearly to zero
                from the mixed-layer base to the surface (Danabasoglu et al. 2008;
                NEMO ldfslp). Off by default; adds "_mltaper" to the run name.
  WIND_VELOCITY Set to "true" to use absolute wind (Δu = u_atm) in the bulk
                formula instead of the OMIP-2 default relative wind
                (Δu = u_atm − u_ocean). For isolating ACC-current feedback.
  NZ            Number of vertical levels. Default: per config (70 for orca, 100 for
                quarter/twelfth degree). Adds "_nz<N>" to the run name when overridden.

  DZ_TOP        Target thickness of the top (surface) cell in meters. If set,
                the ExponentialDiscretization scale is found by bisection so
                that Δz of the surface level matches DZ_TOP within ~0.1%.
                Must satisfy 0 < DZ_TOP < depth/Nz. Default: unset (scale=1300).

Equatorial-MLD tuning knobs (closure parameters; configuration switches):
  NORMALIZE_SALINITY "true" (default) applies the conservative, salt-conserving
                surface-salinity restoring (zero global mean). Set to "false" to use
                the raw un-normalized restoring (the old, non-conserving behavior),
                e.g. for A/B comparison. Default: true.
  NORMALIZE_FRESHWATER Removes the global mean of the atmospheric surface freshwater
                flux, holding the global ocean volume fixed (standard OMIP-2 practice;
                the sea-ice exchange is excluded and the freshwater heat content is
                corrected alongside so the carried heat flux is unchanged). Options:
                  none      no correction (default; "false" also accepted)
                  timestep  remove the instantaneous global mean; pins volume exactly
                            but also removes the seasonal land-water-storage cycle,
                            which is comparable in size ("true" also accepted)
                  annual    remove a running mean relaxed over a year; keeps the
                            seasonal cycle, spins up over the first ~2 years
  ISOPYCNAL     Which closure carries GM/Redi. "standard" (default) is
                IsopycnalSkewSymmetricDiffusivity; "triad" is
                TriadIsopycnalSkewSymmetricDiffusivity, which evaluates the tensor on the
                Griffies triad stencil. The triad closure carries no skew flux formulation,
                so ISOPYCNAL=triad requires SKEW_FORMULATION=diffusive (the default).
                Adds "_triad" to the run name.
  SKEW_FORMULATION How the GM skew transport is applied. Ignored when KSKEW=0.
                  diffusive  add it to the tracer flux (default)
                  advective  build the eddy-induced velocity and advect with it,
                             which also exposes the bolus transport as a model
                             field. Equivalent continuously, not discretely, so
                             the two give different answers run-to-run.
                  boundary_value  the Ferrari, Griffies, Nurser & Vallis (2010)
                             transport: each column solves
                             (c d²/dz² - N²) Y = -N² Y_GM with Y = 0 at the
                             surface and the bottom, low-passing the baroclinic
                             modes. Satisfies the boundary conditions without
                             tapering and interpolates through weakly stratified
                             layers with no floor on N². Applied advectively;
                             requires a depth-independent KSKEW (a number or
                             "nemo" -- "cesm"/"hybrid" are rejected).
                Adds "_gmadv" or "_gmbvp" to the run name.
  BVP_MODE      Mode number M in the speed c = max(BVP_CMIN, (M pi)^-1 int N dz) that
                weights the second-order operator. Larger M filters less and gives a
                larger transport; M=1 is the first baroclinic mode, whose amplitude is
                about half the truncated GM transport. Default: 2.
                Only used when SKEW_FORMULATION=boundary_value.
  BVP_CMIN      Floor c_min on that speed, in m/s. Keeps the transport bounded in
                weakly stratified columns. Default: 0.1.
  RESTORING_UNDER_ICE Set to "false" to stop the surface-salinity restoring acting under
                sea ice (weighted by the open-water fraction 1-ℵ, with the zero-mean
                correction spread over open water only, so no net salt is injected).
                WOA is poorly constrained beneath ice and the restoring there fights the
                ice-ocean salt flux. Requires NORMALIZE_SALINITY=true. Adds "_noicerest"
                to the run name. Default: true (OMIP-2 convention).
  RIVER_SPREAD  Radius in degrees over which each river/iceberg mouth's discharge is
                divided equally among the surrounding wet cells. A geographic radius
                keeps the freshwater flux per unit area resolution-independent; raise it
                if a refined grid drives coastal salinity to zero. Set to "none" for the
                historical cell-count footprint instead. Per-config default: 1.2 for the
                refined grids, "none" for orca/test (pinned to their old behaviour).
  RIVER_SPREAD_CELLS  Cells in that footprint, nearest first: a cap when RIVER_SPREAD is
                a radius, the footprint itself when it is "none". Per-config default:
                unset (uncapped) for the refined grids, 8 for orca/test.
  RIVER_DIVERSION  Fraction (0-1) of the river and iceberg discharge landing in the Atlantic
                that is delivered to the Pacific instead, at the latitude it was diverted from.
                The global freshwater input is unchanged; only the basin receiving it is.
                Adds "_divert<value>" to the run name. Default: 0.
  RIVER_MIXING  Set to "false" to disable the extra vertical diffusivity applied over
                the spread footprint. Default: true.
  RIVER_MIXING_K      That diffusivity in m^2/s (default: 0.1).
  RIVER_MIXING_DEPTH  Depth in m over which it is applied (default: 10).
  CATKE_CWUSTAR `Cᵂu★` of CATKEEquation: surface shear-driven TKE flux
                coefficient. Higher → more wind-injected TKE → deeper
                equatorial ML. Default (Oceananigans): 3.179.
  BACKGROUND_K  Interior background tracer diffusivity κ added underneath CATKE (or RBVD).
                Either a number in m^2/s (uniform), or "bryan_lewis" for the Bryan & Lewis
                (1979) depth profile, 3e-5 in the upper ocean rising across ~2500 m to 1.3e-4
                in the abyss (deep upwelling without diffusing the thermocline), or
                "abyssal_henyey" for Henyey in the thermocline with the same arctangent
                enhancement added beneath it, reaching +5e-5 at 5000 m — the abyssal upwelling
                without the upper-ocean value that sets the drift.
                Default: unset, i.e. the Henyey et al. (1986) latitudinal internal-wave
                scaling κ = max(2e-6, 1e-5 |sin φ|), 2e-6 at the equator to 1e-5 at the poles.
                Raising it strengthens the diapycnal upwelling that closes the AMOC lower limb
                (Bryan 1987: AMOC ~ κ^(2/3)). Adds "_bgk<value>" to the run name.
  CHLOROPHYLL   Optics for the penetrating shortwave. "seawifs" (default) uses the SeaWiFS monthly
                climatology, cycled annually, so the decay scale varies in space and season. A number
                gives globally uniform optics with that chlorophyll in mg/m^3; 0.147 reproduces the
                Jerlov Type I 23 m scale that was the default before the climatology landed.
                "none" removes the penetrating scheme, so the whole shortwave is applied to the
                surface heat flux and the vertical closure's surface buoyancy flux sees it.
                Adds "_chl<value>" to the run name unless it is "seawifs".
  BACKGROUND_NU Uniform interior background viscosity ν in m^2/s. Unset resolves to
                `default_background_viscosity = 1e-4` for every closure that takes one
                (omip_simulation.jl:1104-1106), so a run without "_bgnu" in its name carries
                1e-4, NOT the 3e-5 the production runs set explicitly. Adds "_bgnu<value>".
                Both are rejected by the closures that carry their own background
                (simple/nori/kpp/nemo_tke).
  IMEX_DRAG     Whether the bottom and immersed quadratic drag are applied semi-implicitly
                (J = lambda * u_b with lambda = -mu |u| in the vertical solver's diagonal).
                Default: true. "false" applies both explicitly, the treatment this replaced,
                and adds "_explicitdrag" to the run name. The two differ most in shallow fast
                cells such as narrow straits, where mu |u| dt / dz is order ten. This is the
                change bisected to commit 3c60be3b. The surface momentum flux is always
                semi-implicit and is not affected.
  DRAG_UB       Unresolved velocity u_b in m/s added in quadrature to the resolved speed in the
                quadratic bottom drag, tau = mu u sqrt(u_b^2 + |u|^2), standing for tides and
                other motions the grid does not carry. Default: 0, the purely resolved drag.
                Adds "_ub<value>" to the run name. It bites where the flow is slow: at |u| =
                1 cm/s and u_b = 0.1 the stress is ten times the resolved one, at 50 cm/s it is
                within a percent of it. GFDL's OM4 uses 0.1.
  BAROTROPIC_SUBSTEPS
                *** TAGGED "_substeps<val>" SINCE 2026-09-03 when it differs from the config default
                (orca 300, quarterdegree/twelfthdegree 200, else 100). Before that it had NO tag, so
                runs launched earlier record their substep count nowhere on disk -- read it from
                /proc/<pid>/cmdline while the job is alive. orca's dt = 5400 s NEEDS 300: at the
                generic 100 the gravity-wave Courant number is 1.16 and the free surface NaNs on
                step 1.
                Split-explicit substeps the free surface takes per DT. Default: unset, i.e. 200
                for quarterdegree/twelfthdegree and 100 otherwise. The barotropic gravity wave
                must stay inside a substep, so raising DT above the config default needs
                proportionally more; too few blows the free surface up on the first step. The
                run warns with the count the grid needs.

Environment variables (I/O & runtime):
  BACKEND_SIZE  Number of JRA55 time indices kept in memory (default: 240,
                i.e. 30 days of 3-hourly data ≈ 2 GB RAM for 11 variables)
  FORCING_DIR   Path to JRA55 forcing data (default: ${DATA}forcing_data)
  STAGING_DIR   Base directory for JRA55 staging (default: ./staged_data).
                A per-run subdirectory (STAGING_DIR/<run_name>) is created
                with symlinks from FORCING_DIR; files are progressively
                copied ahead of each simulated year.
                Keeps ~50 GB per run (current + next year).
  OUTPUT_DIR    Base directory for simulation output (default: ".").
                A per-run subdirectory (OUTPUT_DIR/<run_name>_run) is created.
  NODE          Pin job to a specific node (default: 2904)
  THREADS       Number of Julia threads / CPUs per task (default: 4)
  PROFILE       Set to "true" for nsys profiling. Also:
                  - disables the OMIP diagnostic output writers
                    (Average / JLD2 dumps / checkpoint / KE spectrum)
                    — they add per-iteration I/O and compute! overhead
                    that contaminates the trace;
                  - drops `pickup` from `run!` so the simulation always
                    starts from scratch.

Examples:
  ./launch.sh orca
  NCAR=true ./launch.sh orca
  NCAR=true SNOW=true ./launch.sh orca
  CORRECTED=true SNOW=true ./launch.sh orca
  CB=0.1 NCAR=true ./launch.sh orca
  KSKEW=1000 KSYMM=500 ./launch.sh orca
  KSKEW=0 ./launch.sh orca                    # disable eddy closure
  BIHARMONIC=5days ./launch.sh orca           # custom biharmonic timescale
  BIHARMONIC=nothing ./launch.sh orca         # disable biharmonic viscosity
  BIHVISC=1e12 ./launch.sh orca               # constant biharmonic viscosity ν=1e12 m^4/s
  DZ_TOP=2 ./launch.sh orca                   # 2 m top cell (scale chosen by bisection)
  IC_CONDITIONS=blended ./launch.sh orca      # January WOA Monthly blend + summer ice both hemispheres
  CATKE_CWUSTAR=5.0 ./launch.sh orca          # stronger surface TKE injection in CATKE
  BACKGROUND_K=3e-5 ./launch.sh orca          # uniform background κ = 3e-5 m^2/s (AMOC sensitivity)
  BACKGROUND_K=1e-4 ./launch.sh orca          # uniform background κ = 1e-4 m^2/s
  BACKGROUND_K=bryan_lewis ./launch.sh orca   # Bryan-Lewis depth profile, 3e-5 → 1.3e-4
  SKEW_FORMULATION=boundary_value ./launch.sh orca            # Ferrari et al. (2010) GM transport
  SKEW_FORMULATION=boundary_value KSKEW=nemo ./launch.sh orca # ... with the Treguier coefficient
  FORCING_DIR=/other/path/forcing_data STAGING_DIR=/scratch/staged ./launch.sh orca
  PROFILE=true ./launch.sh orca
USAGE
}

CONFIG="${1:-}"
if [[ -z "$CONFIG" ]]; then
    usage
    exit 1
fi
shift || true

case "$CONFIG" in
    halfdegree|half_degree)
        CONFIG="halfdegree"
        ;;
    quarterdegree|quarter_degree)
        CONFIG="quarterdegree"
        ;;
    orca|twelfthdegree) ;;
    -h|--help)
        usage
        exit 0
        ;;
    *)
        echo "Error: unknown configuration '$CONFIG'" >&2
        usage
        exit 1
        ;;
esac

# ── Per-config defaults ───────────────────────────────────────────────
#                     KSKEW  KSYMM  NZ   DT          BIHARMONIC  ARCH                                             GPUS  EXTRA_USING                              FILE_SPLIT  RUN_CMD
case "$CONFIG" in
    halfdegree)
        DEFAULT_KSKEW=250;  DEFAULT_KSYMM=100; DEFAULT_NZ=70;  DEFAULT_DT="30minutes"; DEFAULT_DZ_TOP="2.0"
        DEFAULT_BIHARMONIC="40days"; ARCH="GPU()"; GPUS_PER_NODE=1
        EXTRA_USING=""; FILE_SPLIT=""
        RUN_CMD="sim.stop_time = 300 * 365days
run!(sim, pickup=:latest)"
        ;;
    quarterdegree)
        DEFAULT_KSKEW=0;    DEFAULT_KSYMM=0;   DEFAULT_NZ=100; DEFAULT_DT="20minutes"; DEFAULT_DZ_TOP="1.5"
        DEFAULT_BIHARMONIC="nothing"; ARCH="GPU()"; GPUS_PER_NODE=1
        EXTRA_USING="using Oceananigans.DistributedComputations"
        FILE_SPLIT=""
        RUN_CMD="sim.stop_time = 300 * 365days
run!(sim, pickup =:latest)"
        ;;
    orca)
        DEFAULT_KSKEW=800;  DEFAULT_KSYMM=800; DEFAULT_NZ=70;  DEFAULT_DT=5400;        DEFAULT_DZ_TOP="1.5"
        DEFAULT_BIHARMONIC="50days"; ARCH="GPU()"; GPUS_PER_NODE=1
        EXTRA_USING=""; FILE_SPLIT=""
        RUN_CMD="sim.stop_time = 300 * 365days
run!(sim; pickup = :latest)
"
        ;;
    twelfthdegree)
        DEFAULT_KSKEW=0;    DEFAULT_KSYMM=0;   DEFAULT_NZ=100; DEFAULT_DT="6minutes";  DEFAULT_DZ_TOP="1.5"
        DEFAULT_BIHARMONIC="nothing"; ARCH="Distributed(GPU(), partition=Partition(1, 4))"; GPUS_PER_NODE=4
        EXTRA_USING="using Oceananigans.DistributedComputations"
        FILE_SPLIT="file_splitting_interval = 180days,"
        RUN_CMD="sim.stop_time = 181days
run!(sim)

sim.Δt = 15minutes
sim.stop_time = 300 * 365days
run!(sim; pickup = true)"
        ;;
esac

# Profile mode: drop the `pickup=...` from `run!` so the simulation
# always starts from scratch. Resuming a checkpoint would skip the
# initialization phase the profiler needs to see, and the choice of
# checkpoint is irrelevant for kernel-cost measurements.
if [[ "${PROFILE:-false}" == "true" ]]; then
  RUN_CMD="sim.stop_iteration = 200; run!(sim)"
fi

# 0 means "no eddy closure" (maps to Julia `nothing`)
case "$CONFIG" in
  orca)                        DEFAULT_SUBSTEPS=300 ;;
  quarterdegree|twelfthdegree) DEFAULT_SUBSTEPS=200 ;;
  *)                           DEFAULT_SUBSTEPS=100 ;;
esac

export KSKEW="${KSKEW:-$DEFAULT_KSKEW}"
export KSYMM="${KSYMM:-$DEFAULT_KSYMM}"
export DT="${DT:-$DEFAULT_DT}"
export DZ_TOP="${DZ_TOP:-$DEFAULT_DZ_TOP}"
export NZ="${NZ:-$DEFAULT_NZ}"
export BIHARMONIC="${BIHARMONIC:-$DEFAULT_BIHARMONIC}"
KSKEW_JULIA="$KSKEW"; [[ "$KSKEW" == "0" ]] && KSKEW_JULIA="nothing"
KSYMM_JULIA="$KSYMM"; [[ "$KSYMM" == "0" ]] && KSYMM_JULIA="nothing"
[[ "$KSKEW" == "nemo" ]] && KSKEW_JULIA=":nemo"
[[ "$KSYMM" == "nemo" ]] && KSYMM_JULIA=":nemo"
[[ "$KSKEW" == "cesm" ]] && KSKEW_JULIA=":cesm"
[[ "$KSYMM" == "cesm" ]] && KSYMM_JULIA=":cesm"
[[ "$KSKEW" == "hybrid" ]] && KSKEW_JULIA=":hybrid"
[[ "$KSYMM" == "hybrid" ]] && KSYMM_JULIA=":hybrid"
export KSKEW_JULIA KSYMM_JULIA
export NZ DT ARCH EXTRA_USING FILE_SPLIT RUN_CMD

# ── Boundary reconstruction ───────────────────────────────────────────
# BUFFER_ORDER=1 is the pre-CWENOZ spelling of BOUNDARY_SCHEME=upwind; both name the same run.
TRACER_ORDER="${TRACER_ORDER:-7}"
BUFFER_ORDER="${BUFFER_ORDER:-3}"
BOUNDARY_SCHEME="${BOUNDARY_SCHEME:-default}"

case "$BUFFER_ORDER" in
  3) ;;
  1) if [[ "$BOUNDARY_SCHEME" != "default" && "$BOUNDARY_SCHEME" != "upwind" ]]; then
         echo "BUFFER_ORDER=1 and BOUNDARY_SCHEME=$BOUNDARY_SCHEME ask for two different boundary reconstructions" >&2
         exit 1
     fi
     BOUNDARY_SCHEME="upwind" ;;
  *) echo "BUFFER_ORDER must be 1 or 3; a chain truncated above third order is not expressible as a boundary_scheme" >&2
     exit 1 ;;
esac

TRACER_BOUNDARY_SCHEME="${TRACER_BOUNDARY_SCHEME:-$BOUNDARY_SCHEME}"
MOMENTUM_BOUNDARY_SCHEME="${MOMENTUM_BOUNDARY_SCHEME:-$BOUNDARY_SCHEME}"

for scheme_name in BOUNDARY_SCHEME TRACER_BOUNDARY_SCHEME MOMENTUM_BOUNDARY_SCHEME; do
  case "${!scheme_name}" in
    default|upwind|cwenoz) ;;
    *) echo "$scheme_name must be default, upwind or cwenoz, got '${!scheme_name}'" >&2; exit 1 ;;
  esac
done
export TRACER_ORDER BUFFER_ORDER BOUNDARY_SCHEME TRACER_BOUNDARY_SCHEME MOMENTUM_BOUNDARY_SCHEME

# ── Initial-condition preset ──────────────────────────────────────────
IC_CONDITIONS="${IC_CONDITIONS:-default}"
case "$IC_CONDITIONS" in
  default|summerice) ;;
  blended) IC_BLEND="${IC_BLEND:-500}" ;;
  *) echo "Unknown IC_CONDITIONS: '$IC_CONDITIONS' (expected 'default', 'summerice' or 'blended')" >&2; exit 1 ;;
esac
export IC_CONDITIONS IC_BLEND

# ── Build run name from config + options ──────────────────────────────
RUN_NAME="$CONFIG"
[[ "${CORRECTED:-true}" != "true" ]]           && RUN_NAME="${RUN_NAME}_rawflux"
[[ "${NCAR:-false}" == "true" ]]               && RUN_NAME="${RUN_NAME}_ncar"
[[ "${SNOW:-true}" != "true" ]]                && RUN_NAME="${RUN_NAME}_nosnow"
[[ "${ICE_DYNAMICS:-true}" == "false" ]]       && RUN_NAME="${RUN_NAME}_noicedyn"
[[ "${ICE_LATERAL:-no_slip}" != "no_slip" ]]   && RUN_NAME="${RUN_NAME}_freeslip"
[[ "${ICE_BASAL:-true}" != "true" ]]           && RUN_NAME="${RUN_NAME}_nolandfast"
[[ "${ICE_DRAG:-5.5e-3}" != "5.5e-3" ]]        && RUN_NAME="${RUN_NAME}_cio${ICE_DRAG}"
[[ "${ICE_HEAT_TRANSFER:-0.0057}" != "0.0057" ]] && RUN_NAME="${RUN_NAME}_ah${ICE_HEAT_TRANSFER}"
[[ -n "${ICE_PSTAR:-}" ]]                        && RUN_NAME="${RUN_NAME}_pstar${ICE_PSTAR}"
[[ -n "${ICE_SALINITY:-}" ]]                     && RUN_NAME="${RUN_NAME}_sice${ICE_SALINITY}"
[[ "${ICE_DRAGREF:-6}" != "6" ]]                 && RUN_NAME="${RUN_NAME}_dragref${ICE_DRAGREF}"
[[ "${ICE_LIQUIDUS:-teos10}" != "teos10" ]]      && RUN_NAME="${RUN_NAME}_liq${ICE_LIQUIDUS}"
[[ "${ICE_Z0:-5e-4}" != "5e-4" ]]                && RUN_NAME="${RUN_NAME}_icez0${ICE_Z0}"
[[ -n "${SNOW_CATEGORIES:-}" && "${SNOW_CATEGORIES}" != "${ICE_CATEGORIES:-4}" ]] \
                                                 && RUN_NAME="${RUN_NAME}_snowcat${SNOW_CATEGORIES}"
[[ -n "${ICE_ITD_SHAPE:-}" ]]                    && RUN_NAME="${RUN_NAME}_itd${ICE_ITD_SHAPE//,/-}"
[[ "${TRACER_ORDER:-7}" != "7" ]]                && RUN_NAME="${RUN_NAME}_tracer${TRACER_ORDER}"
# A scheme shared by tracers and momentum keeps the undecorated tag, so BUFFER_ORDER=1 still names "_buford1".
if [[ "$TRACER_BOUNDARY_SCHEME" == "$MOMENTUM_BOUNDARY_SCHEME" ]]; then
  [[ "$TRACER_BOUNDARY_SCHEME" == "upwind" ]]    && RUN_NAME="${RUN_NAME}_buford1"
  [[ "$TRACER_BOUNDARY_SCHEME" == "cwenoz" ]]    && RUN_NAME="${RUN_NAME}_cwenoz"
else
  [[ "$TRACER_BOUNDARY_SCHEME" != "default" ]]   && RUN_NAME="${RUN_NAME}_tr${TRACER_BOUNDARY_SCHEME}"
  [[ "$MOMENTUM_BOUNDARY_SCHEME" != "default" ]] && RUN_NAME="${RUN_NAME}_mom${MOMENTUM_BOUNDARY_SCHEME}"
fi
[[ -n "${TEMPERATURE_VARIATION:-}" ]]            && RUN_NAME="${RUN_NAME}_vT${TEMPERATURE_VARIATION}"
[[ -n "${SALINITY_VARIATION:-}" ]]               && RUN_NAME="${RUN_NAME}_vS${SALINITY_VARIATION}"
[[ -n "${MOMENTUM_VARIATION:-}" ]]               && RUN_NAME="${RUN_NAME}_vu${MOMENTUM_VARIATION}"
[[ "${ICE_TILT:-false}" == "true" ]]             && RUN_NAME="${RUN_NAME}_icetilt"
[[ -n "${IC_BLEND:-}" ]]                         && RUN_NAME="${RUN_NAME}_icblend${IC_BLEND}"
[[ "$IC_CONDITIONS" != "default" ]]              && RUN_NAME="${RUN_NAME}_summerice"
[[ "${ICE_CATEGORIES:-4}" != "4" ]]              && RUN_NAME="${RUN_NAME}_ncat${ICE_CATEGORIES}"
[[ "${CLOSURE:-catke}" == "simple"   ]]        && RUN_NAME="${RUN_NAME}_simple"
[[ "${CLOSURE:-catke}" == "nori"     ]]        && RUN_NAME="${RUN_NAME}_nori"
[[ "${CLOSURE:-catke}" == "rbvd"     ]]        && RUN_NAME="${RUN_NAME}_rbvd"
[[ "${CLOSURE:-catke}" == "kpp"      ]]        && RUN_NAME="${RUN_NAME}_kpp"
[[ "${CLOSURE:-catke}" == "nemo_tke" ]]        && RUN_NAME="${RUN_NAME}_nemotke"
[[ "${WIND_VELOCITY:-false}" == "true" ]]      && RUN_NAME="${RUN_NAME}_wind"
[[ "${ML_TAPER:-false}" == "true" ]]           && RUN_NAME="${RUN_NAME}_mltaper"
[[ "${ISOPYCNAL:-standard}" == "triad" ]]      && RUN_NAME="${RUN_NAME}_triad"
[[ "${PARTIAL_CELLS:-false}" == "true" ]]      && RUN_NAME="${RUN_NAME}_pcells"
[[ -n "${BBL_KAPPA:-}" ]]                      && RUN_NAME="${RUN_NAME}_bbl${BBL_KAPPA}"
[[ -n "${BBL_GAMMA:-}" ]]                      && RUN_NAME="${RUN_NAME}_cg${BBL_GAMMA}"
[[ -n "${OVERFLOW_RESTORE:-}" ]]               && RUN_NAME="${RUN_NAME}_dsow${OVERFLOW_RESTORE}"
[[ -n "${LAB_RESTORE:-}" ]]                     && RUN_NAME="${RUN_NAME}_labrest${LAB_RESTORE}"
[[ "${NORMALIZE_SALINITY:-true}" == "false" ]] && RUN_NAME="${RUN_NAME}_rawsalt"
[[ "${RESTORING_UNDER_ICE:-true}" == "false" ]] && RUN_NAME="${RUN_NAME}_noicerest"
case "${NORMALIZE_FRESHWATER:-timestep}" in
  none|false)    RUN_NAME="${RUN_NAME}_fwnone" ;;
  annual)        RUN_NAME="${RUN_NAME}_fwnormann" ;;
esac
[[ -n "${RIVER_DIVERSION:-}" ]]                && RUN_NAME="${RUN_NAME}_divert${RIVER_DIVERSION}"
[[ "${SKEW_FORMULATION:-diffusive}" == "advective" ]]      && RUN_NAME="${RUN_NAME}_gmadv"
[[ "${SKEW_FORMULATION:-diffusive}" == "boundary_value" ]] && RUN_NAME="${RUN_NAME}_gmbvp"
[[ -n "${BVP_MODE:-}" ]]                       && RUN_NAME="${RUN_NAME}_bvpm${BVP_MODE}"
[[ -n "${BVP_CMIN:-}" ]]                       && RUN_NAME="${RUN_NAME}_bvpc${BVP_CMIN}"
[[ "${CB:-0.01}" != "0.01" ]]                  && RUN_NAME="${RUN_NAME}_cb${CB}"
[[ -n "${CUNB:-}" ]] && [[ "${CUNB}" != "0.0" ]] && RUN_NAME="${RUN_NAME}_cunb${CUNB}"
[[ -n "${CF:-}" ]]   && [[ "${CF}" != "1.0" ]]   && RUN_NAME="${RUN_NAME}_cf${CF}"
[[ -n "${CF0:-}" ]]  && [[ "${CF0}" != "1e9" ]]  && RUN_NAME="${RUN_NAME}_cf0${CF0}"
[[ -n "${CFD:-}" ]]  && [[ "${CFD}" != "0.75" ]] && RUN_NAME="${RUN_NAME}_cfd${CFD}"
[[ -n "${CP:-}" ]]                             && RUN_NAME="${RUN_NAME}_cp${CP}"
[[ -n "${ICE_FW:-}" ]]                         && RUN_NAME="${RUN_NAME}_icefw${ICE_FW}"
[[ "${ICE_VSF:-false}" == "true" ]]            && RUN_NAME="${RUN_NAME}_icevsf"
[[ "${ICE_MELTWATER_TB:-true}" != "true" ]]          && RUN_NAME="${RUN_NAME}_nomeltTb"
[[ "${ICE_MELT_MIX:-false}" == "true" ]]       && RUN_NAME="${RUN_NAME}_icemix"
[[ -n "${ICE_MELT_K:-}" ]]                     && RUN_NAME="${RUN_NAME}k${ICE_MELT_K}"
[[ -n "${ICE_MELT_DEPTH:-}" ]]                 && RUN_NAME="${RUN_NAME}d${ICE_MELT_DEPTH}"
[[ -n "${UNDER_ICE_NU:-}" ]]                   && RUN_NAME="${RUN_NAME}_icenu${UNDER_ICE_NU}"
[[ -n "${UNDER_ICE_NU_DEPTH:-}" ]]             && RUN_NAME="${RUN_NAME}d${UNDER_ICE_NU_DEPTH}"
case "${ICE_ARCH:-}" in
  "")     ;;
  nares)  RUN_NAME="${RUN_NAME}_naresarch${ICE_ARCH_STRESS:-}" ;;
  davis)  RUN_NAME="${RUN_NAME}_davisarch${ICE_ARCH_STRESS:-}" ;;
  *)      RUN_NAME="${RUN_NAME}_arch$(echo "$ICE_ARCH" | tr -d ' ' | tr ',' '_')${ICE_ARCH_STRESS:-}" ;;
esac
# The arch season is physics: an untagged ICE_ARCH_MONTHS would give a seasonal arch the same run name as
# the year-round one, and a relaunch into that directory resumes its checkpoint instead of starting.
[[ -n "${ICE_ARCH:-}" && -n "${ICE_ARCH_MONTHS:-}" ]] && \
  RUN_NAME="${RUN_NAME}m$(echo "$ICE_ARCH_MONTHS" | tr -d ' ' | tr ',' '-')"
[[ "$KSKEW" != "$DEFAULT_KSKEW" ]]             && RUN_NAME="${RUN_NAME}_kskew${KSKEW}"
[[ "$KSYMM" != "$DEFAULT_KSYMM" ]]             && RUN_NAME="${RUN_NAME}_ksymm${KSYMM}"
[[ "$BIHARMONIC" != "$DEFAULT_BIHARMONIC" ]]   && RUN_NAME="${RUN_NAME}_bih${BIHARMONIC}"
[[ "$DT" != "$DEFAULT_DT" ]]                   && RUN_NAME="${RUN_NAME}_dt${DT}"
[[ "${BAROTROPIC_SUBSTEPS:-$DEFAULT_SUBSTEPS}" != "$DEFAULT_SUBSTEPS" ]] && RUN_NAME="${RUN_NAME}_substeps${BAROTROPIC_SUBSTEPS}"
[[ -n "${BIHVISC:-}" ]]                        && RUN_NAME="${RUN_NAME}_bihvisc${BIHVISC}"
[[ "$DZ_TOP" != "$DEFAULT_DZ_TOP" ]]           && RUN_NAME="${RUN_NAME}_dz${DZ_TOP}"
[[ "$NZ" != "$DEFAULT_NZ" ]]                    && RUN_NAME="${RUN_NAME}_nz${NZ}"
[[ -n "${CATKE_CWUSTAR:-}" ]]                  && RUN_NAME="${RUN_NAME}_cwu${CATKE_CWUSTAR}"
[[ -n "${BACKGROUND_K:-}" ]]                   && RUN_NAME="${RUN_NAME}_bgk${BACKGROUND_K}"
[[ "${BACKGROUND_NU:-3e-5}" != "3e-5" ]]       && RUN_NAME="${RUN_NAME}_bgnu${BACKGROUND_NU}"
[[ "${CHLOROPHYLL:-seawifs}" != "seawifs" ]]   && RUN_NAME="${RUN_NAME}_chl${CHLOROPHYLL}"
[[ "${IMEX_DRAG:-true}" == "false" ]]            && RUN_NAME="${RUN_NAME}_explicitdrag"
[[ -n "${DRAG_UB:-}" && "${DRAG_UB:-0}" != "0" ]] && RUN_NAME="${RUN_NAME}_ub${DRAG_UB}"
[[ "${PVEL:-0.254}" != "0.254" ]]              && RUN_NAME="${RUN_NAME}_pvel${PVEL}"
[[ -n "${RUN_TAG:-}" ]]                        && RUN_NAME="${RUN_NAME}_${RUN_TAG}"

# A run launched before the 2026-09-03 rename carries tags for what are now defaults, so the name built
# above cannot match its directory and `pickup` would start it from zero. RUN_NAME_OVERRIDE resumes it.
[[ -n "${RUN_NAME_OVERRIDE:-}" ]] && RUN_NAME="$RUN_NAME_OVERRIDE"

REPORT_NAME="${REPORT_NAME:-${RUN_NAME}_report}"
JOB_NAME="${JOB_NAME:-$RUN_NAME}"

SBATCH_ARGS=()
PARTITION="${PARTITION:-pi_raffaele}"

NODE="${NODE:-2904}"
if [[ "${PARTITION}" != "default" && -n "${NODE}" ]]; then
    SBATCH_ARGS+=(-w "node${NODE}")
fi
SBATCH_ARGS+=(--gres="gpu:${GPUS_PER_NODE}")
# Override the heredoc default `--ntasks-per-node=1` so distributed configs
# (twelfthdegree: 1×4 partition) actually launch GPUS_PER_NODE MPI ranks.
SBATCH_ARGS+=(--ntasks-per-node="${GPUS_PER_NODE}")

export THREADS="${THREADS:-8}"
SBATCH_ARGS+=(--cpus-per-task="${THREADS}")

if [[ "${PARTITION}" != "default" ]]; then
    SBATCH_ARGS+=(--partition="${PARTITION}")
fi

if [[ "${PARTITION}" == "default" ]]; then
    TIME="${TIME:-05:00:00}"
else
    TIME="${TIME:-120:00:00}"
fi
SBATCH_ARGS+=(--time="${TIME}")

MEM="${MEM:-150GB}"
SBATCH_ARGS+=(--mem="${MEM}")

if [[ "${PROFILE:-false}" == "true" ]]; then
    SBATCH_ARGS+=(-o "${RUN_NAME}_profile.out")
    SBATCH_ARGS+=(-e "${RUN_NAME}_profile.err")
    SBATCH_ARGS+=(-J "${JOB_NAME}_profile")
    SBATCH_ARGS+=(--export="ALL,PROFILE=true,REPORT_NAME=${REPORT_NAME},CONFIG=${CONFIG},RUN_NAME=${RUN_NAME}")
else
    SBATCH_ARGS+=(-o "${RUN_NAME}.out")
    SBATCH_ARGS+=(-e "${RUN_NAME}.err")
    SBATCH_ARGS+=(-J "$JOB_NAME")
    SBATCH_ARGS+=(--export="ALL,CONFIG=${CONFIG},RUN_NAME=${RUN_NAME}")
fi

sbatch "${SBATCH_ARGS[@]}" "$@" <<'EOF'
#!/bin/bash
#SBATCH -N 1

source /etc/profile.d/modules.sh
module load nvhpc

# nvhpc only puts CUDA/NCCL/NVSHMEM on LD_LIBRARY_PATH — not HPC-X's OpenMPI
# tree. Export it ourselves so libmpi.so's dlopen of its UCX/PMIx/UCC
# neighbours resolves. Must match the libmpi path baked into MPIPreferences.
# (UCX ships under both ucx/lib and ucx/prof/lib; include both.)
HPCX_DIR="/orcd/software/core/001/pkg/nvhpc/26.1/Linux_x86_64/26.1/comm_libs/13.1/hpcx/hpcx-2.25.1"
export LD_LIBRARY_PATH="$HPCX_DIR/ompi/lib:$HPCX_DIR/ucx/lib:$HPCX_DIR/ucx/prof/lib:$HPCX_DIR/ucc/lib:${LD_LIBRARY_PATH:-}"

# Activate HPC-X so mpirun and friends are on PATH (hpcx_load only adjusts
# PATH on this build — relocation is handled by mpirun for its children).
source "$HPCX_DIR/hpcx-init-ompi.sh"
hpcx_load

# CUDA-aware MPI knobs (HPC-X / UCX).
export OMPI_MCA_opal_cuda_support=1
export UCX_TLS=cuda_copy,cuda_ipc,rc,sm,self
export UCX_MEMTYPE_CACHE=n          # avoids known UCX+CUDA memtype-cache bug
# Disable UCX's CUDA IPC handle cache. With JULIA_CUDA_MEMORY_POOL=none, CUDA
# recycles freed virtual addresses quickly; UCX's IPC cache then serves a
# stale handle and cuIpcOpenMemHandle() fails with "resource already mapped"
# on the importing rank a few minutes into the run (cuda_ipc_cache.c:212).
export UCX_CUDA_IPC_CACHE=n

# Disable CUDA's stream-ordered memory pool. With the pool, allocations come
# from a per-stream cache whose layout depends on prior call order, so two
# runs from the same checkpoint can land kernels on differently-aligned
# memory and pick different FMA paths → last-bit drift that compounds in
# unstable regions. `none` falls back to plain cudaMalloc/Free, which gives
# bit-reproducible kernel inputs at the cost of allocator overhead.
export JULIA_CUDA_MEMORY_POOL=none

JULIA="${JULIA:-$HOME/julia-1.12.5/bin/julia}"

# ── Shared environment ────────────────────────────────────────────────
FORCING_DIR="${FORCING_DIR:-${DATA}forcing_data}"
STAGING_DIR="${STAGING_DIR:-./staged_data}"
CB="${CB:-0.01}"
CUNB="${CUNB:-}"
CF="${CF:-}"
CF0="${CF0:-}"
CFD="${CFD:-}"
CP="${CP:-}"
ICE_FW="${ICE_FW:-}"
ICE_VSF="${ICE_VSF:-false}"
ICE_MELTWATER_TB="${ICE_MELTWATER_TB:-true}"
ICE_MELT_MIX="${ICE_MELT_MIX:-false}"
ICE_MELT_K="${ICE_MELT_K:-}"
ICE_MELT_DEPTH="${ICE_MELT_DEPTH:-}"
ICE_MELT_THRESHOLD="${ICE_MELT_THRESHOLD:-}"
UNDER_ICE_NU="${UNDER_ICE_NU:-}"
UNDER_ICE_NU_DEPTH="${UNDER_ICE_NU_DEPTH:-}"
ICE_ARCH="${ICE_ARCH:-}"
ICE_ARCH_STRESS="${ICE_ARCH_STRESS:-}"
ICE_ARCH_MONTHS="${ICE_ARCH_MONTHS:-}"
IC_CONDITIONS="${IC_CONDITIONS:-default}"
BIHVISC="${BIHVISC:-}"
DZ_TOP="${DZ_TOP:-}"
CATKE_CWUSTAR="${CATKE_CWUSTAR:-}"
BACKGROUND_K="${BACKGROUND_K:-}"
BACKGROUND_NU="${BACKGROUND_NU:-3e-5}"
PVEL="${PVEL:-0.254}"
BAROTROPIC_SUBSTEPS="${BAROTROPIC_SUBSTEPS:-}"
IMEX_DRAG="${IMEX_DRAG:-true}"
DRAG_UB="${DRAG_UB:-}"
CHLOROPHYLL="${CHLOROPHYLL:-seawifs}"
BACKEND_SIZE="${BACKEND_SIZE:-}"
NCAR="${NCAR:-false}"
CORRECTED="${CORRECTED:-true}"
SNOW="${SNOW:-true}"
ICE_DYNAMICS="${ICE_DYNAMICS:-true}"
OUTPUT_DIR="${OUTPUT_DIR:-.}"

# ── Build optional kwargs strings ─────────────────────────────────────

# Per-run staging subdirectory to avoid conflicts between concurrent jobs
STAGING_KWARG=""
if [[ -n "$STAGING_DIR" ]]; then
    RUN_STAGING_DIR="${STAGING_DIR}/${RUN_NAME}"
    STAGING_KWARG="staging_dir = \"${RUN_STAGING_DIR}\","
fi

CB_KWARG=""
[[ -n "$CB" ]] && CB_KWARG="Cᵇ = ${CB},"

CUNB_KWARG=""
[[ -n "$CUNB" ]] && CUNB_KWARG="Cᵘⁿᵇ = ${CUNB},"

CF_KWARG=""
[[ -n "$CF" ]]  && CF_KWARG="Cᶠ = ${CF},"
[[ -n "$CF0" ]] && CF_KWARG="${CF_KWARG} Cᶠ⁰ = ${CF0},"
[[ -n "$CFD" ]] && CF_KWARG="${CF_KWARG} Cᶠᵟ = ${CFD},"

CP_KWARG=""
[[ -n "$CP" ]] && CP_KWARG="Cᵉc = ${CP},"

ICE_FW_KWARG=""
[[ -n "$ICE_FW" ]] && ICE_FW_KWARG="ice_freshwater_fraction = ${ICE_FW},"
[[ "$ICE_VSF" == "true" ]] && ICE_FW_KWARG="${ICE_FW_KWARG}ice_virtual_salt_flux = true,"
[[ "$ICE_MELTWATER_TB" != "true" ]] && ICE_FW_KWARG="${ICE_FW_KWARG}ice_meltwater_at_interface_temperature = false,"

ICE_MELT_KWARG=""
[[ "$ICE_MELT_MIX" == "true" ]]     && ICE_MELT_KWARG="ice_melt_mixing = true,"
[[ -n "$ICE_MELT_K" ]]              && ICE_MELT_KWARG="${ICE_MELT_KWARG}ice_melt_mixing_κ = ${ICE_MELT_K},"
[[ -n "$ICE_MELT_DEPTH" ]]          && ICE_MELT_KWARG="${ICE_MELT_KWARG}ice_melt_mixing_depth = ${ICE_MELT_DEPTH},"
[[ -n "$ICE_MELT_THRESHOLD" ]]      && ICE_MELT_KWARG="${ICE_MELT_KWARG}ice_melt_mixing_threshold = ${ICE_MELT_THRESHOLD},"

UNDER_ICE_NU_KWARG=""
[[ -n "$UNDER_ICE_NU" ]]            && UNDER_ICE_NU_KWARG="under_ice_viscosity = ${UNDER_ICE_NU},"
[[ -n "$UNDER_ICE_NU_DEPTH" ]]      && UNDER_ICE_NU_KWARG="${UNDER_ICE_NU_KWARG}under_ice_viscosity_depth = ${UNDER_ICE_NU_DEPTH},"

ICE_ARCH_KWARG=""
case "$ICE_ARCH" in
  "")     ;;
  nares)  ICE_ARCH_KWARG="ice_arch_region = (-78, -58, 77.5, 82.5)," ; ICE_ARCH_MONTHS_DEFAULT="(12, 7)" ;;
  davis)  ICE_ARCH_KWARG="ice_arch_region = (-80, -48, 65, 80),"     ; ICE_ARCH_MONTHS_DEFAULT="(1, 12)" ;;
  *)      ICE_ARCH_KWARG="ice_arch_region = (${ICE_ARCH}),"          ; ICE_ARCH_MONTHS_DEFAULT="(12, 7)" ;;
esac
if [[ -n "$ICE_ARCH_KWARG" ]]; then
  if [[ -n "$ICE_ARCH_MONTHS" ]]; then
    ICE_ARCH_KWARG="${ICE_ARCH_KWARG}ice_arch_months = (${ICE_ARCH_MONTHS}),"
  else
    ICE_ARCH_KWARG="${ICE_ARCH_KWARG}ice_arch_months = ${ICE_ARCH_MONTHS_DEFAULT},"
  fi
  [[ -n "$ICE_ARCH_STRESS" ]] && ICE_ARCH_KWARG="${ICE_ARCH_KWARG}ice_arch_stress = ${ICE_ARCH_STRESS},"
fi

BIHVISC_KWARG=""
[[ -n "$BIHVISC" ]] && BIHVISC_KWARG="biharmonic_viscosity = ${BIHVISC},"

DZ_TOP_KWARG=""
[[ -n "$DZ_TOP" ]] && DZ_TOP_KWARG="Δz_top = ${DZ_TOP},"

CATKE_CWUSTAR_KWARG=""
[[ -n "$CATKE_CWUSTAR" ]] && CATKE_CWUSTAR_KWARG="Cᵂu★ = ${CATKE_CWUSTAR},"

# A named profile is passed as a Julia Symbol, a number verbatim.
BACKGROUND_K_KWARG=""
case "$BACKGROUND_K" in
    "")                     ;;
    henyey|bryan_lewis|abyssal_henyey)
                            BACKGROUND_K_KWARG="background_vertical_diffusivity = :${BACKGROUND_K}," ;;
    *)                      BACKGROUND_K_KWARG="background_vertical_diffusivity = ${BACKGROUND_K}," ;;
esac

BACKGROUND_NU_KWARG=""
[[ -n "$BACKGROUND_NU" ]] && BACKGROUND_NU_KWARG="background_vertical_viscosity = ${BACKGROUND_NU},"

# The split-explicit free surface needs the barotropic gravity wave to stay inside a substep. The
# per-config default is sized for the config default DT; raising DT needs proportionally more.
IMEX_DRAG_KWARG=""
[[ "$IMEX_DRAG" == "false" ]] && IMEX_DRAG_KWARG="implicit_bottom_drag = false,"

DRAG_UB_KWARG=""
[[ -n "$DRAG_UB" ]] && DRAG_UB_KWARG="bottom_drag_background_velocity = $DRAG_UB,"

BAROTROPIC_SUBSTEPS_KWARG=""
[[ -n "$BAROTROPIC_SUBSTEPS" ]] && BAROTROPIC_SUBSTEPS_KWARG="barotropic_substeps = ${BAROTROPIC_SUBSTEPS},"

CHLOROPHYLL_KWARG="chlorophyll = :seawifs,"
[[ "$CHLOROPHYLL" != "seawifs" ]] && CHLOROPHYLL_KWARG="chlorophyll = ${CHLOROPHYLL},"
[[ "$CHLOROPHYLL" == "none" ]]    && CHLOROPHYLL_KWARG="chlorophyll = :none,"

# Pass the value explicitly (default true = conservative restoring) so the Julia-side default
# never silently overrides a "false" request.
NORMALIZE_SALINITY="${NORMALIZE_SALINITY:-true}"
case "$NORMALIZE_SALINITY" in
    true|false) ;;
    *) echo "NORMALIZE_SALINITY must be 'true' or 'false', got '$NORMALIZE_SALINITY'" >&2; exit 1 ;;
esac
NORMALIZE_SALINITY_KWARG="normalize_salinity = ${NORMALIZE_SALINITY},"

RESTORING_UNDER_ICE="${RESTORING_UNDER_ICE:-true}"
case "$RESTORING_UNDER_ICE" in
    true|false) ;;
    *) echo "RESTORING_UNDER_ICE must be 'true' or 'false', got '$RESTORING_UNDER_ICE'" >&2; exit 1 ;;
esac
RESTORING_UNDER_ICE_KWARG=""
[[ "$RESTORING_UNDER_ICE" == "false" ]] && RESTORING_UNDER_ICE_KWARG="restoring_under_sea_ice = false,"

NORMALIZE_FRESHWATER="${NORMALIZE_FRESHWATER:-timestep}"
case "$NORMALIZE_FRESHWATER" in
    none|false)    NORMALIZE_FRESHWATER_JULIA=":none" ;;
    timestep|true) NORMALIZE_FRESHWATER_JULIA=":timestep" ;;
    annual)        NORMALIZE_FRESHWATER_JULIA=":annual" ;;
    *) echo "NORMALIZE_FRESHWATER must be none|timestep|annual, got '$NORMALIZE_FRESHWATER'" >&2; exit 1 ;;
esac
NORMALIZE_FRESHWATER_KWARG="normalize_freshwater = ${NORMALIZE_FRESHWATER_JULIA},"

RIVER_KWARG=""
if [[ -n "${RIVER_SPREAD:-}" ]]; then
    RIVER_SPREAD_JULIA="$RIVER_SPREAD"
    [[ "$RIVER_SPREAD" == "none" ]] && RIVER_SPREAD_JULIA="nothing"
    RIVER_KWARG="${RIVER_KWARG}river_spread_radius = ${RIVER_SPREAD_JULIA},"
fi
[[ -n "${RIVER_SPREAD_CELLS:-}" ]]  && RIVER_KWARG="${RIVER_KWARG}river_spread_cells = ${RIVER_SPREAD_CELLS},"
[[ -n "${RIVER_MIXING_K:-}" ]]      && RIVER_KWARG="${RIVER_KWARG}river_mixing_κ = ${RIVER_MIXING_K},"
[[ -n "${RIVER_MIXING_DEPTH:-}" ]]  && RIVER_KWARG="${RIVER_KWARG}river_mixing_depth = ${RIVER_MIXING_DEPTH},"
[[ "${RIVER_MIXING:-true}" == "false" ]] && RIVER_KWARG="${RIVER_KWARG}river_mixing = false,"
[[ -n "${RIVER_DIVERSION:-}" ]]     && RIVER_KWARG="${RIVER_KWARG}atlantic_runoff_diversion = ${RIVER_DIVERSION},"

SKEW_FORMULATION="${SKEW_FORMULATION:-diffusive}"
case "$SKEW_FORMULATION" in
    diffusive|advective|boundary_value) ;;
    *) echo "SKEW_FORMULATION must be diffusive|advective|boundary_value, got '$SKEW_FORMULATION'" >&2; exit 1 ;;
esac
SKEW_FORMULATION_KWARG="skew_flux_formulation = :${SKEW_FORMULATION},"

BVP_KWARG=""
[[ -n "${BVP_MODE:-}" ]] && BVP_KWARG="${BVP_KWARG}boundary_value_mode_number = ${BVP_MODE},"
[[ -n "${BVP_CMIN:-}" ]] && BVP_KWARG="${BVP_KWARG}boundary_value_minimum_speed = ${BVP_CMIN},"

BACKEND_KWARG=""
[[ -n "$BACKEND_SIZE" ]] && BACKEND_KWARG="backend_size = ${BACKEND_SIZE},"

FLUX_KWARG=""
[[ "$NCAR" == "true" ]]        && FLUX_KWARG="flux_configuration = :ncar,"
[[ "$CORRECTED" == "true" ]]   && FLUX_KWARG="flux_configuration = :corrected,"

PVELKWARG=""
[[ -n "$PVEL" ]] && PVELKWARG="piston_velocity = ${PVEL},"
IC_BLEND="${IC_BLEND:-}"
IC_BLEND_KWARG=""
[[ -n "$IC_BLEND" ]] && IC_BLEND_KWARG="initial_condition_blend_depth = ${IC_BLEND},"

CLOSURE_KWARG=""
[[ "${CLOSURE:-catke}" == "simple"   ]] && CLOSURE_KWARG="vertical_closure = :simple,"
[[ "${CLOSURE:-catke}" == "nori"     ]] && CLOSURE_KWARG="vertical_closure = :nori,"
[[ "${CLOSURE:-catke}" == "rbvd"     ]] && CLOSURE_KWARG="vertical_closure = :rbvd,"
[[ "${CLOSURE:-catke}" == "kpp"      ]] && CLOSURE_KWARG="vertical_closure = :kpp,"
[[ "${CLOSURE:-catke}" == "nemo_tke" ]] && CLOSURE_KWARG="vertical_closure = :nemo_tke,"

VELOCITY_KWARG=""
[[ "${WIND_VELOCITY:-false}" == "true" ]] && VELOCITY_KWARG="velocity_formulation = :wind,"

ML_TAPER_KWARG=""
[[ "${ML_TAPER:-false}" == "true" ]] && ML_TAPER_KWARG="mixed_layer_tapering = true,"

PARTIAL_CELLS_KWARG=""
[[ "${PARTIAL_CELLS:-false}" == "true" ]] && PARTIAL_CELLS_KWARG="partial_cell_bathymetry = true,"

BBL_KWARG=""
[[ -n "${BBL_KAPPA:-}" ]] && BBL_KWARG="bbl_diffusivity = ${BBL_KAPPA},"
[[ -n "${BBL_GAMMA:-}" ]] && BBL_KWARG="${BBL_KWARG}bbl_transport_coefficient = ${BBL_GAMMA},"
[[ -n "${OVERFLOW_RESTORE:-}" ]] && BBL_KWARG="${BBL_KWARG}overflow_restoring_timescale = ${OVERFLOW_RESTORE}days,"
[[ -n "${LAB_RESTORE:-}" ]] && BBL_KWARG="${BBL_KWARG}labrador_restoring_timescale = ${LAB_RESTORE}days,"

SNOW_KWARG=""
[[ "$SNOW" == "true" ]] && SNOW_KWARG="with_snow = true,"

ICE_DYNAMICS_KWARG=""
[[ "$ICE_DYNAMICS" == "false" ]] && ICE_DYNAMICS_KWARG="with_ice_dynamics = false,"

ICE_LATERAL="${ICE_LATERAL:-no_slip}"
ICE_BASAL="${ICE_BASAL:-true}"
ICE_DRAG="${ICE_DRAG:-5.5e-3}"
ICE_HEAT_TRANSFER="${ICE_HEAT_TRANSFER:-0.0057}"
ICE_PSTAR="${ICE_PSTAR:-}"
ICE_CATEGORIES="${ICE_CATEGORIES:-4}"
SEA_ICE_KWARG="sea_ice_lateral_boundary_condition = :${ICE_LATERAL},"
SEA_ICE_KWARG="${SEA_ICE_KWARG}sea_ice_ocean_drag_coefficient = ${ICE_DRAG},"
SEA_ICE_KWARG="${SEA_ICE_KWARG}sea_ice_ocean_heat_transfer_coefficient = ${ICE_HEAT_TRANSFER},"
[[ "$ICE_BASAL" == "false" ]] && SEA_ICE_KWARG="${SEA_ICE_KWARG}with_landfast_basal_stress = false,"
SEA_ICE_KWARG="${SEA_ICE_KWARG}thickness_categories = ${ICE_CATEGORIES},"
[[ -n "${SNOW_CATEGORIES:-}" ]] && SEA_ICE_KWARG="${SEA_ICE_KWARG}snow_thickness_categories = ${SNOW_CATEGORIES},"
if [[ -n "${ICE_ITD_SHAPE:-}" ]]; then
  IFS=',' read -r ITD_SMIN ITD_SMAX ITD_HSTAR <<< "$ICE_ITD_SHAPE"
  if [[ -z "$ITD_SMIN" || -z "$ITD_SMAX" || -z "$ITD_HSTAR" ]]; then
    echo "ICE_ITD_SHAPE must be <minimum_shape>,<maximum_shape>,<transition_thickness>, got '$ICE_ITD_SHAPE'" >&2
    exit 1
  fi
  SEA_ICE_KWARG="${SEA_ICE_KWARG}itd_shape = ThicknessDependentConductivity(minimum_shape = ${ITD_SMIN}, maximum_shape = ${ITD_SMAX}, transition_thickness = ${ITD_HSTAR}),"
fi
[[ -n "$ICE_PSTAR" ]] && SEA_ICE_KWARG="${SEA_ICE_KWARG}ice_compressive_strength = ${ICE_PSTAR},"
ICE_SALINITY="${ICE_SALINITY:-}"
[[ -n "$ICE_SALINITY" ]] && SEA_ICE_KWARG="${SEA_ICE_KWARG}ice_salinity = ${ICE_SALINITY},"
ISOPYCNAL="${ISOPYCNAL:-standard}"
case "$ISOPYCNAL" in
    standard|triad) ;;
    *) echo "ISOPYCNAL must be 'standard' or 'triad', got '$ISOPYCNAL'" >&2; exit 1 ;;
esac
ISOPYCNAL_KWARG=""
[[ "$ISOPYCNAL" == "triad" ]] && ISOPYCNAL_KWARG="isopycnal_formulation = :triad,"

ADVECTION_KWARG=""
[[ "$TRACER_ORDER" != "7" ]] && ADVECTION_KWARG="${ADVECTION_KWARG}tracer_advection_order = ${TRACER_ORDER},"
[[ "$TRACER_BOUNDARY_SCHEME" != "default" ]]   && ADVECTION_KWARG="${ADVECTION_KWARG}tracer_boundary_scheme = :${TRACER_BOUNDARY_SCHEME},"
[[ "$MOMENTUM_BOUNDARY_SCHEME" != "default" ]] && ADVECTION_KWARG="${ADVECTION_KWARG}momentum_boundary_scheme = :${MOMENTUM_BOUNDARY_SCHEME},"
[[ -n "${TEMPERATURE_VARIATION:-}" ]]           && ADVECTION_KWARG="${ADVECTION_KWARG}temperature_reference_variation = ${TEMPERATURE_VARIATION},"
[[ -n "${SALINITY_VARIATION:-}" ]]              && ADVECTION_KWARG="${ADVECTION_KWARG}salinity_reference_variation = ${SALINITY_VARIATION},"
[[ -n "${MOMENTUM_VARIATION:-}" ]]              && ADVECTION_KWARG="${ADVECTION_KWARG}momentum_reference_variation = ${MOMENTUM_VARIATION},"
ICE_LIQUIDUS="${ICE_LIQUIDUS:-teos10}"
[[ "$ICE_LIQUIDUS" != "teos10" ]] && SEA_ICE_KWARG="${SEA_ICE_KWARG}sea_ice_liquidus = :${ICE_LIQUIDUS},"
ICE_DRAGREF="${ICE_DRAGREF:-6}"
if [[ "$ICE_DRAGREF" == "none" ]]; then
  SEA_ICE_KWARG="${SEA_ICE_KWARG}sea_ice_ocean_drag_reference_depth = nothing,"
else
  SEA_ICE_KWARG="${SEA_ICE_KWARG}sea_ice_ocean_drag_reference_depth = ${ICE_DRAGREF},"
fi
ICE_Z0="${ICE_Z0:-5e-4}"
[[ "$ICE_Z0" != "5e-4" ]] && SEA_ICE_KWARG="${SEA_ICE_KWARG}sea_ice_momentum_roughness_length = ${ICE_Z0},"
ICE_TILT="${ICE_TILT:-false}"
[[ "$ICE_TILT" == "true" ]] && SEA_ICE_KWARG="${SEA_ICE_KWARG}with_ocean_surface_tilt = true,"
[[ "$IC_CONDITIONS" != "default" ]] && SEA_ICE_KWARG="${SEA_ICE_KWARG}northern_sea_ice_initial_date = DateTime(1993, 9, 1),"

# Profile runs disable the OMIP diagnostic output writers (Average,
# JLD2 dumps, checkpoint, KE spectrum). They add per-iteration I/O and
# `compute!` overhead that pollutes nsys traces. The simulation core
# is unchanged, so kernel timings reflect pure stepping cost.
DIAGNOSTICS_KWARG=""
[[ "${PROFILE:-false}" == "true" ]] && DIAGNOSTICS_KWARG="diagnostics = false,"

# ── Build and run Julia expression ────────────────────────────────────
JULIA_EXPR="using ClimaOcean
using Oceananigans
using Oceananigans.Units
using CUDA
using Dates
${EXTRA_USING}

sim = omip_simulation(:${CONFIG};
                      arch = ${ARCH},
                      Nz = ${NZ},
                      depth = 5500,
                      ${DZ_TOP_KWARG}
                      κ_skew = ${KSKEW_JULIA},
                      κ_symmetric = ${KSYMM_JULIA},
                      biharmonic_timescale = ${BIHARMONIC},
                      ${BIHVISC_KWARG}
                      ${CB_KWARG}
                      ${CUNB_KWARG}
                      ${CF_KWARG}
                      ${CP_KWARG}
                      ${ICE_FW_KWARG}
                      ${ICE_MELT_KWARG}
                      ${UNDER_ICE_NU_KWARG}
                      ${ICE_ARCH_KWARG}
                      ${FLUX_KWARG}
                      ${CLOSURE_KWARG}
                      ${VELOCITY_KWARG}
                      ${ML_TAPER_KWARG}
                      ${PARTIAL_CELLS_KWARG}
                      ${BBL_KWARG}
                      ${SNOW_KWARG}
                      ${ICE_DYNAMICS_KWARG}
                      ${SEA_ICE_KWARG}
                      ${DIAGNOSTICS_KWARG}
                      ${PVELKWARG}
                      ${IC_BLEND_KWARG}
                      ${CATKE_CWUSTAR_KWARG}
                      ${BACKGROUND_K_KWARG}
                      ${BACKGROUND_NU_KWARG}
                      ${BAROTROPIC_SUBSTEPS_KWARG}
                      ${IMEX_DRAG_KWARG}
                      ${DRAG_UB_KWARG}
                      ${CHLOROPHYLL_KWARG}
                      ${NORMALIZE_SALINITY_KWARG}
                      ${RESTORING_UNDER_ICE_KWARG}
                      ${NORMALIZE_FRESHWATER_KWARG}
                      ${RIVER_KWARG}
                      ${ADVECTION_KWARG}
                      ${SKEW_FORMULATION_KWARG}
                      ${ISOPYCNAL_KWARG}
                      ${BVP_KWARG}
                      Δt = ${DT},
                      forcing_dir = \"${FORCING_DIR}\",
                      ${STAGING_KWARG}
                      ${BACKEND_KWARG}
                      ${FILE_SPLIT}
                      output_dir = \"${OUTPUT_DIR}/${RUN_NAME}_run\",
                      filename_prefix = \"${RUN_NAME}\")

${RUN_CMD}"

THREADS="${THREADS:-8}"

# Launch via HPC-X's mpirun, NOT srun. SLURM 25.05 on this cluster has only
# pmi2/cray_shasta MPI plugins (no PMIx), and HPC-X's OpenMPI 4.1.9a1 ships
# only the pmix3x client (no native PMI-1/PMI-2 wire support). srun and HPC-X
# therefore cannot handshake. mpirun runs its own PMIx-3 server, queries SLURM
# only for the allocation (hostlist + nranks), and spawns ranks itself.
# --bind-to none avoids fighting SLURM's cgroup over CPU affinity.
NRANKS=$((SLURM_NNODES * SLURM_NTASKS_PER_NODE))
if [[ "${PROFILE:-false}" == "true" ]]; then
    echo "Profiling ${RUN_NAME} -> ${REPORT_NAME}"
    mpirun -np "$NRANKS" --bind-to none \
                      nsys profile --trace=cuda \
                      --output="$REPORT_NAME" \
                      --force-overwrite true \
                      "$JULIA" --project=.. --check-bounds=no -t "${THREADS}" -e "$JULIA_EXPR"
else
    mpirun -np "$NRANKS" --bind-to none \
        "$JULIA" --project=.. --check-bounds=no -t "${THREADS}" -e "$JULIA_EXPR"
fi
EOF
