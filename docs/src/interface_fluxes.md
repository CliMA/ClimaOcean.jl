# Interface fluxes

`ClimaOcean`'s `OceanSeaIceModel` has essentially two goals:
    1. Manage time-stepping multiple component models forward simulatneously,
    2. Compute and/or pass fluxes between the component models.

This page describes how OceanSeaIceModel computes fluxes at the interfaces between
component models.


## Component models

`OceanSeaIceModel`