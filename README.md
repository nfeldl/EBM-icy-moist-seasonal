# EBM-icy-moist-seasonal
Code and data accompanying the paper: Feldl, N., & Merlis, T. M. (2021). Polar amplification in idealized climates: The role of ice, moisture, and seasons. Geophysical Research Letters, 48, e2021GL094130. https://doi.org/10.1029/2021GL094130

## EBM

The EBM can either be executed as a script or imported as a module.

To execute as a script:
```
python sea_ice_MEBM_ghost.py
```

To import as a module:
```
import sea_ice_MEBM_ghost as ebm
ebm.model(grid, Ti, F, moist = moist, albT = albT, seas = seas, thermo = thermo)
```

## GCM

The netcdf file `GCM_results.nc` contains the results of the simulations presented in the paper. The tar file contains the source code for the idealized GCM with thermodynamic ice described in the paper.
