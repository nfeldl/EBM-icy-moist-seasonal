# EBM-icy-moist-seasonal
Code and data accompanying the paper "Polar amplification in idealized climates: the role of ice, moisture, and seasons" by Feldl and Merlis.

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

The netcdf file `GCM_results.nc` contains the results of the simulations presented in the paper. The tar file contains the source code for the idealized GCM with thermodynamic ice used in this study and described in the paper.
