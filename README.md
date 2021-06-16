# EBM-icy-moist-seasonal
Code and data accompanying the paper "Polar amplification in idealized climates: the role of ice, moisture, and seasons" by Feldl and Merlis.

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
