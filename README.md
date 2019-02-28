# boggy-retina

## Description

This repository is a place to collaborate and build out simulations to describe unwanted spread of activation observed in the retina.

## Quick Start

The generate a plot of the spread of extracellular current in a normal retina under stimulation with one electrode, run the following from the [Simulations directory](./Simulations/):

```
AnyElecConfig_RotatedNeurite_4L_Depth
BuildFigures
```

To see the effect of different NFL thicknesses, change the line `simulation = 17;` in `AnyElecConfig_RotatedNeurite_4L_Depth` to `simulation = 18;` or `simulation = 19;`.

In order to modify the conductivity of the NFL, modify parameters `b` and `d` in [Utilities/NTESparams.m](Utilities/NTESparams.m) as demonstrated in the comments of that file. Changing these parameters alters how much _space_ there is between fibers in the NFL.

To modify the conductivity of the GCL, change the conductivity value defined by `sigma_L = 0.1;` in [CellComp4Layer_Ve_f_Plane.m](CellComp4Layer_Ve_f_Plane.m)
