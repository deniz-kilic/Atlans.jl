# Atlantis Subsidence (Atlans)

Atlans is package to calculate 1D soil subsidence using voxel based geological and lithological models, such as [GeoTOP](https://basisregistratieondergrond.nl/inhoud-bro/registratieobjecten/modellen/geotop-gtm/) ([Stafleu et al., 2011](https://doi.org/10.1017/S0016774600000597)). The amount of subsidence is calculated over stressperiods which are typically related to groundwater management cycles.

## Processes which contribute to subsidence
Currently, Atlans calculates the amount of subsidence based on the sum of three contributing processes:

1. **consolidation**
2. **oxidation**
3. **shrinkage** 

The calculated subsidence and the contributions of each process for every stressperiod are stored as Netcdf output. A user can choose to ignore individual processes (for example, do not include shrinkage in the calculations). 

See [Atlans docs]()(in development) for a detailed explanation of the usage and the calculation methods for every process.
