# Atlantis Subsidence (Atlans)

Atlans is package to calculate 1D soil subsidence using voxel based geological and lithological models, such as [GeoTOP](https://basisregistratieondergrond.nl/inhoud-bro/registratieobjecten/modellen/geotop-gtm/) ([Stafleu et al., 2011](https://doi.org/10.1017/S0016774600000597)). The amount of subsidence is calculated over stressperiods which are typically related to groundwater management cycles. The results for every stressperiod are stored as Netcdf output.

## Processes which contribute to subsidence
Currently, Atlans calculates the amount of subsidence based on the sum of three contributing processes: **consolidation**, **oxidation** and **shrinkage** and stores the calculated subsidence and the contributions of each process for every stressperiod in the output Netcdf. A user can choose to ignore individual processes (for example, do not include shrinkage in the calculations). See [Atlans docs]() for a detailed explanation of the usage and the calculation methods for every process.

