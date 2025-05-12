# Atlantis Subsidence (Atlans.jl)

[![CI](https://github.com/Deltares-research/Atlans.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/Deltares-research/Atlans.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/Deltares-research/Atlans.jl/graph/badge.svg)](https://codecov.io/gh/Deltares-research/Atlans.jl)

Atlans is package to calculate 1D soil subsidence using voxel based geological and lithological models, such as [GeoTOP](https://basisregistratieondergrond.nl/inhoud-bro/registratieobjecten/modellen/geotop-gtm/) ([Stafleu et al., 2011](https://doi.org/10.1017/S0016774600000597)). The amount of subsidence is calculated over stressperiods which are typically related to groundwater management cycles.

## Processes which contribute to subsidence
Currently, Atlans calculates the amount of subsidence based on the sum of three contributing processes:

1. **consolidation**
2. **oxidation**
3. **shrinkage** 

The calculated subsidence and the contributions of each process for every stressperiod are stored as Netcdf output. A user can choose to ignore individual processes (for example, do not include shrinkage in the calculations). 

See [Atlans docs](https://deltares-research.github.io/Atlans.jl/)(in development) for a detailed explanation of the usage and the calculation methods for every process.

## Installing as Julia package
In the Julia REPL:
```julia-repl
julia> using Pkg; Pkg.add("Atlans")
```

or, enter the package manager by pressing `]` and type:
```julia-repl
pkg> add Atlans
```

To install the latest development version from the GitHub repository directly, enter the package manager and use:
```julia-repl
pkg> add https://github.com/Deltares-research/Atlans.jl.git
```
