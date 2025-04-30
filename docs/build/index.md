
<a id='API-Reference'></a>

<a id='API-Reference-1'></a>

# API Reference


*This is the private internal documentation of the Atlans API.*

- [API Reference](index.md#API-Reference)
    - [Modules](index.md#Modules)
    - [Types](index.md#Types)
    - [Functions](index.md#Functions)
    - [Constants](index.md#Constants)
    - [Macros](index.md#Macros)
    - [Index](index.md#Index)


<a id='Modules'></a>

<a id='Modules-1'></a>

## Modules


<a id='Types'></a>

<a id='Types-1'></a>

## Types

<a id='Atlans.AdaptiveCellsize' href='#Atlans.AdaptiveCellsize'>#</a>
**`Atlans.AdaptiveCellsize`** &mdash; *Type*.



```julia
AdaptiveCellsize(Δzmax::Float, split_tolerance::Float)
```

Logic for splitting cells in a column to accomodate for a moving phreatic level in combination with organic matter stores. Handles how the thickness of thick voxels  (>Δzmax) should be discretized and determines when splitting occurs. If the thickness of a cell above, or below, the groundwater table is lower than the tolerance, no splitting occurs.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.Clock' href='#Atlans.Clock'>#</a>
**`Atlans.Clock`** &mdash; *Type*.



```julia
Clock(time::Vector{DateTime}, iteration::int, stop_time::DateTime)
```

Object to keep track of the stress periods, number of iterations and stop time of an Atlantis Simulation.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.ExponentialTimeStepper' href='#Atlans.ExponentialTimeStepper'>#</a>
**`Atlans.ExponentialTimeStepper`** &mdash; *Type*.



```julia
ExponentialTimestepper(start::Float, multiplier::T)
```

Struct to discretize time steps (in days) within each stress period.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.Model-Tuple{Type, Type, Type, Type, Type, Vararg{Any, 4}}' href='#Atlans.Model-Tuple{Type, Type, Type, Type, Type, Vararg{Any, 4}}'>#</a>
**`Atlans.Model`** &mdash; *Method*.



```julia
Model(
	groundwater::Type,
	consolidation::Type,
	oxidation::Type,
	preconsolidation::Type,
	shrinkage::Type,
	adaptive_cellsize,
	timestepper,
	path_subsoil,
	path_lookup
)
```

Initialize a model with specified groundwater, consolidation, oxidation and shrinkage processes from a netCDF file and CSV lookup table describing the subsurface parameters, appropriate for the chosen processes.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.ShrinkageColumn' href='#Atlans.ShrinkageColumn'>#</a>
**`Atlans.ShrinkageColumn`** &mdash; *Type*.



```julia
ShrinkageColumn{S}(cells, z, Δz, result, Hv0)
```

Collection of SimpleShrinkage cells to compute shrinkage for. 

**Arguments:**

  * `cells::Vector{S}`: Collection of cells containing the shrinkage process.
  * `z::Vector{Float}`: Depth of the cells.
  * `Δz::Vector{Float}`: Thickness of the cells.
  * `result::Vector{Float}`: Computed shrinkage of each cell.
  * `Hv0::Float`: Absolute depth above phreatic level to compute shrinkage for in cells.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.SimpleShrinkage' href='#Atlans.SimpleShrinkage'>#</a>
**`Atlans.SimpleShrinkage`** &mdash; *Type*.



```julia
SimpleShrinkage(Δz, n, τ, r, shrinkage)
```

Simple voxel with attributes to compute shrinkage for.

#Arguments

  * `Δz::Float`: Thickness of the voxel. [m]
  * `n::Float`: Shrinkage factor of the voxel. [-]
  * `L::Float`: Mass fraction of lutum.
  * `H::Float`: Mass fraction of organic.
  * `τ::Float`: Time dependent factor for shrinkage process. [days]
  * `r::Float`: Direction of shrinkage, r is 3 indicates isoptropic. [-]
  * `sf::Float`: TODO: look-up in document [-]
  * `shrinkage::Float`: Computed shrinkage or elevation change over time. [m]


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.Simulation-Tuple{Atlans.Model, String, Dates.DateTime}' href='#Atlans.Simulation-Tuple{Atlans.Model, String, Dates.DateTime}'>#</a>
**`Atlans.Simulation`** &mdash; *Method*.



```julia
Simulation(model, path_output, stop_time, forcings, additional_times)
```

Setup a simulation from an initialized model.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.SoilColumn' href='#Atlans.SoilColumn'>#</a>
**`Atlans.SoilColumn`** &mdash; *Type*.



x, y, z are all midpoints.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.VerticalDomain' href='#Atlans.VerticalDomain'>#</a>
**`Atlans.VerticalDomain`** &mdash; *Type*.



Temporary structure used to create SoilColumns.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>


<a id='Functions'></a>

<a id='Functions-1'></a>

## Functions

<a id='Atlans.Qcreep_derivative-Tuple{Atlans.AbcIsotache, Float64, Float64}' href='#Atlans.Qcreep_derivative-Tuple{Atlans.AbcIsotache, Float64, Float64}'>#</a>
**`Atlans.Qcreep_derivative`** &mdash; *Method*.



Derivative of Qcreep with respect to head.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.U-Tuple{Atlans.ConsolidationProcess, Float64}' href='#Atlans.U-Tuple{Atlans.ConsolidationProcess, Float64}'>#</a>
**`Atlans.U`** &mdash; *Method*.



Terzaghi, degree of consolidation


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.add_time-Tuple{Any, Any}' href='#Atlans.add_time-Tuple{Any, Any}'>#</a>
**`Atlans.add_time`** &mdash; *Method*.



Add a new time to the unlimited time dimension, and return the index


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.advance!-Tuple{Any}' href='#Atlans.advance!-Tuple{Any}'>#</a>
**`Atlans.advance!`** &mdash; *Method*.



Advances the clock by one iteration.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.advance_forcingperiod!-Tuple{Any, Any}' href='#Atlans.advance_forcingperiod!-Tuple{Any, Any}'>#</a>
**`Atlans.advance_forcingperiod!`** &mdash; *Method*.



Advance a single stress period for all columns.

Timesteps are determined by the total duration and the chosen timestepper.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.advance_forcingperiod!-Tuple{Any}' href='#Atlans.advance_forcingperiod!-Tuple{Any}'>#</a>
**`Atlans.advance_forcingperiod!`** &mdash; *Method*.



Advance the simulation by a single forcing period. Reads new input, computes, writes output.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.advance_forcingperiod!-Tuple{Atlans.SoilColumn, Vector{Float64}}' href='#Atlans.advance_forcingperiod!-Tuple{Atlans.SoilColumn, Vector{Float64}}'>#</a>
**`Atlans.advance_forcingperiod!`** &mdash; *Method*.



```julia
advance_forcingperiod!(column, timesteps)
```

Advances a prepared column by a number of timesteps.

Note, the correct order of execution is:

  * prepare a forcing period: compute pre-load stress
  * apply forcing: change load
  * advance forcing period


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.advance_timestep!-Tuple{Atlans.SoilColumn, Float64}' href='#Atlans.advance_timestep!-Tuple{Atlans.SoilColumn, Float64}'>#</a>
**`Atlans.advance_timestep!`** &mdash; *Method*.



```julia
advance_timestep!(column, Δt)
```

Advance a single timestep.

During a timestep the following states are computed:

```
* head
* pore pressure
* total stress
* effective stress
```

Then, consolidation and oxidation are computed.

Finally, thickness and elevation are updated.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.cellsplit!-Tuple{Atlans.OxidationColumn{Atlans.CarbonStore}, Vararg{Any, 4}}' href='#Atlans.cellsplit!-Tuple{Atlans.OxidationColumn{Atlans.CarbonStore}, Vararg{Any, 4}}'>#</a>
**`Atlans.cellsplit!`** &mdash; *Method*.



For CarbonStore, organic and mineral mass should be split according to cell height (Δz).


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.cellsplit!-Tuple{Union{Atlans.ConsolidationColumn{Atlans.NullConsolidation, Atlans.OverConsolidationRatio}, Atlans.OxidationColumn{Atlans.NullOxidation}, Atlans.ShrinkageColumn{Atlans.NullShrinkage}}, Vararg{Any, 4}}' href='#Atlans.cellsplit!-Tuple{Union{Atlans.ConsolidationColumn{Atlans.NullConsolidation, Atlans.OverConsolidationRatio}, Atlans.OxidationColumn{Atlans.NullOxidation}, Atlans.ShrinkageColumn{Atlans.NullShrinkage}}, Vararg{Any, 4}}'>#</a>
**`Atlans.cellsplit!`** &mdash; *Method*.



```julia
cellsplit!(column, _, newlength, _, _,)
```

Logic for cellsplit! if one of the processes is ignored (e.g. NullConsolidation).


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.compress_γ_dry-Tuple{Atlans.ConsolidationProcess, Float64}' href='#Atlans.compress_γ_dry-Tuple{Atlans.ConsolidationProcess, Float64}'>#</a>
**`Atlans.compress_γ_dry`** &mdash; *Method*.



Consolidation reduces pore space, pushes out the air.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.compress_γ_wet-Tuple{Atlans.ConsolidationProcess, Float64}' href='#Atlans.compress_γ_wet-Tuple{Atlans.ConsolidationProcess, Float64}'>#</a>
**`Atlans.compress_γ_wet`** &mdash; *Method*.



Consolidation reduces pore space, pushes out the water.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.consolidate-Tuple{Atlans.DrainingAbcIsotache, Any, Any}' href='#Atlans.consolidate-Tuple{Atlans.DrainingAbcIsotache, Any, Any}'>#</a>
**`Atlans.consolidate`** &mdash; *Method*.



```julia
consolidate(abc::DrainingAbcIsotache, σ′, Δt)
```

Compute consolidation for a single cell.

The cell contains a state U for Terzaghi's degree of consolidation. This state is updated every consolidate step. During every Δt, a new U is computed.  The increase in U can be directly related to the pore pressure, and so the increase in effective stress is equal to the effective stress prior to loading (abc.σ′) and the new effective stress (σ′) reached when U becomes 1.0.

The degree of consolidation plays only one role: it distributes the load (final σ′ - initial σ′) over time and as the column might be submerging, the increase in σ′ may become lower, providing a negative feedback mechanism.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.create_timesteps-Tuple{Atlans.ExponentialTimeStepper, Any}' href='#Atlans.create_timesteps-Tuple{Atlans.ExponentialTimeStepper, Any}'>#</a>
**`Atlans.create_timesteps`** &mdash; *Method*.



Based on the duration and the timestepper, create the required timesteps.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.currenttime-Tuple{Any}' href='#Atlans.currenttime-Tuple{Any}'>#</a>
**`Atlans.currenttime`** &mdash; *Method*.



Return current time.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.discretize-Tuple{Any, Float64}' href='#Atlans.discretize-Tuple{Any, Float64}'>#</a>
**`Atlans.discretize`** &mdash; *Method*.



In how many equal parts should a thick cell be divided?


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.draining_abc_isotache_column-NTuple{8, Any}' href='#Atlans.draining_abc_isotache_column-NTuple{8, Any}'>#</a>
**`Atlans.draining_abc_isotache_column`** &mdash; *Method*.



Turn a collection of vectors into a collection of DrainingAbcIsotache cells.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.effective_stress!-Tuple{Atlans.AbstractConsolidationColumn}' href='#Atlans.effective_stress!-Tuple{Atlans.AbstractConsolidationColumn}'>#</a>
**`Atlans.effective_stress!`** &mdash; *Method*.



Compute effective stress for entire column


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.formulate-Tuple{Atlans.AbcIsotache, Vararg{Float64, 4}}' href='#Atlans.formulate-Tuple{Atlans.AbcIsotache, Vararg{Float64, 4}}'>#</a>
**`Atlans.formulate`** &mdash; *Method*.



First term of Taylor expansion

f(x) ~= f(a) + f′(a) * (x - a)


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.initialize-Tuple{Type{Atlans.CarbonStore}, Any, Any, Any}' href='#Atlans.initialize-Tuple{Type{Atlans.CarbonStore}, Any, Any, Any}'>#</a>
**`Atlans.initialize`** &mdash; *Method*.



```julia
initialize(::Type{CarbonStore}, domain, subsoil, I)
```

Initialize a OxidationColumn for a domain at location I based subsurface input.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.initialize-Tuple{Type{Atlans.CarbonStore}, Atlans.VerticalDomain, Dict}' href='#Atlans.initialize-Tuple{Type{Atlans.CarbonStore}, Atlans.VerticalDomain, Dict}'>#</a>
**`Atlans.initialize`** &mdash; *Method*.



```julia
initialize(::Type{CarbonStore}, domain::VerticalDomain, lookup_table::Dict)
```

Initialize a OxidationSurcharge column that can be added to an OxidationColumn when Surcharge is applied as a forcing during a forcingperiod.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.initialize-Tuple{Type{Atlans.DrainingAbcIsotache}, Type, Any, Any, Any}' href='#Atlans.initialize-Tuple{Type{Atlans.DrainingAbcIsotache}, Type, Any, Any, Any}'>#</a>
**`Atlans.initialize`** &mdash; *Method*.



```julia
initialize(::Type{DrainingAbcIsotache}, domain, subsoil, I)
```

Initialize a ConsolidationColumn for a domain at location I based subsurface input.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.initialize-Tuple{Type{Atlans.DrainingAbcIsotache}, Type, Atlans.VerticalDomain, Dict}' href='#Atlans.initialize-Tuple{Type{Atlans.DrainingAbcIsotache}, Type, Atlans.VerticalDomain, Dict}'>#</a>
**`Atlans.initialize`** &mdash; *Method*.



```julia
initialize(
    ::Type{DrainingAbcIsotache},
    preconsolidation::Type,
    domain::VerticalDomain,
    lookup_table::Dict
)
```

Initialize a ConsolidationSurcharge column that can be added to a ConsolidationColumn when Surcharge is applied as a forcing during a forcingperiod.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.initialize-Tuple{Type{Atlans.HydrostaticGroundwater}, Atlans.Phreatic, Atlans.VerticalDomain}' href='#Atlans.initialize-Tuple{Type{Atlans.HydrostaticGroundwater}, Atlans.Phreatic, Atlans.VerticalDomain}'>#</a>
**`Atlans.initialize`** &mdash; *Method*.



```julia
initialize(::Type{HydrostaticGroundwater}, phreatic::Phreatic, domain::VerticalDomain)
```

Initialize a GroundwaterSurcharge column that can be added to a HydrostaticGroundwater column when Surcharge is applied as a forcing during a forcingperiod.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.initialize-Tuple{Type{Atlans.NullConsolidation}, Type, Any, Any}' href='#Atlans.initialize-Tuple{Type{Atlans.NullConsolidation}, Type, Any, Any}'>#</a>
**`Atlans.initialize`** &mdash; *Method*.



```julia
initialize(::Type{NullConsolidation}, preconsolidation::Type, domain, _)
```

Initialize an empty ConsolidationSurcharge column when the consolidation process is ignored.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.initialize-Tuple{Type{Atlans.NullConsolidation}, Vararg{Any, 4}}' href='#Atlans.initialize-Tuple{Type{Atlans.NullConsolidation}, Vararg{Any, 4}}'>#</a>
**`Atlans.initialize`** &mdash; *Method*.



```julia
initialize(::Type{NullConsolidation}, domain, subsoil, I)
```

Initialize an empty ConsolidationColumn (i.e. consolidation is ignored) at location I.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.initialize-Tuple{Type{Atlans.NullOxidation}, Any, Any, Any}' href='#Atlans.initialize-Tuple{Type{Atlans.NullOxidation}, Any, Any, Any}'>#</a>
**`Atlans.initialize`** &mdash; *Method*.



```julia
initialize(::Type{SimpleShrinkage}, domain, subsoil, I)
```

Initialize an empty OxidationColumn (i.e. oxidation is ignored) at location I.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.initialize-Tuple{Type{Atlans.NullOxidation}, Any, Any}' href='#Atlans.initialize-Tuple{Type{Atlans.NullOxidation}, Any, Any}'>#</a>
**`Atlans.initialize`** &mdash; *Method*.



```julia
initialize(::Type{NullOxidation}, domain, _)
```

Initialize an empty OxidationSurcharge column when the oxidation process is ignored.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.initialize-Tuple{Type{Atlans.NullShrinkage}, Any, Any, Any}' href='#Atlans.initialize-Tuple{Type{Atlans.NullShrinkage}, Any, Any, Any}'>#</a>
**`Atlans.initialize`** &mdash; *Method*.



```julia
initialize(::Type{SimpleShrinkage}, domain, subsoil, I)
```

Initialize an empty ShrinkageColumn (i.e. shrinkage is ignored) at location I.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.initialize-Tuple{Type{Atlans.NullShrinkage}, Any, Any}' href='#Atlans.initialize-Tuple{Type{Atlans.NullShrinkage}, Any, Any}'>#</a>
**`Atlans.initialize`** &mdash; *Method*.



```julia
initialize(::Type{NullShrinkage}, domain, _)
```

Initialize an empty ShrinkageSurcharge column when the shrinkage process is ignored.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.initialize-Tuple{Type{Atlans.SimpleShrinkage}, Any, Any, Any}' href='#Atlans.initialize-Tuple{Type{Atlans.SimpleShrinkage}, Any, Any, Any}'>#</a>
**`Atlans.initialize`** &mdash; *Method*.



```julia
initialize(::Type{SimpleShrinkage}, domain, subsoil, I)
```

Initialize a ShrinkageColumn for a domain at location I based subsurface input.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.initialize-Tuple{Type{Atlans.SimpleShrinkage}, Atlans.VerticalDomain, Dict}' href='#Atlans.initialize-Tuple{Type{Atlans.SimpleShrinkage}, Atlans.VerticalDomain, Dict}'>#</a>
**`Atlans.initialize`** &mdash; *Method*.



```julia
initialize(::Type{SimpleShrinkage}, domain::VerticalDomain, lookup_table::Dict)
```

Initialize a ShrinkageSurcharge column that can be added to an ShrinkageColumn when Surcharge is applied as a forcing during a forcingperiod.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.parse_loglevel-Tuple{AbstractString}' href='#Atlans.parse_loglevel-Tuple{AbstractString}'>#</a>
**`Atlans.parse_loglevel`** &mdash; *Method*.



```julia
parse_loglevel(input_level::AbstractString)::LogLevel
parse_loglevel(input_level::Integer)::LogLevel
```

Parse a log level from either an integer or string.

**Examples**

```
parse_loglevel("info") -> Logging.Info
parse_loglevel(0) -> LogLevel(0) (== Logging.Info)
```


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.periodduration-Tuple{Any}' href='#Atlans.periodduration-Tuple{Any}'>#</a>
**`Atlans.periodduration`** &mdash; *Method*.



Compute duration of forcing period from current time.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.pow-Tuple{Any, Any}' href='#Atlans.pow-Tuple{Any, Any}'>#</a>
**`Atlans.pow`** &mdash; *Method*.



Faster method for exponentiation


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.prepare_domain-NTuple{7, Any}' href='#Atlans.prepare_domain-NTuple{7, Any}'>#</a>
**`Atlans.prepare_domain`** &mdash; *Method*.



Temporary structure used to create SoilColumns.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.prepare_forcingperiod!' href='#Atlans.prepare_forcingperiod!'>#</a>
**`Atlans.prepare_forcingperiod!`** &mdash; *Function*.



```julia
prepare_forcingperiod!(column, split_tolerance)
```

Prepare a single forcing period.

This:

```
* splits the soil column at maximum oxidation depth (if necessary)
* computes pore pressure, total stress, effective stress prior to this stress periods loading
* sets the effective stress in every consolidation cell
* Reset U and t for DrainingConsolidation processes.
```

Note that splitting requires knowing where the phreatic level ends up after applying all forcings. This means that changes to phreatic level and deep subsidence must be accounted for. Furthermore, the split must be applied before applying the changes to compute the pre-loading effective stress.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.prepare_forcingperiod!-Tuple{Atlans.ConsolidationColumn{Atlans.DrainingAbcIsotache, P} where P<:Atlans.Preconsolidation}' href='#Atlans.prepare_forcingperiod!-Tuple{Atlans.ConsolidationColumn{Atlans.DrainingAbcIsotache, P} where P<:Atlans.Preconsolidation}'>#</a>
**`Atlans.prepare_forcingperiod!`** &mdash; *Method*.



Reset degree of consolidation and time.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.prepare_surcharge_column-Tuple{Atlans.Surcharge, Atlans.SoilColumn, CartesianIndex}' href='#Atlans.prepare_surcharge_column-Tuple{Atlans.Surcharge, Atlans.SoilColumn, CartesianIndex}'>#</a>
**`Atlans.prepare_surcharge_column`** &mdash; *Method*.



```julia
prepare_surcharge_column(sur::Surcharge, column::SoilColumn, I::CartesianIndex)
```

Create a SurchargeColumn with correct groundwater, consolidation, oxidation and shrinkage Surcharge columns. The correct Atlantis processes (e.g. DrainingAbcIsotache) that each of the Surcharge columns are built-up from, are derived from the input SoilColumn. The CartesianIndex reads to lithology and thickness to build the Surcharge column from at the correct location in the Surcharge forcing input.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.prepare_timestep!-Tuple{Atlans.SoilColumn, Any}' href='#Atlans.prepare_timestep!-Tuple{Atlans.SoilColumn, Any}'>#</a>
**`Atlans.prepare_timestep!`** &mdash; *Method*.



```julia
prepare_timestep!(column)
```

Prepare a single timestep. This updates stresses in the column and accounts for e.g. drowning of the column.

This computes:

```
* Pore pressure
* Total stress
* Effective stress
```


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.prepare_timestep!-Tuple{Atlans.SurchargeColumn, Any}' href='#Atlans.prepare_timestep!-Tuple{Atlans.SurchargeColumn, Any}'>#</a>
**`Atlans.prepare_timestep!`** &mdash; *Method*.



```julia
prepare_timestep!(column::SurchargeColumn, Δt)
```

Set the initial stresses for a SurchargeColumn.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.relative_oxidation_rate' href='#Atlans.relative_oxidation_rate'>#</a>
**`Atlans.relative_oxidation_rate`** &mdash; *Function*.



```julia
relative_oxidation_rate(T::Float)
```

Empirical relation between the decay rate of organic soils and air temperature in Celsius from Hendriks and Vermeulen (1997). Relation is valid between 0° and 19.5° Celsius.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.repeat_elements-Tuple{Any, Any}' href='#Atlans.repeat_elements-Tuple{Any, Any}'>#</a>
**`Atlans.repeat_elements`** &mdash; *Method*.



Repeat function not in base / stdlib


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.run!-Tuple{Any}' href='#Atlans.run!-Tuple{Any}'>#</a>
**`Atlans.run!`** &mdash; *Method*.



Run all forcing periods of the simulation.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.set_periods!-Tuple{Any, Any}' href='#Atlans.set_periods!-Tuple{Any, Any}'>#</a>
**`Atlans.set_periods!`** &mdash; *Method*.



Collect the period boundaries from the forcing input.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.set_surcharge!-Tuple{Atlans.SoilColumn, Atlans.SurchargeColumn}' href='#Atlans.set_surcharge!-Tuple{Atlans.SoilColumn, Atlans.SurchargeColumn}'>#</a>
**`Atlans.set_surcharge!`** &mdash; *Method*.



```julia
set_surcharge!(column::SoilColumn, surcharge::SurchargeColumn)
```

Combine a SoilColumn and SurchargeColumn when Surcharge is added as a forcing during a forcingperiod.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.shrink-Tuple{Atlans.SimpleShrinkage, Float64}' href='#Atlans.shrink-Tuple{Atlans.SimpleShrinkage, Float64}'>#</a>
**`Atlans.shrink`** &mdash; *Method*.



```julia
shrink(voxel, Δt)
```

Shrink a voxel for given time interval.

#Arguments

  * `voxel::SimpleShrinkage`: Voxel to shrink.
  * `Δt::Float`: Time interval. [days]


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.subside!-Tuple{Atlans.SoilColumn}' href='#Atlans.subside!-Tuple{Atlans.SoilColumn}'>#</a>
**`Atlans.subside!`** &mdash; *Method*.



```julia
subside!(column)
```

Apply consolidation, oxidation and shrinkage to thickness


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.total_stress!-Tuple{Atlans.AbstractConsolidationColumn, Any}' href='#Atlans.total_stress!-Tuple{Atlans.AbstractConsolidationColumn, Any}'>#</a>
**`Atlans.total_stress!`** &mdash; *Method*.



Compute total stress for entire column


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.transfer_stress!-Tuple{Atlans.AbstractConsolidationColumn}' href='#Atlans.transfer_stress!-Tuple{Atlans.AbstractConsolidationColumn}'>#</a>
**`Atlans.transfer_stress!`** &mdash; *Method*.



Transfer computed stress to the cells of the ConsolidationColumn.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.update_alpha-Tuple{Atlans.CarbonStore, Float64}' href='#Atlans.update_alpha-Tuple{Atlans.CarbonStore, Float64}'>#</a>
**`Atlans.update_alpha`** &mdash; *Method*.



```julia
update_alpha(cell::CarbonStore, T::Float)
```

Return a new CarbonStore cell with an oxidation rate (α) that is corrected for air temperature.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.update_z!-Tuple{Atlans.SoilColumn}' href='#Atlans.update_z!-Tuple{Atlans.SoilColumn}'>#</a>
**`Atlans.update_z!`** &mdash; *Method*.



```julia
function update_z!(column)
```

Compute new midpoints and surface level.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.volume_organic-Tuple{Any, Any}' href='#Atlans.volume_organic-Tuple{Any, Any}'>#</a>
**`Atlans.volume_organic`** &mdash; *Method*.



Empirical equation to compute specific volume of organic material.

As the organic matter in a soil breaks down, density increases. The "airiest" parts are the first to go.


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>

<a id='Atlans.weight-NTuple{5, Float64}' href='#Atlans.weight-NTuple{5, Float64}'>#</a>
**`Atlans.weight`** &mdash; *Method*.



Weight of (part of) a single cell


<a target='_blank' href='https://gitlab.com/deltares/subsidence/atlans.jl.git' class='documenter-source'>source</a><br>


<a id='Constants'></a>

<a id='Constants-1'></a>

## Constants


<a id='Macros'></a>

<a id='Macros-1'></a>

## Macros


<a id='Index'></a>

<a id='Index-1'></a>

## Index

- [`Atlans.AdaptiveCellsize`](index.md#Atlans.AdaptiveCellsize)
- [`Atlans.Clock`](index.md#Atlans.Clock)
- [`Atlans.ExponentialTimeStepper`](index.md#Atlans.ExponentialTimeStepper)
- [`Atlans.Model`](index.md#Atlans.Model-Tuple{Type, Type, Type, Type, Type, Vararg{Any, 4}})
- [`Atlans.ShrinkageColumn`](index.md#Atlans.ShrinkageColumn)
- [`Atlans.SimpleShrinkage`](index.md#Atlans.SimpleShrinkage)
- [`Atlans.Simulation`](index.md#Atlans.Simulation-Tuple{Atlans.Model, String, Dates.DateTime})
- [`Atlans.SoilColumn`](index.md#Atlans.SoilColumn)
- [`Atlans.VerticalDomain`](index.md#Atlans.VerticalDomain)
- [`Atlans.Qcreep_derivative`](index.md#Atlans.Qcreep_derivative-Tuple{Atlans.AbcIsotache, Float64, Float64})
- [`Atlans.U`](index.md#Atlans.U-Tuple{Atlans.ConsolidationProcess, Float64})
- [`Atlans.add_time`](index.md#Atlans.add_time-Tuple{Any, Any})
- [`Atlans.advance!`](index.md#Atlans.advance!-Tuple{Any})
- [`Atlans.advance_forcingperiod!`](index.md#Atlans.advance_forcingperiod!-Tuple{Any})
- [`Atlans.advance_forcingperiod!`](index.md#Atlans.advance_forcingperiod!-Tuple{Any, Any})
- [`Atlans.advance_forcingperiod!`](index.md#Atlans.advance_forcingperiod!-Tuple{Atlans.SoilColumn, Vector{Float64}})
- [`Atlans.advance_timestep!`](index.md#Atlans.advance_timestep!-Tuple{Atlans.SoilColumn, Float64})
- [`Atlans.cellsplit!`](index.md#Atlans.cellsplit!-Tuple{Union{Atlans.ConsolidationColumn{Atlans.NullConsolidation, Atlans.OverConsolidationRatio}, Atlans.OxidationColumn{Atlans.NullOxidation}, Atlans.ShrinkageColumn{Atlans.NullShrinkage}}, Vararg{Any, 4}})
- [`Atlans.cellsplit!`](index.md#Atlans.cellsplit!-Tuple{Atlans.OxidationColumn{Atlans.CarbonStore}, Vararg{Any, 4}})
- [`Atlans.compress_γ_dry`](index.md#Atlans.compress_γ_dry-Tuple{Atlans.ConsolidationProcess, Float64})
- [`Atlans.compress_γ_wet`](index.md#Atlans.compress_γ_wet-Tuple{Atlans.ConsolidationProcess, Float64})
- [`Atlans.consolidate`](index.md#Atlans.consolidate-Tuple{Atlans.DrainingAbcIsotache, Any, Any})
- [`Atlans.create_timesteps`](index.md#Atlans.create_timesteps-Tuple{Atlans.ExponentialTimeStepper, Any})
- [`Atlans.currenttime`](index.md#Atlans.currenttime-Tuple{Any})
- [`Atlans.discretize`](index.md#Atlans.discretize-Tuple{Any, Float64})
- [`Atlans.draining_abc_isotache_column`](index.md#Atlans.draining_abc_isotache_column-NTuple{8, Any})
- [`Atlans.effective_stress!`](index.md#Atlans.effective_stress!-Tuple{Atlans.AbstractConsolidationColumn})
- [`Atlans.formulate`](index.md#Atlans.formulate-Tuple{Atlans.AbcIsotache, Vararg{Float64, 4}})
- [`Atlans.initialize`](index.md#Atlans.initialize-Tuple{Type{Atlans.NullOxidation}, Any, Any, Any})
- [`Atlans.initialize`](index.md#Atlans.initialize-Tuple{Type{Atlans.DrainingAbcIsotache}, Type, Atlans.VerticalDomain, Dict})
- [`Atlans.initialize`](index.md#Atlans.initialize-Tuple{Type{Atlans.NullConsolidation}, Type, Any, Any})
- [`Atlans.initialize`](index.md#Atlans.initialize-Tuple{Type{Atlans.SimpleShrinkage}, Any, Any, Any})
- [`Atlans.initialize`](index.md#Atlans.initialize-Tuple{Type{Atlans.NullShrinkage}, Any, Any, Any})
- [`Atlans.initialize`](index.md#Atlans.initialize-Tuple{Type{Atlans.DrainingAbcIsotache}, Type, Any, Any, Any})
- [`Atlans.initialize`](index.md#Atlans.initialize-Tuple{Type{Atlans.NullShrinkage}, Any, Any})
- [`Atlans.initialize`](index.md#Atlans.initialize-Tuple{Type{Atlans.NullOxidation}, Any, Any})
- [`Atlans.initialize`](index.md#Atlans.initialize-Tuple{Type{Atlans.CarbonStore}, Any, Any, Any})
- [`Atlans.initialize`](index.md#Atlans.initialize-Tuple{Type{Atlans.CarbonStore}, Atlans.VerticalDomain, Dict})
- [`Atlans.initialize`](index.md#Atlans.initialize-Tuple{Type{Atlans.SimpleShrinkage}, Atlans.VerticalDomain, Dict})
- [`Atlans.initialize`](index.md#Atlans.initialize-Tuple{Type{Atlans.HydrostaticGroundwater}, Atlans.Phreatic, Atlans.VerticalDomain})
- [`Atlans.initialize`](index.md#Atlans.initialize-Tuple{Type{Atlans.NullConsolidation}, Vararg{Any, 4}})
- [`Atlans.parse_loglevel`](index.md#Atlans.parse_loglevel-Tuple{AbstractString})
- [`Atlans.periodduration`](index.md#Atlans.periodduration-Tuple{Any})
- [`Atlans.pow`](index.md#Atlans.pow-Tuple{Any, Any})
- [`Atlans.prepare_domain`](index.md#Atlans.prepare_domain-NTuple{7, Any})
- [`Atlans.prepare_forcingperiod!`](index.md#Atlans.prepare_forcingperiod!)
- [`Atlans.prepare_forcingperiod!`](index.md#Atlans.prepare_forcingperiod!-Tuple{Atlans.ConsolidationColumn{Atlans.DrainingAbcIsotache, P} where P<:Atlans.Preconsolidation})
- [`Atlans.prepare_surcharge_column`](index.md#Atlans.prepare_surcharge_column-Tuple{Atlans.Surcharge, Atlans.SoilColumn, CartesianIndex})
- [`Atlans.prepare_timestep!`](index.md#Atlans.prepare_timestep!-Tuple{Atlans.SurchargeColumn, Any})
- [`Atlans.prepare_timestep!`](index.md#Atlans.prepare_timestep!-Tuple{Atlans.SoilColumn, Any})
- [`Atlans.relative_oxidation_rate`](index.md#Atlans.relative_oxidation_rate)
- [`Atlans.repeat_elements`](index.md#Atlans.repeat_elements-Tuple{Any, Any})
- [`Atlans.run!`](index.md#Atlans.run!-Tuple{Any})
- [`Atlans.set_periods!`](index.md#Atlans.set_periods!-Tuple{Any, Any})
- [`Atlans.set_surcharge!`](index.md#Atlans.set_surcharge!-Tuple{Atlans.SoilColumn, Atlans.SurchargeColumn})
- [`Atlans.shrink`](index.md#Atlans.shrink-Tuple{Atlans.SimpleShrinkage, Float64})
- [`Atlans.subside!`](index.md#Atlans.subside!-Tuple{Atlans.SoilColumn})
- [`Atlans.total_stress!`](index.md#Atlans.total_stress!-Tuple{Atlans.AbstractConsolidationColumn, Any})
- [`Atlans.transfer_stress!`](index.md#Atlans.transfer_stress!-Tuple{Atlans.AbstractConsolidationColumn})
- [`Atlans.update_alpha`](index.md#Atlans.update_alpha-Tuple{Atlans.CarbonStore, Float64})
- [`Atlans.update_z!`](index.md#Atlans.update_z!-Tuple{Atlans.SoilColumn})
- [`Atlans.volume_organic`](index.md#Atlans.volume_organic-Tuple{Any, Any})
- [`Atlans.weight`](index.md#Atlans.weight-NTuple{5, Float64})

