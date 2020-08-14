# PhysPropsRules

[![Build Status](https://travis-ci.com/longemen3000/PhysPropsRules.jl.svg?branch=master)](https://travis-ci.com/longemen3000/PhysPropsRules.jl)

Basic block for interfacing with thermodynamic models.
The idea is to have a common interface of the form:

```julia
property(model,specs,options)
```

## Basics
this package proposes the following conventions for naming properties:

 - `mol_$PROPERTY` is a property of units U/mol (molar volume, molar Helmholtz energy, molar enthalpy, etc)

 - `mass_$PROPERTY` is a property of units U/kg.

 - `total_$PROPERTY` is a property of units U. (total Helmholtz energy has joule units)

At this moment, the following properties are defined:
 - `volume` (`total_volume`,`mol_volume`,`mass_volume`)
 - `density` (`mol_density`,`mass_density`)
 - `pressure`
 - `temperature`
 - `enthalpy` (`total_enthalpy`,`mol_entalphy`,`mass_enthalpy`)
 - `entropy` (`total_entropy`,`mol_entropy`,`mass_entropy`)
 - `gibbs` (`total_gibbs`,`mol_gibbs`,`mass_gibbs`) 
 - `helmholtz` (`total_helmholtz`,`mol_helmholtz`,`mass_helmholtz`) 
 - `internal_energy` (`total_internal_energy`,`mol_internal_energy`,`mass_internal_energy`) 
 - `moles`
 - `mass`

For defining property specifications, the package defines the `AbstractSpec` type and the `Spec` container. an individual specification can be defined by:

```julia
using Unitful, PhysPropsRules
h0 = spec(mol_h="kg/mol")
```
you can create various specifications with the `specs` function:

```julia
props = specs(v=990u"m3/mol",T=350u"K",mass=3u"mg")
```
and extract values for the use of property models:

```julia
mass_rho0=mass_density(FromSpecs(),props,"kg/L",18.0u"g/mol")
```
so, instead of defining multiple functions for each combination of specifications, you just dispatch your new model on:
```julia
function property(model::MyModel,specs::Specs;kwargs...)
    v = mol_volume(FromSpecs(),specs,unit,molecular_weight(model)
    T = mol_volume(FromSpecs(),specs,u"°C",molecular_weight(model)) #default: SI units (Kelvin)
    return property_impl(v,T)
end
```
## Specification object (`Spec`)

a `Spec` is just a tagged value. it can be constructed by two ways:
- type-value constructor: 
  ```julia
   t0 = spec(temperature,300.15)
   t0 = spec(temperature,30u"°C") #normalized
   t0 = spec(temperature,30u"°C",false) #not normalized
  ```
  by default, unitful quantities are unit stripped and normalized to SI units, you can use the argument `normalize_units` to change that default.
  there is a slight thing to mention: `temperature` is really a function, binded to the `Temperature <: AbstractSpec` type. this is for convenience, as we will see later.
- keyword constructor: 
    ```julia
   t0 = spec(t=300.15)
   t0 = spec(T=30u"°C") 
   t0 = spec(T = 30u"°C",normalize_units= false) 
  ```
all keyword arguments are stored in the `KW_TO_SPEC` dict keys.

The result of those operations is a `Spec` struct:
```julia-repl
julia> t0 = spec(T = 30.0u"°C")     
Temperature : 303.15 K

julia> typeof(t0)
Spec{PhysPropsRules.Types.Temperature,Float64}

julia> val_t0 = value(t0) #extracting value
303.15

julia> spec_t0 = specification(t0) #specification of t0 
::Temperature
```
an example of a enthalpy specification:
```julia-repl
julia> h0 = spec(mol_h = 3000.0)
Molar enthalpy : 3000.0 J mol^-1   

julia> specification(h0)
::Molar enthalpy

julia> typeof(h0)
Spec{PhysPropsRules.Types.Enthalpy{PhysPropsRules.Types.MOLAR},Float64} ##?
```
in this case, the specification type is a parametric singleton struct: `Enthalpy{MOLAR}`. the `MOLAR` parameter is known as a _spec modifier_ , and is used for dispatch in unit transformations and conversions (from molar to mass units, for example)

There are two special cases with two parameters: volume amounts (molar, total and specific volume, mass and molar density) and material compound proportions (mass numbers, mol numbers, mass fractions and mol fractions). volume amounts are tagged with the specification `VolumeAmount{<:SpecModifier,<:SpecModifier}` and material compounds tagged with the specification `MaterialCompounds{<:SpecModifier,<:SpecModifier}`. lets see some examples:

```julia-repl
julia> using PhysPropsRules.Types #importing the spec types for better printing


julia> typeof(specification(spec(mol_v = 0.005)))
VolumeAmount{MOLAR,VOLUME}

julia> typeof(specification(spec(mass_rho = 875.2)))
VolumeAmount{TOTAL,DENSITY}

julia> typeof(specification(spec(xn = [0.5,0.5])))
MaterialCompounds{MOLAR,FRACTION}

julia> typeof(specification(spec(m = [40.2,35.3])))
MaterialCompounds{MASS,TOTAL}
```


This package is experimental and many features could (and will) change, please write any sugerences on the issues!


