# PhysPropsRules

[![Build Status](https://travis-ci.com/longemen3000/PhysPropsRules.jl.svg?branch=master)](https://travis-ci.com/longemen3000/PhysPropsRules.jl)

Basic block for interfacing with thermodynamic models.
The idea is yo have a common interface of the form:

```julia
property(model,specs,options)
```
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

For defining property specifications, the package defines the `AbstractSpec` type and the `Spec` container. an individual specification can be defined by:

```julia
using Unitful, PhysPropsRules
h0 = Spec(Enthalpy(),30.0u"kg/mol")
```
you can create various specifications with the `specs` function:

```julia
props = specs(v=990u"m3/mol",T=350u"K",mass=3u"mg")
```
and extract values for the use of property models. calling the `specs` functions return a `Specs` struct, containing the collection of Spec's, and checking if the Gibbs phase rule is consistent.
Stracting a specification from the `Specs` struct can be done as if the specifications themselves are a thermodinamic model: 

```julia
mass_rho=mass_density(FromSpecs(),props,mw=18.0u"g/mol",units=)
```
so, instead of defining multiple functions for each combination of specifications, you just dispatch your new model on:
```julia
function property(model::MyModel,specs::Specs;kwargs...)
    v = mol_volume(FromSpecs(),specs,mw=molecular_weight(model),units = u"mol/cm3")
    T = mol_volume(FromSpecs(),specs) #default: SI units (Kelvin)
    return property_impl(v,T)
end
```
This package is experimental and many features could (and will) change, please write any sugerences on the issues!


