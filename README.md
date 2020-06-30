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

For defining property specifications, the package defines the `AbstractSpec` type and the `Spec` container. an individual specification can be defined by:

```julia
using Unitful, PhysPropsRules
h0 = Spec(Enthalpy(),30.0u"kg/mol")
```
you can create various specifications with the `specs` function:




