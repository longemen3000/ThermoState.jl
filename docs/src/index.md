```@meta
CurrentModule = ThermoState
```

# ThermoState
ThermoState.jl is a basic block for interfacing and specifying thermodynamic models. this package proposes and provides tools to create a common property interface of the form: 

```julia
property(model,state::ThermodynamicState,unit)::ThermodynamicState
```

This package proposes the following conventions for naming properties:

 - `mol_$PROPERTY` is a property of units U/mol (molar volume, molar Helmholtz energy, molar enthalpy, etc)

 - `mass_$PROPERTY` is a property of units U/kg.

 - `total_$PROPERTY` is a property of units U. (total Helmholtz energy has joule units)



