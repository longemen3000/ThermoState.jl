## Contents

```@contents
Pages = ["api.md"]
```

## Index

```@index
Pages = ["api.md"]
```
## Specification and Thermodynamic State
```@docs
ThermoState.spec
ThermoState.state
```

## Utilities

```@docs
ThermoState.@to_units
ThermoState.@spec_str
ThermoState.normalize_units
ThermoState.convert_unit
ThermoState.default_units
```

## Defined Properties

### PVT Properties
```@docs
ThermoState.temperature
ThermoState.pressure
ThermoState.mol_volume
ThermoState.mass_volume
ThermoState.total_volume
ThermoState.mol_density
ThermoState.mass_density
```

### Energy Properties
```@docs
ThermoState.mol_helmholtz
ThermoState.mass_helmholtz
ThermoState.total_helmholtz
ThermoState.mol_gibbs
ThermoState.mass_gibbs
ThermoState.total_gibbs
ThermoState.mol_enthalpy
ThermoState.mass_enthalpy
ThermoState.total_enthalpy
ThermoState.mol_entropy
ThermoState.mass_entropy
ThermoState.total_entropy
ThermoState.mol_internal_energy
ThermoState.mass_internal_energy
ThermoState.total_internal_energy
```

### Mass and Molar Properties
ThermoState.mol_fraction
ThermoState.mass_fraction
ThermoState.mol_number
ThermoState.mass_number
ThermoState.moles
ThermoState.mass
ThermoState.molar_mass

### Other Properties
ThermoState.options
ThermoState.phase
ThermoState.quality
ThermoState.molecular_weight
```