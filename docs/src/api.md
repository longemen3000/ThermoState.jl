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
ThermoState.get_spec
ThermoState.has_spec

ThermoState.state_type

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

The following properties are accepted by an state and have an accesor function:

| Property             |Units     |Accessor function      |Keywords for `spec` and `state` |
|----------------------|----------|-----------------------|--------------------------------|
| Molar volume         |m3/mol    |`mol_volume`           | `v`, `mol_v`                   |
| Mass volume          |m3/kg     |`mass_volume`          | `mass_v`                       |
| Total volume         |m3        |`total_volume`         | `total_v`, `V`                 |
| Molar density        |mol/m3    |`mol_density`          | `ρ`, `mol_ρ`, `rho`, `mol_rho` |
| Mass density         |mol/kg    |`mass_density`         | `mass_ρ`, `mass_rho`           |
| Temperature          |K         |`temperature`          | `t`, `T`                       |
| Pressure             |Pa        |`pressure`             | `p`, `P`                       |
| Molar enthalpy       |J/mol     |`mol_enthalpy`         | `h`, `mol_h`                   |
| Mass enthalpy        |J/kg      |`mass_enthalpy`        | `mass_h`                       |
| Total enthalpy       |J         |`total_enthalpy`       | `total_h`                      |
| Molar Gibbs          |J/mol     |`mol_gibbs`            | `g`, `mol_g`                   |
| Mass Gibbs           |J/kg      |`mass_gibbs`           | `mass_g`                       |
| Total Gibbs          |J         |`total_gibbs`          | `total_g`                      |
| Molar Helmholtz      |J/mol     |`mol_helmholtz`        | `a`, `mol_a`                   |
| Mass Helmholtz       |J/kg      |`mass_helmholtz`       | `mass_a`                       |
| Total Helmholtz      |J         |`total_helmholtz`      | `total_a`                      |
| Molar Internal Energy|J/mol     |`mol_internal_energy`  | `a`, `mol_a`                   |
| Mass Internal Energy |J/kg      |`mass_internal_energy` | `mass_a`                       |
| Total Internal Energy|J         |`total_internal_energy`| `total_a`                      |
| Molar entropy        |J/(mol K) |`mol_entropy`          | `s`, `mol_s`                   |
| Mass entropy         |J/(kg K)  |`mass_entropy`         | `mass_s`                       |
| Total entropy        |J/K       |`total_entropy`        | `total_s`                      |
| Amount of Moles      |mol       |`moles`                | `moles`                        |
| Amount of Mass       |kg        |`mass`                 | `mass`                         |
| Molar Fraction       |no units  |`mol_fraction`         | `xn`                           |
| Mass Fraction        |no units  |`mass_fraction`        | `xm`                           |
| Molar Numbers        |mol       |`mol_number`           | `n`                            |
| Mass Numbers         |kg        |`mass_number`          | `m`                            |
| Molar Vapor Fraction |no units  |`mol_vapor_fraction`   | `vfrac`                        |
| Mass Vapor Fraction  |no units  |`mass_vapor_fraction`  | `quality`                      |
| Current Phase        |no units  |`phase`                | `phase`, `sat`                 |
| Aditional Options    |no units  |`options`              | `options`                      |
---

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
```@docs
ThermoState.mol_fraction
ThermoState.mass_fraction
ThermoState.mol_number
ThermoState.mass_number
ThermoState.moles
ThermoState.mass
ThermoState.molar_mass
```

### Other Properties
```@docs
ThermoState.options
ThermoState.phase
ThermoState.quality
ThermoState.molecular_weight
```