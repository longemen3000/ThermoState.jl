
For defining property specifications, the package defines the `AbstractSpec` types and the `Spec` container. an individual specification can be defined by:

```julia
using Unitful, ThermoState
h0 = spec(mol_h="kg/mol")
```
You can create various specifications with the `state` function:

```julia
st = state(v=990u"m^3/mol",T=350u"K",mass=3u"mg")
```
And extract values for the use of property models:

```julia
mass_rho0=mass_density(FromState(),st,"kg/L",18.0u"g/mol")
```

## Specification Object (`Spec`)

A `Spec` is just a tagged value. it can be constructed by two ways:
- type-value constructor: 
  ```julia
   t0 = spec(Types.Temperature(),300.15)
   t0 = spec(Types.Temperature(),30u"°C") #normalized
   t0 = spec(Types.Temperature(),30u"°C",false) #not normalized
  ```

- type-value constructor (using `@spec_str` macro): 
  ```julia
   t0 = spec(spec"t",300.15)
   t0 = spec(spec"t",30u"°C") #normalized
   t0 = spec(spec"t",30u"°C",false) #not normalized
  ```
  by default, unitful quantities are unit stripped and normalized to SI units, you can use the argument `normalize_units` to change that default.

- keyword constructor: 
    ```julia
   t0 = spec(t= 300.15)
   t0 = spec(T= 30u"°C") 
   t0 = spec(T = 30u"°C",normalize_units= false) 
  ```

All keyword arguments are stored in the `KW_TO_SPEC` dict keys.
The Main difference between the type-value constructor and the keyword constructor is that the first can be resolved at compile time, where as the second has a runtime cost.


The result of those operations is a `Spec` struct:
```julia-repl
julia> t0 = spec(T = 30.0u"°C")     
spec(t = 303.15[K])

julia> typeof(t0)
Spec{ThermoState.Types.Temperature,Float64}

julia> val_t0 = value(t0) #extracting value
303.15

julia> spec_t0 = specification(t0) #specification of t0 
Temperature()
```
an example of a enthalpy specification:
```julia-repl
julia> h0 = spec(mol_h = 3000.0)
spec(mol_h = 3000.0[J mol^-1])

julia> specification(h0)
Enthalpy{MOLAR}()

julia> typeof(h0)
Spec{ThermoState.Types.Enthalpy{ThermoState.Types.MOLAR},Float64} ##?
```
In this case, the specification type is a parametric singleton struct: `Enthalpy{MOLAR}`. the `MOLAR` parameter is known as a _spec modifier_ , and is used for dispatch in unit transformations and conversions (from molar to mass units, for example).

There are two special cases with two parameters: volume amounts (molar, total and specific volume, mass and molar density) and material compound proportions (mass numbers, mol numbers, mass fractions and mol fractions). volume amounts are tagged with the specification `VolumeAmount{<:SpecModifier,<:SpecModifier}` and material compounds tagged with the specification `MaterialCompounds{<:SpecModifier,<:SpecModifier}`. lets see some examples:

```julia-repl
julia> using ThermoState.Types #importing the spec types for shorter printing


julia> typeof(specification(spec(mol_v = 0.005)))
VolumeAmount{MOLAR,VOLUME}

julia> typeof(spec"mol_v")
VolumeAmount{MOLAR,VOLUME}

julia> typeof(specification(spec(mass_rho = 875.2)))
VolumeAmount{TOTAL,DENSITY}

julia> typeof(specification(spec(xn = [0.5,0.5])))
MaterialCompounds{MOLAR,FRACTION}

julia> typeof(specification(spec(m = [40.2,35.3])))
MaterialCompounds{MASS,TOTAL_AMOUNT}
```

There are two special types: `PhaseTag` (keyword = `phase`) and `Options` (keyword = `options`). the first is useful to signal an underlying model to calculate just one phase (for example, in cubic equations, the root calculation gives the vapor and liquid volumes). `Options` accepts a named tuple that can be passed to any underlying model to specify numerical options (choice of differentiation method, maximum iterations,etc).
