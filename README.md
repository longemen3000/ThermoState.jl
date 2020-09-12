# ThermoState

[![Build Status](https://travis-ci.com/longemen3000/ThermoState.jl.svg?branch=master)](https://travis-ci.com/longemen3000/ThermoState.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://longemen3000.github.io/ThermoState.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://longemen3000.github.io/ThermoState.jl/dev)

Basic block for interfacing with thermodynamic models.
The idea is to have a common interface of the form:

```julia
property(model,state,unit,options)
```

## Basics
This package proposes the following conventions for naming properties:

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
So, instead of defining multiple functions for each combination of specifications, you just dispatch your new model on:
```julia
function property(model::MyModel,state::ThermodynamicState;kwargs...)
    v = mol_volume(FromSpecs(),state,unit,molecular_weight(model)
    T = mol_volume(FromSpecs(),state,u"°C",molecular_weight(model)) #default: SI units (Kelvin)
    return property_impl(v,T)
end
```
## Specification object (`Spec`)

A `Spec` is just a tagged value. it can be constructed by two ways:
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
All keyword arguments are stored in the `KW_TO_SPEC` dict keys.

The result of those operations is a `Spec` struct:
```julia-repl
julia> t0 = spec(T = 30.0u"°C")     
Temperature : 303.15 K

julia> typeof(t0)
Spec{ThermoState.Types.Temperature,Float64}

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
Spec{ThermoState.Types.Enthalpy{ThermoState.Types.MOLAR},Float64} ##?
```
In this case, the specification type is a parametric singleton struct: `Enthalpy{MOLAR}`. the `MOLAR` parameter is known as a _spec modifier_ , and is used for dispatch in unit transformations and conversions (from molar to mass units, for example).

There are two special cases with two parameters: volume amounts (molar, total and specific volume, mass and molar density) and material compound proportions (mass numbers, mol numbers, mass fractions and mol fractions). volume amounts are tagged with the specification `VolumeAmount{<:SpecModifier,<:SpecModifier}` and material compounds tagged with the specification `MaterialCompounds{<:SpecModifier,<:SpecModifier}`. lets see some examples:

```julia-repl
julia> using ThermoState.Types #importing the spec types for better printing


julia> typeof(specification(spec(mol_v = 0.005)))
VolumeAmount{MOLAR,VOLUME}

julia> typeof(specification(spec(mass_rho = 875.2)))
VolumeAmount{TOTAL,DENSITY}

julia> typeof(specification(spec(xn = [0.5,0.5])))
MaterialCompounds{MOLAR,FRACTION}

julia> typeof(specification(spec(m = [40.2,35.3])))
MaterialCompounds{MASS,TOTAL_AMOUNT}
```

There are two special types: `PhaseTag` (keyword = `phase`) and `Options` (keyword = `options`). the first is useful to signal an underlying model to calculate just one phase (for example, in cubic equations, the root calculation gives the vapor and liquid volumes). `Options` accepts a named tuple that can be passed to any underlying model to specify numerical options (choice of differentiation method, maximum iterations,etc).

## `ThermodynamicState`

A specification is not sufficient to define a system. The [Gibbs' Phase Rule](https://en.wikipedia.org/wiki/Phase_rule) specifies that a system in equilibrium has `F` degrees of freedom, where `F = NumberOfComponents - NumberOfPhases + 2 `. 

The `ThermodynamicState` struct is a collection of `Spec`s. When creating this object, the Gibbs' Phase Rule is evaluated on the specification arguments to check its validity. lets see one example: 

```julia-repl
julia> a = state(t=300.0,ρ=5u"mol/m^3")
2 Property Specifications:
 Temperature : 300.0 K
 Molar density : 5 mol m^-3
 ```
In this case, neither a phase nor any amount of compounds was specified. In those cases, the function assumes when the system has one mol, one phase, and/or a single component.

Another way to build a `ThermodynamicState` struct is by directly passing `Spec`s as arguments:
```julia-repl
 julia>
 h0 = spec(mol_h = 3000.0)
 t0 = spec(t=401.0)
 st = state(h0,t0)
 2 Property Specifications:
 Molar enthalpy : 3000.0 J mol^-1
 Temperature : 401.0 K
```
You can skip the check of the Gibbs' Phase Rule using the `check = false` keyword, or, in the case of building the state using keywords, decide to not normalize units via the 
`normalize_units= false` keyword.

## Obtaining properties from a `ThermodynamicState` struct

as said in the Basics section, you can query properties from the created `ThermodynamicState` struct. the way of obtaining those properties is by calling a property function. The interface proposed by this package is the following:

```julia
prop = property(model::MyModel,state::ThermodynamicState,unit::Unitful.unit=[default],args...)

```
dispatching on the model type, you can calculate properties from the specifications contained in `state`.
This package exports one single model: `FromState`, that doesn't calculate (almost) anything, just checks if the selected property is in the argument, and if it is, it returns its numerical (unit stripped) value.
By default, the units of the number returned correspond to SI units, you can obtain a number with appropiate units by passing a corresponding `unit` argument.
the `FromState` model, given a corresponding molecular weight `mw` argument, can calculate derived properties, for example:

```julia-repl
a = state(t=300.0,ρ=5.0u"mol/L")
v = mass_volume(FromSpecs(),a,u"cm^3/g",50)
4.0
```
The `mw` argument can be, depending of the situation, a number of vector of number. if not units are provided, a default units of `g/mol` are assumed. this conversion can be done on any unit that accepts molar, total and mass specifications, mass and moles themselves, molar and mass fractions,and molar and mass numbers.

## Variable `ThermodynamicState`s

sometimes is needed a more direct approach to evaluation of properties.For example, you may want to create a function that accepts only temperature to pass it to an ODE system or an optimization system. for this purpose, the Singleton `VariableSpec` is provided. if you pass it to a `Spec` and create a `ThermodynamicState` object (or assign via keywords) the resulting state will be callable:

```julia-repl
a_t = state(t=VariableSpec(),ρ=5.0u"mol/L")
a_t(303)
```


## Implementing a model using the `ThermoState` interface

using this package, we can implement a basic ideal gas model that only calculates the pressure, given a temperature and molar volume, using the relation `Pv=RT`:


```julia
using ThermoState, Unitful
struct MyIdealGas
    mw::Float64
end
function pressure(model::MyIdealGas,state::ThermodynamicState,unit=u"Pa")
    v = mol_volume(FromSpecs(),state,u"m^3/mol",model.mw)
    t = temperature(FromSpecs(),state) #Kelvin as default, doesnt require mw
    p =  (8.314*t/v)u"Pa"
    return Unitful.ustrip(Unitful.uconvert(unit,p))
end

a = state(mass = 3u"kg",total_v = 30u"m^3",t=30u"°C")
model = MyIdealGas(18.01)
p = pressure(model,a)

```
Using variable state:

```julia
tx = state(mass = 3u"kg",total_v = 30u"m^3",t=var) #one free variable
p = pressure(model,tx(30u"°C"))
p_list = map(t-> pressure(model,tx(t)),273.0:373.0)
```

## State of this package

At the moment of writing this, this package is a  experimental state and many features could (and will) change, please write any sugerences on the issues!, pull requests are very appreciated!


