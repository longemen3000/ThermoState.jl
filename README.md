# ThermoState.jl

[![Build Status](https://travis-ci.com/longemen3000/ThermoState.jl.svg?branch=master)](https://travis-ci.com/longemen3000/ThermoState.jl)
[![Build Status](https://github.com/longemen3000/ThermoState.jl/workflows/CI/badge.svg)](https://github.com/longemen3000/ThermoState.jl/actions)
[![Codecov](https://codecov.io/gh/longemen3000/ThermoState.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/longemen3000/ThermoState.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://longemen3000.github.io/ThermoState.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://longemen3000.github.io/ThermoState.jl/dev)

ThermoState.jl is a basic block for interfacing and specifying thermodynamic models. this package proposes and provides tools to create a common property interface of the form: 

```julia
property(model,state,unit)
```

## Basics
This package proposes the following conventions for naming properties:

 - `mol_$PROPERTY` is a property of units U/mol (molar volume, molar Helmholtz energy, molar enthalpy, etc)

 - `mass_$PROPERTY` is a property of units U/kg.

 - `total_$PROPERTY` is a property of units U. (total Helmholtz energy has joule units)

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


You can view the keywords in the `KW_TO_SPEC` constant.

The Following functions aren't accepted as a keyword to `state` or `spec`, but are defined on the package:


| Property             |Units     |Accessor function      |Notes                                          |
|----------------------|----------|-----------------------|-----------------------------------------------|
| Molar Mass           |kg/mol    |`molar_mass`           |molecular weight ponderated by material amounts|
| Molar Cₚ              |J/(mol K) |`mol_cp`               |no impl, for interop use                      |
| Mass Cₚ               |J/(kg K)  |`mass_cp`              |no impl, for interop use                      |
| Molar Cᵥ             |J/(mol K) |`mol_cv`               |no impl, for interop use                       |
| Mass Cᵥ              |J/(kg K)  |`mass_cv`              |no impl, for interop use                       |
| Speed of sound       |m/s       |`sound_speed`          |no impl, for interop use                       |
---


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

## Specification object (`Spec`)

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

## `ThermodynamicState`

A specification is not sufficient to define a system. The [Gibbs' Phase Rule](https://en.wikipedia.org/wiki/Phase_rule) specifies that a system in equilibrium has `F` degrees of freedom, where `F = NumberOfComponents - NumberOfPhases + 2 `. 

The `ThermodynamicState` struct is a collection of `Spec`s. When creating this object, the Gibbs' Phase Rule is evaluated on the specification arguments to check its validity. lets see one example: 

```julia-repl
julia> a = state(t=300.0,ρ=5u"mol/m^3")
ThermodynamicState with 2 properties:
  Temperature : 300.0[K]
  Molar density : 5[mol m^-3]
 ```
In this case, neither a phase nor any amount of compounds was specified. In those cases, the function assumes when the system has one mol, one phase, and/or a single component.

Another way to build a `ThermodynamicState` struct is by directly passing `Spec`s as arguments:
```julia-repl
 julia>
 h0 = spec(mol_h = 3000.0)
 t0 = spec(t=401.0)
 st = state(h0,t0)
julia>  st = state(h0,t0)
ThermodynamicState with 2 properties:
  Molar enthalpy : 3000.0[J mol^-1]
  Temperature : 401.0[K]
```
You can skip the check of the Gibbs' Phase Rule using the `check = false` keyword, or, in the case of building the state using keywords, decide to not normalize units via the 
`normalize_units= false` keyword.

## Obtaining properties from a `ThermodynamicState` struct

As said in the Basics section, you can query properties from the created `ThermodynamicState` struct. the way of obtaining those properties is by calling a property function. The interface proposed by this package is the following:

```julia
prop = property(model::MyModel,state::ThermodynamicState,unit::Unitful.unit=[default],args...)

```
Dispatching on the model type, you can calculate properties from the specifications contained in `state`.
This package exports one single model: `FromState`, that doesn't calculate (almost) anything, just checks if the selected property is in the argument, and if it is, it returns its numerical (unit stripped) value.
By default, the units of the number returned correspond to SI units, you can obtain a number with appropiate units by passing a corresponding `unit` argument.
the `FromState` model, given a corresponding molecular weight `mw` argument, can calculate derived properties, for example:

```julia-repl
mw = 50
a = state(t=300.0,ρ=5.0u"mol/L")
v = mass_volume(FromSpecs(),a,u"cm^3/g",mw)
4.0
```
The `mw` argument can be, depending of the situation, a number of vector of number. if not units are provided, a default units of `g/mol` are assumed. this conversion can be done on any unit that accepts molar, total and mass specifications, mass and moles themselves, molar and mass fractions,and molar and mass numbers.

## Variable `ThermodynamicState`s
Sometimes a more direct approach is needed when properties.For example, you may want to create a function that accepts only temperature to pass it to an ODE system or an optimization system. for this purpose, the Singleton `VariableSpec` is provided. if you pass it to a `Spec` or create a `ThermodynamicState` object, the resulting `state` or `spec` will be callable:

```julia-repl
julia> a_t = state(t=VariableSpec(),ρ=5.0u"mol/L")
ThermodynamicState(x₁) with 2 properties:
  Temperature : x₁
  Molar density : 5000.0[mol m^-3]
 
julia> a_t(300.15)
ThermodynamicState with 2 properties:
  Temperature : 300.15[K]
  Molar density : 5000.0[mol m^-3]

```
up to 3 `VariableSpec` can be added to each `ThermodynamicState`. a spec can only have one `VariableSpec` 

Its important to note that the functor created accepts the arguments in the same order as the arguments passed to the state function. for example,  `st = state(t=VariableSpec(),P=VariableSpec())` will return a functor of the form `st(t,p)` whereas `st = state(P=VariableSpec(),t=VariableSpec())` will create a functor with reversed order: `st(p,t)`

## Dispatching on the state type with `state_type`

Good. we now have a struct designed to store thermodynamic properties. now we can create functions that dispatch on a specific combination of thermodynamic specifications, using the function `state_type(st::ThermodynamicState)`:

```julia-repl
julia> st = state(t=373.15,ρ=5.0u"mol/L");state_type(st)
(VolumeAmount{MOLAR,DENSITY}(), Temperature(), SingleComponent())        
```
This tuple of thermodynamic specifications are an ordered representation of the specifications contained:
```julia-repl
st1 = state(t=373.15,ρ=5.0u"mol/L")
st2 = state(ρ=5.0u"mol/m^3",t=373.15,normalize_units=false)
state_type(s1) === state_type(s2) #true
```
Some abstract tuple types are saved on `ThermoState.QuickStates`. the tuple types exported are:

- `SinglePT,MultiPT`
- `SingleVT,MultiVT`
- `SinglePV,MultiPV`
- `SinglePS,MultiPS`
- `SinglePH,MultiPH`
- `SingleSatT,MultiSatT` (Two phase equilibrium)
- `SingleSatP,MultiSatP`
- `SingleΦT,MultiΦT`  (vapor fraction, general)
- `SingleΦP,MultiΦP`
- `SingleΦmT,MultiΦmT` (mass vapor fraction, or quality)
- `SingleΦmP,MultiΦmP`
- `SingleΦnT,MultiΦnT`(molar vapor fraction)
- `SingleΦnP,MultiΦnP`

## Exported utilities

## `@to_units`

Sometimes, `Unitful` quantities are needed. by default (and convention) all property accessors return a number without units, in SI system. You can automatically prefix the correct unit to the accessor function adding the `@to_units` macro at the start of the expression

```julia
st = state(t=1.0u"K",p=2.0u"Pa")
t0 = temperature(FromState(),st) #returns 1.0
t1 = @to_units temperature(FromState(),st) #returns 1.0 K
t1 = @to_units temperature(FromState(),st)
```

## `@spec_str`

Calling `spec(;key = value,normalize_units=true)` is simple, but it has a runtime cost. on the other part, `spec(sp::AbstractSpec,value,normalize_units::Bool=true)` can be determined at compile time. A problem with this interface is that writing the correct type can be cumbersome, for example, for some molar numbers:

```julia
n0 = spec(MaterialCompounds{MOLAR,TOTAL_AMOUNT}(),[0.1,0.3]) #very long type declaration
```

The `@spec_str` helps in this situation, creating the `AbstractSpec` type corresponding to the input keyword:

```julia
n1 = spec(spec"n",[0.1,0.3]) #shorter
n0 == n1 #true
```
The statement `spec(spec"n",[0.1,0.3])` can be defined at compile time.

### `normalize_units(val)`

On normal numbers, it is the identity, but on numbers or vectors of `Unitful.Quantity`,it converts the unit to an equivalent SI unit and strips the unit information.

```julia
x = 0.0u"°C"
normalize_units(x) #273.15
```

### `convert_unit(from,to,val)`

Converts an unit from the unit stored in `from` to the unit stored in `to`. when both units are equal, it justs returns `val`. if `val` itself is an `unit`, then it convert the from the unit in `val` to the unit in `to`. 

```julia
convert_unit(u"Pa",u"kPa",1000.0) #1.0
convert_unit(u"Pa",u"kPa",1u"atm") #101.325
```

###  `default_units(val)`

Returns the SI unit of a thermodynamic specification type or a function name corresponding to those types:

```julia
#from a function
default_units(mol_density) #u"mol/m^3"

#from a thermodynamic specification type
default_units(Pressure()) #u"Pa"
```

## Implementing a model using the `ThermoState` interface

Using this package, we can implement a basic ideal gas model that only calculates the pressure, given a temperature and molar volume, using the relation `pvₙ=Rt`:

```julia
using ThermoState, Unitful, 
using ThermoState.QuickStates #provides SingleVT for dispatch on state_type
import ThermoState: pressure #import all functions to overload, if you have a custom property, this is not necessary.

import ThermoState: molecular_weight #you can use your own name in your files, but it is recommended to use this to better interop between packages.

struct MyIdealGas
    mw::Float64
end
#t
molecular_weight(model::MyIdealGas) = model.mw

#your implementation of pressure, with v and t. 
#pass a state type to dispatch on the available properties.
function pressure_impl(mt::SingleVT,model::MyIdealGas,v,t)
    return 8.314*t/v
end

function pressure(model::MyIdealGas,st::ThermodynamicState,unit=u"Pa")
return pressure(state_type(st),model,st,unit)
end

function pressure(mt::SingleVT,model::MyIdealGas,st::ThermodynamicState,unit)
    v = mol_volume(FromState(),st,u"m^3/mol",mw)
    t = temperature(FromState(),st,u"K") 
    val = pressure_impl(mt,model,v,t)
    return convert_unit(u"Pa",unit,val)
end

a = state(mass = 3u"kg",total_v = 30u"m^3",t=30u"°C")
model = MyIdealGas(18.01)
p = pressure(model,a)

```
Using a variable state:

```julia
tx = state(mass = 3u"kg",total_v = 30u"m^3",t=var) #one free variable
p = pressure(model,tx(30u"°C"))
p_list = map(t-> pressure(model,tx(t)),273.0:373.0)
```

Here the actual function that does the work has the form `$property_impl(mt,model,args...)`, that accepts positional arguments of the form indicated by the result obtained by `state_type(st)`. in this example, it doesnt seem too useful, as there is only one posible implementation: `SingleVT`. However, this  helps when there is more than one posible pair of input to the function. For example, with this dispatch, you can do this:

```julia

struct WaterModel end

...

#Enthalpy from  Pressure - Temperature 
function mol_enthalpy_impl(::SinglePT,::WaterModel,p,t)
...
end

#Enthalpy from  Volume - Temperature 
function mol_enthalpy_impl(::SingleVT,::WaterModel,v,t)
...
end

#Enthalpy from  Entropy - Temperature 
function mol_enthalpy_impl(::SingleST,::WaterModel,s,t)
...
end


#Enthalpy from  Pressure - Entropy 
function mol_enthalpy_impl(::SinglePS,::WaterModel,p,s)
...
end

#ThermoState interface, only needed once
function mol_enthalpy(model::WaterModel,st::ThermodynamicState,unit="J/mol")
  return mol_enthalpy(state_type(st),model,st,unit)
end

#dispatch on the apropiate implementation, extracting the args from st
function mol_enthalpy(mt::SingleVT,model::WaterModel,st::ThermodynamicState,unit)
    v = mol_volume(FromState(),st,u"m^3/mol",mw)
    t = temperature(FromState(),st,u"K") 
    val = mol_enthalpy(mt,model,v,t)
    return convert_unit(u"J/mol",unit,val)
end

```

You can also provide automatic conversion to mass an total units, defining the following functions:

```julia

function mass_enthalpy(model::WaterModel,st::ThermodynamicState,unit="J/kg")
  mol_h = mol_enthalpy(model,st) #we dont care about state_type, just the result
  #to obtain mass_h, we need to divide by kg/mol
  kg_per_mol = molar_mass(FromState(),st,u"kg/mol",molecular_weight(model)) #we can do this or overload molar_mass(WaterModel,st)

  mass_h = mol_h/kg_per_mol
  return convert_unit(u"J/kg",unit,mass_h)
end

function total_enthalpy(model::WaterModel,st::ThermodynamicState,unit="J")
  mol_h = mol_enthalpy(model,st) 
  #to obtain total_h, we need to divide by moles
  mol = moles(FromState(),st,u"mol",molecular_weight(model)) #we can do this or overload moles(WaterModel,st)
  total_h = mol_h*mol
  return convert_unit(u"J/kg",unit,total_h)
end
```
As seen here, this is an easy, but repetitive and boring task. Thankfully, julia metaprogramming helps a lot here. if you have a lot of properties, you can evaluate all those functions at once this `@eval`. This is an actual piece of code used on WaterIF97.jl (not published yet). the implementation functions where defined before hand, whereas the `ThermoState` interface is mostly defined in this `@eval` loop

```julia
for op in [:helmholtz, :gibbs, :internal_energy, :enthalpy,:cp,:cv,:volume,:entropy]
        mol_op_impl = Symbol(:mol_,op,:_impl)
        mass_op_impl = Symbol(:mass_,op,:_impl)
        total_op_impl = Symbol(:total_,op,:_impl)
        mol_op = Symbol(:mol_,op)
        mass_op = Symbol(:mass_,op)
        total_op = Symbol(:total_,op)
        if op == :volume
            _unit = u"m^3/kg"
            mol_unit = u"m^3/mol"
            total_unit = u"m^3"

        elseif op in (:cv,:cp,:entropy)
            _unit = u"J/(kg*K)"
            mol_unit = u"J/(mol*K)"
            total_unit = u"J/(K)"
        else
            _unit = u"J/(kg)"
            mol_unit = u"J/(mol)"
            total_unit = u"J"
        end
     
        @eval begin

            #dispatch basic mass impl to the correct region, P,T
            function $mass_op_impl(mt::SinglePT,model::IndustrialWater,p,t)
                id = region_id(mt,model,p,t)
                return $mass_op_impl(mt,IF97Region{id}(),p,t)
            end


            function $mass_op(model::IndustrialWater,st::ThermodynamicState,unit=$_unit)
                return $mass_op(state_type(st),model,st,unit)
            end
            # P T impl
            function $mass_op(mt::SinglePT,model::IndustrialWater,st::ThermodynamicState,unit)
                p = pressure(FromState(),st)
                t = temperature(FromState(),st)
                res = $mass_op_impl(mt,model,p,t)
                return convert_unit($_unit,unit,res)
            end

            #mol op
            function $mol_op(model::IndustrialWater,st::ThermodynamicState,unit=$mol_unit)
                prod = molar_mass(FromState(),st,u"kg/mol",molecular_weight(model))
                res =  $mass_op(state_type(st),model,st,$_unit)*prod
                return convert_unit($mol_unit,unit,res)
            end 
        end

        if !(op in (:cv,:cp))
            #total ops, cv and cp dont have total operations
            @eval begin
                function $total_op(model::IndustrialWater,st::ThermodynamicState,unit=$total_unit)
                    prod = mass(FromState(),st,u"kg",molecular_weight(model))
                    res =  $mass_op(mt,model,st,unit)*prod
                    return convert_unit($total_unit,unit,res)
                end
            end
        end

    if op != :enthalpy
        @eval begin
            #if not enthalpy, define PH impl
            function $mass_op_impl(mt::SinglePH,model::IndustrialWater,p,h)
                id = region_id(mt,model,p,h)
                return $mass_op_impl(mt,IF97Region{id}(),p,h)
            end

            function $mass_op(mt::SinglePH,model::IndustrialWater,st::ThermodynamicState,unit)
                p = pressure(FromState(),st)
                h = mass_enthalpy(FromState(),st)
                res = $mass_op_impl(mt,model,p,h)
                return convert_unit($_unit,unit,res)
            end
        end   
    end

    if op != :entropy
        @eval begin
            #if not entropy, define PS impl
            function $mass_op_impl(mt::SinglePS,model::IndustrialWater,p,s)
                id = region_id(mt,model,p,h)
                t = temperature_impl(mt,IF97Region{id}(),p,s) #transform to SinglePT
                _mt = QuickStates.pt()
                return $mass_op_impl(_mt,model,p,t)
            end

            function $mass_op(mt::SinglePS,model::IndustrialWater,st::ThermodynamicState,unit)
                p = pressure(FromState(),st)
                s = mass_entropy(FromState(),st)
                res = $mass_op_impl(mt,model,p,s)
                return convert_unit($_unit,unit,res)
            end
        end   
    end
end
```
With this loop, we defined at 22 property accessor functions, in molar, mass an total units, with 3 implementations each (except enthalpy and entropy, with 2 implementations).







