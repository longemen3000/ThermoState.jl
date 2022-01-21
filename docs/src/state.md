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
- `SingleT,MultiT` (a state that has a temperature)
