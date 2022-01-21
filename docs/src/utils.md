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
