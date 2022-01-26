## Quick Start

Lets suppose a mix of water (18.01 `g/mol`) and ethanol (46.07 `g/mol`)

```julia
using Unitful,ThermoState
#m = mass amount 
#n = molar amount
#xn = molar fraction
#xm = mass fraction
mw = [18.01,46.07]
st = state(t = 300u"°C",mass_rho=805u"g/L",n =[5.1,4.2]u"mol")
```
- What's is the molecular weight of the mix?

```julia-repl
julia> molar_mass(FromState(),st,mw)
0.03068225806451613
```
The result is in `kg/mol`, we could use `g/mol` instead:
```julia-repl
julia> molar_mass(FromState(),st,u"g/mol",mw)
30.68225806451613
```


- How many moles do we have in our mix?

```julia-repl
julia> moles(FromState(),st,mw)
9.3
```

For this specific combination of properties and state, we don't need a molecular weight vector:
```julia-repl
julia> moles(FromState(),st)
9.3
```

- What's the volume of the mix in `inch³`
```julia-repl
julia> total_volume(FromState(),st,u"inch^3",mw)
21.63083261951725
```


