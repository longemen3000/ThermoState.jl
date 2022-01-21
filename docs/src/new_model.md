# Implementing a model using the `ThermoState` interface

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
