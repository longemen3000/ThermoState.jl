"""
    normalize_units(val)
On normal numbers, it is the identity, but on numbers or vectors of `Unitful.Quantity`,it converts the unit to an equivalent SI unit and strips the unit information.

# Examples
```julia-repl
julia> normalize_units(0.0u"°C")
273.15

```
```julia-repl
julia> normalize_units(273.15)
273.15
```
"""
normalize_units(x) =  Unitful.ustrip(Unitful.upreferred(x))
normalize_units(x::AbstractVector) =  Unitful.ustrip.(Unitful.upreferred.(x))


#upreferred, but the standard unit with just numbers is transformed to kg/mol
function mw_mul(x,mw::Unitful.Quantity)
    return upreferred(x*mw)
end

function mw_mul(x,mw)
    return 0.001*x*mw
end

function mw_div(x,mw::Unitful.Quantity)
    return upreferred(x/mw)
end

function mw_div(x,mw)
    return 1000.0*x/mw
end

"""
    convert_unit(from,to,val)
Converts an unit from the unit stored in `from` to the unit stored in `to`. 

When both units are equal, it justs returns `val`. 

If `val` itself is an `unit`, then it converts the from the unit in `val` to the unit in `to`. 

# Examples

```julia-repl
julia> convert_unit(u"Pa",u"kPa",1000.0)      
1.0
```

```julia-repl
julia> convert_unit(u"Pa",u"kPa",1u"atm")     
4053//40
```
"""
function convert_unit end
function convert_unit(from::T,to::T,val::N) where {T,N}
    return val
end

function convert_unit(from::T1,to::T2,val::N) where {T1,T2,N<:Number}
    return Unitful.ustrip(Unitful.uconvert(to,val*from))
end

function convert_unit(from::T1,to::T2,val::N) where {T1,T2,N<:AbstractVector{<: Number}}
    return Unitful.ustrip.(Unitful.uconvert.(to,val*from))
end

function convert_unit(from::T1,to::T2,val::N) where {T1,T2,N<:AbstractVector{<:Unitful.Quantity}}
    return Unitful.ustrip.(Unitful.uconvert.(to,val))
end
