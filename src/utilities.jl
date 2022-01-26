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
function normalize_units(x::Unitful.Quantity{T}) where T
    x::T =  Unitful.ustrip(Unitful.upreferred(x))
end

function normalize_units(x::AbstractVector{T}) where T<:Unitful.Quantity{Q} where Q
    return Unitful.ustrip.(Unitful.upreferred.(x))
end

normalize_units(x::Number) = x
normalize_units(x::AbstractVector) = x


#upreferred, but the standard unit with just numbers is transformed to kg/mol
function mw_mul(x,mw::Unitful.Quantity)
    return normalize_units(x*mw)
end

function mw_mul(x::Unitful.Quantity,mw::Unitful.Quantity)
    return normalize_units(x*mw)
end

function mw_mul(x,mw)
    return 0.001*x*mw
end

function mw_div(x,mw::Unitful.Quantity)
    return normalize_units(x/mw)
end

function mw_div(x::Unitful.Quantity,mw::Unitful.Quantity)
    return normalize_units(x/mw)
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

is_real(x) = false
is_real(x::Real) = true
is_real(x::Bool) = false
is_real(x::AbstractVector{T} where T<:Real) = true 
is_unitful(x) = false
is_unitful(x::Unitful.Quantity) = true
function is_unitful(x::T) where T <: AbstractArray{T2} where T2 <: Unitful.Quantity
    return true 
end
