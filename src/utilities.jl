"""
    normalize_units(x)

For normal numbers, this is the identity function.
 For `Unitful` quantities, it converts to SI units and strips the `Unitful` type.

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

function convert_unit(from::T,to::T,val::N) where {T,N<:Number}
    return val
end

function convert_unit(from::T1,to::T2,val::N) where {T1,T2,N<:Number}
    return Unitful.ustrip(Unitful.uconvert(to,val*from))
end

function convert_unit(from::T1,to::T2,val::N) where {T1,T2,N<:Unitful.Quantity}
    return Unitful.ustrip(Unitful.uconvert(to,val))
end
