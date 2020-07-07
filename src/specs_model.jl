struct FromSpecs end


@inline function spec_in(val,specs)
    for spec in specs
        if specification(spec) == val
            return spec
        end
    end
    throw(error(string(val) * "not found in specifications"))
end


#requires no conversion
_uconvert1(unit,value) = value,true,true 
_uconvert(unit::Unitful.Quantity,value::Unitful.Quantity) = Unitful.ustrip(Unitful.uconvert(unit,value)),true,true

#requires conversion to default unit

#this requires uconvert(unit/default,value)
_uconvert1(unit::Unitful.Quantity,value) = value,false,false
#this requires uconvert(default,value)
_uconvert1(unit,value::Unitful.Quantity) = value,false,true

function uconvert_and_strip(value,unit,default)
    val,done,has_units = _uconvert1(unit,value)
    if done
        return val
    else
        if has_units
            return Unitful.ustrip(Unitful.uconvert(default,value))
        else
            return Unitful.ustrip(Unitful.uconvert(unit/default,value))
        end
    end
end
function pressure(model::FromSpecs,props::Specs;mw=nothing,unit=nothing)
    _p = spec_in(Pressure(),props.specs)
    return uconvert_and_strip(value(_p),unit,u"Pa")
end
