struct FromSpecs end



function _default_units(x::T,is_mol,is_total,inverted) where T <: AbstractSpec

    if !is_mol & is_total & !inverted#total units
        return total_units(x)
    elseif is_mol & !is_total & !inverted
        return mol_units(x)
    elseif !is_mol & is_total & !inverted
        return mass_units(x)
    elseif !is_mol & is_total & inverted #total units
        return 1/total_units(x)
    elseif is_mol & !is_total & inverted
        return 1/mol_units(x)
    elseif !is_mol & is_total & inverted
        return 1/mass_units(x)
    end
end

function _conform1(s::Spec{T,U},is_mol,is_total,inverted) where {T,U}
    val = value(s)
    _spec = specification(s)
    if s.inverted != inverted
        val = 1/val
    end
    if (s.is_mol==is_mol) &(s.is_total==is_total)
        if !(U<:Unitful.Quantity)
            return val,true
        else
            u = _default_units(_spec,false,true,false)
            if unit(u) == unit(val)
                return Unitful.ustrip(val)
            else
                return Unitful.ustrip(Unitful.uconvert(u,val))
            end
        end
    else
        return val,false
    end
end

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

function uconvert_and_strip2(value,unit,default)
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
    spec_p = spec_in(Pressure(),props.specs)
    p,_ = _conform1(spec_p,false,true,false)
    if unit !== nothing
        default_unit = unit/_default_units(Pressure(),false,true,false)
        return Unitful.ustrip(Unitful.uconvert(default_unit,p))
    else
        return p
    end
end

function temperature(model::FromSpecs,props::Specs;mw=nothing,unit=nothing)
    spec_t = spec_in(Temperature(),props.specs)
    t,_ = _conform1(spec_t,false,true,false)
    if unit !== nothing
        default_unit = unit/_default_units(Temperature(),false,true,false)
        return Unitful.ustrip(Unitful.uconvert(default_unit,t))
    else
        return t
    end
end


