struct FromSpecs end

#unified ustrip uconvert
function _ucs(u,x,normalize_units=false)
    if normalize_units
        return Unitful.ustrip(Unitful.uconvert(u,x))
    else
        return x
    end
end

function _ucs(u,x::AbstractVector,normalize_units=false)
    if normalize_units
        return Unitful.ustrip.(Unitful.uconvert.(u,x))
    else
        return x
    end
end

function _default_units(x::T,is_mol,is_total,inverted) where T <: AbstractSpec
    if !is_mol & is_total & !inverted#total units
        return total_units(x)
    elseif is_mol & !is_total & !inverted
        return mol_units(x)
    elseif !is_mol & is_total & !inverted
        return mass_units(x)
    elseif !is_mol & is_total & inverted #total units
        return inv(total_units(x))
    elseif is_mol & !is_total & inverted
        return inv(mol_units(x))
    elseif !is_mol & is_total & inverted
        return inv(mass_units(x))
    end
end

function conform1(s::Spec{T,U},is_mol,is_total,inverted) where {T,U}
    val = value(s)
    _spec = specification(s)
    if s.inverted != inverted
        val = 1/val
    end
    if (s.is_mol==is_mol) &(s.is_total==is_total)
        return _ups(val,true)
    else
        return val,false
    end
end

function throw_get_spec(val,specs)
    res = get_spec(val,specs)
    if res !== nothing
        return res
    else
        throw(error(string(val) * " not found in specifications"))
    end
end



function pressure(model::FromSpecs,props::Specs,unit::T=u"Pa") where T <: Unitful.PressureUnits
    spec_p = throw_get_spec(Pressure(),props.specs)
    p,_ = conform1(spec_p,false,true,false)
    if unit !== u"Pa"
        default_unit = _ucs(unit/_default_units(Pressure(),false,true,false),one(p),true)
        return default_unit*p 
    else
        return p
    end
end

function temperature(model::FromSpecs,props::Specs,unit::T=u"K") where T <: Unitful.TemperatureUnits
    spec_t = throw_get_spec(Temperature(),props.specs)
    t,_ = conform1(spec_t,false,true,false)
    if unit !== u"K"
        default_unit = _ucs(unit/_default_units(Temperature(),false,true,false),one(t),true)
        return default_unit*t
    else
        return t
    end
end
