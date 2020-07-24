##TODO: CHECK UNITS

const MolarEnergyUnits = Unitful.Units{U,(Unitful.𝐍^-1)*Unitful.𝐌*Unitful.𝐋^2/Unitful.𝐓^2,A} where A where U
const MassEnergyUnits = Unitful.Units{U,Unitful.𝐋^2/Unitful.𝐓^2,A} where A where U
const MolUnits = Unitful.Units{U,Unitful.𝐍,A} where A where U
const INTENSIVE_UNITS = (Helmholtz(),Gibbs(),InternalEnergy(),Enthalpy()) 

#Entropy
const EntropyUnits = Unitful.Units{U,(Unitful.𝚯^-1)*Unitful.𝐌*Unitful.𝐋^2/Unitful.𝐓^2,A} where A where U
const MassEntropyUnits = Unitful.Units{U,(Unitful.𝚯^-1)*Unitful.𝐋^2/Unitful.𝐓^2,A} where A where U
const MolEntropyUnits = Unitful.Units{U,(Unitful.𝐍^-1)*(Unitful.𝚯^-1)*Unitful.𝐌*Unitful.𝐋^2/Unitful.𝐓^2,A} where A where U



function pressure(model::FromSpecs,props::Specs,unit::T=u"Pa") where T <: Unitful.PressureUnits
    spec_p = throw_get_spec(Pressure(),props.specs)
    p = conform_pt(spec_p,false,true,false)
    if unit !== u"Pa"
        default_unit = _ucs(unit/_default_units(Pressure(),false,true,false),one(p),true)
        return default_unit*p 
    else
        return p
    end
end



function temperature(model::FromSpecs,props::Specs,unit::T=u"K") where T <: Unitful.TemperatureUnits
    spec_t = throw_get_spec(Temperature(),props.specs)
    t = conform_pt(spec_t,false,true,false)
    if unit !== u"K"
        default_unit = _ups(one(t)*u"K"/unit,true)
        return default_unit*t
    else
        return t
    end
end



function mass(model::FromSpecs,props::Specs,unit::T=u"kg",mw=nothing) where T <: Unitful.MassUnits
    m = mass2(props,mw)
    if unit !== u"kg"
        default_unit = _ups(one(m)*u"kg"/unit,true)
        return default_unit*m
    else
        return m
    end
end



function moles(model::FromSpecs,props::Specs,unit::T=u"mol",mw=nothing) where T <: MolUnits
    m = moles2(props,mw)
    if unit !== u"mol"
        default_unit = _ups(one(m)*u"mol"/unit,true)
        return default_unit*m
    else
        return m
    end
end


for (op,sp) in zip((:mol_helmoltz, :mol_gibbs, :mol_internal_energy, :mol_enthalpy),INTENSIVE_UNITS)
    @eval begin 
        function $op(model::FromSpecs,props::Specs,unit::T=u"J/mol",mw=nothing) where T <: MolarEnergyUnits
            sval = props[$sp]
            val = conform3(sval,props,false,true,false,mw)
            if unit !== u"J/mol"
                default_unit = _ups(one(val)*u"J/mol"/unit,true)
                return default_unit*val
            else
                return val
            end
        end   
    end
end

for (op,sp) in zip((:mass_helmoltz, :mass_gibbs, :mass_internal_energy, :mass_enthalpy),INTENSIVE_UNITS)
    @eval begin 
        function $op(model::FromSpecs,props::Specs,unit::T=u"J/kg",mw=nothing) where T <: MassEnergyUnits
            sval = props[$sp]
            val = conform3(sval,props,false,false,false,mw)
            if unit !== u"J/kg"
                default_unit = _ups(one(val)*u"J/kg"/unit,true)
                return default_unit*val
            else
                return val
            end
        end     
    end
end

for (op,sp) in zip((:total_helmoltz, :total_gibbs, :total_internal_energy, :total_enthalpy),INTENSIVE_UNITS)
    @eval begin 
            function $op(model::FromSpecs,props::Specs,unit::T=u"J",mw=nothing) where T <: Unitful.EnergyUnits
                sval = props[$sp]
                val = conform3(sval,props,true,false,false,mw)
                if unit !== u"J"
                    default_unit = _ups(one(val)*u"J"/unit,true)
                    return default_unit*val
                else
                    return val
                end
            end
    end
end

function mass_entropy(model::FromSpecs,props::Specs,unit::T=u"J/(K*mol)",mw=nothing) where T <: MolEntropyUnits
    sval = props[Entropy()]
    val = conform3(sval,props,false,true,false,mw)
    if unit !== u"J/(K*mol)"
        default_unit = _ups(one(val)*u"J/(K*mol)"/unit,true)
        return default_unit*val
    else
        return val
    end
end

function mass_entropy(model::FromSpecs,props::Specs,unit::T=u"J/(K*kg)",mw=nothing) where T <: MassEntropyUnits
    sval = props[Entropy()]
    val = conform3(sval,props,false,false,false,mw)
    if unit !== u"J/(K*kg)"
        default_unit = _ups(one(val)*u"J/(K*kg)"/unit,true)
        return default_unit*val
    else
        return val
    end
end

function total_entropy(model::FromSpecs,props::Specs,unit::T=u"J/(K)",mw=nothing) where T <: EntropyUnits
    sval = props[Entropy()]
    val = conform3(sval,props,true,false,false,mw)
    if unit !== u"J/K"
        default_unit = _ups(one(val)*u"J/(K)"/unit,true)
        return default_unit*val
    else
        return val
    end
end







