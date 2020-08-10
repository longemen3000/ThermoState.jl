##TODO: CHECK UNITS

const MolarEnergyUnits = Unitful.Units{U,(Unitful.𝐍^-1)*Unitful.𝐌*Unitful.𝐋^2/Unitful.𝐓^2,A} where A where U
const MassEnergyUnits = Unitful.Units{U,Unitful.𝐋^2/Unitful.𝐓^2,A} where A where U
const MolUnits = Unitful.Units{U,Unitful.𝐍,A} where A where U
const INTENSIVE_ENERGY_UNITS = (Helmholtz,Gibbs,InternalEnergy,Enthalpy) 

#Entropy
const EntropyUnits = Unitful.Units{U,(Unitful.𝚯^-1)*Unitful.𝐌*Unitful.𝐋^2/Unitful.𝐓^2,A} where A where U
const MassEntropyUnits = Unitful.Units{U,(Unitful.𝚯^-1)*Unitful.𝐋^2/Unitful.𝐓^2,A} where A where U
const MolEntropyUnits = Unitful.Units{U,(Unitful.𝐍^-1)*(Unitful.𝚯^-1)*Unitful.𝐌*Unitful.𝐋^2/Unitful.𝐓^2,A} where A where U

#volume and density units
const MolDensityUnits = Unitful.MolarityUnits
const MassDensityUnits = Unitful.MolalityUnits
const MassVolumeUnits = Unitful.Units{U,((Unitful.𝐋)^3)/(Unitful.𝐌),A} where A where U
const MolVolumeUnits = Unitful.Units{U,((Unitful.𝐋)^3)/(Unitful.𝐍),A} where A where U


function pressure(model::FromSpecs,props::Specs,unit::T=u"Pa") where T <: Unitful.PressureUnits
    sval = throw_get_spec(Pressure(),props)
    val = to_spec(props,sval,nothing,Pressure())
    if unit !== u"Pa"
        default_unit = _ups(one(val)*u"Pa"/unit,true)
        return default_unit*val 
    else
        return val
    end
end

function temperature(model::FromSpecs,props::Specs,unit::T=u"K") where T <: Unitful.TemperatureUnits
    sval = throw_get_spec(Temperature(),props)
    val = to_spec(props,sval,nothing,Temperature())
    if unit !== u"K"
        default_unit = _ups(one(val)*u"K"/unit,true)
        return default_unit*val 
    else
        return val
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


for (op,sp) in zip((:mol_helmoltz, :mol_gibbs, :mol_internal_energy, :mol_enthalpy),INTENSIVE_ENERGY_UNITS)
    @eval begin 
        function $op(model::FromSpecs,props::Specs,unit::T=u"J/mol",mw=nothing) where T <: MolarEnergyUnits
            sval = throw_get_spec($sp,props)
            val = to_spec(props,sval,nothing,MOLAR())
            if unit !== u"J/mol"
                default_unit = _ups(one(val)*u"J/mol"/unit,true)
                return default_unit*val 
            else
                return val
            end
        end
    end
end

for (op,sp) in zip((:mass_helmoltz, :mass_gibbs, :mass_internal_energy, :mass_enthalpy),INTENSIVE_ENERGY_UNITS)
    @eval begin 
        function $op(model::FromSpecs,props::Specs,unit::T=u"J/kg",mw=nothing) where T <: MassEnergyUnits
            sval = throw_get_spec($sp,props)
            val = to_spec(props,sval,nothing,MASS())
            if unit !== u"J/kg"
                default_unit = _ups(one(val)*u"J/kg"/unit,true)
                return default_unit*val 
            else
                return val
            end
        end    
    end
end

for (op,sp) in zip((:total_helmoltz, :total_gibbs, :total_internal_energy, :total_enthalpy),INTENSIVE_ENERGY_UNITS)
    @eval begin 
            function $op(model::FromSpecs,props::Specs,unit::T=u"J",mw=nothing) where T <: Unitful.EnergyUnits
                sval = throw_get_spec($sp,props)
                val = to_spec(props,sval,nothing,TOTAL())
                if unit !== u"J"
                    default_unit = _ups(one(val)*u"J"/unit,true)
                    return default_unit*val 
                else
                    return val
                end
            end
    end
end

function mol_entropy(model::FromSpecs,props::Specs,unit::T=u"J/(K*mol)",mw=nothing) where T <: MolEntropyUnits
    sval = throw_get_spec(Entropy,props)
    val = to_spec(sval,props,mw,MOLAR())
    if unit !== u"J/(K*mol)"
        default_unit = _ups(one(val)*u"J/(K*mol)"/unit,true)
        return default_unit*val
    else
        return val
    end
end

function mass_entropy(model::FromSpecs,props::Specs,unit::T=u"J/(K*kg)",mw=nothing) where T <: MassEntropyUnits
    sval = throw_get_spec(Entropy,props)
    val = to_spec(sval,props,mw,MASS())
    if unit !== u"J/(K*kg)"
        default_unit = _ups(one(val)*u"J/(K*kg)"/unit,true)
        return default_unit*val
    else
        return val
    end
end

function total_entropy(model::FromSpecs,props::Specs,unit::T=u"J/(K)",mw=nothing) where T <: EntropyUnits
    sval = throw_get_spec(Entropy,props)
    val = to_spec(sval,props,mw,MASS())
    if unit !== u"J/K"
        default_unit = _ups(one(val)*u"J/(K)"/unit,true)
        return default_unit*val
    else
        return val
    end
end

function total_volume(model::FromSpecs,props::Specs,unit::T=u"m^3",mw=nothing) where T <: Unitful.VolumeUnits
    sval = throw_get_spec(VolumeAmount,props)
    val = to_spec(sval,props,mw,VolumeAmount{TOTAL,VOLUME}())
    if unit !== u"m^3"
        default_unit = _ups(one(val)*u"m^3"/unit,true)
        return default_unit*val
    else
        return val
    end
end

function mass_volume(model::FromSpecs,props::Specs,unit::T=u"(m^3)/kg",mw=nothing) where T <: MassVolumeUnits
    sval = throw_get_spec(VolumeAmount,props)
    val = to_spec(sval,props,mw,VolumeAmount{MASS,VOLUME}())
    if unit !== u"(m^3)/kg"
        default_unit = _ups(one(val)*u"(m^3)/kg"/unit,true)
        return default_unit*val
    else
        return val
    end
end

function mol_volume(model::FromSpecs,props::Specs,unit::T=u"(m^3)/mol",mw=nothing) where T <: MolVolumeUnits
    sval = throw_get_spec(VolumeAmount,props)
    val = to_spec(sval,props,mw,VolumeAmount{MOL,VOLUME}())
    if unit !== u"(m^3)/mol"
        default_unit = _ups(one(val)*u"(m^3)/mol"/unit,true)
        return default_unit*val
    else
        return val
    end
end

function mass_density(model::FromSpecs,props::Specs,unit::T=u"kg/m^3",mw=nothing) where T <: MassDensityUnits
    sval = throw_get_spec(VolumeAmount,props)
    val = to_spec(sval,props,mw,VolumeAmount{MASS,DENSITY}())
    if unit !== u"kg/m^3"
        default_unit = _ups(one(val)*u"kg/m^3"/unit,true)
        return default_unit*val
    else
        return val
    end
end

function mol_density(model::FromSpecs,props::Specs,unit::T=u"mol/m^3",mw=nothing) where T <: MolDensityUnits
    sval = throw_get_spec(VolumeAmount,props)
    val = to_spec(sval,props,mw,VolumeAmount{MOLAR,DENSITY}())
    if unit !== u"mol/m^3"
        default_unit = _ups(one(val)*u"mol/m^3"/unit,true)
        return default_unit*val
    else
        return val
    end
end





