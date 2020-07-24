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

function _default_units(x::T,is_total,is_mol,inverted) where T <: AbstractSpec
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

#fast conform  for pressure and temperature
function conform_pt(s::Spec{T,U},is_total,is_mol,inverted) where {T,U}
    val = value(s)
    if (s.is_mol==is_mol) &(s.is_total==is_total)
        if s.inverted != inverted
            val = inv(val)
        end
        return _ups(val,true),true
    else
        throw(error("spec cannot be converted to the correct units"))
    end
end



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

function amount_type(specs::Specs)
    amount = has_spec(MaterialAmount(),specs)
    compounds = has_spec(MaterialCompounds(),specs)
    return amount,compounds
end

function mass2(specs,mw)
    _amount,_compounds = amount_type(specs)
    if _amount #amount is specified
        amount = get_spec(MaterialAmount(),specs)
        if !_compounds #monocomponent
            if amount.is_mol #monocomponent, mol
                return _ups(mw_mul(value(amount),mw),true)
            else #monocomponent, mass
                return _ups(value(amount),true)
            end
        else #multicomponent fractions
            compounds = get_spec(MaterialCompounds(),specs)
            if !compounds.is_total #only fractions
                if amount.is_mol #moles
                    _moles = value(amount)
                    if compounds.is_mol #moles, mol fractions
                        #mw_i * xi = (gc/molc) * molc/molmix = gc/molmix (*molmix)
                        sum_x_mw = mapreduce(mw_mul,+,value(compounds),mw)
                        return _ups(_moles*sum_x_mw,true)
                    else #moles, mass fractions
                        #mw_i / wi = (molc/gc) * gc/gmix = molmix/gmix 
                        #mass = molmix / (molmix/gmix)
                        sum_x_mw = mapreduce(mw_div,+,value(compounds),mw)
                        return _ups(_moles/sum_x_mw,true)
                    end
                else #mass specified
                    return _ups(value(amount),true)
                end
            end
        end
    else 
        if !_compounds #monocomponent, one mol
            return _ups(1.0,true)
        else #vector of compounds
            compounds = get_spec(MaterialCompounds(),specs)
            if compounds.is_total #only option
                if !compounds.is_mol
                    return _ups(sum(value(compounds)),true)
                else
               ###### needs mass
                    ### exists mol numbers
                    #mw_i * ni = (gc/molc) * molc = gc, sum(gc) = gmix
                    sum_n_mw = mapreduce(mw_mul,+,value(compounds),mw)
                    return _ups(sum_n_mw,true)
                end
            end
        end
    end
end

function moles2(specs,mw)
    _amount,_compounds = amount_type(specs)
    if _amount #amount is specified
        amount = get_spec(MaterialAmount(),specs)
        if !_compounds #monocomponent
            if amount.is_mol #mol specified
                return _ups(value(amount),true)
            else #mass specified
                return _ups(mw_div(value(amount),mw),true)
            end
        else #multicomponent system
            compounds = get_spec(MaterialCompounds(),specs)
            if !compounds.is_total #only fractions valid
                
                if amount.is_mol #mol specified
                    return _ups(value(amount),true)
                else #mass specified
                    _mass = value(amount)
                    if !compounds.is_mol #mass fractions and mass
                        #mw_i * xmi = gc/gmix / (molc/gc)= molc/gmix (sum)-> molmix/gmix
                        # molmix = molmix/gmix * (gmix)
                        sum_x_mw = mapreduce(mw_div,+,value(compounds),mw)
                        return _ups(_mass*sum_x_mw,true)
                    else #mol fraction and mass
                        #xni * mwi = mc/mmix * (gc/molc)= gc/molmix (sum)-> gmix/molmix
                        #molmix = gmix / (gmix/molmix)
                        sum_x_mw = mapreduce(mw_mul,+,value(compounds),mw)
                        return _ups(value(amount)/sum_x_mw,true)

                    end
                end
            end
        end
    else 
        if !_compounds #monocomponent, one mol
            return one(_ups(one(mw),true))
        else #vector of compounds
            compounds = get_spec(MaterialCompounds(),specs)
            if compounds.is_total #only option
                if !compounds.is_mol
                    ###### needs mol
                    ### exists mass numbers
                    #mw_i * xmi = gc* (molc/gc)= molc, molmix = sum(molc)
                    sum_n_mw = mapreduce(mw_div,+,value(compounds),mw)
                    return _ups(sum_n_mw,true)                
                else
                    _ups(sum(value(compounds)),true)
                end
            end
        end
    end
end

function kg_per_mol2(specs,mw)
    _amount,_compounds = amount_type(specs)
    if _amount #amount is specified
        amount = get_spec(MaterialAmount(),specs)
        if !_compounds #monocomponent
            return _ups(1.0,true)
        else #multicomponent system
            compounds = get_spec(MaterialCompounds(),specs)
            if !compounds.is_total #only fractions valid
                if compounds.is_mol & amount.is_mol #moles, mol fractions
                    _moles = value(amount)
                    #mw_i * xi = (gc/molc) * molc/molmix = gc/molmix sum_> gmix/molmix
                    sum_x_mw = mapreduce(mw_mul,+,value(compounds),mw)
                    return _ups(sum_x_mw,true)
                elseif (!compounds.is_mol) & amount.is_mol#moles, mass fractions
                    _moles = value(amount)
                    #mw_i / wi = (molc/gc) * gc/gmix sum-> molmix/gmix 
                    sum_x_mw = mapreduce(mw_div,+,value(compounds),mw)
                    return _ups(1.0/sum_x_mw,true)
                elseif  (!amount.is_mol) & (!compounds.is_mol) #mass fractions and mass
                    _mass = value(amount)
                    #mw_i * xmi = gc/gmix / (molc/gc)= molc/gmix (sum)-> molmix/gmix
                    sum_x_mw = mapreduce(mw_div,+,value(compounds),mw)
                    return _ups(1.0/sum_x_mw,true)
                elseif  (!amount.is_mol) & (compounds.is_mol) #mol fraction and mass
                    _mass = value(amount)
                    #xni * mwi = mc/mmix * (gc/molc)= gc/molmix (sum)-> gmix/molmix
                    #molmix = gmix / (gmix/molmix)
                    sum_x_mw = mapreduce(mw_mul,+,value(compounds),mw)
                    return _ups(sum_x_mw,true)
                end
            end
        end
    else #amount not specified
        if !_compounds #monocomponent, one mol
            return _ups(mw_mul(1.0,mw),true)
        else #multicomponent
            compounds = get_spec(MaterialCompounds(),specs)
            if compounds.is_total #only option
                if !compounds.is_mol #mw from mol numbers
                    #molc * gc/molc = gc sum-> gmix
                    gmix = mapreduce(mw_mul,+,value(compounds),mw)
                    molmix = sum(value(compounds))
                    return _ups(gmix/molmix,true)                
                else #mw from mass numbers
                    #gc / (gc/molc) = molc sum-> molmix
                    molmix = mapreduce(mw_div,+,value(compounds),mw)
                    gmix = sum(value(compounds))
                    return _ups(gmix/molmix,true)                
                end
            end
        end
    end
end

#conforms with transformations on mol <-> mass <-> total
function conform2(s::Spec{T,U},_specs,_is_total::Bool,_is_mol::Bool,_inverted::Bool,mw) where {T,U}
    val = value(s)
    if s.inverted #standard form
        val = inv(val)
    end
    
    if _is_mol && !_is_total #required molar units
        if is_total(s)
            return _ups(val/moles2(_specs,mw),true)
        else #mass
            return _ups(val/kg_per_mol2(_specs,mw),true)
        end
    elseif !_is_mol && !_is_total #required mass units
        if is_total(s)
            return _ups(val/mass2(_specs,mw),true)
        else
            return _ups(val*kg_per_mol2(_specs,mw),true)
        end
        
    elseif !_is_mol && _is_total #required total units
        if is_molar(s)
            return _ups(val*moles2(_specs,mw),true)
        else
            return _ups(val/mass2(_specs,mw),true)
        end
    end
end

#conforms with inversions, used in molar densities
function conform3(s::Spec{T,U},_specs,is_total::Bool,is_mol::Bool,inverted::Bool,mw) where {T,U}
    valc2 =  conform2(s,_specs,is_total,is_mol,inverted,mw)
    if inverted #requires inversion
        if !(s.inverted) #not inverted
            valc2 =  inv(valc2)
        end
        return valc2
    else #does not require inversion
        if s.inverted #but is inverted
            valc2 =  inv(valc2)
        end
        return valc2        
    end 
end

#unified conform for specs that require transformations
function conform0(s::Spec{T,U},_specs,is_total::Bool,is_mol::Bool,inverted::Bool,mw) where {T,U}
    val = value(s)
    if (s.is_mol==is_mol) &(s.is_total==is_total)
        if s.inverted != inverted
            val = inv(val)
        end
        return _ups(val,true)
    else
        return conform3(s,_specs,is_total,is_mol,inverted,mw)
    end
end

function mol_frac2(specs,mw)
    _amount,_compounds = amount_type(specs)
    if !_compounds 
        return [1.0]
    else
        elements = get_spec(MaterialCompounds(),specs)
        if elements.!is_total #a fraction
            if elements.is_mol #molar fraction
                return value(elements)
            else# mass fraction
            end
        else #mass or mol number
            if elements.is_mol #molar number
                n = sum(value(elements))
                return _ups(value(elements)/n,true)
            else# mass number
            end 
        end
    end
end

function mass_frac2(specs,mw)
    _amount,_compounds = amount_type(specs)
    if !_compounds 
        return [1.0]
    else
        elements = get_spec(MaterialCompounds(),specs)
        if elements.!is_total #a fraction
            if elements.is_mol #molar fraction
            else# mass fraction
                return value(elements)
            end
        else #mass or mol number
            if elements.is_mol #molar number
            else# mass number
                n = sum(value(elements))
                return _ups(value(elements)/n,true)
            end 
        end
    end
end

function mol_num2(specs,mw)
    _amount,_compounds = amount_type(specs)
    if !_compounds 
        return [1.0]
    else
        elements = get_spec(MaterialCompounds(),specs)
        if elements.!is_total #a fraction
            if elements.is_mol #molar fraction
                x = value(elements)
                n = value(get_spec(MaterialAmount(),specs))
                return _ups(x*n,true)
            else# mass fraction
            end
        else #mass or mol number
            if elements.is_mol #molar number
                return value(elements)
            else# mass number
            end 
        end
    end
end

function mass_num2(specs,mw)
    _amount,_compounds = amount_type(specs)
    if !_compounds 
        return [1.0]
    else
        elements = get_spec(MaterialCompounds(),specs)
        if elements.!is_total #a fraction
            if elements.is_mol #molar fraction
            else# mass fraction
                x = value(elements)
                n = value(get_spec(MaterialAmount(),specs))
                return _ups(x*n,true)
            end
        else #mass or mol number
            if elements.is_mol #molar number
            else# mass number
                return value(elements)
            end 
        end
    end
end