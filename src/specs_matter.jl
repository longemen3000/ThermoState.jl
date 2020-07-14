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
    return 1000*x/mw
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
            if amount.is_mol #mol specified
                return _ups(mw_mul(value(amount),mw),true)
            else #mass specified
                return _ups(value(amount),true)
            end
        else #multicomponent fractions
            compounds = get_spec(MaterialCompounds(),specs)
            if !compounds.is_total #only this case is valid, it cant be total compounds and a amount
                if amount.is_mol #mol specified
                    ###### needs mass
                    ### exists mol fractions, moles
                    #mw_i * xi = (gc/molc) * molc/molmix = gc/molmix (*molmix)
                    sum_x_mw = mapreduce(mw_mul,+,value(compounds),mw)
                    return _ups(value(amount)*sum_x_mw,true)
                else #mass specified
                    return _ups(value(amount),true)
                end
            end
        end
    else 
        if !_compounds #monocomponent, one mol
            return _ups(mw,true)
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
        else #multicomponent fractions
            compounds = get_spec(MaterialCompounds(),specs)
            if !compounds.is_total #only this case is valid, it cant be total compounds and a amount
                if amount.is_mol #mol specified
                    return _ups(value(amount),true)
                else #mass specified
                    ###### needs mol
                    ### exists mass fractions, moles
                    #mw_i * xmi = gc/gmix * (molc/gc)= molc/gmix (*gmix)
                    sum_x_mw = mapreduce(mw_div,+,value(compounds),mw)
                    return _ups(value(amount)*sum_x_mw,true)
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

function mass(model::FromSpecs,props::Specs,unit::T=u"kg",mw) where T <: Unitful.MassUnits
    m = mass2(props,mw)
    if unit !== u"kg"
        default_unit = _ucs(unit/_default_units(Mass(),true,false,false),one(m),true)
        return default_unit*m
    else
        return m
    end
end

function moles(model::FromSpecs,props::Specs,unit::T=u"kg",mw) where T <: Unitful.MassUnits
    m = moles2(props,mw)
    if unit !== u"mol"
        default_unit = _ucs(unit/_default_units(Moles(),true,false,false),one(m),true)
        return default_unit*m
    else
        return m
    end
end