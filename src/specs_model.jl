struct FromSpecs end
const UnitReal = Union{Real,Unitful.Quantity}
const URVec = AbstractVector{T} where T<:UnitReal
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

function amount_type(specs::Tuple)
    return AMOUNT_CONST[_specs_components(specs)]
end

function amount_type(specs::Specs{Nothing}) 
    return amount_type(specs.specs)
end

function amount_type(specs::Specs{T}) where T<:Tuple
    amount,compounds = specs.amount_type
    return amount,compounds
end

moles2(a::Specs,mw) = moles3(amount_type(a),a,mw)
mass2(a::Specs,mw) =  mass3(amount_type(a),a,mw)
kg_per_mol2(a::Specs,mw) =  kg_per_mol3(amount_type(a),a,mw)

##one mol cases

#those cases should rarely be called
function moles3(::Tuple{T,OneMol},specs::Specs,mw::Nothing) where T
    return 1.0
end

function moles3(::Tuple{T,OneMol},specs::Specs,mw::T2) where T where T2 <:UnitReal
    return one(T2)
end

#error?
function moles3(::Tuple{T,OneMol},specs::Specs,mw::T2) where T where T2<:AbstractVector{T3} where T3 <: UnitReal
    return one(T3)
end

function moles3(::Tuple{SingleComponent,MaterialAmount{MOLAR}},specs::Specs,mw::T) where T <:UnitReal
    return _ups(value(get_spec(MaterialAmount{MOLAR}(),specs)),true)
end

function moles3(::Tuple{SingleComponent,MaterialAmount{MASS}},specs::Specs,mw::T) where T <:UnitReal
    # g/(g/mol) = mol 
    x = mw_div(value(get_spec(MaterialAmount{MASS}(),specs)),mw)
    return _ups(x,true)
end


## total ammounts:
function moles3(::Tuple{T,T},specs::Specs,mw) where T<:MaterialCompounds{MOLAR,TOTAL_AMOUNT}
    return _ups(sum(value(get_spec(T(),specs))),true)
end

function moles3(::Tuple{T,T},specs::Specs,mw::T2) where T<:MaterialCompounds{MASS,TOTAL_AMOUNT} where T2<: URVec
    ###### needs mol
    ### exists mass numbers
    #mw_i * xmi = gc* (molc/gc)= molc, molmix = sum(molc)
    compounds = value(get_spec(T(),specs))
    sum_n_mw = mapreduce(mw_div,+,compounds,mw)
    return _ups(sum_n_mw,true)  
end

function moles3(::Tuple{T,MaterialAmount{MOLAR}},specs::Specs,mw) where T
    return _ups(value(get_spec(MaterialAmount{MOLAR}(),specs)),true)
end

#inverse operations with fractions
function moles3(::Tuple{MaterialCompounds{MOLAR,FRACTION},MaterialAmount{MASS}},specs::Specs,mw::T) where T<:URVec
    #mol fraction and mass
    #xni * mwi = mc/mmix * (gc/molc)= gc/molmix (sum)-> gmix/molmix
    #molmix = gmix / (gmix/molmix)
    compounds = get_spec(MaterialCompounds{MOLAR,FRACTION}(),specs)
    amount = get_spec(MaterialAmount{MASS}(),specs)
    sum_x_mw = mapreduce(mw_mul,+,value(compounds),mw)
    return _ups(value(amount)/sum_x_mw,true)  
end

function moles3(::Tuple{MaterialCompounds{MASS,FRACTION},MaterialAmount{MASS}},specs::Specs,mw::T) where T<:URVec
    #mass fractions and mass
    #mw_i * xmi = gc/gmix / (molc/gc)= molc/gmix (sum)-> molmix/gmix
    # molmix = molmix/gmix * (gmix)
    compounds = get_spec(MaterialCompounds{MASS,FRACTION}(),specs)
    amount = get_spec(MaterialAmount{MASS}(),specs)
    sum_x_mw = mapreduce(mw_div,+,value(compounds),mw)
    return _ups(value(amount)*sum_x_mw,true)  
end

function mass3(::Tuple{SingleComponent,OneMol},specs::Specs,mw::T2) where T2 <:UnitReal
    return _ups(mw_mul(one(T2),mw),true)
end

#mass defined
function mass3(::Tuple{T,MaterialAmount{MASS}},specs::Specs,mw) where T
    return _ups(value(get_spec(MaterialAmount{MASS}(),specs)),true)
end

## total ammounts:
function mass3(::Tuple{T,T},specs::Specs,mw) where T<:MaterialCompounds{MASS,TOTAL_AMOUNT}
    return _ups(sum(value(get_spec(T(),specs))),true)
end

function mass3(::Tuple{T,T},specs::Specs,mw::T2) where T<:MaterialCompounds{MOLAR,TOTAL_AMOUNT} where T2<:URVec
    ###### needs mass
    ### exists mol numbers
    #mw_i * ni = (gc/molc) * molc = gc, sum(gc) = gmix
    compounds = get_spec(T(),specs)
    sum_n_mw = mapreduce(mw_mul,+,value(compounds),mw)
    return _ups(sum_n_mw,true)
end

#inverse operations with fractions
function mass3(::Tuple{MaterialCompounds{MOLAR,FRACTION},MaterialAmount{MOLAR}},specs::Specs,mw::T2) where T2<:URVec
    #moles, mol fractions
    #mw_i * xi = (gc/molc) * molc/molmix = gc/molmix (*molmix)
    amount = get_spec(MaterialAmount{MOLAR}(),specs)
    compounds = get_spec(MaterialCompounds{MOLAR,FRACTION}(),specs)
    sum_x_mw = mapreduce(mw_mul,+,value(compounds),mw)
    return _ups(value(amount)*sum_x_mw,true)
end

function mass3(::Tuple{MaterialCompounds{MASS,FRACTION},MaterialAmount{MOLAR}},specs::Specs,mw::T2) where T2<:URVec
    #mass fractions and mass
    #mw_i * xmi = gc/gmix / (molc/gc)= molc/gmix (sum)-> molmix/gmix
    # molmix = molmix/gmix * (gmix)
    compounds = get_spec(MaterialCompounds{MASS,FRACTION}(),specs)
    amount = get_spec(MaterialAmount{MASS}(),specs)
    sum_x_mw = mapreduce(mw_div,+,value(compounds),mw)
    return _ups(value(amount)*sum_x_mw,true)  
end

#single component cases, all the same
function kg_per_mol3(::Tuple{SingleComponent,OneMol},specs::Specs,mw::T2) where T2 <:UnitReal
    return _ups(mw_mul(one(T2),mw),true)
end

function kg_per_mol3(::Tuple{SingleComponent,MaterialAmount},specs::Specs,mw::T2) where {T2 <:UnitReal}
    return _ups(mw_mul(one(T2),mw),true)
end

#total amounts:
function kg_per_mol3(::Tuple{T,T},specs::Specs,mw) where T<:MaterialCompounds{MASS,TOTAL_AMOUNT}
    #mw from mass numbers
    #gc / (gc/molc) = molc sum-> molmix
    compounds = get_spec(T(),specs)
    molmix = mapreduce(mw_div,+,value(compounds),mw)
    gmix = sum(value(compounds))
    return _ups(gmix/molmix,true)  
end

function kg_per_mol3(::Tuple{T,T},specs::Specs,mw) where T<:MaterialCompounds{MOLAR,TOTAL_AMOUNT}   
    #mw from mol numbers
    #molc * gc/molc = gc sum-> gmix
    compounds = get_spec(T(),specs)
    gmix = mapreduce(mw_mul,+,value(compounds),mw)
    molmix = sum(value(compounds))
    return _ups(gmix/molmix,true) 
end

#fraction amounts
function kg_per_mol3(::Tuple{MaterialCompounds{MOLAR,FRACTION},MaterialAmount{T}},specs::Specs,mw) where T
    #mol fraction and mass
    #xni * mwi = mc/mmix * (gc/molc)= gc/molmix (sum)-> gmix/molmix
    compounds = get_spec(MaterialCompounds{MOLAR,FRACTION}(),specs)
    sum_x_mw = mapreduce(mw_mul,+,value(compounds),mw)
    return _ups(one(sum_x_mw)/sum_x_mw,true)
end

function kg_per_mol3(::Tuple{MaterialCompounds{MASS,FRACTION},MaterialAmount{T}},specs::Specs,mw) where T
    #mass fractions and mass
    compounds = get_spec(MaterialCompounds{MASS,FRACTION}(),specs)
    #xmi * mwi = gc/gmix / gc/molc = molc/gmix (sum)-> molmix/gmix
    sum_x_mw = mapreduce(mw_div,+,value(compounds),mw)
    return _ups(one(sum_x_mw)/sum_x_mw,true)
end

#invariant
function to_spec(sps,sp::Spec{SP},mw,::MOLAR) where {SP<:AbstractIntensiveSpec{MOLAR}} 
    return _ups(value(sp),true)
end

function to_spec(sps,sp::Spec{SP},mw,any_value) where {SP<:Union{Pressure,Temperature}} 
    return _ups(value(sp),true)
end

function to_spec(sps,sp::Spec{SP},mw,::MASS) where {SP<:AbstractIntensiveSpec{MASS}} 
    return _ups(value(sp),true)
end

function to_spec(sps,sp::Spec{SP},mw,::TOTAL) where {SP<:AbstractIntensiveSpec{TOTAL}} 
    return _ups(value(sp),true)
end

#to mass

#from molar to mass
function to_spec(sps,sp::Spec{SP},mw,::MASS) where {SP<:AbstractIntensiveSpec{MOLAR}} 
    return _ups(value(sp),true)/kg_per_mol2(sps,mw)
end


function to_spec(sps,sp::Spec{SP},mw,::MASS) where {SP<:AbstractIntensiveSpec{TOTAL}}
    return _ups(value(sp),true)/mass2(sps,mw)
end

#to mol
function to_spec(sps,sp::Spec{SP},mw,::MOLAR) where {SP<:AbstractIntensiveSpec{MASS}}
    return _ups(value(sp),true)*kg_per_mol2(sps,mw)
end

function to_spec(sps,sp::Spec{SP},mw,::MOLAR) where {SP<:AbstractIntensiveSpec{TOTAL}}
    return _ups(value(sp),true)/moles2(sps,mw)
end



#volume and density 
#same value
function to_spec_vol(sps,sp::Spec{T},mw,::T) where {T<:VolumeAmount}
    return _ups(value(sp),true)
end

#inversion

function to_spec_vol(sps,sp::Spec{VolumeAmount{T,VOLUME}},mw,::VolumeAmount{T,DENSITY}) where T
    val = value(sp)
    return _ups(one(val)/val,true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{T,DENSITY}},mw,::VolumeAmount{T,VOLUME}) where T
    val = value(sp)
    return _ups(one(val)/val,true)
end

#same type, volume
function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,VOLUME}},mw,::VolumeAmount{MASS,VOLUME})
    val = value(sp)
    return _ups(val*kg_per_mol2(sp,mw),true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS,VOLUME}},mw,::VolumeAmount{MOLAR,VOLUME})
    val = value(sp)
    return _ups(val/kg_per_mol2(sp,mw),true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,VOLUME}},mw,::VolumeAmount{TOTAL,VOLUME})
    val = value(sp)
    return _ups(val*moles2(sp,mw),true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL,VOLUME}},mw,::VolumeAmount{MOLAR,VOLUME})
    val = value(sp)
    return _ups(val/moles2(sp,mw),true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS,VOLUME}},mw,::VolumeAmount{TOTAL,VOLUME})
    val = value(sp)
    return _ups(val*mass2(sp,mw),true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL,VOLUME}},mw,::VolumeAmount{MASS,VOLUME})
    val = value(sp)
    return _ups(val/mass2(sp,mw),true)
end

#same type, density
function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,DENSITY}},mw,::VolumeAmount{MASS,DENSITY})
    val = value(sp)
    return _ups(val/kg_per_mol2(sp,mw),true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS,DENSITY}},mw,::VolumeAmount{MOLAR,DENSITY})
    val = value(sp)
    return _ups(val*kg_per_mol2(sp,mw),true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,DENSITY}},mw,::VolumeAmount{TOTAL,DENSITY})
    val = value(sp)
    return _ups(val/moles2(sp,mw),true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL,DENSITY}},mw,::VolumeAmount{MOLAR,DENSITY})
    val = value(sp)
    return _ups(val*moles2(sp,mw),true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS,DENSITY}},mw,::VolumeAmount{TOTAL,DENSITY})
    val = value(sp)
    return _ups(val/mass2(sp,mw),true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL,DENSITY}},mw,::VolumeAmount{MASS,DENSITY})
    val = value(sp)
    return _ups(val*mass2(sp,mw),true)
end


#all different combinations
function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,VOLUME }},mw,::VolumeAmount{TOTAL, DENSITY})
    val = value(sp)
    totv = val*moles2(sps,mw)
    return _ups(one(totv)/totv,true)
  end
  
function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,VOLUME }},mw,::VolumeAmount{MASS, DENSITY})
    val = value(sp)
    val2 = val*kg_per_mol2(sps,mw)
    return _ups(one(val2)/val2,true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL, VOLUME }},mw,::VolumeAmount{MOLAR, DENSITY})
    val = value(sp)
    val2 = val/moles2(sps,mw)
    return _ups(one(val2)/val2,true)
end
function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL, VOLUME }},mw,::VolumeAmount{MASS, DENSITY})
    val = value(sp)
    val2 = val/mass2(sps,mw)
    return _ups(one(val2)/val2,true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS, VOLUME }},mw,::VolumeAmount{MOLAR, DENSITY})
    val = value(sp)
    val2 = val/kg_per_mol2(sps,mw)
    return _ups(one(val2)/val2,true)
end
function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS, VOLUME }},mw,::VolumeAmount{TOTAL, DENSITY})
    val = value(sp)
    val2 = val*mass2(sps,mw)
    return _ups(one(val2)/val2,true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL, DENSITY }},mw,::VolumeAmount{MOLAR, VOLUME})
    val = value(sp)
    val2 = val*moles2(sps,mw)
    return _ups(one(val2)/val2,true)
end
function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS, DENSITY }},mw,::VolumeAmount{MOLAR, VOLUME})
    val = value(sp)
    val2 = val/kg_per_mol2(sps,mw)
    return _ups(one(val2)/val2,true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,DENSITY }},mw,::VolumeAmount{TOTAL, VOLUME})
    val = value(sp)
    val2 = val/moles2(sps,mw)
    return _ups(one(val2)/val2,true)
end 
function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS, DENSITY }},mw,::VolumeAmount{TOTAL, VOLUME})
    val = value(sp)
    val2 = val/mass2(sps,mw)
    return _ups(one(val2)/val2,true)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,DENSITY }},mw,::VolumeAmount{MASS, VOLUME})
    val = value(sp)
    val2 = val*kg_per_mol2(sps,mw)
    return _ups(one(val2)/val2,true)
end
function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL, DENSITY }},mw,::VolumeAmount{MASS, VOLUME})
    val = value(sp)
    val2 = val*mass2(sps,mw)
    return _ups(one(val2)/val2,true)
end
  
  
function to_spec(sps,sp::Spec{SP1},mw,x::SP2) where {SP1<:VolumeAmount,SP2<:VolumeAmount}
    return to_spec_vol(sps,sp,mw,x)
end


#material compounds part

#invariant
function to_spec_mat(sps,sp::Spec{T},mw,::T) where {T<:MaterialCompounds}
    return _ups(value(sp),true)
end

#fraction <-> total part, same base
function to_spec_mat(sps,sp::Spec{MaterialCompounds{MOLAR,FRACTION}},mw,::MaterialCompounds{MOLAR,TOTAL_AMOUNT}) 
    return moles2(sps,mw) .* _ups(value(sp),true)
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MASS,FRACTION}},mw,::MaterialCompounds{MASS,TOTAL_AMOUNT}) 
    return mass2(sps,mw) .* _ups(value(sp),true)
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MOLAR,TOTAL_AMOUNT}},mw,::MaterialCompounds{MOLAR,FRACTION}) 
    return _ups(value(sp),true) ./ moles2(sps,mw)
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MASS,TOTAL_AMOUNT}},mw,::MaterialCompounds{MASS,FRACTION}) 
    return _ups(value(sp),true) ./ mass2(sps,mw)
end

#molar total <-> mass total
function to_spec_mat(sps,sp::Spec{MaterialCompounds{MOLAR,TOTAL_AMOUNT}},mw,::MaterialCompounds{MASS,TOTAL_AMOUNT}) 
    return _ups(map(mw_mul,value(sp),mw),true)
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MASS,TOTAL_AMOUNT}},mw,::MaterialCompounds{MOLAR,TOTAL_AMOUNT}) 
    return  _ups(map(mw_div,value(sp),mw),true)
end

#molar fraction <-> mass fraction
function to_spec_mat(sps,sp::Spec{MaterialCompounds{MOLAR,FRACTION}},mw,::MaterialCompounds{MASS,FRACTION}) 
    n =  _ups(value(sp),true) .* moles2(sps,mw)
    m = _ups(map(mw_mul,n,mw),true)
    return m ./ sum(m)
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MASS,FRACTION}},mw,::MaterialCompounds{MOLAR,FRACTION}) 
    m =  _ups(value(sp),true) .* mass2(sps,mw)
    n = _ups(map(mw_div,m,mw),true)
    return n ./ sum(n)
end

#crossed properties:

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MOLAR,FRACTION}},mw,::MaterialCompounds{MASS,TOTAL_AMOUNT}) 
    n =  _ups(value(sp),true) .* moles2(sps,mw)
    m = _ups(map(mw_mul,n,mw),true)
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MASS,FRACTION}},mw,::MaterialCompounds{MOLAR,TOTAL_AMOUNT}) 
    m =  _ups(value(sp),true) .* mass2(sps,mw)
    n = _ups(map(mw_div,m,mw),true)
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MOLAR,TOTAL_AMOUNT}},mw,::MaterialCompounds{MASS,FRACTION}) 
    m =  _map(mw_mul,value(sp),mw)
   return _ups(m ./ sum(m),true)
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MASS,TOTAL_AMOUNT}},mw,::MaterialCompounds{MOLAR,FRACTION}) 
    n =  _map(mw_div,value(sp),mw)
   return _ups(m ./ sum(m),true)
end

function to_spec_compounds(sps::Specs,mw,x::SP2) where {SP2<:MaterialCompounds}
    return to_spec_compounds(amount_type(sps),sps,mw,x)
end

function to_spec_compounds(amount_type::Tuple{SingleComponent,T},sps,mw,x) where {T}
    return single_component_to_spec(sps,mw,x)
end

function to_spec_compounds(amount_type::Tuple{T1,T2},sps,mw,x) where {T1<:MaterialCompounds,T2}
    sp = get_spec(MaterialCompounds,sps)
    return to_spec_mat(sp,sps,mw,x)
end

function single_component_to_spec(sps,mw,x::MaterialCompounds{MOLAR,FRACTION})
    return [1.0]
end

function single_component_to_spec(sps,mw,x::MaterialCompounds{MASS,FRACTION})
    return [1.0]
end

function single_component_to_spec(sps,mw,x::MaterialCompounds{MOLAR,TOTAL_AMOUNT})
    return [moles2(sps,mw)]
end

function single_component_to_spec(sps,mw,x::MaterialCompounds{MASS,TOTAL_AMOUNT})
    return [mass2(sps,mw)]
end

