

"""
    FromState

a basic thermodynamic model, used to extract thermodynamic specifications from a state without using a Equation of State.
"""
struct FromState <: ThermoModel end
const URVec = AbstractVector{T} where T<:Number

moles2(a::ThermodynamicState,mw) = moles3(amount_type(a),a,mw)
mass2(a::ThermodynamicState,mw) =  mass3(amount_type(a),a,mw)
kg_per_mol2(a::ThermodynamicState,mw) =  kg_per_mol3(amount_type(a),a,mw)


function mw_error(::Nothing)
    throw(error("molecular weight information is required for this combination of properties."))
end

function mw_error(mw)
    return nothing
end
##one mol cases

#those cases should rarely be called
function moles3(::Tuple{T,OneMol},st,mw::Nothing) where T
    return one(numtype(st))
end

function moles3(::Tuple{T,OneMol},st,mw::T2) where {T,T2 <:Number}
    return one(promote_type(T2,numtype(st)))
end

#error?
function moles3(::Tuple{T,OneMol},st,mw::T2) where T where T2<:AbstractVector{T3} where T3 <: Number
    return one(promote_type(T3,numtype(st)))
end

function moles3(::Tuple{SingleComponent,MaterialAmount{MOLAR}},st,mw::T) where T <:Number
    return normalize_units(value(get_spec(MaterialAmount{MOLAR}(),st)))*_one(st)
end

function moles3(::Tuple{SingleComponent,MaterialAmount{MASS}},st,mw::T) where T <:Number
    # g/(g/mol) = mol 
    mw_error(mw)
    x = mw_div(value(get_spec(MaterialAmount{MASS}(),st)),mw)
    return normalize_units(x)*_one(st)
end

## total ammounts:
function moles3(::Tuple{T,T},st,mw) where T<:MaterialCompounds{MOLAR,TOTAL_AMOUNT}
    return normalize_units(sum(value(get_spec(T(),st))))*_one(st)
end

function moles3(::Tuple{T,T},st,mw::T2) where T<:MaterialCompounds{MASS,TOTAL_AMOUNT} where T2<: URVec
    ###### needs mol
    ### exists mass numbers
    #mw_i * xmi = gc* (molc/gc)= molc, molmix = sum(molc)
    mw_error(mw)
    compounds::Spec = get_spec(T(),st)
    v = compounds.val
    sum_n_mw = zero(numtype(st))
    for i in 1:length(mw)
        sum_n_mw += mw_div(v[i],mw[i])
    end
    #sum_n_mw = mapreduce(mw_div,+,compounds,mw)
    return normalize_units(sum_n_mw)  
end

function moles3(::Tuple{T,MaterialAmount{MOLAR}},st,mw) where T
    return normalize_units(value(get_spec(MaterialAmount{MOLAR}(),st)))*_one(st)
end

#inverse operations with fractions
function moles3(::Tuple{MaterialCompounds{MOLAR,FRACTION},MaterialAmount{MASS}},st,mw::T) where T<:URVec
    #mol fraction and mass
    #xni * mwi = mc/mmix * (gc/molc)= gc/molmix (sum)-> gmix/molmix
    #molmix = gmix / (gmix/molmix)
    mw_error(mw)
    compounds::Spec = get_spec(MaterialCompounds{MOLAR,FRACTION}(),st)
    amount::Spec = get_spec(MaterialAmount{MASS}(),st)
    v = compounds.val
    #sum_x_mw = mapreduce(mw_mul,+,value(compounds),mw)
    sum_x_mw = zero(numtype(st))
    for i in 1:length(mw)
        sum_x_mw += mw_mul(v[i],mw[i])
    end
    return normalize_units(value(amount)/sum_x_mw) 
end

function moles3(::Tuple{MaterialCompounds{MASS,FRACTION},MaterialAmount{MASS}},st,mw::T) where T<:URVec
    #mass fractions and mass
    #mw_i * xmi = gc/gmix / (molc/gc)= molc/gmix (sum)-> molmix/gmix
    # molmix = molmix/gmix * (gmix)
    mw_error(mw)
    compounds::Spec = get_spec(MaterialCompounds{MASS,FRACTION}(),st)
    amount::Spec  = get_spec(MaterialAmount{MASS}(),st)
    v = compounds.val
    sum_x_mw = zero(numtype(st))
    for i in 1:length(mw)
        sum_x_mw += mw_div(v[i],mw[i])
    end
    #sum_x_mw = mapreduce(mw_div,+,value(compounds),mw)
    return normalize_units(value(amount)*sum_x_mw)  
end

function mass3(::Tuple{SingleComponent,OneMol},st,mw::T2) where T2 <:Number
    mw_error(mw)
    return normalize_units(mw_mul(one(T2),mw))*_one(st)
end

#mass defined
function mass3(::Tuple{T,MaterialAmount{MASS}},st,mw) where T
    return normalize_units(value(get_spec(MaterialAmount{MASS}(),st)))*_one(st)
end

function mass3(::Tuple{SingleComponent,MaterialAmount{MOLAR}},st,mw)
    mw_error(mw)
    mol = value(get_spec(MaterialAmount{MOLAR}(),st))
    return normalize_units(mw_mul(mol,mw))*_one(st)
end

## total ammounts:
function mass3(::Tuple{T,T},st,mw) where T<:MaterialCompounds{MASS,TOTAL_AMOUNT}
    return normalize_units(sum(value(get_spec(T(),st))))*_one(st)
end

function mass3(::Tuple{T,T},st::S,mw::T2) where {S<:ThermodynamicState,T<:MaterialCompounds{MOLAR,TOTAL_AMOUNT},T2<:URVec}
    ###### needs mass
    ### exists mol numbers
    #mw_i * ni = (gc/molc) * molc = gc, sum(gc) = gmix
    mw_error(mw)
    compounds = get_spec(T(),st)::Spec
    v = compounds.val
    sum_n_mw = zero(numtype(st))
    for i in 1:length(mw)
        sum_n_mw += mw_mul(v[i],mw[i])
    end
    return normalize_units(sum_n_mw)
end

#inverse operations with fractions
function mass3(::Tuple{MaterialCompounds{MOLAR,FRACTION},MaterialAmount{MOLAR}},st,mw::T2) where T2<:URVec
    #moles, mol fractions
    #mw_i * xi = (gc/molc) * molc/molmix = gc/molmix (*molmix)
    mw_error(mw)
    amount = get_spec(MaterialAmount{MOLAR}(),st)
    compounds::Spec = get_spec(MaterialCompounds{MOLAR,FRACTION}(),st)
    sum_x_mw = zero(numtype(st))
    v = compounds.val
    for i in 1:length(mw)
        sum_x_mw += mw_mul(v[i],mw[i])
    end    
    return normalize_units(value(amount)*sum_x_mw)
end

function mass3(::Tuple{MaterialCompounds{MASS,FRACTION},MaterialAmount{MOLAR}},st,mw::T2) where T2<:URVec
    #mass fractions and moles
    #sum_x_mw = xmi/mwi = mc/mmix / mc/molc = mc/mmix * molc/mc = molc/mmix -> sum -> molmix/mmix
    #a = molmix/mmix
    #inv(a) = mmix/molmix
    # molmix = molmix/gmix * (gmix)
    mw_error(mw)
    compounds::Spec = get_spec(MaterialCompounds{MASS,FRACTION}(),st)
    v = compounds.val
    amount = get_spec(MaterialAmount{MOLAR}(),st)
    sum_x_mw = zero(numtype(st))
    for i in 1:length(mw)
        sum_x_mw += mw_div(v[i],mw[i])
    end
    return normalize_units(value(amount)/sum_x_mw) 
end

#single component cases, all the same
function kg_per_mol3(::Tuple{SingleComponent,OneMol},st,mw::T2) where T2 <:Number
    return normalize_units(mw_mul(one(T2),mw))*_one(st)
end

function kg_per_mol3(::Tuple{SingleComponent,MaterialAmount},st,mw::T2) where {T2 <:Number}
    return normalize_units(mw_mul(one(T2),mw))*_one(st)
end

#total amounts:
function kg_per_mol3(::Tuple{T,T},st,mw) where T<:MaterialCompounds{MASS,TOTAL_AMOUNT}
    #mw from mass numbers
    #gc / (gc/molc) = molc sum-> molmix
    mw_error(mw)
    compounds::Spec = get_spec(T(),st)
    v = compounds.val
    molmix = zero(numtype(st))
    for i in 1:length(mw)
        molmix += mw_div(v[i],mw[i])
    end
    gmix = sum(v)
    return normalize_units(gmix/molmix) 
end

function kg_per_mol3(::Tuple{T,T},st,mw) where T<:MaterialCompounds{MOLAR,TOTAL_AMOUNT}   
    #mw from mol numbers
    #molc * gc/molc = gc sum-> gmix
    mw_error(mw)
    compounds::Spec = get_spec(T(),st)
    v = compounds.val
    gmix = zero(numtype(st))
    for i in 1:length(mw)
        gmix += mw_mul(v[i],mw[i])
    end
    #gmix = mapreduce(mw_mul,+,value(compounds),mw)
    molmix = sum(value(compounds))
    return normalize_units(gmix/molmix)*_one(st) 
end

#fraction amounts
function kg_per_mol3(::Tuple{MaterialCompounds{MOLAR,FRACTION},MaterialAmount{T}},st,mw) where T
    #mol fraction and mass
    #xni * mwi = mc/mmix * (gc/molc)= gc/molmix (sum)-> gmix/molmix
    mw_error(mw)
    compounds::Spec = get_spec(MaterialCompounds{MOLAR,FRACTION}(),st)
    #sum_x_mw = mapreduce(mw_mul,+,value(compounds),mw)*_one(st)
    v = compounds.val
    sum_x_mw = zero(numtype(st))
    for i in 1:length(mw)
        sum_x_mw += mw_mul(v[i],mw[i])
    end
    return normalize_units(one(sum_x_mw)/sum_x_mw)
end

function kg_per_mol3(::Tuple{MaterialCompounds{MASS,FRACTION},MaterialAmount{T}},st,mw) where T
    mw_error(mw)
    #mass fractions and mass
    compounds::Spec = get_spec(MaterialCompounds{MASS,FRACTION}(),st)
    #xmi * mwi = gc/gmix / gc/molc = molc/gmix (sum)-> molmix/gmix
    #sum_x_mw = mapreduce(mw_div,+,value(compounds),mw)*_one(st)
    v = compounds.val
    sum_x_mw = zero(numtype(st))
    for i in 1:length(mw)
        sum_x_mw += mw_div(v[i],mw[i])
    end
    return normalize_units(one(sum_x_mw)/sum_x_mw)
end

#invariant
function to_spec(st,sp::Spec{SP},mw,::MOLAR) where {SP<:AbstractIntensiveSpec{MOLAR}} 
    return normalize_units(value(sp))*_one(st)
end

function to_spec(st,sp::Spec{SP},mw,any_value) where {SP<:Union{Pressure,Temperature}} 
    return normalize_units(value(sp))*_one(st)
end

function to_spec(st,sp::Spec{SP},mw,::MASS) where {SP<:AbstractIntensiveSpec{MASS}} 
    return normalize_units(value(sp))*_one(st)
end

function to_spec(st,sp::Spec{SP},mw,::TOTAL) where {SP<:AbstractIntensiveSpec{TOTAL}} 
    return normalize_units(value(sp))*_one(st)
end

#to mass

#from molar to mass
function to_spec(st,sp::Spec{SP},mw,::MASS) where {SP<:AbstractIntensiveSpec{MOLAR}} 
    return normalize_units(value(sp))/kg_per_mol2(st,mw)*_one(st)
end


function to_spec(st,sp::Spec{SP},mw,::MASS) where {SP<:AbstractIntensiveSpec{TOTAL}}
    return normalize_units(value(sp))/mass2(st,mw)*_one(st)
end


#to mol
function to_spec(st,sp::Spec{SP},mw,::MOLAR) where {SP<:AbstractIntensiveSpec{MASS}}
    return normalize_units(value(sp))*kg_per_mol2(st,mw)
end

function to_spec(st,sp::Spec{SP},mw,::MOLAR) where {SP<:AbstractIntensiveSpec{TOTAL}}
    return normalize_units(value(sp))/moles2(st,mw)
end

#to total
function to_spec(st,sp::Spec{SP},mw,::TOTAL) where {SP<:AbstractIntensiveSpec{MOLAR}}
    return normalize_units(value(sp))*moles2(st,mw)
end

function to_spec(st,sp::Spec{SP},mw,::TOTAL) where {SP<:AbstractIntensiveSpec{MASS}}
    return normalize_units(value(sp))*mass2(st,mw)
end


#volume and density 
#same value
function to_spec_vol(st,sp::Spec{T},mw,::T) where {T<:VolumeAmount}
    return normalize_units(value(sp))
end

#inversion

function to_spec_vol(st,sp::Spec{VolumeAmount{T,VOLUME}},mw,::VolumeAmount{T,DENSITY}) where T
    val = value(sp)
    return normalize_units(one(val)/val)
end

function to_spec_vol(st,sp::Spec{VolumeAmount{T,DENSITY}},mw,::VolumeAmount{T,VOLUME}) where T
    val = value(sp)
    return normalize_units(one(val)/val)
end

#same type, volume
function to_spec_vol(st,sp::Spec{VolumeAmount{MOLAR,VOLUME}},mw,::VolumeAmount{MASS,VOLUME})
    val = value(sp)
    return normalize_units(val/kg_per_mol2(st,mw))
end

function to_spec_vol(st,sp::Spec{VolumeAmount{MASS,VOLUME}},mw,::VolumeAmount{MOLAR,VOLUME})
    val = value(sp)
    return normalize_units(val*kg_per_mol2(st,mw))
end

function to_spec_vol(st,sp::Spec{VolumeAmount{MOLAR,VOLUME}},mw,::VolumeAmount{TOTAL,VOLUME})
    val = value(sp)
    return normalize_units(val*moles2(st,mw))
end

function to_spec_vol(st,sp::Spec{VolumeAmount{TOTAL,VOLUME}},mw,::VolumeAmount{MOLAR,VOLUME})
    val = value(sp)
    return normalize_units(val/moles2(st,mw))
end

function to_spec_vol(st,sp::Spec{VolumeAmount{MASS,VOLUME}},mw,::VolumeAmount{TOTAL,VOLUME})
    val = value(sp)
    return normalize_units(val*mass2(st,mw))
end

function to_spec_vol(st,sp::Spec{VolumeAmount{TOTAL,VOLUME}},mw,::VolumeAmount{MASS,VOLUME})
    val = value(sp)
    return normalize_units(val/mass2(st,mw))
end

#same type, density
function to_spec_vol(st,sp::Spec{VolumeAmount{MOLAR,DENSITY}},mw,::VolumeAmount{MASS,DENSITY})
    val = value(sp)
    return normalize_units(val*kg_per_mol2(st,mw))
end

function to_spec_vol(st,sp::Spec{VolumeAmount{MASS,DENSITY}},mw,::VolumeAmount{MOLAR,DENSITY})
    val = value(sp)
    return normalize_units(val/kg_per_mol2(st,mw))
end

function to_spec_vol(st,sp::Spec{VolumeAmount{MOLAR,DENSITY}},mw,::VolumeAmount{TOTAL,DENSITY})
    val = value(sp)
    return normalize_units(val/moles2(st,mw))
end

function to_spec_vol(st,sp::Spec{VolumeAmount{TOTAL,DENSITY}},mw,::VolumeAmount{MOLAR,DENSITY})
    val = value(sp)
    return  (val*moles2(st,mw))
end

function to_spec_vol(st,sp::Spec{VolumeAmount{MASS,DENSITY}},mw,::VolumeAmount{TOTAL,DENSITY})
    val = value(sp)
    return normalize_units(val/mass2(st,mw))
end

function to_spec_vol(st,sp::Spec{VolumeAmount{TOTAL,DENSITY}},mw,::VolumeAmount{MASS,DENSITY})
    val = value(sp)
    return normalize_units(val*mass2(st,mw))
end


#all different combinations
function to_spec_vol(st,sp::Spec{VolumeAmount{MOLAR,VOLUME}},mw,::VolumeAmount{TOTAL, DENSITY})
    val = value(sp)
    totv = val*moles2(st,mw)
    return normalize_units(one(totv)/totv)
  end
  
function to_spec_vol(st,sp::Spec{VolumeAmount{MOLAR,VOLUME}},mw,::VolumeAmount{MASS, DENSITY})
    val = value(sp)
    val2 = val/kg_per_mol2(st,mw)
    return normalize_units(one(val2)/val2)
end

function to_spec_vol(st,sp::Spec{VolumeAmount{TOTAL,VOLUME}},mw,::VolumeAmount{MOLAR, DENSITY})
    val = value(sp)
    val2 = val/moles2(st,mw)
    return normalize_units(one(val2)/val2)
end
function to_spec_vol(st,sp::Spec{VolumeAmount{TOTAL,VOLUME}},mw,::VolumeAmount{MASS, DENSITY})
    val = value(sp)
    val2 = val/mass2(st,mw)
    return normalize_units(one(val2)/val2)
end

function to_spec_vol(st,sp::Spec{VolumeAmount{MASS,VOLUME}},mw,::VolumeAmount{MOLAR, DENSITY})
    val = value(sp)
    val2 = val*kg_per_mol2(st,mw)
    return normalize_units(one(val2)/val2)
end
function to_spec_vol(st,sp::Spec{VolumeAmount{MASS,VOLUME}},mw,::VolumeAmount{TOTAL, DENSITY})
    val = value(sp)
    val2 = val*mass2(st,mw)
    return normalize_units(one(val2)/val2)
end

function to_spec_vol(st,sp::Spec{VolumeAmount{TOTAL,DENSITY}},mw,::VolumeAmount{MOLAR, VOLUME})
    val = value(sp)
    val2 = val*moles2(st,mw)
    return normalize_units(one(val2)/val2)
end
function to_spec_vol(st,sp::Spec{VolumeAmount{MASS,DENSITY}},mw,::VolumeAmount{MOLAR, VOLUME})
    val = value(sp)
    val2 = val/kg_per_mol2(st,mw)
    return normalize_units(one(val2)/val2)
end

function to_spec_vol(st,sp::Spec{VolumeAmount{MOLAR,DENSITY}},mw,::VolumeAmount{TOTAL, VOLUME})
    val = value(sp)
    val2 = val/moles2(st,mw)
    return normalize_units(one(val2)/val2)
end 
function to_spec_vol(st,sp::Spec{VolumeAmount{MASS,DENSITY}},mw,::VolumeAmount{TOTAL, VOLUME})
    val = value(sp)
    val2 = val/mass2(st,mw)
    return normalize_units(one(val2)/val2)
end

function to_spec_vol(st,sp::Spec{VolumeAmount{MOLAR,DENSITY }},mw,::VolumeAmount{MASS, VOLUME})
    val = value(sp)
    val2 = val*kg_per_mol2(st,mw)
    return normalize_units(one(val2)/val2)
end
function to_spec_vol(st,sp::Spec{VolumeAmount{TOTAL, DENSITY }},mw,::VolumeAmount{MASS, VOLUME})
    val = value(sp)
    val2 = val*mass2(st,mw)
    return normalize_units(one(val2)/val2)
end
  
  
function to_spec(st,sp::Spec{SP1},mw,x::SP2) where {SP1<:VolumeAmount,SP2<:VolumeAmount}
    return to_spec_vol(st,sp,mw,x)
end


#material compounds part



#fraction <-> total part, same base
function to_spec_mat(st,sp::Spec{MaterialCompounds{MOLAR,FRACTION}},mw,::MaterialCompounds{MOLAR,TOTAL_AMOUNT}) 
    return moles2(st,mw) .* normalize_units(value(sp))
end

function to_spec_mat(st,sp::Spec{MaterialCompounds{MASS,FRACTION}},mw,::MaterialCompounds{MASS,TOTAL_AMOUNT}) 
    return mass2(st,mw) .* normalize_units(value(sp))
end

function to_spec_mat(st,sp::Spec{MaterialCompounds{MOLAR,TOTAL_AMOUNT}},mw,::MaterialCompounds{MOLAR,FRACTION}) 
    return normalize_units(value(sp)) ./ moles2(st,mw)
end

function to_spec_mat(st,sp::Spec{MaterialCompounds{MASS,TOTAL_AMOUNT}},mw,::MaterialCompounds{MASS,FRACTION}) 
    return normalize_units(value(sp)) ./ mass2(st,mw)
end

#molar total <-> mass total
function to_spec_mat(st,sp::Spec{MaterialCompounds{MOLAR,TOTAL_AMOUNT}},mw,::MaterialCompounds{MASS,TOTAL_AMOUNT}) 
    return normalize_units(map(mw_mul,value(sp),mw))
end

function to_spec_mat(st,sp::Spec{MaterialCompounds{MASS,TOTAL_AMOUNT}},mw,::MaterialCompounds{MOLAR,TOTAL_AMOUNT}) 
    return  normalize_units(map(mw_div,value(sp),mw))
end

#molar fraction <-> mass fraction
function to_spec_mat(st,sp::Spec{MaterialCompounds{MOLAR,FRACTION}},mw,::MaterialCompounds{MASS,FRACTION}) 
    n =  normalize_units(value(sp)) .* moles2(st,mw)
    m = normalize_units(map(mw_mul,n,mw))
    return m ./ sum(m)
end

function to_spec_mat(st,sp::Spec{MaterialCompounds{MASS,FRACTION}},mw,::MaterialCompounds{MOLAR,FRACTION}) 
    m =  normalize_units(value(sp)) .* mass2(st,mw)
    n = map(mw_div,m,mw)
    return n ./ sum(n)
end

#crossed properties:

function to_spec_mat(st,sp::Spec{MaterialCompounds{MOLAR,FRACTION}},mw,::MaterialCompounds{MASS,TOTAL_AMOUNT}) 
    n =  value(sp).* moles2(st,mw)
    m = normalize_units(map(mw_mul,n,mw))
end

function to_spec_mat(st,sp::Spec{MaterialCompounds{MASS,FRACTION}},mw,::MaterialCompounds{MOLAR,TOTAL_AMOUNT}) 
    m =  value(sp) .* mass2(st,mw)
    n = normalize_units(map(mw_div,m,mw))
end

function to_spec_mat(st,sp::Spec{MaterialCompounds{MOLAR,TOTAL_AMOUNT}},mw,::MaterialCompounds{MASS,FRACTION}) 
    mw_error(mw)
    m =  map(mw_mul,value(sp),mw)
    return normalize_units(m ./ sum(m))
end

function to_spec_mat(st,sp::Spec{MaterialCompounds{MASS,TOTAL_AMOUNT}},mw,::MaterialCompounds{MOLAR,FRACTION}) 
    mw_error(mw)
    n =  map(mw_div,value(sp),mw)
    return normalize_units(n ./ sum(n))
end

function to_spec_compounds(st,mw,x::SP2) where {SP2<:MaterialCompounds}
    return to_spec_compounds(amount_type(st),st,mw,x)
end

function to_spec_compounds(amount_type::Tuple{SingleComponent,T},st,mw,x) where {T}
    return single_component_to_spec(st,mw,x)
end


function to_spec_compounds(amount_type::Tuple{T1,T2},st,mw,x::T1) where {T1<:MaterialCompounds,T2}
    sp = get_spec(T1,st)
    return normalize_units(value(sp))
end

function to_spec_compounds(amount_type::Tuple{T1,T2},st,mw,x) where {T1<:MaterialCompounds,T2}
    sp = get_spec(T1,st)
    return to_spec_mat(st,sp,mw,x)
end

function single_component_to_spec(st,mw,x::MaterialCompounds{MOLAR,FRACTION})
    return [_one(st)]
end

function single_component_to_spec(st,mw,x::MaterialCompounds{MASS,FRACTION})
    return [_one(st)]
end

function single_component_to_spec(st,mw,x::MaterialCompounds{MOLAR,TOTAL_AMOUNT})
    return [moles2(st,mw)]
end

function single_component_to_spec(st,mw,x::MaterialCompounds{MASS,TOTAL_AMOUNT})
    return [mass2(st,mw)]
end

