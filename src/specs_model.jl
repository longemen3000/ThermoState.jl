struct FromSpecs end
const UnitReal = Union{Real,Unitful.Quantity}
const URVec = AbstractVector{T} where T<:UnitReal

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
    return normalize_units(value(get_spec(MaterialAmount{MOLAR}(),specs)))
end

function moles3(::Tuple{SingleComponent,MaterialAmount{MASS}},specs::Specs,mw::T) where T <:UnitReal
    # g/(g/mol) = mol 
    x = mw_div(value(get_spec(MaterialAmount{MASS}(),specs)),mw)
    return normalize_units(x)
end


## total ammounts:
function moles3(::Tuple{T,T},specs::Specs,mw) where T<:MaterialCompounds{MOLAR,TOTAL_AMOUNT}
    return normalize_units(sum(value(get_spec(T(),specs))))
end

function moles3(::Tuple{T,T},specs::Specs,mw::T2) where T<:MaterialCompounds{MASS,TOTAL_AMOUNT} where T2<: URVec
    ###### needs mol
    ### exists mass numbers
    #mw_i * xmi = gc* (molc/gc)= molc, molmix = sum(molc)
    compounds = value(get_spec(T(),specs))
    sum_n_mw = mapreduce(mw_div,+,compounds,mw)
    return normalize_units(sum_n_mw)  
end

function moles3(::Tuple{T,MaterialAmount{MOLAR}},specs::Specs,mw) where T
    return normalize_units(value(get_spec(MaterialAmount{MOLAR}(),specs)))
end

#inverse operations with fractions
function moles3(::Tuple{MaterialCompounds{MOLAR,FRACTION},MaterialAmount{MASS}},specs::Specs,mw::T) where T<:URVec
    #mol fraction and mass
    #xni * mwi = mc/mmix * (gc/molc)= gc/molmix (sum)-> gmix/molmix
    #molmix = gmix / (gmix/molmix)
    compounds = get_spec(MaterialCompounds{MOLAR,FRACTION}(),specs)
    amount = get_spec(MaterialAmount{MASS}(),specs)
    sum_x_mw = mapreduce(mw_mul,+,value(compounds),mw)
    return normalize_units(value(amount)/sum_x_mw)  
end

function moles3(::Tuple{MaterialCompounds{MASS,FRACTION},MaterialAmount{MASS}},specs::Specs,mw::T) where T<:URVec
    #mass fractions and mass
    #mw_i * xmi = gc/gmix / (molc/gc)= molc/gmix (sum)-> molmix/gmix
    # molmix = molmix/gmix * (gmix)
    compounds = get_spec(MaterialCompounds{MASS,FRACTION}(),specs)
    amount = get_spec(MaterialAmount{MASS}(),specs)
    sum_x_mw = mapreduce(mw_div,+,value(compounds),mw)
    return normalize_units(value(amount)*sum_x_mw)  
end

function mass3(::Tuple{SingleComponent,OneMol},specs::Specs,mw::T2) where T2 <:UnitReal
    return normalize_units(mw_mul(one(T2),mw))
end

#mass defined
function mass3(::Tuple{T,MaterialAmount{MASS}},specs::Specs,mw) where T
    return normalize_units(value(get_spec(MaterialAmount{MASS}(),specs)))
end

## total ammounts:
function mass3(::Tuple{T,T},specs::Specs,mw) where T<:MaterialCompounds{MASS,TOTAL_AMOUNT}
    return normalize_units(sum(value(get_spec(T(),specs))))
end

function mass3(::Tuple{T,T},specs::Specs,mw::T2) where T<:MaterialCompounds{MOLAR,TOTAL_AMOUNT} where T2<:URVec
    ###### needs mass
    ### exists mol numbers
    #mw_i * ni = (gc/molc) * molc = gc, sum(gc) = gmix
    compounds = get_spec(T(),specs)
    sum_n_mw = mapreduce(mw_mul,+,value(compounds),mw)
    return normalize_units(sum_n_mw)
end

#inverse operations with fractions
function mass3(::Tuple{MaterialCompounds{MOLAR,FRACTION},MaterialAmount{MOLAR}},specs::Specs,mw::T2) where T2<:URVec
    #moles, mol fractions
    #mw_i * xi = (gc/molc) * molc/molmix = gc/molmix (*molmix)
    amount = get_spec(MaterialAmount{MOLAR}(),specs)
    compounds = get_spec(MaterialCompounds{MOLAR,FRACTION}(),specs)
    sum_x_mw = mapreduce(mw_mul,+,value(compounds),mw)
    return normalize_units(value(amount)*sum_x_mw)
end

function mass3(::Tuple{MaterialCompounds{MASS,FRACTION},MaterialAmount{MOLAR}},specs::Specs,mw::T2) where T2<:URVec
    #mass fractions and mass
    #mw_i * xmi = gc/gmix / (molc/gc)= molc/gmix (sum)-> molmix/gmix
    # molmix = molmix/gmix * (gmix)
    compounds = get_spec(MaterialCompounds{MASS,FRACTION}(),specs)
    amount = get_spec(MaterialAmount{MASS}(),specs)
    sum_x_mw = mapreduce(mw_div,+,value(compounds),mw)
    return normalize_units(value(amount)*sum_x_mw)  
end

#single component cases, all the same
function kg_per_mol3(::Tuple{SingleComponent,OneMol},specs::Specs,mw::T2) where T2 <:UnitReal
    return normalize_units(mw_mul(one(T2),mw))
end

function kg_per_mol3(::Tuple{SingleComponent,MaterialAmount},specs::Specs,mw::T2) where {T2 <:UnitReal}
    return normalize_units(mw_mul(one(T2),mw))
end

#total amounts:
function kg_per_mol3(::Tuple{T,T},specs::Specs,mw) where T<:MaterialCompounds{MASS,TOTAL_AMOUNT}
    #mw from mass numbers
    #gc / (gc/molc) = molc sum-> molmix
    compounds = get_spec(T(),specs)
    molmix = mapreduce(mw_div,+,value(compounds),mw)
    gmix = sum(value(compounds))
    return normalize_units(gmix/molmix)  
end

function kg_per_mol3(::Tuple{T,T},specs::Specs,mw) where T<:MaterialCompounds{MOLAR,TOTAL_AMOUNT}   
    #mw from mol numbers
    #molc * gc/molc = gc sum-> gmix
    compounds = get_spec(T(),specs)
    gmix = mapreduce(mw_mul,+,value(compounds),mw)
    molmix = sum(value(compounds))
    return normalize_units(gmix/molmix) 
end

#fraction amounts
function kg_per_mol3(::Tuple{MaterialCompounds{MOLAR,FRACTION},MaterialAmount{T}},specs::Specs,mw) where T
    #mol fraction and mass
    #xni * mwi = mc/mmix * (gc/molc)= gc/molmix (sum)-> gmix/molmix
    compounds = get_spec(MaterialCompounds{MOLAR,FRACTION}(),specs)
    sum_x_mw = mapreduce(mw_mul,+,value(compounds),mw)
    return normalize_units(one(sum_x_mw)/sum_x_mw)
end

function kg_per_mol3(::Tuple{MaterialCompounds{MASS,FRACTION},MaterialAmount{T}},specs::Specs,mw) where T
    #mass fractions and mass
    compounds = get_spec(MaterialCompounds{MASS,FRACTION}(),specs)
    #xmi * mwi = gc/gmix / gc/molc = molc/gmix (sum)-> molmix/gmix
    sum_x_mw = mapreduce(mw_div,+,value(compounds),mw)
    return normalize_units(one(sum_x_mw)/sum_x_mw)
end

#invariant
function to_spec(sps,sp::Spec{SP},mw,::MOLAR) where {SP<:AbstractIntensiveSpec{MOLAR}} 
    return normalize_units(value(sp))
end

function to_spec(sps,sp::Spec{SP},mw,any_value) where {SP<:Union{Pressure,Temperature}} 
    return normalize_units(value(sp))
end

function to_spec(sps,sp::Spec{SP},mw,::MASS) where {SP<:AbstractIntensiveSpec{MASS}} 
    return normalize_units(value(sp))
end

function to_spec(sps,sp::Spec{SP},mw,::TOTAL) where {SP<:AbstractIntensiveSpec{TOTAL}} 
    return normalize_units(value(sp))
end

#to mass

#from molar to mass
function to_spec(sps,sp::Spec{SP},mw,::MASS) where {SP<:AbstractIntensiveSpec{MOLAR}} 
    return normalize_units(value(sp))/kg_per_mol2(sps,mw)
end


function to_spec(sps,sp::Spec{SP},mw,::MASS) where {SP<:AbstractIntensiveSpec{TOTAL}}
    return normalize_units(value(sp))/mass2(sps,mw)
end

#to mol
function to_spec(sps,sp::Spec{SP},mw,::MOLAR) where {SP<:AbstractIntensiveSpec{MASS}}
    return normalize_units(value(sp))*kg_per_mol2(sps,mw)
end

function to_spec(sps,sp::Spec{SP},mw,::MOLAR) where {SP<:AbstractIntensiveSpec{TOTAL}}
    return normalize_units(value(sp))/moles2(sps,mw)
end



#volume and density 
#same value
function to_spec_vol(sps,sp::Spec{T},mw,::T) where {T<:VolumeAmount}
    return normalize_units(value(sp))
end

#inversion

function to_spec_vol(sps,sp::Spec{VolumeAmount{T,VOLUME}},mw,::VolumeAmount{T,DENSITY}) where T
    val = value(sp)
    return normalize_units(one(val)/val)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{T,DENSITY}},mw,::VolumeAmount{T,VOLUME}) where T
    val = value(sp)
    return normalize_units(one(val)/val)
end

#same type, volume
function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,VOLUME}},mw,::VolumeAmount{MASS,VOLUME})
    val = value(sp)
    return normalize_units(val*kg_per_mol2(sps,mw))
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS,VOLUME}},mw,::VolumeAmount{MOLAR,VOLUME})
    val = value(sp)
    return normalize_units(val/kg_per_mol2(sps,mw))
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,VOLUME}},mw,::VolumeAmount{TOTAL,VOLUME})
    val = value(sp)
    return normalize_units(val*moles2(sps,mw))
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL,VOLUME}},mw,::VolumeAmount{MOLAR,VOLUME})
    val = value(sp)
    return normalize_units(val/moles2(sps,mw))
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS,VOLUME}},mw,::VolumeAmount{TOTAL,VOLUME})
    val = value(sp)
    return normalize_units(val*mass2(sps,mw))
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL,VOLUME}},mw,::VolumeAmount{MASS,VOLUME})
    val = value(sp)
    return normalize_units(val/mass2(sps,mw))
end

#same type, density
function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,DENSITY}},mw,::VolumeAmount{MASS,DENSITY})
    val = value(sp)
    return normalize_units(val/kg_per_mol2(sps,mw))
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS,DENSITY}},mw,::VolumeAmount{MOLAR,DENSITY})
    val = value(sp)
    return normalize_units(val*kg_per_mol2(sps,mw))
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,DENSITY}},mw,::VolumeAmount{TOTAL,DENSITY})
    val = value(sp)
    return normalize_units(val/moles2(sps,mw))
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL,DENSITY}},mw,::VolumeAmount{MOLAR,DENSITY})
    val = value(sp)
    return normalize_units(val*moles2(sps,mw))
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS,DENSITY}},mw,::VolumeAmount{TOTAL,DENSITY})
    val = value(sp)
    return normalize_units(val/mass2(sps,mw))
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL,DENSITY}},mw,::VolumeAmount{MASS,DENSITY})
    val = value(sp)
    return normalize_units(val*mass2(sps,mw))
end


#all different combinations
function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,VOLUME }},mw,::VolumeAmount{TOTAL, DENSITY})
    val = value(sp)
    totv = val*moles2(sps,mw)
    return normalize_units(one(totv)/totv)
  end
  
function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,VOLUME }},mw,::VolumeAmount{MASS, DENSITY})
    val = value(sp)
    val2 = val*kg_per_mol2(sps,mw)
    return normalize_units(one(val2)/val2)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL, VOLUME }},mw,::VolumeAmount{MOLAR, DENSITY})
    val = value(sp)
    val2 = val/moles2(sps,mw)
    return normalize_units(one(val2)/val2)
end
function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL, VOLUME }},mw,::VolumeAmount{MASS, DENSITY})
    val = value(sp)
    val2 = val/mass2(sps,mw)
    return normalize_units(one(val2)/val2)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS, VOLUME }},mw,::VolumeAmount{MOLAR, DENSITY})
    val = value(sp)
    val2 = val/kg_per_mol2(sps,mw)
    return normalize_units(one(val2)/val2)
end
function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS, VOLUME }},mw,::VolumeAmount{TOTAL, DENSITY})
    val = value(sp)
    val2 = val*mass2(sps,mw)
    return normalize_units(one(val2)/val2)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL, DENSITY }},mw,::VolumeAmount{MOLAR, VOLUME})
    val = value(sp)
    val2 = val*moles2(sps,mw)
    return normalize_units(one(val2)/val2)
end
function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS, DENSITY }},mw,::VolumeAmount{MOLAR, VOLUME})
    val = value(sp)
    val2 = val/kg_per_mol2(sps,mw)
    return normalize_units(one(val2)/val2)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,DENSITY }},mw,::VolumeAmount{TOTAL, VOLUME})
    val = value(sp)
    val2 = val/moles2(sps,mw)
    return normalize_units(one(val2)/val2)
end 
function to_spec_vol(sps,sp::Spec{VolumeAmount{MASS, DENSITY }},mw,::VolumeAmount{TOTAL, VOLUME})
    val = value(sp)
    val2 = val/mass2(sps,mw)
    return normalize_units(one(val2)/val2)
end

function to_spec_vol(sps,sp::Spec{VolumeAmount{MOLAR,DENSITY }},mw,::VolumeAmount{MASS, VOLUME})
    val = value(sp)
    val2 = val*kg_per_mol2(sps,mw)
    return normalize_units(one(val2)/val2)
end
function to_spec_vol(sps,sp::Spec{VolumeAmount{TOTAL, DENSITY }},mw,::VolumeAmount{MASS, VOLUME})
    val = value(sp)
    val2 = val*mass2(sps,mw)
    return normalize_units(one(val2)/val2)
end
  
  
function to_spec(sps,sp::Spec{SP1},mw,x::SP2) where {SP1<:VolumeAmount,SP2<:VolumeAmount}
    return to_spec_vol(sps,sp,mw,x)
end


#material compounds part

#invariant
function to_spec_mat(sps,sp::Spec{T},mw,::T) where {T<:MaterialCompounds}
    return normalize_units(value(sp))
end

#fraction <-> total part, same base
function to_spec_mat(sps,sp::Spec{MaterialCompounds{MOLAR,FRACTION}},mw,::MaterialCompounds{MOLAR,TOTAL_AMOUNT}) 
    return moles2(sps,mw) .* normalize_units(value(sp))
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MASS,FRACTION}},mw,::MaterialCompounds{MASS,TOTAL_AMOUNT}) 
    return mass2(sps,mw) .* normalize_units(value(sp))
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MOLAR,TOTAL_AMOUNT}},mw,::MaterialCompounds{MOLAR,FRACTION}) 
    return normalize_units(value(sp)) ./ moles2(sps,mw)
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MASS,TOTAL_AMOUNT}},mw,::MaterialCompounds{MASS,FRACTION}) 
    return normalize_units(value(sp)) ./ mass2(sps,mw)
end

#molar total <-> mass total
function to_spec_mat(sps,sp::Spec{MaterialCompounds{MOLAR,TOTAL_AMOUNT}},mw,::MaterialCompounds{MASS,TOTAL_AMOUNT}) 
    return normalize_units(map(mw_mul,value(sp),mw))
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MASS,TOTAL_AMOUNT}},mw,::MaterialCompounds{MOLAR,TOTAL_AMOUNT}) 
    return  normalize_units(map(mw_div,value(sp),mw))
end

#molar fraction <-> mass fraction
function to_spec_mat(sps,sp::Spec{MaterialCompounds{MOLAR,FRACTION}},mw,::MaterialCompounds{MASS,FRACTION}) 
    n =  normalize_units(value(sp)) .* moles2(sps,mw)
    m = normalize_units(map(mw_mul,n,mw))
    return m ./ sum(m)
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MASS,FRACTION}},mw,::MaterialCompounds{MOLAR,FRACTION}) 
    m =  normalize_units(value(sp)) .* mass2(sps,mw)
    n = map(mw_div,m,mw)
    return n ./ sum(n)
end

#crossed properties:

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MOLAR,FRACTION}},mw,::MaterialCompounds{MASS,TOTAL_AMOUNT}) 
    n =  value(sp).* moles2(sps,mw)
    m = normalize_units(map(mw_mul,n,mw))
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MASS,FRACTION}},mw,::MaterialCompounds{MOLAR,TOTAL_AMOUNT}) 
    m =  value(sp) .* mass2(sps,mw)
    n = normalize_units(map(mw_div,m,mw))
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MOLAR,TOTAL_AMOUNT}},mw,::MaterialCompounds{MASS,FRACTION}) 
    m =  _map(mw_mul,value(sp),mw)
   return normalize_units(m ./ sum(m))
end

function to_spec_mat(sps,sp::Spec{MaterialCompounds{MASS,TOTAL_AMOUNT}},mw,::MaterialCompounds{MOLAR,FRACTION}) 
    n =  _map(mw_div,value(sp),mw)
   return normalize_units(m ./ sum(m))
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

