
 #returns the specification symbol in multicomponent
#returns :singlecomponent if there are no specifications
#throws error if more than one specification is found
const AMOUNT_CONST_FIRST =Dict{Int,Any}(
     000 => SingleComponent()
    ,110 => MaterialCompounds{MOLAR,FRACTION}()
    ,120 => MaterialCompounds{MASS,FRACTION}()
    ,210 => MaterialCompounds{MOLAR,TOTAL_AMOUNT}()
    ,220 => MaterialCompounds{MASS,TOTAL_AMOUNT}()
    ,310 => HumiditySpec{WetBulbTemperature}()
    ,320 => HumiditySpec{HumidityRatio}()
    ,330 => HumiditySpec{MolarHumidity}()
    ,340 => HumiditySpec{MassHumidity}()
    ,350 => HumiditySpec{RelativeHumidity}()
    ,360 => HumiditySpec{HumidityDewPoint}()
)

const AMOUNT_CONST_SECOND =Dict{Int,Any}(
     000 => OneMol()
    ,001 => MaterialAmount{MOLAR}()
    ,002 => MaterialAmount{MASS}()
)


#the following functions are essencial in the ordering
#do not change the relative order!!

#equilibria properties are always first
_stateorder(x::PhaseFractions) = 1
_stateorder(x::TwoPhaseEquilibrium) = 2
_stateorder(x::VaporFraction) = 3

#pressure is always first in the scalar properties
_stateorder(x::Pressure) = 4
_stateorder(x::Helmholtz) = 5
_stateorder(x::Gibbs) = 6
_stateorder(x::InternalEnergy) = 7
_stateorder(x::VolumeAmount) = 8
_stateorder(x::Enthalpy) = 9
_stateorder(x::Entropy) = 10
#temperature always for last of the scalar properties
_stateorder(x::Temperature) = 11


#From here, it doesnt matter the order, as those are treated specially

_stateorder(x::HumiditySpec) = 1000

_stateorder(x::MaterialCompounds) = 10002

_stateorder(x::MaterialAmount ) = 10001
_stateorder(x::PhaseTag) = 10000
_stateorder(x::Options) = 10004
_stateorder(x::MolecularWeight) = 10005
_stateorder(x::VariableSpec) = 10006



function _static_specs(x::Type{ThermodynamicState{S,C}}) where {S,C}
    tvec = S.parameters
    res1 = [tvec[i].parameters[1]() for i in 1:length(tvec)]
    if C === Tuple{}
        return res1
    else
        tvec2 = C.parameters
        res2 = [tvec2.parameters[i]() for i in 1:length(tvec)]
        append!(res1,res2)
        return res1
    end
end

function _static_specs_types(x::Type{ThermodynamicState{S,C}}) where {M,C,S}
    tvec = S.parameters
    res1 = [tvec[i].parameters[1] for i in 1:length(tvec)]
    if C === Tuple{}
        return res1
    else
        res2 = collect(C.parameters)
        append!(res1,res2)
        return res1
    end
end

@generated function static_amount_type(x::T) where T<: ThermodynamicState
    sps = _static_specs_types(x)
    mass_int = mapreduce(_reduce_check_mass,+,sps)
    val = mass_int
    modval = mod(val,10)
    if val in (210,220)
        res =  (AMOUNT_CONST_FIRST[val],AMOUNT_CONST_FIRST[val])
    else
        res =  (AMOUNT_CONST_FIRST[val-modval],AMOUNT_CONST_SECOND[modval])
    end
    return :($res)
end


function _state_from_type(x::Type{ThermodynamicState{S,C}}) where {S,C}
    state_specs =  _static_specs(x)
    has_compound = false
    for spec in state_specs
        if spec isa MaterialCompounds
            has_compound = true
            comp_spec = spec 
            break
        elseif spec isa HumiditySpec
            has_compound = true
            comp_spec = spec 
            break
        end
    end
    if !has_compound
        comp_spec = SingleComponent()
    end
    sorted_specs = sort(state_specs,by= _stateorder)
    res = (sorted_specs[1],sorted_specs[2],comp_spec)
end

function _static_get_spec_idx(x::Type{ThermodynamicState{S,C}},sp::Type{S2})::Int where {S,C,S2}
    state_specs =  _static_specs(x)
    for (i,spec) in pairs(state_specs)
        if spec isa sp
            return i
        end
    end
    return 0
end

"""
    state_type(x::ThermodynamicState)

Returns a sorted tuple of thermodynamic specifications of the input state,
 the tuple has 3 elements, corresponding to two general specifications and a last amount specification,
  indicating if the state is single or multicomponent.


"""
@generated function state_type(x::ThermodynamicState)
    res = _state_from_type(x)
    return :($res)
end

function get_spec(val::T,st::Tuple) where T<:AbstractSpec
    for spec in st
        if specification(spec) isa T
            return spec
        end
    end
    return nothing
end

function get_spec(val::Type{T},st::Tuple) where T<: AbstractSpec
    for spec in st
        if specification(spec) isa T
            return spec
        end
    end
    return nothing
end

get_spec(val,st::ThermodynamicState) = get_spec(val,st.specs)
 
function throw_get_spec(val,st::ThermodynamicState)
    return static_get_spec(st,val)
end



function has_spec(val::Type{T},st::ThermodynamicState) where T<: AbstractSpec
    static_has_spec(st,val)
end


function amount_type(st::T) where T<:ThermodynamicState
    return static_amount_type(st)
end

@generated function static_get_spec(st::ThermodynamicState,val::Type{T}) where T<:AbstractSpec
    idx =_static_get_spec_idx(st,T)
    if iszero(idx)
        return :(nothing)
    else
        return :(st.specs[$idx])
    end
end

@generated function static_has_spec(st::ThermodynamicState,val::Type{T}) where T<:AbstractSpec
    idx =_static_get_spec_idx(st,T)
    res = !iszero(idx)
    return :($res)
end
