
 #returns the specification symbol in multicomponent
#returns :singlecomponent if there are no specifications
#throws error if more than one specification is found


const AMOUNT_CONST =Dict{Int,Any}(
     0000 => OneMol()
    ,0001 => MaterialAmount{MOLAR}()
    ,0002 => MaterialAmount{MASS}() 
    ,0099 => SingleComponent()
    ,1100 => MaterialCompounds{MOLAR,FRACTION}()
    ,1200 => MaterialCompounds{MASS,FRACTION}()
    ,2100 => MaterialCompounds{MOLAR,TOTAL_AMOUNT}()
    ,2200 => MaterialCompounds{MASS,TOTAL_AMOUNT}()
    ,3100 => HumiditySpec{WetBulbTemperature}()
    ,3200 => HumiditySpec{HumidityRatio}()
    ,3300 => HumiditySpec{MolarHumidity}()
    ,3400 => HumiditySpec{MassHumidity}()
    ,3500 => HumiditySpec{RelativeHumidity}()
    ,3600 => HumiditySpec{HumidityDewPoint}()
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


#using specs here allows to delete all reduce_check_mass(::Type)
function _static_specs(x::Type{ThermodynamicState{S,C}}) where {S,C}
    @nospecialize S,C
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

@generated function static_amount_type(x::T) where T<: ThermodynamicState
    sps = _static_specs(x)
    val = sum(_reduce_check_mass(i) for i in sps)
    modval = mod(val,MATERIAL_SINGLE_MAX)
    if val in (2100,2200)
        _1 = val
        _2 = val
    else
        _1 = max(val-modval,99)
        _2 = modval
    end
    res =  (AMOUNT_CONST[_1],AMOUNT_CONST[_2])

    return :($res)
end


function _state_from_type(x::Type{T}) where T<: ThermodynamicState
    @nospecialize T
    state_specs =  _static_specs(T)
    has_compound = false
    for spec in state_specs
        if spec isa MaterialCompounds
            has_compound = true
            comp_spec = spec 
            break
        end
    end

    #double iteration, used in moist air to create an invalid
    #state where compound and humidity spec are present
    #in this case, the resulting type is MultiPT instead of HumPT
    for spec in state_specs
        if (spec isa HumiditySpec) & !has_compound
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
    return res
end

__unwrap(x::Type{T}) where T = x.parameters[1]
__unwrap(x) = x
function _dynamic_get_spec_idx(@nospecialize(x::Type{ThermodynamicState{S,C}}),@nospecialize(sp::Type{S2}))::Int where {S,C,S2}
    state_specs =  _static_specs(x)
    S3 = __unwrap(sp)
    for (i,spec) in pairs(state_specs)
        if typeof(spec) <: S3
            return i
        end
    end
    return 0
end

@generated function _static_get_spec_idx(@nospecialize(x::ThermodynamicState{S,C}),@nospecialize(sp::Type{S2}))::Int where {S,C,S2}
    res = _dynamic_get_spec_idx(x,sp)
    return :($res)
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

"""
    get_spec(val::AbstractSpec,st::ThermodynamicState)
    get_spec(val::Type{AbstractSpec},st::ThermodynamicState)

Returns the first specification of the type `val` present in the state `st`. returns `nothing` if not found.

## Example

```julia-repl
julia> a1 = state(n = [1.1,2.2,3.3],T=273.15,mol_h = 101.342)
ThermodynamicState with 3 properties:
  Molar amounts : [1.1, 2.2, 3.3][mol]
  Temperature : 273.15[K]
  Molar enthalpy : 101.342[J mol^-1]

julia> get_spec(Enthalpy,a1)
  spec(mol_h = 101.342[J mol^-1])

julia> get_spec(VolumeAmount,a1) == nothing
  true
```

"""
function get_spec end
get_spec(val::T,st::ThermodynamicState) where T<:AbstractSpec = get_spec(T,st)


function get_spec(val::Type{T},st::ThermodynamicState) where T<: AbstractSpec
    i = _static_get_spec_idx(st,T)
    i === 0 && return nothing
    return st.specs[i]
    #return nothing
end
 
function throw_get_spec(val,st::ThermodynamicState)::Spec
    res = get_spec(val,st)
    res === nothing && throw(error("$st does not have $val"))
    return res
end

"""
    has_spec(val::AbstractSpec,st::ThermodynamicState)::Bool
    has_spec(val::Type{AbstractSpec},st::ThermodynamicState)::Bool

checks if the specification of the type `val` present in the state `st`.

## Example

```julia-repl
julia> a1 = state(n = [1.1,2.2,3.3],T=273.15,mol_h = 101.342)
ThermodynamicState with 3 properties:
  Molar amounts : [1.1, 2.2, 3.3][mol]
  Temperature : 273.15[K]
  Molar enthalpy : 101.342[J mol^-1]

julia> has_spec(Enthalpy,a1)
  true

julia> has_spec(VolumeAmount,a1) == nothing
  false
```

"""
function has_spec end
function has_spec(val::Type{T},st::ThermodynamicState) where T<: AbstractSpec
    static_has_spec(st,val)
end

"""
    amount_type(st::ThermodynamicState)

given a state, returns a tuple with information about how is the amount of matter specified.

The first element of the tuple is related to the number of components present in the state. if it has only one component. it will return `SingleComponent`. if a multicomponent specification (`MaterialCompounds`)
is defined, it will return that type instead.

The second element is related to the type of material specified (molar or mass). if there isn't any specification, it will assume one mol, via the type `OneMol`. It will return a `MaterialAmount`  
if it is specified.

## Examples

```julia-repl
julia> tp = state(t=1,p=2)
ThermodynamicState with 2 properties:
  Temperature : 1[K]
  Pressure : 2[Pa]

julia> amount_type(tp)
(SingleComponent(), OneMol())

julia> tpn = state(t=1,p=2,mass=3)
ThermodynamicState with 3 properties:
  Temperature : 1[K]
  Pressure : 2[Pa]
  Mass : 3[kg]

julia> amount_type(tpn)
(SingleComponent(), MaterialAmount{MASS}())

julia> tpxm = state(t=1,p=2,mass=3,xn = [0.1,0.9])      
ThermodynamicState with 4 properties:
  Temperature : 1[K]
  Pressure : 2[Pa]
  Mass : 3[kg]
  Molar fraction : [0.1, 0.9]

julia> amount_type(tpxm)
(MaterialCompounds{MOLAR, FRACTION}(), MaterialAmount{MASS}())

julia> tpm = state(t=1,p=2,m = [1,2])
ThermodynamicState with 3 properties:
  Temperature : 1[K]
  Pressure : 2[Pa]
  Mass amounts : [1, 2][kg]

julia> amount_type(tpm)
(MaterialCompounds{MASS, TOTAL_AMOUNT}(), MaterialCompounds{MASS, TOTAL_AMOUNT}())
```
"""
function amount_type(st::T) where T<:ThermodynamicState
    return static_amount_type(st)
end

@generated function static_has_spec(st::ThermodynamicState,val::Type{T}) where T<:AbstractSpec
    idx =_static_get_spec_idx(st,T)
    res = !iszero(idx)
    return :($res)
end

reduce_numtype(::T) where T <: Spec = reduce_numtype(T)
reduce_numtype(::Type{Spec{V,T}}) where {V,T<:Unitful.Quantity{W}} where W = W  
reduce_numtype(::Type{Spec{V,T}}) where {V,T<:AbstractVector{<:Unitful.Quantity{W}}} where W = W   
reduce_numtype(::Type{Spec{V,T}}) where {V,T<:AbstractVector{<:W}} where W = W   
reduce_numtype(::Type{Spec{V,T}}) where {V,T} = Bool  
reduce_numtype(::Type{Spec{V,T}}) where {V,T<:Number} = T
reduce_numtype(::Type{T}) where T<:Number = T
reduce_numtype(A,B) = promote_type(reduce_numtype(A),reduce_numtype(B))
function numtype(st::T) where T <: ThermodynamicState
    res = foldl(reduce_numtype,st.specs)
    return res
end
_one(st::T) where T <: ThermodynamicState = one(numtype(st))