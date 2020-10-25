"""
    default_units(::Union{Type{AbstractSpec},Type{property_function}})

returns a default unit of a spec

# Examples
```julia-repl
default_units(pressure) #function name
Pa

default_units(Gibbs{MOLAR}) #spec type
J/mol
```
"""
function default_units end

# total units, function necessary to defining a new spec
default_units(x::AbstractSpec) = default_units(typeof(x))
default_units(::Type{MOLAR}) = u"mol"
default_units(::Type{MASS}) = u"kg"
default_units(::Type{TOTAL}) = Unitful.NoUnits

function default_units(::Type{T1}) where T1
    return Unitful.NoUnits
end


function default_units(::Type{T1}) where T1 <: AbstractEnergySpec{T2} where T2
    return u"J"/default_units(T2)
end

function default_units(::Type{Entropy{T}}) where T
    return u"J/K"/default_units(T)
end

default_units(::Type{Pressure}) = u"Pa"
default_units(::Type{Temperature}) = u"K"

function default_units(::Type{VolumeAmount{T,VOLUME}}) where T
    return u"m^3"/default_units(T)
end

function default_units(::Type{VolumeAmount{T,DENSITY}}) where T
    return default_units(T)/u"m^3"
end

function default_units(::Type{MaterialAmount{T1}}) where T1
    return default_units(T1)
end

function default_units(::Type{SingleComponent})
    return Unitful.NoUnits
end

function default_units(::Type{MaterialCompounds{T,FRACTION}}) where T
    return Unitful.NoUnits
end

function default_units(x::Type{MaterialCompounds{T,TOTAL_AMOUNT}}) where T
    return default_units(T)
end

default_units(x::Type{PhaseFractions}) = Unitful.NoUnits
default_units(x::Type{VaporFraction}) = Unitful.NoUnits

default_units(::Type{HumiditySpec{HumidityDewPoint}}) = u"K"
default_units(::Type{HumiditySpec{HumidityRatio}}) = Unitful.NoUnits
default_units(::Type{HumiditySpec{WetBulbTemperature}}) = u"K"
default_units(::Type{HumiditySpec{RelativeHumidity}}) = Unitful.NoUnits
default_units(::Type{HumiditySpec{MolarHumidity}}) = Unitful.NoUnits
default_units(::Type{HumiditySpec{MassHumidity}}) = Unitful.NoUnits


struct Spec{T <: AbstractSpec,U}
    type::T
    val::U
end

value(s::Spec) = s.val
specification(s::Spec) = s.type

#preferred stripped
function _ups(x,normalize_units=true)
    if normalize_units
        return Unitful.ustrip.(Unitful.upreferred.(x))
    else
        return x
    end
end

#check units and normalizes
function check_and_norm(::SP,val::U,normalize_units::Bool=true) where {SP<:AbstractSpec,U}
    is_real(val) && return val
    if dimension(default_units(SP)) == dimension(eltype(val))
        return _ups(val,normalize_units)
    else
        throw(ArgumentError("the input value is not a type of " * string(SP)))
    end
end


#special cases, VariableSpec and Misssing
function spec(sp::AbstractSpec, val::VariableSpec,normalize_units::Bool=true)
    return Spec(sp,val)
end

function spec(sp::AbstractSpec, val::Missing,normalize_units::Bool=true)
    return Spec(sp,missing)
end

function spec(sp::AbstractIntensiveSpec,  val::Number,normalize_units::Bool=true)
    val = check_and_norm(sp,val,normalize_units)
    return Spec(sp,val)
end

function spec(sp::Union{Pressure,Temperature},  val::Number,normalize_units::Bool=true)
    val = check_and_norm(sp,val,normalize_units)
    return Spec(sp,val)
end

function spec(sp::HumiditySpec{T},  val::Number,normalize_units::Bool=true) where T <: Union{HumidityDewPoint,WetBulbTemperature}
    val = check_and_norm(sp,val,normalize_units)
    return Spec(sp,val)
end

function spec(sp::HumiditySpec,  val::Number,normalize_units::Bool=true)
    if _is_frac(val)
        return Spec(sp, float(val))
    else
        throw(ArgumentError("the value " * string(val) * " is not between 0 and 1."))
    end
end


function spec(sp::Pressure, val::Nothing,normalize_units::Bool=true)
    return Spec(sp,val)
end

function spec(sp::MaterialAmount, val::Number,normalize_units::Bool=true)
    val = check_and_norm(sp,val,normalize_units)
    return Spec(sp,val)
end

function spec(sp::MaterialCompounds{T1,TOTAL_AMOUNT}, val::V,normalize_units::Bool=true) where {T1,V<:AbstractVector}
    if isone(length(val))
        return ArgumentError("invalid material compounds specification, provide a vector.")
    else
        val = check_and_norm(sp,val,normalize_units)
        return Spec(sp,val)
    end
end

_is_frac(x::Real) = (zero(x) <= x <= one(x))
function _is_frac_vec(u::T1) where T1<:AbstractVector{<:Real}
    return all(x -> x ≥ zero(x), u) && isapprox(sum(u), one(eltype(u)))
end

function spec(sp::MaterialCompounds{T1,FRACTION}, val::V,normalize_units::Bool=true) where {T1,V<:AbstractVector}
    if _is_frac_vec(val)
        return Spec(sp,val)
    else
        throw(ArgumentError("the the vector of values is not a valid fraction"))
    end
end

function spec(t::AbstractFractionSpec, u::AbstractVector{Real},normalize_units::Bool=true)
    if _is_frac_vec(u)
        return Spec(t, u)
    else
        return ArgumentError("the the vector of values is not a valid fraction")
    end
end

function spec(t::AbstractFractionSpec, u::Real,normalize_units::Bool=true)
    if _is_frac(u)
        return Spec(t, float(u))
    else
        return ArgumentError("the value " * string(u) * " is not between 0 and 1.")
    end
end

#catch all for non numerical specs
function spec(t::T, u,normalize_units::Bool=true) where T<:CategoricalSpec
    return Spec(t,u)
end

function spec(;kwargs...)
    if hasproperty(kwargs.data,:normalize_units)
        norm_units = getproperty(kwargs.data,:normalize_units)
        if length(kwargs.data)==2
            kw1,kw2 = keys(kwargs)
            if kw1 == :normalize_units
                return spec(KW_TO_SPEC[kw2],last(Base.values(kwargs.data)),norm_units)
            else
                return spec(KW_TO_SPEC[kw1],first(Base.values(kwargs.data)),norm_units)
            end
        else
            throw("invalid keyword combination.")
        end
    else
        if length(kwargs.data)==1
            return spec(KW_TO_SPEC[first(keys(kwargs.data))],first(Base.values(kwargs.data)))
        else
            throw("invalid keyword combination.")
        end
    end
end
"""
    struct ThermodynamicState{S,C}

Struct containing thermodynamic specifications, constituting a state of an specific amount of matter.

"""
struct ThermodynamicState{S,C}
    specs::S
    callables::C
end

ThermodynamicState(st::Tuple) = ThermodynamicState(st,())

Base.values(s::ThermodynamicState) = s.specs

#reduce-based mass transform, to check properties

_reduce_check_mass(x::Spec) = 0
_reduce_check_mass(x::Spec{MaterialCompounds{MOLAR,FRACTION}}) = 110
_reduce_check_mass(x::Spec{MaterialCompounds{MASS,FRACTION}}) = 120
_reduce_check_mass(x::Spec{MaterialCompounds{MOLAR,TOTAL_AMOUNT}}) = 210
_reduce_check_mass(x::Spec{MaterialCompounds{MASS,TOTAL_AMOUNT}}) = 220
_reduce_check_mass(x::Spec{MaterialAmount{MOLAR}}) = 1
_reduce_check_mass(x::Spec{MaterialAmount{MASS}}) = 2

_reduce_check_mass(x::Type{MaterialCompounds{MOLAR,FRACTION}}) = 110
_reduce_check_mass(x::Type{MaterialCompounds{MASS,FRACTION}}) = 120
_reduce_check_mass(x::Type{MaterialCompounds{MOLAR,TOTAL_AMOUNT}}) = 210
_reduce_check_mass(x::Type{MaterialCompounds{MASS,TOTAL_AMOUNT}}) = 220

_reduce_check_mass(x::Type{HumiditySpec{WetBulbTemperature}})  = 310
_reduce_check_mass(x::Type{HumiditySpec{HumidityRatio}}) = 320
_reduce_check_mass(x::Type{HumiditySpec{MolarHumidity}}) = 330
_reduce_check_mass(x::Type{HumiditySpec{MassHumidity}}) = 340
_reduce_check_mass(x::Type{HumiditySpec{RelativeHumidity}}) = 350
_reduce_check_mass(x::Type{HumiditySpec{HumidityDewPoint}}) = 360


#mass modifiers
_reduce_check_mass(x::Type{MaterialAmount{MOLAR}}) = 1
_reduce_check_mass(x::Type{MaterialAmount{MASS}}) = 2
_reduce_check_mass(x::Type) = 0


_reduce_check_mass(x::Symbol) = _reduce_check_mass(Val(x))
_reduce_check_mass(::Val{:xn})=110
_reduce_check_mass(::Val{:xm})=120
_reduce_check_mass(::Val{:n})=210
_reduce_check_mass(::Val{:m})=220
_reduce_check_mass(::Val{:moles})=1
_reduce_check_mass(::Val{:mass})=2
_reduce_check_mass(::Val{T} where T)=0

_reduce_check_mass(::Val{:hum_wetbulb}) = 310
_reduce_check_mass(::Val{:hum_ratio}) = 320
_reduce_check_mass(::Val{:hum_molfrac}) = 330
_reduce_check_mass(::Val{:hum_massfrac}) = 340
_reduce_check_mass(::Val{:rel_hum}) = 350
_reduce_check_mass(::Val{:hum_dewpoint}) = 360

#=

=#
#reduce-based phase transform, to check properties

_reduce_check_phase(x::Spec) = 0
_reduce_check_phase(x::Spec{TwoPhaseEquilibrium}) = 1
_reduce_check_phase(x::Spec{VaporFraction}) = 1
_reduce_check_phase(x::Spec{PhaseFractions}) = 10

_reduce_check_phase(x::Symbol) = _reduce_check_phase(Val(x))
_reduce_check_phase(::Val{:vle})=1
_reduce_check_phase(::Val{:lle})=1
_reduce_check_phase(::Val{:sat})=1
_reduce_check_phase(::Val{:quality})=1
_reduce_check_phase(::Val{:vfrac})=1
_reduce_check_phase(::Val{:phase_fracs})=10
_reduce_check_phase(::Val{T} where T)=0

#function to correctly dispatch on the function terms
keys_or_tuple(x::Tuple) = x
keys_or_tuple(x::NamedTuple) = keys(x)

function _specs_components(args)::Int64
    return mapreduce(_reduce_check_mass,+,keys_or_tuple(args))
end

#this needs future work to specify
function _specs_phase_basis(args)::Int64
    return mapreduce(_reduce_check_phase,+,keys_or_tuple(args))
end

_reduce_mass_spec_c(x::Spec) = 0
_reduce_mass_spec_c(x::Spec{MaterialCompounds}) = length(value(x))

_specs_C(kwargs::NamedTuple,kw::Int64) = 1
_specs_C(tup::Tuple,kw::Int64)::Int64 = 1

_reduce_mass_spec_p(x::Spec) = 0
_reduce_mass_spec_p(x::Spec{PhaseFractions}) = length(value(x))

function _specs_P(kw::Int64)::Int64
    if iszero(kw)
        return 1
    elseif isone(kw)
        return 2
    else
        return 0
    end
end

#optional values not counted during degrees of freedom calc

_reduce_check_opt(x::Spec) = 0
_reduce_check_opt(x::Spec{PhaseTag}) = 1
_reduce_check_opt(x::Spec{Options}) = 1

_reduce_check_opt(x::Symbol) = _reduce_check_opt(Val(x))
_reduce_check_opt(::Val{:phase})=1
_reduce_check_opt(::Val{:options})=1
_reduce_check_opt(::Val{T} where T)=0

function _specs_F(args,mass_basis::Int64,phase_basis::Int64)::Int64
    F = length(args)
    opt = mapreduce(_reduce_check_opt,+,keys_or_tuple(args))
    F = F - opt
    mass_basis_mod = mod(mass_basis,10)
    if iszero(mass_basis)
    elseif (mass_basis in (1,2)) | iszero(mass_basis_mod) #one_specified
        F = F -1
    elseif mass_basis_mod in (1,2) #two specified
        F = F-2
    else
        throw(error("incorrect mass specification."))
    end

    if phase_basis == 0
    elseif phase_basis in (1,10)
    F = F - 1
    else
        throw(error("incorrect phase specification."))
    end
    return F
end


function spec_equal(x1::Spec{T},x2::Spec{T})::Bool where {T}
    return true
end

function spec_equal(x1::Spec{T1},x2::Spec{T2})::Bool where {T1,T2}
    return false
end
function spec_tuple_unique(a)::Bool
    len = length(a)
    if len == 1
        return true
    elseif len == 2
        return !spec_equal(a[1],a[2])
    elseif len == 3
        return !(spec_equal(a[1],a[2]) | spec_equal(a[3],a[2]) | spec_equal(a[3],a[2]))
    else
        @inbounds for i = 1:len
            for j = (i+1):len
                if spec_equal(a[i],a[j])
                    return false
                end
            end
        end
    end
    return true
end

function check_spec(args)
    mass_basis = _specs_components(args)
    phase_basis = _specs_phase_basis(args)
    C = _specs_C(args,mass_basis)
    P = _specs_P(phase_basis)
    F = _specs_F(args,mass_basis,phase_basis)
    DF = C - P + 2 - F#behold, the gibbs phase rule!
    if DF<0
        throw(error("the variables are overspecified by " * string(abs(DF)) * " degrees of freedom."))
    elseif DF>0
        throw(error("the variables are underspecified by "* string(abs(DF)) * " degrees of freedom."))
    else
        return mass_basis
    end
end
_is_variable_spec(x) = false
_is_variable_spec(x::VariableSpec) = true
_is_variable_spec(x::Spec{T,VariableSpec} where T) =true


"""
    state(;normalize_units=true,check=true,kwargs...)

Creates a `ThermodynamicState` from the arguments passed.
If one or more of the arguments is the value `VariableSpec()`, the state created will be a callable,
 returning a complete state when evaluated.

`Unitful` quantities are normalized to SI units and unit-stripped.
 This can be disabled with `normalize_units=false`

When creating a thermodynamic state, the input arguments are checked for consistency with the gibbs rule.
 This check can be skipped with `check=false`

"""
function state(;check=true,normalize_units=true,kwargs...)
    if !any(_is_variable_spec,values(kwargs.data))
        f0 = k -> spec(KW_TO_SPEC[k],getproperty(kwargs.data,k),normalize_units)
        tup = map(f0,keys(kwargs.data))
        callables = ()
    else
        tup = ((spec(KW_TO_SPEC[k],v,normalize_units) for (k, v) in kwargs if !_is_variable_spec(v))...,)
        callables = ((KW_TO_SPEC[k] for (k,v) in kwargs if _is_variable_spec(v))...,)
    end
    if check
        mass_basis = check_spec(kwargs.data)
        return ThermodynamicState(tup,callables)
    else
        return ThermodynamicState(tup,callables)
    end
end

function state(args::Vararg{Spec};check=true)
    if !any(_is_variable_spec,args)
        callables = ()
        tup = args
    else
        tup = ((arg for arg in args if !_is_variable_spec(arg))...,)
        callables  = ((arg.type for arg in args if _is_variable_spec(arg))...,)
    end
    if check
        all_unique = spec_tuple_unique(args) #this function is the main speed bottleneck
        if !all_unique
            throw(error("specifications are not unique."))
        end
        mass_basis = check_spec(args)

        return ThermodynamicState(tup,callables)
    else
        return ThermodynamicState(tup,callables)
    end
end

function (f::ThermodynamicState{T,Tuple{S1}})(x1::T1;normalize_units=true) where {T,S1,T1}
    return ThermodynamicState(
    (spec(only(f.callables),x1,normalize_units),
    f.specs...))
end

function (f::ThermodynamicState{S,Tuple{S1,S2}})(x1::T1,x2::T2;normalize_units=true) where {S,S1,S2,T1,T2}
    return ThermodynamicState(
            (spec(first(f.callables),x1),
            spec(last(f.callables),x2),
            f.specs...))
end

function (f::ThermodynamicState{S,Tuple{S1,S2,S3}})(x1::T1,x2::T2,x3::T3;normalize_units=true) where {S,S1,S2,S3,T1,T2,T3}
    return ThermodynamicState(
            (spec(first(f.callables),x1,normalize_units),
            spec(f.callables[2],x2,normalize_units),
            spec(last(f.callables),x3,normalize_units),
            f.specs...))
end

function (f::Spec{T,VariableSpec})(x;normalize_units = true) where T
    return spec(specification(f),x,normalize_units)
end

function Base.Dict(st::ThermodynamicState)
    return Dict(SPEC_TO_KW[specification(sp)] => value(sp) for sp in st.specs)
end

module StatePoints

#state points
    abstract type StatePoint end

    struct CriticalPoint <: StatePoint end
    struct NormalBoilingPoint <: StatePoint end
    struct TriplePoint <: StatePoint end
    struct NormalConditions <: StatePoint end
    struct StandardConditions <: StatePoint end

    export CriticalPoint,NormalBoilingPoint,TriplePoint,NormalConditions,StandardConditions
end


function state(x::StatePoints.StandardConditions)
    _t = Spec(Temperature(),273.15)
    _p = Spec(Pressure(),100000.0)
    return ThermodynamicState((_p,_t))
end

function state(x::StatePoints.NormalConditions)
    _t = Spec(Temperature(),293.15)
    _p = Spec(Pressure(),101325.0)
    return ThermodynamicState((_p,_t))
end

function state(x::StatePoints.NormalBoilingPoint)
    _p = Spec(Pressure(),101325.0)
    _sat = Spec(TwoPhaseEquilibrium(),true)

    return ThermodynamicState((_sat,_p))
end

pressure(model,x::StatePoints.NormalBoilingPoint,unit=u"Pa") = pressure(model,state(x),unit)
pressure(model,x::StatePoints.NormalConditions,unit=u"Pa") = pressure(model,state(x),unit)
pressure(model,x::StatePoints.StandardConditions,unit=u"Pa") = pressure(model,state(x),unit)
temperature(model,x::StatePoints.NormalConditions,unit=u"K") = temperature(model,state(x),unit)
temperature(model,x::StatePoints.StandardConditions,unit=u"K") = temperature(model,state(x),unit)
