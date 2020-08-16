module Types
    abstract type AbstractSpec end
    abstract type SpecModifier <: AbstractSpec  end

    struct SingleComponent <: AbstractSpec end
    struct OneMol <: AbstractSpec end

    abstract type MassBasisModifier <: SpecModifier end
    struct MOLAR <: MassBasisModifier end
    struct MASS <: MassBasisModifier end
    struct TOTAL <: MassBasisModifier end

    abstract type CompoundModifier <: SpecModifier end
    struct FRACTION <: CompoundModifier end
    struct TOTAL_AMOUNT <: CompoundModifier end

    abstract type VolumeModifier <: SpecModifier end
    struct DENSITY <: VolumeModifier end
    struct VOLUME <: VolumeModifier end

    abstract type AbstractIntensiveSpec{T1<:MassBasisModifier} <: AbstractSpec  end
    abstract type AbstractTotalSpec <: AbstractSpec  end
    abstract type AbstractFractionSpec <: AbstractSpec  end
    abstract type CategoricalSpec <: AbstractSpec  end




    abstract type AbstractEnergySpec{T1} <: AbstractIntensiveSpec{T1} end

    struct Enthalpy{T} <: AbstractEnergySpec{T} end
    struct InternalEnergy{T} <: AbstractEnergySpec{T} end
    struct Gibbs{T} <: AbstractEnergySpec{T} end
    struct Helmholtz{T} <: AbstractEnergySpec{T} end

    struct Entropy{T} <: AbstractIntensiveSpec{T} end

    struct VolumeAmount{T1,T2<:VolumeModifier} <: AbstractIntensiveSpec{T1} end

    struct Pressure <: AbstractSpec end
    struct Temperature <: AbstractSpec end

    struct Mass <: AbstractTotalSpec end
    struct Moles <: AbstractTotalSpec end

    # those are vectors, 

    struct MaterialCompounds{T1<:MassBasisModifier,T2<:CompoundModifier} <: AbstractTotalSpec end

    struct MaterialAmount{T1<:MassBasisModifier} <: AbstractTotalSpec end


    struct PhaseFractions <: AbstractFractionSpec end
    struct VaporFraction <: AbstractFractionSpec end

    struct PhaseTag <: CategoricalSpec end
    struct TwoPhaseEquilibrium <: CategoricalSpec end

    struct Options <: CategoricalSpec end

    #not defined for now
    struct MolecularWeight <: AbstractSpec end

    export AbstractSpec 
    export SpecModifier 
    export AbstractIntensiveSpec 
    export AbstractTotalSpec 
    export AbstractFractionSpec 
    export CategoricalSpec 
    export SingleComponent 
    export OneMol  
    export MOLAR  
    export MASS  
    export TOTAL  
    export FRACTION 
    export TOTAL_AMOUNT 
    export DENSITY  
    export VOLUME  
    export AbstractEnergySpec  
    export Enthalpy 
    export InternalEnergy
    export Gibbs  
    export Helmholtz  
    export Entropy  
    export VolumeAmount 
    export Pressure 
    export Temperature 
    export Mass 
    export Moles 
    export MaterialCompounds 
    export MaterialAmount 
    export PhaseFractions 
    export VaporFraction 
    export PhaseTag 
    export TwoPhaseEquilibrium 
    export Options 
    export MolecularWeight

end

using .Types

const KW_TO_SPEC = IdDict{Symbol,Any}(
:h =>  Enthalpy{MOLAR}()
,:g =>  Gibbs{MOLAR}()
,:a =>  Helmholtz{MOLAR}()
,:u =>  InternalEnergy{MOLAR}()

,:mol_h =>  Enthalpy{MOLAR}()
,:mol_g =>  Gibbs{MOLAR}()
,:mol_a =>  Helmholtz{MOLAR}()
,:mol_u =>  InternalEnergy{MOLAR}()

,:mass_h =>  Enthalpy{MASS}()
,:mass_g =>  Gibbs{MASS}()
,:mass_a =>  Helmholtz{MASS}()
,:mass_u =>  InternalEnergy{MASS}()

,:total_h =>  Enthalpy{TOTAL}()
,:total_g =>  Gibbs{TOTAL}()
,:total_a =>  Helmholtz{TOTAL}()
,:total_u =>  InternalEnergy{TOTAL}()

,:s =>  Entropy{MOLAR}()
,:mol_s =>  Entropy{MOLAR}()
,:mass_s =>  Entropy{MASS}()
,:total_s =>  Entropy{TOTAL}()

,:p =>  Pressure()
,:P =>  Pressure()
,:t =>  Temperature()
,:T => Temperature()

,:v =>  VolumeAmount{MOLAR,VOLUME}()
,:mol_v =>  VolumeAmount{MOLAR,VOLUME}()
,:mass_v =>  VolumeAmount{MASS,VOLUME}()
,:total_v =>  VolumeAmount{TOTAL,VOLUME}()

,:rho =>  VolumeAmount{MOLAR,DENSITY}() 
,:mol_rho =>  VolumeAmount{MOLAR,DENSITY}() 
,:mass_rho =>  VolumeAmount{MASS,DENSITY}() 

,:ρ =>  VolumeAmount{MOLAR,DENSITY}()
,:mol_ρ =>  VolumeAmount{MOLAR,DENSITY}()
,:mass_ρ =>  VolumeAmount{MASS,DENSITY}() 
,:mass =>  MaterialAmount{MASS}()
,:moles =>  MaterialAmount{MOLAR}()
,:xn =>  MaterialCompounds{MOLAR,FRACTION}()
,:xm =>   MaterialCompounds{MASS,FRACTION}()
,:n =>  MaterialCompounds{MOLAR,TOTAL_AMOUNT}()
,:m =>   MaterialCompounds{MASS,TOTAL_AMOUNT}()
,:mw =>  MolecularWeight()

,:vfrac =>  VaporFraction() #looking for better name
,:phase_fracs =>  PhaseFractions() #looking for better name

,:phase =>PhaseTag()

,:sat => TwoPhaseEquilibrium()
,:vle => TwoPhaseEquilibrium()
,:lle => TwoPhaseEquilibrium()

,:single_component => SingleComponent()
,:one_mol => OneMol()
,:options => Options()
)


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

function default_units(::Type{MaterialCompounds{T,FRACTION}}) where T
    return Unitful.NoUnits
end

function default_units(x::Type{MaterialCompounds{T,TOTAL_AMOUNT}}) where T
    return default_units(T)
end



default_units(x::Type{PhaseFractions}) = Unitful.NoUnits
default_units(x::Type{VaporFraction}) = Unitful.NoUnits

struct Spec{T <: AbstractSpec,U}
    type::T
    val::U
end

value(s::Spec) = s.val
specification(s::Spec) = s.type



#preferred stripped
function _ups(x,normalize_units=true)
    if normalize_units
        return Unitful.ustrip(Unitful.upreferred(x))
    else
        return x
    end
end

function _ups(x::AbstractVector,normalize_units=true)
    if normalize_units
        return Unitful.ustrip.(Unitful.upreferred.(x))
    else
        return x
    end
end

#check units and normalizes
function check_and_norm(::SP,val::U,normalize_units::Bool=true) where {SP<:AbstractSpec,U<:Unitful.Quantity}
    if dimension(default_units(SP)) == dimension(val)
    return _ups(val,normalize_units)
    else
        throw(ArgumentError("the input value is not a type of " * string(SP)))
    end
end
function check_and_norm_vec(::SP,val::U,normalize_units::Bool=true) where {SP<:AbstractSpec,U<:AbstractVector{<: Unitful.Quantity}}
    if dimension(default_units(SP)) == dimension(eltype(val))
    return _ups(val,normalize_units)
    else
        throw(ArgumentError("the input value is not a type of " * string(SP)))
    end
end


function check_and_norm_vec(::SP,val::U,normalize_units::Bool=true) where {SP<:AbstractSpec,U<:Real}
    throw(ArgumentError("invalid material compounds specification, provide a vector."))
end

function check_and_norm_vec(::SP,val::U,normalize_units::Bool=true) where {SP<:AbstractSpec,U<:AbstractVector{<:Real}}
    return _ups(val,normalize_units)
end

function check_and_norm(::SP,val::U,normalize_units::Bool=true) where {SP<:AbstractSpec,U<:Real}
    return _ups(val,normalize_units)

end


function spec(sp::AbstractIntensiveSpec, val,normalize_units::Bool=true)
    val = check_and_norm(sp,val,normalize_units)
    return Spec(sp,val)
end

function spec(sp::Union{Pressure,Temperature}, val,normalize_units::Bool=true)
    val = check_and_norm(sp,val,normalize_units)
    return Spec(sp,val)
end

function spec(sp::MaterialAmount, val,normalize_units::Bool=true)
    val = check_and_norm(sp,val,normalize_units)
    return Spec(sp,val)
end

function spec(sp::MaterialCompounds{T1,TOTAL_AMOUNT}, val::V,normalize_units::Bool=true) where {T1,V<:AbstractVector} 
    val = check_and_norm_vec(sp,val,normalize_units)
    return Spec(sp,val)
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
        throw(ArgumentError("the the vector of values is not a valid fraction"))
    end
end

function spec(t::AbstractFractionSpec, u::Real,normalize_units::Bool=true)
    if _is_frac(u)
        return Spec(t, float(u))
    else
        throw(ArgumentError("the value " * string(u) * " is not between 0 and 1."))
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

struct Specs{T1,T2}
    amount_type::T1
    specs::T2
    checked::Bool
end

Base.values(s::Specs) = s.specs

function get_spec(val,specs::Tuple)
    for spec in specs
        if typeof(specification(spec)) <: typeof(val) 
            return spec
        end
    end
    return nothing
end

function get_spec(val::Type{T},specs::Tuple) where T<: AbstractSpec
    for spec in specs
        if typeof(specification(spec)) <: T
            return spec
        end
    end
    return nothing
end

get_spec(val,spec::Specs) = get_spec(val,spec.specs)
 
function throw_get_spec(val,specs::Specs)
    res = get_spec(val,specs)
    if res !== nothing
        return res
    else
        throw(error(string(val) * " not found in specifications"))
    end
end

@inline function Base.getindex(specs::Specs,val::Int)
    return specs.specs[val]
end

#linear search on values
@inline function Base.getindex(specs::Specs,val::AbstractSpec)
    return throw_get_spec(val,specs)
end

function has_spec(val,tup::Tuple)::Bool
    @inbounds for i in 1:length(tup)
        if typeof(specification(tup[i])) <: typeof(val) 
            return true
        end
    end
    return false
end

function has_spec(val::Type{T},tup::Tuple) where T<: AbstractSpec
    @inbounds for i in 1:length(tup)
        if typeof(specification(tup[i])) <: T 
            return true
        end
    end
    return false
end

has_spec(val,sp::Specs) = has_spec(val,sp.specs)

#returns the specification symbol in multicomponent
#returns :singlecomponent if there are no specifications
#throws error if more than one specification is found

const AMOUNT_CONST =Dict{Int,Any}(
    0 => (SingleComponent(),OneMol())
    ,1 => (SingleComponent(),MaterialAmount{MOLAR}())
    ,2 => (SingleComponent(),MaterialAmount{MASS}())
    ,110 => (MaterialCompounds{MOLAR,FRACTION}(),OneMol())
    ,111 => (MaterialCompounds{MOLAR,FRACTION}(),MaterialAmount{MOLAR}())
    ,112 => (MaterialCompounds{MOLAR,FRACTION}(),MaterialAmount{MASS}())
    ,120 => (MaterialCompounds{MASS,FRACTION}(),OneMol())
    ,121 => (MaterialCompounds{MASS,FRACTION}(),MaterialAmount{MOLAR}())
    ,122 => (MaterialCompounds{MASS,FRACTION}(),MaterialAmount{MASS}())
    ,210 => (MaterialCompounds{MOLAR,TOTAL_AMOUNT}(),MaterialCompounds{MOLAR,TOTAL_AMOUNT}())
    ,220 => (MaterialCompounds{MASS,TOTAL_AMOUNT}(),MaterialCompounds{MASS,TOTAL_AMOUNT}())
    )

    AMOUNT_CONST2 = Dict{Int64,Symbol}(
        110 => :xn
        ,111 => :xn
        ,112 => :xn
        ,120 => :xm
        ,121 => :xm
        ,122 => :xm
        ,210 => :n
        ,220 => :m
    )

#reduce-based mass transform, to check properties


_reduce_check_mass(x::Spec) = 0
_reduce_check_mass(x::Spec{MaterialCompounds{MOLAR,FRACTION}}) = 110
_reduce_check_mass(x::Spec{MaterialCompounds{MASS,FRACTION}}) = 120
_reduce_check_mass(x::Spec{MaterialCompounds{MOLAR,TOTAL_AMOUNT}}) = 210
_reduce_check_mass(x::Spec{MaterialCompounds{MASS,TOTAL_AMOUNT}}) = 220
_reduce_check_mass(x::Spec{MaterialAmount{MOLAR}}) = 1
_reduce_check_mass(x::Spec{MaterialAmount{MASS}}) = 2

#for symbols
_reduce_check_mass(x::Symbol) = _reduce_check_mass(Val(x))
_reduce_check_mass(::Val{:xn})=110
_reduce_check_mass(::Val{:xm})=120
_reduce_check_mass(::Val{:n})=210
_reduce_check_mass(::Val{:m})=220
_reduce_check_mass(::Val{:moles})=1
_reduce_check_mass(::Val{:mass})=2
_reduce_check_mass(::Val{:n})=210
_reduce_check_mass(::Val{T} where T)=0

#reduce-based phase transform, to check properties

_reduce_check_phase(x::Spec) = 0
_reduce_check_phase(x::Spec{TwoPhaseEquilibrium}) = 1
_reduce_check_phase(x::Spec{VaporFraction}) = 1
_reduce_check_phase(x::Spec{PhaseFractions}) = 10

_reduce_check_phase(x::Symbol) = _reduce_check_phase(Val(x))
_reduce_check_phase(::Val{:vle})=1
_reduce_check_phase(::Val{:lle})=1
_reduce_check_phase(::Val{:sat})=1
_reduce_check_phase(::Val{:vfrac})=1
_reduce_check_phase(::Val{:phase_fracs})=10
_reduce_check_phase(::Val{T} where T)=0


function _specs_components(kwargs::NamedTuple)::Int64
    return mapreduce(_reduce_check_mass,+,keys(kwargs))
end
#gives the appropiate symbol to extract amount of matter


#this needs future work to specify
function _specs_phase_basis(kwargs::NamedTuple)::Int64
    return mapreduce(_reduce_check_phase,+,keys(kwargs))
end

function _specs_C(kwargs::NamedTuple,kw::Int64)::Int64
    if kw < 10
        return 1
    else
        return length(getproperty(kwargs,AMOUNT_CONST2[kw]))
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

function _specs_P(kwargs::NamedTuple,kw::Int64)::Int64
    if (kw == 0)
        return 1
    elseif kw == 1
        return 2
    else
        return length(kwargs.phase_fracs)
    end
end



function _specs_F(kwargs::NamedTuple,mass_basis::Int64,phase_basis::Int64)::Int64
    F = length(kwargs)
    opt = mapreduce(_reduce_check_opt,+,keys(kwargs))
    F = F - opt
    
    if mass_basis in (1,2,110,120,210,220) #one_specified
        F = F -1
    elseif mass_basis in (111,112,121,122) #two specified
        F = F-2
    elseif mass_basis == 0
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

function specs(;check=true,normalize_units=true,kwargs...)
    f0 = k -> spec(KW_TO_SPEC[k],getproperty(kwargs.data,k),normalize_units)
    tup = map(f0,keys(kwargs.data))
    if check == true
        mass_basis = _specs_components(kwargs.data)
        phase_basis = _specs_phase_basis(kwargs.data)
        C = _specs_C(kwargs.data,mass_basis)
        P = _specs_P(kwargs.data,phase_basis)
        F = _specs_F(kwargs.data,mass_basis,phase_basis)
        DF = C - P + 2 - F#behold, the gibbs phase rule!
        if DF<0
            throw(error("the variables are overspecified."))
        elseif DF>0
            throw(error("the variables are underspecified."))
        else
            return Specs(AMOUNT_CONST[mass_basis],tup,true)
        end
    else
        return Specs(nothing,tup,false)
    end
end

function specs_grid(;check=true,normalize_units=true,kwargs...)
    kw = NamedTuple{keys(kwargs.data)}
    function f0(v)
        _kwargs = kw(v)
        specs(;check=check,normalize_units=normalize_units,_kwargs...)
    end
    return map(f0,Iterators.product(values(kwargs.data)...))
end

function _specs_components(tup::Tuple)::Int64
  return mapreduce(_reduce_check_mass,+,tup)
end

function _specs_phase_basis(tup::Tuple)::Int64
    return mapreduce(_reduce_check_phase,+,tup)
end


_reduce_mass_spec_c(x::Spec) = 0
_reduce_mass_spec_c(x::Spec{MaterialCompounds}) = length(value(x))

function _specs_C(tup::Tuple,kw::Int64)::Int64
    if kw < 10
        return 1
    else
        return mapreduce(_reduce_mass_spec_c,+,tup)
    end
end

_reduce_mass_spec_p(x::Spec) = 0
_reduce_mass_spec_p(x::Spec{PhaseFractions}) = length(value(x))

function _specs_P(tup::Tuple,kw::Int64)::Int64
    if (kw == 0)
        return 1
    elseif kw == 1
        return 2
    else
        return mapreduce(_reduce_mass_spec_p,+,tup)
    end
end

function _specs_F(tup::Tuple,mass_basis::Int64,phase_basis::Int64)::Int64
    F = length(tup)
    opt = mapreduce(_reduce_check_opt,+,tup)
    F = F - opt

    
    if mass_basis in (1,2,110,120,210,220) #one_specified
        F = F -1
    elseif mass_basis in (111,112,121,122) #two specified
        F = F-2
    elseif mass_basis == 0
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

function specs(args::Vararg{Spec};check=true)
    
    if check == true
        all_unique = spec_tuple_unique(args) #this function is the main speed bottleneck
        if !all_unique
            throw(error("specifications are not unique."))
        end
        mass_basis = _specs_components(args)
        phase_basis = _specs_phase_basis(args)
        C = _specs_C(args,mass_basis)
        P = _specs_P(args,phase_basis)
        F = _specs_F(args,mass_basis,phase_basis)
        DF = C - P + 2 - F#behold, the gibbs phase rule!
        if DF<0
            throw(error("the variables are overspecified."))
        elseif DF>0
            throw(error("the variables are underspecified."))
        else
            return Specs(AMOUNT_CONST[mass_basis],args,true)
        end
    else
        return Specs(nothing,args,false)
    end
end
