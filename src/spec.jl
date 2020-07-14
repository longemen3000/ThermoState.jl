abstract type AbstractSpec end
abstract type AbstractIntensiveSpec <: AbstractSpec  end
abstract type AbstractTotalSpec <: AbstractSpec  end
abstract type AbstractFractionSpec <: AbstractSpec  end
abstract type CategoricalSpec <: AbstractSpec  end

struct Enthalpy <: AbstractIntensiveSpec end
struct Entropy <: AbstractIntensiveSpec end
struct InternalEnergy <: AbstractIntensiveSpec end
struct Gibbs <: AbstractIntensiveSpec end
struct Helmholtz <: AbstractIntensiveSpec end
struct VolumeAmount <: AbstractIntensiveSpec end

abstract type AbstractVolumeSpec <: AbstractIntensiveSpec end
struct MolVolume <: AbstractVolumeSpec end
struct MolDensity <: AbstractVolumeSpec end
struct MassVolume <: AbstractVolumeSpec end
struct MassDensity <: AbstractVolumeSpec end
struct TotalVolume <: AbstractVolumeSpec end


struct Pressure <: AbstractIntensiveSpec end
struct Temperature <: AbstractIntensiveSpec end

struct Mass <: AbstractTotalSpec end
struct Moles <: AbstractTotalSpec end

# those are vectors, 
struct MassNumbers <: AbstractTotalSpec end
struct MolNumbers <: AbstractTotalSpec end
struct MaterialCompounds <: AbstractTotalSpec end
struct MaterialAmount <: AbstractTotalSpec end


struct MassFractions <: AbstractFractionSpec end
struct MolFractions <: AbstractFractionSpec end
struct PhaseFractions <: AbstractFractionSpec end
struct VaporFraction <: AbstractFractionSpec end
# for now,in a different category
struct MolecularWeight <: AbstractSpec end
struct PurePhase <: CategoricalSpec end

# total units, function necessary to defining a new spec
total_units(x::Enthalpy) = u"J"
total_units(x::Entropy) = u"J/K"
total_units(x::InternalEnergy) = u"J"
total_units(x::Gibbs) = u"J"
total_units(x::Helmholtz) = u"J"
total_units(x::TotalVolume) = u"m^3"

total_units(x::Pressure) = u"Pa"
total_units(x::Temperature) = u"K"
total_units(x::Mass) = u"kg"
total_units(x::Moles) = u"mol"
total_units(x::MassNumbers) = u"kg"
total_units(x::MolNumbers) = u"mol"

total_units(x::MassFractions) = Unitful.NoUnits
total_units(x::MolFractions) = Unitful.NoUnits
total_units(x::PhaseFractions) = Unitful.NoUnits
total_units(x::VaporFraction) = Unitful.NoUnits

function total_units(x::T,inverted::Bool) where T<:AbstractSpec
    if inverted
        return inv(total_units(x))
    else
        return total_units(x)
    end
end
function mol_units(x::T) where T <: AbstractIntensiveSpec
    return total_units(x) / u"mol"
end

function mol_units(x::T,inverted::Bool) where T<:AbstractIntensiveSpec
    if inverted
        return inv(mol_units(x))
    else
        return mol_units(x)
    end
end

function mass_units(x::T) where T <: AbstractIntensiveSpec
    return total_units(x) / u"kg"
end
function mass_units(x::T,inverted::Bool) where T<:AbstractSpec
    if inverted
        return inv(mass_units(x))
    else
        return mass_units(x)
    end
end

struct Spec{T <: AbstractSpec,U}
    type::T
    val::U
    is_mol::Bool
    is_total::Bool
    inverted::Bool
end

value(s::Spec) = s.val
specification(s::Spec) = s.type

function Base.show(io::IO, x::Spec{T}) where T
    specname = "Spec{" * string(T) * "}(" * string(value(x)) * ")"
    print(io, specname)
end

#preferred stripped
function _ups(x,normalize_units=false)
    if normalize_units
        return Unitful.ustrip(Unitful.upreferred(x))
    else
        return x
    end
end

function _ups(u,x::AbstractVector,normalize_units=false)
    if normalize_units
        return Unitful.ustrip.(Unitful.upreferred.(x))
    else
        return x
    end
end

function _udim(u::Unitful.Quantity)
    return dimension(u)
end

function _udim(u::AbstractVector{Unitful.Quantity})
    return dimension(eltype(u))
end

function check_spec_units(t::T, u::U,perm::Bool=true,normalize_units::Bool=false) where {T<:AbstractSpec,U<:Unitful.Quantity}    
    if perm == true
        if total_units(t)== _udim(u)
            is_total = true
            is_mol = false
            inverted = false
            _u =  _ups(u,normalize_units)
        elseif Unitful.dimension(mol_units(t)) == _udim(u)
            is_total = false
            is_mol = true
            inverted = false
            _u = _ups(u,normalize_units)
        elseif Unitful.dimension(mass_units(t)) == _udim(u) 
            is_total = false
            is_mol = false
            inverted = true
            _u = _ups(u,normalize_units)
        elseif Unitful.dimension(total_units(t,true)) == _udim(u)
            is_total = true
            is_mol = false
            inverted = true
            _u = _ups(u,normalize_units)
        elseif Unitful.dimension(mol_units(t,true)) == _udim(u)
            is_total = false
            is_mol = true
            inverted = true
            _u = _ups(u,normalize_units)
        elseif Unitful.dimension(mass_units(t,true)) == _udim(u)  
            is_total = false
            is_mol = false
            inverted = true
            _u = _ups(u,normalize_units)
        else
            throw(ArgumentError("the INPUT value is not a type of " * string(typeof(t))))
        end
        return _u,is_total, is_mol, inverted
    else
        if Unitful.dimension(total_units(t)) == _udim(u)
            is_total = true
            is_mol = false
            inverted = false
            _u = _ups(u,normalize_units)
            return _u, is_total, is_mol, inverted
        else
            throw(ArgumentError("the input value is not a type of " * string(typeof(t))))
        end
        
    end
end


_is_inverted(::AbstractSpec) = false
_is_inverted(::MolDensity) = true
_is_inverted(::MassDensity) = true



_is_molar(::MassVolume) = false
_is_molar(::MassDensity) = false

_is_total(::TotalVolume) = true


_is_molar(::AbstractIntensiveSpec) = true
_is_total(::AbstractIntensiveSpec) = false

_is_molar(::Union{Pressure,Temperature}) = false
_is_total(::Union{Pressure,Temperature}) = true

_is_total(::AbstractTotalSpec) = true

_is_total(::MolFractions) = false
_is_total(::MassFractions) = false

_is_molar(::MolNumbers) = true
_is_molar(::MassNumbers) = false
_is_molar(::MolFractions) = true
_is_molar(::MassFractions) = false
_is_molar(::Mass) = false
_is_molar(::Moles) = true


#special case for just real numbers, the tags are calculated on a type basis
function check_spec_units(t::T, u::U,perm::Bool=true,normalize_units::Bool=false) where {T<:AbstractSpec,U<:Real}    
    return u,_is_total(t),_is_molar(t),_is_inverted(t)
end

function Spec(t::AbstractIntensiveSpec, u,normalize_units::Bool=false)
    _u,is_total, is_mol, inverted = check_spec_units(t, u,true,normalize_units)
    return Spec(t, _u, is_total, is_mol, inverted)
end

function Spec(t::AbstractVolumeSpec, u,normalize_units::Bool=false)
    _u,is_total, is_mol, inverted = check_spec_units(t, u,true,normalize_units)
    return Spec(VolumeAmount(), _u, is_total, is_mol, inverted)
end

function Spec(t::Union{Pressure,Temperature}, u,normalize_units::Bool=false)
    _u,is_total, is_mol, inverted = check_spec_units(t, u,false,normalize_units)
    return Spec(t, _u, is_total, is_mol, inverted)
end

function Spec(t::AbstractTotalSpec, u::Unitful.Quantity,normalize_units::Bool=false)
    _u,is_total, is_mol, inverted = check_spec_units(t, u,false,normalize_units)
    return Spec(t, _u, is_total, is_mol, inverted)
end


function Spec(t::AbstractTotalSpec, u,normalize_units::Bool=false)
    return Spec(t, u, false, false, false) # total unit as default
end

function Spec(t::Union{MassNumbers,MolNumbers}, u,normalize_units::Bool=false)
    _u,is_total, is_mol, inverted = check_spec_units(t, u,false,normalize_units)
    return Spec(MaterialCompounds(), _u,is_total, is_mol, inverted)
end

function Spec(t::Union{Mass,Moles}, u,normalize_units::Bool=false)
    _u,is_total, is_mol, inverted = check_spec_units(t, u,false,normalize_units)
    return Spec(MaterialAmount(), _u,is_total, is_mol, inverted)
end


_is_frac(x::Real) = (zero(x) <= x <= one(x))
_is_frac(u::AbstractVector{Real}) = all(x -> x ≥ zero(x), u) && isapprox(sum(u), one(eltype(u)))

function Spec(t::Union{MolFractions,MassFractions} ,u,normalize_units::Bool=false)
    if _is_frac(u)
        return Spec(MaterialCompounds(), u, false, _is_molar(t), false) 
    else
        throw(ArgumentError("the the vector of values is not a valid fraction"))
    end
end

function Spec(t::AbstractFractionSpec, u::AbstractVector{Real})
    if _is_frac(u)
        return Spec(t, u, false, false, false) 
    else
        throw(ArgumentError("the the vector of values is not a valid fraction"))
    end
end

function Spec(t::AbstractFractionSpec, u::Real)
    if _is_frac(u)
        return Spec(t, u, false, false, false)
    else
        throw(ArgumentError("the value " * string(u) * " is not between 0 and 1."))
    end 
end

function Spec(t::T, u::Symbol,normalize_units::Bool=false) where T<:CategoricalSpec
    return Spec(t, u, false, false, false)
end

struct Specs{T}
    specs::T
    checked::Bool
end

values(s::Specs) = s.specs

function get_spec(val,specs)
    for spec in specs.specs
        if specification(spec) == val
            return spec
        end
    end
    return nothing
end

function has_spec(val,specs)::Bool
    for spec in specs.specs
        if specification(spec) == val
            return true
        end
    end
    return false
end

#returns the specification symbol in multicomponent
#returns :singlecomponent if there are no specifications
#throws error if more than one specification is found
function _specs_components(kwargs)
   
    xn = hasproperty(kwargs,:xn)
   xm = hasproperty(kwargs,:xm)
   n = hasproperty(kwargs,:n)
   m = hasproperty(kwargs,:m)
   mass = hasproperty(kwargs,:mass)
   moles = hasproperty(kwargs,:moles)
#if no specification
    amount_basis = (if ! (moles | mass) #one mol basis default
        :one_mol
    elseif moles & !(mass)
        :moles
    else mass & !(moles)
        :mass 
    end)

    if !(xn | xm | n | m )
        return :singlecomponent,amount_basis
    elseif (xn & !(xm | n | m))
        return :xn,amount_basis
    elseif (xm & !(xn | n | m))
        return :xm,amount_basis
    elseif (n & !(xm | xn | m))
        return :n,:n
    elseif (m & !(xm | xn | n))
        if ! (moles | mass) 
            return :m,:m
        end
    else
        throw(error("incorrect mass or molar specifications."))
    end
end
#gives the appropiate symbol to extract amount of matter


#this needs future work to specify
function _specs_phase_basis(kwargs)::Symbol
    pure_phase = hasproperty(kwargs,:pure_phase)
    vfrac = hasproperty(kwargs,:vfrac)
    phase_fracs = hasproperty(kwargs,:phase_fracs)
    if !(vfrac | phase_fracs | pure_phase) #assumes one phase default
        return :one_phase
    elseif pure_phase & !(vfrac | phase_fracs)
        return :pure_phase
    elseif vfrac & !(phase_fracs | pure_phase)
        return :two_phase
    elseif phase_fracs & !(vfrac | pure_phase)
        return :multi_phase
    else
        throw(error("there are more than one phase specification"))
    end
end

const KW_TO_SPEC = Dict{Symbol,AbstractSpec}(
:h =>  Enthalpy()
,:s =>  Entropy()
,:u =>  InternalEnergy()
,:p =>  Pressure()
,:t =>  Temperature()
,:v =>  MolVolume()
,:V =>  TotalVolume()
,:rho =>  MolDensity() #for flexibility
,:ρ =>  MolDensity() #also for flexibility
,:vm =>  MassVolume()
,:rhom =>  MassDensity() #for flexibility
,:ρm =>  MassDensity() #also for flexibility
,:mass =>  Mass()
,:moles =>  Moles()
,:xn =>  MolFractions()
,:xm =>  MassFractions()
,:n =>  MolNumbers()
,:m =>  MassNumbers()
,:mw =>  MolecularWeight()
,:vfrac =>  VaporFraction() #looking for better name
,:phase_fracs =>  PhaseFractions() #looking for better name
,:pure_phase =>PurePhase()
)


function _specs_C(kwargs,kw)::Int64
    if kw == :singlecomponent
        return 1
    else
        return length(getproperty(kwargs,kw))
    end
end

function _specs_P(kwargs,kw::Symbol)::Int64
    if (kw == :one_phase) |  (kw == :pure_phase)
        return 1
    elseif kw == :two_phase
        return 2
    else
        return length(kwargs.phase_fracs)
    end
end

function specs(;check=true,normalize_units=true,kwargs...)
    f0 = k -> Spec(KW_TO_SPEC[k],getproperty(kwargs.data,k),normalize_units)
    if check == true
        component_basis,amount_basis = _specs_components(kwargs.data)
        phase_basis = _specs_phase_basis(kwargs.data)
        C = _specs_C(kwargs.data,component_basis)
        P = _specs_P(kwargs.data,phase_basis)
        F = C - P + 2 #behold, the gibbs phase rule!
        F1 = length(kwargs.data)
        (component_basis != :singlecomponent) && (F1 = F1-1)

        if (amount_basis != :one_mol) & (amount_basis != :m) & (amount_basis != :n)
            F1 = F1-1
        end
        (phase_basis != :one_phase) && (F1 = F1-1)
        if F1 > F
            throw(error("the variables are overspecified."))
        elseif F1 > F
            throw(error("the variables are underspecified."))
        else
        
            tup = map(f0,keys(kwargs.data))
            return Specs(tup,true)
        end
    else
        tup = map(f0,keys(kwargs.data))
        return Specs(tup,false)
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