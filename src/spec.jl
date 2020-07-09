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
struct Volume <: AbstractIntensiveSpec end

struct Pressure <: AbstractIntensiveSpec end
struct Temperature <: AbstractIntensiveSpec end

struct Mass <: AbstractTotalSpec end
struct Moles <: AbstractTotalSpec end

# those are vectors, 
struct MassNumbers <: AbstractTotalSpec end
struct MolNumbers <: AbstractTotalSpec end


struct MassFractions <: AbstractFractionSpec end
struct MolFractions <: AbstractFractionSpec end
struct PhaseFractions <: AbstractFractionSpec end
struct VaporFraction <: AbstractFractionSpec end
# for now,in a different category
struct MolecularWeight <: AbstractSpec end

# total units, function necessary to defining a new spec
total_units(x::Enthalpy) = u"J"
total_units(x::Entropy) = u"J/K"
total_units(x::InternalEnergy) = u"J"
total_units(x::Gibbs) = u"J"
total_units(x::Helmholtz) = u"J"
total_units(x::Volume) = u"m^3"

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

mol_units(x::AbstractIntensiveSpec) = total_units(x) / u"mol"
mass_units(x::AbstractIntensiveSpec) = total_units(x) / u"kg"


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

# extensive specs are chosen as molar specs, in the sense that they can be extensive or
# intensive units (u/kg) or (u/mol).
#ustrip and convert
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

function _udim(u::Unitful.Quantity)
    return dimension(u)
end

function _udim(u::AbstractVector{Unitful.Quantity})
    return dimension(eltype(u))
end

function check_spec_units(t::T, u::U,perm::Bool=true,normalize_units::Bool=false) where {T<:AbstractSpec,U<:Unitful.Quantity}    
    x = total_units(t) 
    if perm == true
        if _udim(u) == Unitful.dimension(x)
            is_total = true
            is_mol = false
            inverted = false
            _u =  _ucs(x,u,normalize_units)
        elseif _udim(u) == Unitful.dimension(x / u"mol")
            is_total = false
            is_mol = true
            inverted = false
            _u = _ucs(x / u"mol",u,normalize_units)
        elseif _udim(u) == Unitful.dimension(x / u"kg")   
            is_total = false
            is_mol = false
            inverted = true
            _u = _ucs(x / u"kg",u,normalize_units)
        elseif _udim(u) == Unitful.dimension(1 / x)
            is_total = true
            is_mol = false
            inverted = true
            _u = _ucs(1 / x,u,normalize_units)
        elseif _udim(u) == Unitful.dimension(u"mol" / x)
            is_total = false
            is_mol = true
            inverted = true
            _u = _ucs(u"mol" / x,u,normalize_units)
        elseif _udim(u) == Unitful.dimension(u"kg" / x)   
            is_total = false
            is_mol = false
            inverted = true
            _u = _ucs(u"kg" / x,u,normalize_units)
        else
            throw(ArgumentError("the input value is not a type of " * string(typeof(t))))
        end
        return _u,is_total, is_mol, inverted
    else
        if _udim(u) == Unitful.dimension(x)
            is_total = true
            is_mol = false
            inverted = false
            _u = _ucs(x,u,normalize_units)
            return _u, is_total, is_mol, inverted
        else
            throw(ArgumentError("the input value is not a type of " * string(typeof(t))))
        end
        
    end
end

function Spec(t::Type{T},u;normalize_units::Bool=false)  where T <: AbstractSpec
    return Spec(t(),u,normalize_units=normalize_units)
end

function Spec(t::AbstractIntensiveSpec, u::Unitful.Quantity;normalize_units::Bool=false)
    _u,is_total, is_mol, inverted = check_spec_units(t, u,true,normalize_units)
    return Spec(t, _u, is_total, is_mol, inverted)
end

function Spec(t::Union{Pressure,Temperature}, u::T;normalize_units::Bool=false) where T<:Number
    return Spec(t, u, false, true, false) 
end

function Spec(t::Union{Pressure,Temperature}, u::Unitful.Quantity;normalize_units::Bool=false)
    _u,is_total, is_mol, inverted = check_spec_units(t, u,false,normalize_units)
    return Spec(t, _u, is_total, is_mol, inverted)
end

function Spec(t::AbstractIntensiveSpec, u::T;normalize_units::Bool=false) where T<:Number
    return Spec(t, u, true, false, false) # molar unit as default
end

function Spec(t::AbstractTotalSpec, u::Unitful.Quantity;normalize_units::Bool=false)
    _u,is_total, is_mol, inverted = check_spec_units(t, u,false,normalize_units)
    return Spec(t, _u, is_total, is_mol, inverted)
end

function Spec(t::AbstractTotalSpec, u;normalize_units::Bool=false)
    return Spec(t, u, false, false, false) # total unit as default
end

function Spec(t::Union{MassNumbers,MolNumbers}, u::AbstractVector{Unitful.Quantity};normalize_units::Bool=false)
    _u,is_total, is_mol, inverted = check_spec_units(t, u,false,normalize_units)
    return Spec(t, _u, is_total, is_mol, inverted)
end

function Spec(t::Union{MassFractions,MolFractions}, u::AbstractVector{Unitful.Quantity};normalize_units::Bool=false)
    if isapprox(sum(u), 1.0) & (any(x -> x < zero(x), u))
        return Spec(t, u, false, true, false) 
    else
        throw(ArgumentError("the the vector of values is not a valid fraction"))
    end
end

function Spec(t::AbstractFractionSpec, u::AbstractVector)
    if isapprox(sum(u), 1.0) & (any(x -> x < zero(x), u))
        return Spec(t, u, false, true, false) 
    else
        throw(ArgumentError("the the vector of values is not a valid fraction"))
    end
end

function Spec(t::AbstractFractionSpec, u::Real)
    if 0 <= u <= 1
        return Spec(t, u, false, true, false)
    else
        throw(ArgumentError("the value " * string(u) * " is not between 0 and 1."))
    end 
end

struct Specs{T}
    specs::T
    checked::Bool
end

#returns the specification symbol in multicomponent
#returns :singlecomponent if there are no specifications
#throws error if more than one specification is found
function _specs_components(kwargs)::Symbol
   xn = hasproperty(kwargs,:xn)
   xm = hasproperty(kwargs,:xm)
   n = hasproperty(kwargs,:n)
   m = hasproperty(kwargs,:m)

#if no specification
    if (xn | xm | n | m) == false
        return :singlecomponent
    elseif (xn & !(xm | n | m))
        return :xn
    elseif (xm & !(xn | n | m))
        return :xm
    elseif (n & !(xm | xn | m))
        return :n
    elseif (m & !(xm | xn | n))
        return :m
    else
        throw(error("there are more than one mass or molar specifications."))
    end
end
#gives the appropiate symbol to extract amount of matter
function _specs_amount_basis(kwargs,component::Symbol)::Symbol
    mass = hasproperty(kwargs,:mass)
    moles = hasproperty(kwargs,:moles)
    
    if (component == :singlecomponent) | (component == :xn) | (component == :xm)
        if ! (moles | mass) #one mol basis default
            return :one_mol
        elseif moles & !(mass)
            return :moles
        elseif mass & !(moles)
            return :mass 
        else
            throw(error("there are more than one mass or molar specifications."))
        end
    elseif component == :n 
        if !(moles | mass) 
            return :n
        end
    elseif component == :m
        if ! (moles | mass) 
            return :m
        end
    else
        throw(error("invalid component specification"))
    end
end

#this needs future work to specify
function _specs_phase_basis(kwargs)::Symbol
    vfrac = hasproperty(kwargs,:vfrac)
    phase_fracs = hasproperty(kwargs,:phase_fracs)
    if !(vfrac | phase_fracs) #assumes one phase default
        return :one_phase
    elseif vfrac & !(phase_fracs)
        return :two_phase
    elseif phase_fracs & !(vfrac)
        return :multi_phase
    else
        throw(error("there are more than one phase specification"))
    end
end

const kw_to_spec = Dict{Symbol,AbstractSpec}(
:h =>  Enthalpy()
,:s =>  Entropy()
,:u =>  InternalEnergy()
,:p =>  Pressure()
,:t =>  Temperature()
,:v =>  Volume()
,:rho =>  Volume() #for flexibility
,:ρ =>  Volume() #also for flexibility
,:mass =>  Mass()
,:moles =>  Moles()
,:xn =>  MolFractions()
,:xm =>  MassFractions()
,:n =>  MolNumbers()
,:m =>  MassNumbers()
,:mw =>  MolecularWeight()
,:vfrac =>  VaporFraction() #looking for better name
,:phase_fracs =>  PhaseFractions() #looking for better name
)
function _specs_C(kwargs,kw)::Int64
    if kw == :singlecomponent
        return 1
    else
        return length(getproperty(kwargs,kw))
    end
end

function _specs_P(kwargs,kw)::Int64
    if kw == :one_phase
        return 1
    elseif kw == :two_phase
        return 2
    else
        return length(kwargs.phase_fracs)
    end
end

function specs(;check=true,normalize_units=true,kwargs...)
    f0 = k -> Spec(kw_to_spec[k],getproperty(kwargs.data,k);normalize_units=normalize_units)
    if check == true
        component_basis = _specs_components(kwargs.data)
        amount_basis = _specs_amount_basis(kwargs.data,component_basis)
        phase_basis = _specs_phase_basis(kwargs.data)
        C = _specs_C(kwargs.data,component_basis)
        P = _specs_P(kwargs.data,phase_basis)
        F = C - P + 2 #behold, the gibbs phase rule!
        F1 = length(kwargs.data)
        (component_basis != :singlecomponent) && (F1 = F1-1)
        !in(amount_basis,(:one_mol,:m,:n)) && (F1 = F1-1)  
        (component_basis != :one_phase) && (F1 = F1-1)
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



