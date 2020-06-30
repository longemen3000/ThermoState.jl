abstract type AbstractSpec end
abstract type AbstractExtensiveSpec <: AbstractSpec  end
abstract type AbstractTotalSpec <: AbstractSpec  end
abstract type AbstractFractionSpec <: AbstractSpec  end

struct Enthalpy <: AbstractExtensiveSpec end
struct Entropy <: AbstractExtensiveSpec end
struct InternalEnergy <: AbstractExtensiveSpec end
struct Gibbs <: AbstractExtensiveSpec end
struct Helmholtz <: AbstractExtensiveSpec end
struct Volume <: AbstractExtensiveSpec end

struct Pressure <: AbstractTotalSpec end
struct Temperature <: AbstractTotalSpec end
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

# for checking if there is a amount of mass specified:
const TOTAL_MASS_SPECS = Union{Mass,Moles,MassNumbers,MolNumbers}


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
total_units(x::MassFractions) = u"kg"
total_units(x::MolFractions) = u"mol"

total_units(x::MassFractions) = Unitful.NoUnits
total_units(x::MolFractions) = Unitful.NoUnits
total_units(x::PhaseFractions) = Unitful.NoUnits
total_units(x::VaporFraction) = Unitful.NoUnits

mol_units(x::AbstractExtensiveSpec) = total_units(x) / u"mol"
mass_units(x::AbstractExtensiveSpec) = total_units(x) / u"kg"

struct Spec{T <: AbstractSpec,U}
    type::T
    val::U
    is_mol::Bool
    is_total::Bool
    inverted::Bool
end

function Base.show(io::IO, x::Spec{T}) where T
    specname = string(nameof(T)) * " specification: " * string(x.val)
    println(io, specname)
end

# extensive specs are chosen as molar specs, in the sense that they can be extensive or
# intensive units (u/kg) or (u/mol).


function check_spec_units(t::AbstractSpec, u::Unitful.Quantity)
    x = total_units(t) 
    is_total = true
    is_mol = true
    inverted = true
    if Unitful.dimension(u) == Unitful.dimension(x)
        is_total = true
        is_mol = false
        inverted = false
    elseif Unitful.dimension(u) == Unitful.dimension(x / u"mol")
        is_total = false
        is_mol = true
        inverted = false
    elseif Unitful.dimension(u) == Unitful.dimension(x / u"kg")   
        is_total = false
        is_mol = false
        inverted = true
    elseif Unitful.dimension(u) == Unitful.dimension(1 / x)
        is_total = true
        is_mol = false
        inverted = true
    elseif Unitful.dimension(u) == Unitful.dimension(u"mol" / x)
        is_total = false
        is_mol = true
        inverted = true
    elseif Unitful.dimension(u) == Unitful.dimension(u"kg" / x)   
        is_total = false
        is_mol = false
        inverted = true
    elseif Unitful.dimension(u) == Unitful.dimension(1 / x)
        is_total = true
        is_mol = true
        inverted = true
    else
        throw(ArgumentError("the input value is not a type of " * string(typeof(t))))
    end
    return is_total, is_mol, inverted
end

function Spec(t::AbstractSpec, u::Unitful.Quantity)
    is_total, is_mol, inverted = check_spec_units(t, u)
    return Spec(t, u, is_total, is_mol, inverted)
end

function Spec(t::AbstractSpec, u::AbstractVector{Unitful.Quantity})
    is_total, is_mol, inverted = check_spec_units(t, eltype(u))
    return Spec(t, u, is_total, is_mol, inverted)
end


# for total specs, a molar or weight value does not make sense.
function Spec(t::AbstractTotalSpec, u::Unitful.Quantity)
    x = total_units(t)
    if Unitful.dimension(u) == Unitful.dimension(x)
        return Spec(t, u, false, false, false) 
    else
        throw(ArgumentError("the input value is not a type of " * string(typeof(t))))
    end
end

function Spec(t::AbstractTotalSpec, u::AbstractVector{Unitful.Quantity})
    x = total_units(t)
    if Unitful.dimension(eltype(u)) == Unitful.dimension(x)
        return Spec(t, u, false, false, false) 
    else
        throw(ArgumentError("the input value is not a type of " * string(typeof(t))))
    end
end

function Spec(t::AbstractExtensiveSpec, u)
    return Spec(t, u, false, true, false) # molar unit as default
end

function Spec(t::AbstractTotalSpec, u)
    return Spec(t, u, false, false, false) # total unit as default
end

function Spec(t::AbstractFractionSpec, u::AbstractVector)
    if isapprox(sum(u), 1.0) & (any(x -> x < zero(x), 1.0))
        return Spec(t, u, false, true, false) # molar unit as default
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


# get_unit(mol_entropy,model,specs,unit = u"kJ"/kg)
function check_validity(specs)
    # check if there is any total spec
    any_total = any(x -> x.is_total, specs)
    any_mass =  any(x -> (<:)(typeof(x.type), TOTAL_MASS_SPECS), specs)

    
    # need to check if there is any mass value to extract:
    if any_total && any_mass
        return true
    else

    end


    function specs(;p=nothing,t=nothing,v=nothing,h=nothing,u=nothing,s=nothing,rho=nothing,phi=nothing,xn=nothing,xm=nothing,n=nothing,m=nothing)
        if !isnothing(v) && !isnothing(rho)
            throw(error("volume and density can't be specified at the same time"))
        end
        return nothing
    end
