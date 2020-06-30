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

#those are vectors, 
struct MassNumbers <: AbstractTotalSpec end
struct MolNumbers <: AbstractTotalSpec end


struct MassFractions <: AbstractFractionSpec end
struct MolFractions <: AbstractFractionSpec end
struct PhaseFractions <: AbstractFractionSpec end
struct VaporFraction <: AbstractFractionSpec end
#for now,in a different category
struct MolecularWeight <: AbstractSpec end

#for checking if there is a amount of mass specified:
const TOTAL_MASS_SPECS = Union{Mass,Moles,MassNumbers,MolNumbers}


#total units, function necessary to defining a new spec
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

mol_units(x::AbstractExtensiveSpec) = total_units(x)/u"mol"
mass_units(x::AbstractExtensiveSpec) = total_units(x)/u"kg"

struct Spec{T<:AbstractSpec,U}
    type::T
    val::U
    is_mol::Bool
    is_total::Bool
    inverted::Bool
end

function Base.show(io::IO, x::Spec{T}) where T
    specname = string(nameof(T)) * " specification: " * string(x.val)
    println(io , specname)
end

#extensive specs are chosen as molar specs, in the sense that they can be extensive or
#intensive units (u/kg) or (u/mol).


function check_spec_units(t::AbstractSpec,u::Unitful.Quantity)
    x = total_units(t) 
    is_total = true
    is_mol = true
    inverted = true
    if Unitful.dimension(u) == Unitful.dimension(x)
        is_total = true
        is_mol = false
        inverted = false
    elseif Unitful.dimension(u) == Unitful.dimension(x/u"mol")
        is_total = false
        is_mol = true
        inverted = false
    elseif Unitful.dimension(u) == Unitful.dimension(x/u"kg")   
        is_total = false
        is_mol = false
        inverted = true
    elseif Unitful.dimension(u) == Unitful.dimension(1/x)
        is_total = true
        is_mol = false
        inverted = true
    elseif Unitful.dimension(u) == Unitful.dimension(u"mol"/x)
        is_total = false
        is_mol = true
        inverted = true
    elseif Unitful.dimension(u) == Unitful.dimension(u"kg"/x)   
        is_total = false
        is_mol = false
        inverted = true
    elseif Unitful.dimension(u) == Unitful.dimension(1/x)
        is_total = true
        is_mol = true
        inverted = true
    else
        throw(ArgumentError("the input value is not a type of " * string(typeof(t))))
    end
    return is_total,is_mol,inverted
end

function Spec(t::AbstractSpec,u::Unitful.Quantity)
    is_total,is_mol,inverted = check_spec_units(t,u)
    return Spec(t,u,is_total,is_mol,inverted)
end

function Spec(t::AbstractSpec,u::AbstractVector{Unitful.Quantity})
    is_total,is_mol,inverted = check_spec_units(t,eltype(u))
    return Spec(t,u,is_total,is_mol,inverted)
end


#for total specs, a molar or weight value does not make sense.
function Spec(t::AbstractTotalSpec,u::Unitful.Quantity)
    x = total_units(t)
    if Unitful.dimension(u) == Unitful.dimension(x)
        return Spec(t,u,false,false,false) 
    else
        throw(ArgumentError("the input value is not a type of " * string(typeof(t))))
    end
end

function Spec(t::AbstractTotalSpec,u::AbstractVector{Unitful.Quantity})
    x = total_units(t)
    if Unitful.dimension(eltype(u)) == Unitful.dimension(x)
        return Spec(t,u,false,false,false) 
    else
        throw(ArgumentError("the input value is not a type of " * string(typeof(t))))
    end
end

function Spec(t::AbstractExtensiveSpec,u)
    return Spec(t,u,false,true,false) #molar unit as default
end

function Spec(t::AbstractTotalSpec,u)
    return Spec(t,u,false,false,false) #total unit as default
end

function Spec(t::AbstractFractionSpec,u::AbstractVector)
    if isapprox(sum(u),1.0) & (any(x->x<zero(x),1.0))
    return Spec(t,u,false,true,false) #molar unit as default
    else
        throw(ArgumentError("the the vector of values is not a valid fraction"))
    end
end

function Spec(t::AbstractFractionSpec,u::Real)
    if 0 <= u <= 1
        return Spec(t,u,false,true,false)
    else
        throw(ArgumentError("the value " * string(u) * " is not between 0 and 1."))
    end 
end
#get_unit(mol_entropy,model,specs,unit = u"kJ"/kg)
function check_validity(specs)
    #check if there is any total spec
    any_total = any(x->x.is_total,specs)
    any_mass =  any(x->(<:)(typeof(x.type),TOTAL_MASS_SPECS),specs)

    
    #need to check if there is any mass value to extract:
    if any_total && any_mass
        return true
    else

end
function get_unit(value::typeof(mol_entropy),model,specs;unit = nothing)
    if isnothing(unit)
    else

    end
end

fastlog2(x::Float64) = reinterpret(Int64,x)/4503599627370496 - 1023.0

function aprox_sqrt1(z::Float64)
    val_int = reinterpret(Int64,z)
    val_int -= 1 << 52
    val_int >>= 1
    val_int += 1 << 61
    return reinterpret(Float64,val_int)
end


@inline function approx_fisqrt1(x::Float64)
    i = reinterpret(Int64,x)
    ix = i + 9218868437227405312
    i >>= 1
    ii = 6902882268410593524 - i
    i =  6902882268410593524 - i
    y = reinterpret(Float64,i)
    yy = reinterpret(Float64,ii)
    y = yy*(4.7642669737958503 - x*y*y)
    mhalf = reinterpret(Float64,ix)
    t = muladd(mhalf,y*y,0.50000031699508796)
    y = muladd(y,t,y)
    t = muladd(mhalf, y*y, 0.50000000000007538)
    y = muladd(y,t,y)
    return y
end




@inline function invsqrt_47d(x::Float64)
    i = reinterpret(Int64,x)
    ix = i + 9218868437227405312
    i = i >> 1
    ii = 6902882269428073933 - i
    i =  6902882269428073933 - i
    y = reinterpret(Float64,i)
    yy = reinterpret(Float64,ii)
    y = yy*(4.7642670025852993 - x*y*y)
    mhalf = reinterpret(Float64,ix)
    t = muladd(mhalf,y*y,0.50000031697852854)
    y = muladd(y,t,y)
    c = x*y
    r =muladd(y,c,-1.0)
    c = muladd(0.375,r,-0.5)
    y = muladd(r*y,c,y)
    return y
end

function fisr(x::Float64)
    x2 = 0.5*x
    i = reinterpret(Int64,x)
    i = 6910469410427058089 - ( i>>1 )
    y = reinterpret(Float64,i)
    return y
    #return y * ( 1.5 - ( x2 * y * y ) )
end

function mysqrt(S::Float64)
    x2 = 0.5*S
    i = reinterpret(UInt64,S)
    i = 0x5FE6EB50C7B537A9 - ( i>>1 )
    yn = reinterpret(Float64,i)
    yn = yn*muladd(yn*yn,-x2,1.5)
    #yn = yn * ( 1.5 - ( x2 * yn * yn ) )
    xn = S*yn
    hn = yn*0.5
    rn = 0.0
    #iter 1
        rn = muladd(-xn,hn,0.5)
        xn = muladd(xn,rn,xn)
        hn = muladd(hn,rn,hn)
    #iter 2
    rn = muladd(-xn,hn,0.5)
    xn = muladd(xn,rn,xn)
    hn = muladd(hn,rn,hn)
    #iter 3
    rn = muladd(-xn,hn,0.5)
    xn = muladd(xn,rn,xn)
    hn = muladd(hn,rn,hn)
    #iter 4
    rn = muladd(-xn,hn,0.5)
    xn = muladd(xn,rn,xn)
    hn = muladd(hn,rn,hn)
    #iter 5
    rn = muladd(-xn,hn,0.5)
    xn = muladd(xn,rn,xn)
    return xn
end

function test_sqrt(f)
    N = 10^4
    v = filter(!isnan, abs.(reinterpret.(Float64, rand(UInt64, N))))
    my_v = f.(v)
    base_v = sqrt.(v)
    std_v = sqrt.(big.(v))
    my_error = abs.((my_v-std_v) ./ std_v)
    base_error = abs.((base_v-std_v) ./ std_v)
    @show maximum(my_error)
    println()
    @show maximum(base_error)
    println()
    @show sum(my_error)/N
    println()
    @show sum(base_error)/N
    println()
end