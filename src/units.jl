
const LAVOISIER_DEFAULT_UNITS = (
    no_units = Unitful.NoUnits,
    mass = u"kg",
    time = u"s",
    mol = u"mol",
    length = u"m",
    temperature = u"K",
    area = u"m^2",
    volume = u"m^3",

    mol_volume = u"(m^3)/mol",
    mass_volume = u"(m^3)/kg",
    mol_density = u"mol/(m^3)",
    mass_density = u"kg/(m^3)",

    energy = u"J",
    potency = u"W",
    molecular_weight = u"g/mol",
    pressure = u"Pa",
    heat_capacity = u"J/K", #entalphy
    energy_density = u"J/m^3",
    mass_specific_energy = u"J/kg",
    mol_specific_energy = u"J/mol", #helmholtz, gibbs

    mass_specific_heat_capacity = u"J/(kg*K)",
    mol_specific_heat_capacity = u"J/(kg*mol)"
)

#ideal gas constant, from unitful
const R_IDEAL_GAS = Unitful.ustrip(Unitful.R) 


"""
    units(property,val,[preferred_unit=nothing])    

given a property and a value, units return the proper measurement unit for that property
the keyword `preferred_unit` allows a custom unit, converts from the default unit.
the default unit of a property can be accessed by default_units(property)


# Examples
```julia-repl
julia> units(pressure,100000,u"bar")
1 bar

julia> units(pressure,100000)
100000 Pa
```
"""
function units(property,val,preferred_unit=nothing)
    unit = default_units(property)
    if isnothing(preferred_unit)
        return unit*val
    else 
        if Unitful.dimension(unit) == Unitful.dimension(preferred_unit)
            return Unitful.uconvert(preferred_unit,unit*val)
        else
            return error(string(preferred_unit) * "cannot be converted from " * string(unit) * "." )
        end
    end
end


"""
    default_units(::Type{property})    

returns a default unit, uses by default the LAVOISIER_DEFAULT_UNITS named tuple for the units.
it is recomended to use the @units macro

# Examples
```julia-repl
default_units(pressure)
bar
```
"""
default_units(property) = Unitful.NoUnits
default_units(property::typeof(mol_helmholtz)) = u"J/mol"
default_units(property::typeof(residual_helmholtz)) = u"J/mol"
default_units(property::typeof(pressure)) = u"Pa"
default_units(property::typeof(compressibility_factor)) = Unitful.NoUnits
default_units(property::typeof(mol_entropy)) = u"J/(mol*K)"
default_units(property::typeof(mol_enthalpy)) = u"J/mol"
default_units(property::typeof(mol_internal_energy)) = u"J/mol"
default_units(property::typeof(mol_isochoric_heat_capacity)) = u"J/(mol*K)"
default_units(property::typeof(mol_isobaric_heat_capacity)) = u"J/(mol*K)"
 

default_units(property::typeof(mass_helmholtz)) = u"J/kg"
default_units(property::typeof(mass_entropy)) = u"J/(kg*K)"
default_units(property::typeof(mass_enthalpy)) = u"J/kg"
default_units(property::typeof(mass_internal_energy)) = u"J/kg"
default_units(property::typeof(mass_isochoric_heat_capacity)) = u"J/(kg*K)"
default_units(property::typeof(mass_isobaric_heat_capacity)) = u"J/(kg*K)"

default_units(property::typeof(sound_speed)) = u"m/s"
default_units(property::typeof(mol_fraction)) = Unitful.NoUnits
default_units(property::typeof(mol_number)) = u"mol"
default_units(property::typeof(mass_fraction)) = Unitful.NoUnits
default_units(property::typeof(mass_number)) = u"kg"
default_units(property::typeof(mol_density)) = u"mol/(m^3)"
default_units(property::typeof(mass_volume)) = u"(m^3)/kg"
default_units(property::typeof(mol_volume)) =u"(m^3)/mol"
default_units(property::typeof(mass_density)) =u"kg/(m^3)"

default_units(property::typeof(moles)) = u"mol"
default_units(property::typeof(mass)) = u"kg"
default_units(property::typeof(temperature)) = u"K"
default_units(property::typeof(molecular_weight)) = u"g/mol"
default_units(property::typeof(compounds_number)) = Unitful.NoUnits
default_units(property::typeof(covolumes)) =u"(m^3)/mol"

#the combinatorics of volume mass and moles
default_units(property::typeof(critical_mol_volume)) =u"(m^3)/mol"
default_units(property::typeof(critical_mass_volume)) = u"(m^3)/kg"
default_units(property::typeof(critical_mol_density)) = u"mol/(m^3)"
default_units(property::typeof(critical_mass_density)) =u"kg/(m^3)"
default_units(property::typeof(critical_compressibility_factor)) =u"kg/(m^3)"
default_units(property::typeof(critical_temperature)) = u"K"
default_units(property::typeof(critical_pressure)) = u"Pa"
default_units(property::typeof(acentric_factor)) = Unitful.NoUnits



"""
    @units f(x) unit

return `f(x)` with Unitful units, and converts to the provided unit.
the measurement unit is determined by `default_units(a::typeof(f))`

#example:
```julia-repl
julia> critical_pressure(water) bar
220.64 bar
```

"""
macro units(fx,unit)
    f = fx.args[1]
    unit_str = string(unit)
    quote
        units($f,$fx,Unitful.@u_str($unit_str))
    end
end


"""
    @units f(x)

return `f(x)` with Unitful units.
the measurement unit is determined by `default_units(a::typeof(f))`

#example:
```julia-repl
julia> critical_pressure(water)
2.2064e7 Pa
```

"""
macro units(fx)
    f = fx.args[1]
    #args = tuple($(fx.args[2:end]...))
    quote
        units($f,$fx)
    end
end




"""
    kg_per_gmol(model,x=nothing)

return how many kg of compound are in one gmol of a mix made of multiple compounds with molar fractions of x:
returns a numeric value with units: kg(mix)/gmol(mix). If x is not specified, the model is assumed to be of only one compound

"""
function kg_per_gmol(model,x=nothing)
    mw = molecular_weight(model)
    if !isnothing(x)
        res = zero(eltype(x))
        for i in eachindex(x)
            res += x[i]/mw[i]  #mol(compound)/mol(mix) * (mol(mix)/g(mix)) = mol(compound)/g(mix)
        end
            #sum(mol(compound_i)/g(mix)) = mol(mix)/g(mix) 
        g = one(eltype(x))/res #g(mix)/mol(mix)
        return 0.001*g # kg/mol
    else
        return 0.001*only(mw)
    end
end


"""
    gmol_per_kg(model,x=nothing)

return how many gmol of compound are in one kg of a mix made of multiple compounds with mass fractions of x:
returns a numeric value with units: gmol(mix)/kg(mix). If x is not specified, the model is assumed to be of only one compound

"""
function gmol_per_kg(model,x=nothing)
    mw = molecular_weight(model)
    if !isnothing(x) 
    res = zero(eltype(x))
    for i in eachindex(x)
        res += x[i]*mw[i]  #g(compound)/g(mix) * (g(mix)*mol(mix)) = g(compound)/mol(mix)
    end
        #sum( g(compound)/mol(mix) = g(mix)/mol(mix)
    g = one(eltype(x))/res #gmol/g
    return 1000.0*g # gmol/kg
    else
        return 1000.0*only(mw)
    end
end
