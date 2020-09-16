# ideal gas constant, from unitful
default_units(property::typeof(molar_mass)) = u"kg/mol"
default_units(property::typeof(mol_helmholtz)) = u"J/mol"
default_units(property::typeof(mass_helmholtz)) = u"J/kg"
default_units(property::typeof(total_helmholtz)) = u"J"

default_units(property::typeof(mol_gibbs)) = u"J/mol"
default_units(property::typeof(mass_gibbs)) = u"J/kg"
default_units(property::typeof(total_gibbs)) = u"J"

default_units(property::typeof(mol_internal_energy)) = u"J/mol"
default_units(property::typeof(mass_internal_energy)) = u"J/kg"
default_units(property::typeof(total_internal_energy)) = u"J"

default_units(property::typeof(mol_enthalpy)) = u"J/mol"
default_units(property::typeof(mass_enthalpy)) = u"J/kg"
default_units(property::typeof(total_enthalpy)) = u"J"

default_units(property::typeof(mol_entropy)) = u"J/(mol*K)"
default_units(property::typeof(mass_entropy)) = u"J/(kg*K)"
default_units(property::typeof(total_entropy)) = u"J/(K)"

default_units(property::typeof(mol_volume)) = u"(m^3)/mol"
default_units(property::typeof(mass_volume)) = u"(m^3)/kg"
default_units(property::typeof(total_volume)) = u"(m^3)"

default_units(property::typeof(mol_density)) = u"mol/(m^3)"
default_units(property::typeof(mass_density)) = u"kg/(m^3)"

default_units(property::typeof(moles)) = u"mol"
default_units(property::typeof(mass)) = u"kg"
default_units(property::typeof(temperature)) = u"K"
default_units(property::typeof(pressure)) = u"Pa"

default_units(property::typeof(mol_fraction)) = Unitful.NoUnits
default_units(property::typeof(mol_number)) = u"mol"
default_units(property::typeof(mass_fraction)) = Unitful.NoUnits
default_units(property::typeof(mass_number)) = u"kg"





function to_spec end

to_spec(property::typeof(mol_helmholtz)) = Helmholtz{MOLAR}()
to_spec(property::typeof(mass_helmholtz)) = Helmholtz{MASS}()
to_spec(property::typeof(total_helmholtz)) = Helmholtz{TOTAL}()

to_spec(property::typeof(mol_gibbs)) = Gibbs{MOLAR}()
to_spec(property::typeof(mass_gibbs)) =Gibbs{MASS}()
to_spec(property::typeof(total_gibbs)) = Gibbs{TOTAL}()

to_spec(property::typeof(mol_internal_energy)) =InternalEnergy{MOLAR}()
to_spec(property::typeof(mass_internal_energy)) = InternalEnergy{MASS}()
to_spec(property::typeof(total_internal_energy)) = InternalEnergy{TOTAL}()

to_spec(property::typeof(mol_enthalpy)) = Enthalpy{MOLAR}()
to_spec(property::typeof(mass_enthalpy)) = Enthalpy{MASS}()
to_spec(property::typeof(total_enthalpy)) = Enthalpy{TOTAL}()

to_spec(property::typeof(mol_entropy)) = Entropy{MOLAR}()
to_spec(property::typeof(mass_entropy)) = Entropy{MASS}()
to_spec(property::typeof(total_entropy)) = Entropy{TOTAL}()

to_spec(property::typeof(mol_volume)) = VolumeAmount{MOLAR,VOLUME}()
to_spec(property::typeof(mass_volume)) = VolumeAmount{MASS,VOLUME}()
to_spec(property::typeof(total_volume)) = VolumeAmount{TOTAL,VOLUME}()

to_spec(property::typeof(mol_density)) = VolumeAmount{MOLAR,VOLUME}()
to_spec(property::typeof(mass_density)) = VolumeAmount{MASS,DENSITY}()


to_spec(property::typeof(moles)) = MaterialAmount{MOLAR}()
to_spec(property::typeof(mass)) = MaterialAmount{MASS}()
to_spec(property::typeof(temperature)) = Temperature()
to_spec(property::typeof(pressure)) = Pressure()

to_spec(property::typeof(mol_fraction)) = MaterialCompounds{MOLAR,FRACTION}()
to_spec(property::typeof(mol_number)) = MaterialCompounds{MOLAR,TOTAL_AMOUNT}()
to_spec(property::typeof(mass_fraction)) = MaterialCompounds{MASS,FRACTION}()
to_spec(property::typeof(mass_number)) = MaterialCompounds{MASS,TOTAL_AMOUNT}()

#to_spec(property::typeof(sound_speed)) = u"m/s"
#to_spec(property::typeof(compressibility_factor)) = Unitful.NoUnits



#to_spec(property::typeof(molecular_weight)) = u"g/mol"
#to_spec(property::typeof(compounds_number)) = Unitful.NoUnits
#to_spec(property::typeof(covolumes)) = u"(m^3)/mol"

# the combinatorics of volume mass and moles
#to_spec(property::typeof(critical_mol_volume)) = u"(m^3)/mol"
#to_spec(property::typeof(critical_mass_volume)) = u"(m^3)/kg"
#to_spec(property::typeof(critical_mol_density)) = u"mol/(m^3)"
#to_spec(property::typeof(critical_mass_density)) = u"kg/(m^3)"
#to_spec(property::typeof(critical_compressibility_factor)) = u"kg/(m^3)"
#to_spec(property::typeof(critical_temperature)) = u"K"
#to_spec(property::typeof(critical_pressure)) = u"Pa"
#to_spec(property::typeof(acentric_factor)) = Unitful.NoUnits

==#





"""
    @to_units f(x)

return `f(x)` with Unitful units.
the measurement unit is determined by `default_units(a::typeof(f))`

#example:
```julia-repl
julia> critical_pressure(water)
2.2064e7 Pa
```

"""
macro to_units(fx)
    f = fx.args[1]
    if length(fx.args) <=3
        return quote
            @__MODULE__().default_units($f)*$fx

        end |> esc
    else
        return quote
            $fx*$(fx.args[4])
        end |> esc
    end
end


macro to_spec(fx)
    f = fx.args[1]
    model = fx.args[2]
    specs = fx.args[3]
    if length(fx.args) <=3
        return quote
            begin
                pp =  @__MODULE__()
                if typeof($model) === FromSpecs
                    pp.get_spec(pp.Types.Temperature(),$specs) #change here
                end
            end
        end |> esc
    else
        return quote
            begin
                pp =  @__MODULE__()
                if typeof($model) === FromSpecs
                    pp.get_spec(pp.Types.Temperature(),$specs) #change here
                else
                    val = pp.@to_units $fx
                    pp.spec(Temperature(),val)
                end
            end
        end |> esc
    end
end

#generic transfommation to get spec from defined functions
function spec(f::F,val,normalize_units=true) where F<:Function
    return spec(to_spec(f)val,normalize_units)
end

#directly calling from constructed specs:

(x::Helmholtz{MOLAR})(args...)=mol_helmholtz(args...)
(x::Helmholtz{MASS})(args...)=mass_helmholtz(args...)
(x::Helmholtz{TOTAL})(args...)=total_helmholtz(args...)

(x::Gibbs{MOLAR})(args...)=mol_helmholtz(args...)
(x::Gibbs{MASS})(args...)=mass_helmholtz(args...)
(x::Gibbs{TOTAL})(args...)=total_helmholtz(args...)

(x::InternalEnergy{MOLAR})(args...)=mol_internal_energy(args...)
(x::InternalEnergy{MASS})(args...)=mass_internal_energy(args...)
(x::InternalEnergy{TOTAL})(args...)=total_internal_energy(args...)

(x::Enthalpy{MOLAR})(args...)=mol_enthalpy(args...)
(x::Enthalpy{MASS})(args...)=mass_enthalpy(args...)
(x::Enthalpy{TOTAL})(args...)=total_enthalpy(args...)

(x::Entropy{MOLAR})(args...)=mol_entropy(args...)
(x::Entropy{MASS})(args...)=mol_entropy(args...)
(x::Entropy{TOTAL})(args...)=mol_entropy(args...)


(x::VolumeAmount{MOLAR,VOLUME})(args...)=mol_volume(args...)
(x::VolumeAmount{MASS,VOLUME})(args...)=mass_volume(args...)
(x::VolumeAmount{TOTAL,VOLUME})(args...)=total_volume(args...)

(x::VolumeAmount{MOLAR,DENSITY})(args...)=mol_density(args...)
(x::VolumeAmount{MASS,DENSITY})(args...)=mass_density(args...)

(x::Temperature)(args...)=temperature(args...)
(x::Pressure)(args...)=pressure(args...)

(x::MaterialAmount{MOLAR})(args...)=moles(args...)
(x::MaterialAmount{MASS})(args...)=mass(args...)


(x::MaterialCompounds{MOLAR,FRACTION})(args...)=mol_fraction(args...)
(x::MaterialCompounds{MOLAR,TOTAL_AMOUNT})(args...)=mol_number(args...)
(x::MaterialCompounds{MASS,FRACTION})(args...)=mass_fraction(args...)
(x::MaterialCompounds{MASS,TOTAL_AMOUNT})(args...)=mass_number(args...)
(x::SingleComponent)(args...) = nothing


_units_tuple(st) = (default_units(typeof(st[1])),
default_units(typeof(st[2])),
default_units(typeof(st[3])))

function properties(st::T,
        	state::ThermodynamicState,
            mw = nothing,
            units=_units_tuple(st)) where T <: Tuple{AbstractSpec,AbstractSpec,AbstractSpec}
        
    f(val,unit) = val(FromState(),state,unit,mw)
    
    return (f(st[1],units[1]),f(st[2],units[2]),f(st[3],units[3]))
end
