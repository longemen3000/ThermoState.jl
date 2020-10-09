
"""
    mol_helmholtz(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:mol_a`, `:a`


Compute the molar helmholtz energy (J/mol), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mol_helmholtz end

"""
    mass_helmholtz(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:mass_a`


Compute the specific (mass) helmholtz energy (J/kg), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mass_helmholtz end

"""
    total_helmholtz(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:total_a`


Compute the total helmholtz energy (J), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function total_helmholtz end


#gibbs energy

"""
    mol_gibbs(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:mol_g`, `:g`


Compute the molar gibbs energy (J/mol), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mol_gibbs end

"""
    mass_gibbs(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:mass_g`


Compute the specific (mass) gibbs energy (J/kg), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mass_gibbs end

"""
    total_gibbs(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:total_g`


Compute the total gibbs energy (J), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function total_gibbs end

#enthalpy

"""
    mol_enthalpy(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:mol_h`, `:h`


Compute the molar enthalpy (J/mol), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mol_enthalpy end

"""
    mass_enthalpy(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:mass_h`


Compute the specific (mass) enthalpy (J/kg), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mass_enthalpy end

"""
    total_enthalpy(model,state::ThermodynamicState,unit,[options...])::Real


Keyword symbols: `:total_h`


Compute the total enthalpy (J), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function total_enthalpy end

#internal energy


"""
    mol_internal_energy(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:mol_u`, `:u`


Compute the molar internal energy (J/mol), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mol_internal_energy end

"""
    mass_internal_energy(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:mass_u`


Compute the specific (mass) internal energy (J/kg), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mass_internal_energy end

"""
    total_internal_energy(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:total_u`


Compute the total internal energy (J), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function total_internal_energy end


#entropy


"""
    mol_entropy(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:mol_s`, `:s`


Compute the molar entropy (J/(mol K)), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mol_entropy end

"""
    mass_entropy(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:mass_s`


Compute the specific (mass) entropy (J/(kg K)), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mass_entropy end

"""
    total_entropy(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:total_s`


Compute the total entropy (J/K), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function total_entropy end

#temperature an pressure
"""
    temperature(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:t`, `:T`


Compute the temperature (K), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function temperature end

"""
pressure(model,state::ThermodynamicState,unit,[options...])::Real

    Keyword symbols: `:p`, `:P`


Compute the pressure (Pa), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function pressure end

#volume amounts

"""
    mol_volume(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:mol_v`, `:v`


Compute the molar volume (m³/mol), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mol_volume end

"""
    mass_volume(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:mass_v`


Compute the specific (mass) volume (m³/kg), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mass_volume end

"""
    total_volume(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:total_v`, `:V`


Compute the total volume (m³), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function total_volume end

"""
    mol_density(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:mol_rho`, `:rho`, `:mol_ρ`, `:ρ`


Compute the molar density (mol/m³), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mol_density end

"""
    mass_density(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:mass_rho`, `:mass_ρ`


Compute the specific (mass) density (kg/m³), from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mass_density end


"""
    mol_fraction(model,state::ThermodynamicState,unit,[options...])::AbstractVector

Keyword symbols: `:xn`


Compute the molar fraction of each present compound present in the provided system, from or using the specifications as base.


The unit argument is not used here. but its present for consistency.
"""
function mol_fraction end

"""
    mass_fraction(model,state::ThermodynamicState,unit,[options...])::AbstractVector

Keyword symbols: `:xm`


Compute the mass fraction of each present compound present in the provided system, from or using the specifications as base.


The unit argument is not used here. but its present for consistency.
"""
function mass_fraction end

"""
    mol_number(model,state::ThermodynamicState,unit,[options...])::AbstractVector

Keyword symbols: `:n`


Compute the amount of moles (in mol) of each present compound present in the provided system, from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mol_number end

"""
    mass_number(model,state::ThermodynamicState,unit,[options...])::AbstractVector

Keyword symbols: `:m`


Compute the amount of mass (in kg) of each present compound present in the provided system, from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mass_number end

"""
    moles(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:moles`


Compute the total amount of moles (in mol) of the provided system, from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function moles end

"""
    mass(model,state::ThermodynamicState,unit,[options...])::Real

Keyword symbols: `:mass`


Compute the total amount of mass (in kg) of the provided system, from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function mass end

"""
    molar_mass(model,state::ThermodynamicState,unit,[options...])::Real


Compute the amount of mass per unit of moles (in kg/mol) of the provided state, from or using the specifications as base.


If you need unitful units, use the macro `@to_units` before the function call

"""
function molar_mass end


"""
    options(model,state::ThermodynamicState)::NamedTuple

Keyword symbols: `:options`


Returns any options related to the model, stored in the state or the model.

Returns an empty named tuple otherwise.

"""
function options end

"""
    phase(model,state::ThermodynamicState)::NamedTuple

Keyword symbols: `:phase`


Returns a symbol indicating the phase specified in the state or the model.

Returns `:unspecified` otherwise.

"""
function phase end

"""
    quality(model,state::ThermodynamicState)::NamedTuple

Keyword symbols: `:quality`


Returns the vapor quality stored in the state

Returns NaN otherwise.

This function is an alias for `mass_vapor_fraction`.

"""
function quality end

function mass_cp end
function mol_cp end
function mass_cv end
function mol_cv end
function sound_speed end

