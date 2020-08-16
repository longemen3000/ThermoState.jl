"""
    mol_helmholtz(model,specs::Specs,unit,[options...])::Real

Compute the molar helmholtz energy (J/mol), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function mol_helmholtz end

"""
    mass_helmholtz(model,specs::Specs,unit,[options...])::Real

Compute the specific (mass) helmholtz energy (J/kg), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function mass_helmholtz end

"""
    total_helmholtz(model,specs::Specs,unit,[options...])::Real

Compute the total helmholtz energy (J), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function total_helmholtz end


#gibbs energy

"""
    mol_gibbs(model,specs::Specs,unit,[options...])::Real

Compute the molar gibbs energy (J/mol), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function mol_gibbs end

"""
    mass_gibbs(model,specs::Specs,unit,[options...])::Real

Compute the specific (mass) gibbs energy (J/kg), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function mass_gibbs end

"""
    total_gibbs(model,specs::Specs,unit,[options...])::Real

Compute the total gibbs energy (J), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function total_gibbs end

#enthalpy

"""
    mol_enthalpy(model,specs::Specs,unit,[options...])::Real

Compute the molar enthalpy (J/mol), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function mol_enthalpy end

"""
    mass_enthalpy(model,specs::Specs,unit,[options...])::Real

Compute the specific (mass) enthalpy (J/kg), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function mass_enthalpy end

"""
    total_enthalpy(model,specs::Specs,unit,[options...])::Real

Compute the total enthalpy (J), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function total_enthalpy end

#internal energy


"""
    mol_internal_energy(model,specs::Specs,unit,[options...])::Real

Compute the molar internal energy (J/mol), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function mol_internal_energy end

"""
    mass_internal_energy(model,specs::Specs,unit,[options...])::Real

Compute the specific (mass) internal energy (J/kg), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function mass_internal_energy end

"""
    total_internal_energy(model,specs::Specs,unit,[options...])::Real

Compute the total internal energy (J), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function total_internal_energy end


#entropy


"""
    mol_entropy(model,specs::Specs,unit,[options...])::Real

Compute the molar entropy (J/(mol K)), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function mol_entropy end

"""
    mass_entropy(model,specs::Specs,unit,[options...])::Real

Compute the specific (mass) entropy (J/(kg K)), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function mass_entropy end

"""
    total_entropy(model,specs::Specs,unit,[options...])::Real

Compute the total entropy (J/K), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function total_entropy end

#temperature an pressure
"""
    temperature(model,specs::Specs,unit,[options...])::Real

Compute the temperature (K), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function temperature end

"""
    pressure(model,specs::Specs,unit,[options...])::Real

Compute the pressure (Pa), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function pressure end

#volume amounts

"""
    mol_volume(model,specs::Specs,unit,[options...])::Real

Compute the molar volume (m³/mol), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function mol_volume end

"""
    mass_volume(model,specs::Specs,unit,[options...])::Real

Compute the specific (mass) volume (m³/kg), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function mass_volume end

"""
    total_volume(model,specs::Specs,unit,[options...])::Real

Compute the total volume (m³), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function total_volume end

"""
    mol_density(model,specs::Specs,unit,[options...])::Real

Compute the molar density (mol/m³), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function mol_density end

"""
    mass_density(model,specs::Specs,unit,[options...])::Real

Compute the specific (mass) density (kg/m³), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function mass_density end
