#functions for calculated properties, they require a model and a property (for example Cp(T))
"""
    mol_helmholtz(model,specs::Specs,unit,[options...])::Real

Compute the molar helmholtz energy (J/mol).
If you need unitful units, use the macro @to_units before the function call

"""
function mol_helmholtz end

"""
    residual_helmholtz(model,args...)::Real
    residual_helmholtz(model,specs::Specs,unit::Unitful.Units)::Real

Compute the molar helmholtz energy evaluated in the arguments.
It shall return a numeric value with units: J/mol.
used in equilibria calculation
"""
function residual_helmholtz end   

"""
    mol_gibbs(model,specs::Specs,unit,[options...])::Real

Compute the molar gibbs energy (J/mol), from or using the specifications as base.
If you need unitful units, use the macro @to_units before the function call

"""
function mol_gibbs end


"""
    pressure(model,specs::Specs,[options])::Real

Compute the pressure of a model evaluated in the arguments.
If the model  uses the pressure as an independent property, it return the pressure stored in the phase.
It shall return a numeric value with units: Pa.

"""

function pressure  end


"""
    compressibility_factor(model,specs::Specs,options)::Real

Compute the  compressibility factor of a model evaluated in the arguments.
It shall return a numeric value with no units.

"""
function compressibility_factor  end

"""
    mol_entropy(model,specs::Specs,options)::Real

Compute the  molar entropy of a model evaluated in the arguments.
It shall return a numeric value with units: J/(K mol).

"""
function mol_entropy end  

"""
    mol_enthalpy(model,args...)::Real
    mol_enthalpy(model,specs::Specs,options)::Real

Compute the molar enthalpy of a model evaluated in the arguments.
It shall return a numeric value with units: J/mol.

""" 
function mol_enthalpy end   

"""

    mol_internal_energy(model,specs::Specs,options)::Real

Compute the molar internal energy of a model evaluated in the arguments.
It shall return a numeric value with units: J/mol
.
""" 
function mol_internal_energy end  

"""
    mol_isochoric_heat_capacity(model,specs::Specs,options)::Real

Compute the molar specific heat capacity at constant volume (Cv) of a model evaluated in the arguments.
It shall return a numeric value with units: J/(K mol).
it can also be called by `mol_cv`, but it must be defined in terms of `mol_isochoric_heat_capacity`

""" 
function mol_isochoric_heat_capacity end 
const mol_cv = mol_isochoric_heat_capacity

"""
    mol_isobaric_heat_capacity(model,specs::Specs,options)::Real

Compute the molar specific heat capacity at constant pressure (Cp) of a model evaluated in the arguments.
It shall return a numeric value with units: J/(K mol).
it can also be called by `mol_cp`, but it must be defined in terms of `mol_isobaric_heat_capacity`

""" 
function mol_isobaric_heat_capacity end   
const mol_cp = mol_isobaric_heat_capacity

"""
    sound_speed(model,specs::Specs,options)::Real

Compute the sound of speed of a model evaluated in the arguments.
It shall return a numeric value with units: m/s.

""" 
function sound_speed end  

"""
    mol_fraction(spec)::AbstractVector{<:Real}

Return the molar fraction stored in the phase.
It shall return a vector of numeric values corresponding to the molar fractions.

""" 
function mol_fraction end   

"""
    mass_helmholtz(model,specs::Specs,options)::Real

Compute the specific mass helmholtz energy. (ideal helmholtz + residual helmholtz), evaluated in the arguments.
It shall return a numeric value with units: J/kg.

"""
function mass_helmholtz end

"""
    mass_gibbs(model,specs::Specs,options)::Real

Compute the mass gibbs energy, evaluated in the arguments.
It shall return a numeric value with units: J/kg.

"""
function mass_gibbs end

"""
    mass_entropy(model,specs::Specs,options)::Real

Compute the mass entropy, evaluated in the arguments.
It shall return a numeric value with units: J/(K*kg).

"""
function mass_entropy end 

"""
    mass_enthalpy(model,specs::Specs,options)::Real

Compute the mass gibbs energy, evaluated in the arguments.
It shall return a numeric value with units: J/kg.

"""
function mass_enthalpy end

"""
    mass_internal_energy(model,specs::Specs,options)::Real

Compute the mass internal energy, evaluated in the arguments.
It shall return a numeric value with units: J/kg.

"""
function mass_internal_energy end


"""
    mass_isochoric_heat_capacity(model,specs::Specs,options)::Real

Compute the mass specific heat capacity at constant volume (Cv) of a model evaluated in the arguments.
It shall return a numeric value with units: J/(K kg).
it can also be called by `mass_cv`, calculated in terms of `mol_isochoric_heat_capacity`

"""
function mass_isochoric_heat_capacity end
const mass_cv = mass_isochoric_heat_capacity

"""
    mass_isobaric_heat_capacity(model,specs::Specs,options)::Real

Compute the mass specific heat capacity at constant pressure (Cp) of a model evaluated in the arguments.
It shall return a numeric value with units: J/(K kg).
it can also be called by `mass_cp`, calculated in terms of `mol_isobaric_heat_capacity`

"""
function mass_isobaric_heat_capacity end
const mass_cp = mass_isobaric_heat_capacity



"""
    mol_number(Spec)::AbstractVector{<:Real}

Return the moles of each compound stored in the phase.
It shall return a vector of numeric values corresponding to the amount of moles of each compound.

""" 
function mol_number end   

"""
    mass_fraction(Spec)::AbstractVector{<:Real}

Return the mass fraction stored in the phase.
It shall return a vector of numeric values corresponding to the mass fractions.

""" 

function mass_fraction end

"""
    mass_number(Spec)::AbstractVector{<:Real}

Return the mass of each compound stored in the phase.
It shall return a vector of numeric values corresponding to the mass of each compound.

""" 
function mass_number end 

"""
    mol_density(model,args...)::Real
    mol_density(model,specs::Specs,options)::Real

Return the molar density in the phase.
If the model  uses a volume as an independent property, it return the molar density stored in the phase.
It shall return a numeric value with units: mol/m3

""" 
function mol_density end  

"""
    mass_volume(model,args...)::Real
    mass_volume(model,specs::Specs,options)::Real

Compute the mass volume of a model evaluated in the arguments.
If the model uses a volume as an independent property, it return the mass volume stored in the phase.
It shall return a numeric value with units: m3/kg

"""  
function mass_volume end   

"""
    mass_density(model,args...)::Real
    mass_density(model,specs::Specs,options)::Real

Compute the mass density of a model evaluated in the arguments.
If the model uses a volume as an independent property, it return the mass density stored in the phase.
It shall return a numeric value with units: kg/m3

"""  
function mass_density end   

"""
    mol_volume(model,args...)::Real
    mol_volume(model,specs::Specs,options)::Real

Return the molar density in the phase.
If the model  uses a volume as an independent property, it return the molar density stored in the phase.
It shall return a numeric value with units: m3/mol

"""
function mol_volume end

"""
    moles(Spec)::Real

Return the total amount of moles present in the phase.
It shall return a numeric value with units: mol

""" 
function moles end   

"""
    mass(Spec)::Real

Return the total mass present in the phase.
It shall return a numeric value with units: kg

"""
function mass end   

"""
    temperature(model,args...)::Real
    temperature(model,specs::Specs,options)::Real

Compute the temperature of a model evaluated in the arguments.
If the model uses a temperature as an independent property, it return the temperature stored in the phase.
It shall return a numeric value with units: K

"""  
function temperature end


"""
    molecular_weight(model)::AbstractVector{<:Real}

Return the total molecular weight of the compounds present in the model.
It shall return a vector of numeric values with units: g/mol

"""
function molecular_weight end

"""
    compounds_number(model)::Int64

Return the number of the compounds present in the model.


"""
function compounds_number end

"""
    covolumes(model,args...)::AbstractVector{<:Real}
    covolumes(model,phase)::AbstractVector{<:Real}

Returns the minimum molar volume that the model can handle.
It shall return a vector numeric value with units: m3/mol

"""
function covolumes end

#the combinatorics of volume mass and moles

function critical_mol_volume end
function critical_mass_volume end
function critical_mol_density end
function critical_mass_density end

function critical_compressibility_factor end
function critical_temperature end
function critical_pressure end
function acentric_factor end
