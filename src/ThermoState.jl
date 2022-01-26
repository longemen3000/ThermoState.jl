module ThermoState


using Unitful
import Unitful: @u_str

include("utilities.jl")
include("docs.jl")
include("types.jl")
include("spec.jl")
include("state_type.jl")
include("units.jl")
include("show.jl")
include("specs_model.jl")
include("fromspecs_props.jl")

#Spec types
export Spec,ThermodynamicState,VariableSpec

#spec functions
export spec,state

#spec utilities
export specification,get_spec,has_spec,amount_type
#FromState Model
export FromState

#unitful utilities
export normalize_units,default_units,convert_unit

#property functions
export pressure,temperature,mass,moles,molar_mass
export mass_volume, mol_volume, total_volume
export mass_density, mol_density
export mass_enthalpy, mol_enthalpy,total_enthalpy
export mass_gibbs, mol_gibbs, total_gibbs
export mass_helmholtz, mol_helmholtz, total_helmholtz
export mass_internal_energy, mol_internal_energy, total_internal_energy
export mass_entropy, mol_entropy, total_entropy
export mass_fraction, mol_fraction
export mass_number, mol_number
export options, phase, quality, mass_vapor_fraction,mol_vapor_fraction
export molecular_weight

#macros
export @to_units, @spec_str

#state type for dispach
export state_type

#modules
export Types,QuickStates,StatePoints

end # module
