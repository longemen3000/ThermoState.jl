module ThermoState


using Unitful
import Unitful: @u_str

#include("utilities.jl")
include("docs.jl")
include("types.jl")
include("spec.jl")
include("units.jl")
include("show.jl")
include("specs_model.jl")
include("fromspecs_props.jl")
include("state_type.jl")



#Spec Dispatch Types
#Spec types
export Spec,ThermodynamicState,VariableSpec

#spec functions
export spec,state

#spec utilities
export value,specification,values,get_spec,normalize_units
#FromSpecs Model
export FromState

#property functions
export pressure,temperature,mass,moles
export mass_volume, mol_volume, total_volume
export mass_density, mol_density
export mass_enthalpy, mol_enthalpy,total_enthalpy
export mass_gibbs, mol_gibbs, total_gibbs
export mass_helmholtz, mol_helmholtz, total_helmholtz
export mass_internal_energy, mol_internal_energy, total_internal_energy
export mass_entropy, mol_entropy, total_entropy
export mass_fraction, mol_fraction
export mass_number, mol_number
export options, phase, quality

export @to_units
export state_type

export SinglePT,MultiPT

end # module
