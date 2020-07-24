module PhysPropsRules


using Unitful
import Unitful: @u_str

include("utilities.jl")
include("properties.jl")
include("units.jl")
include("spec.jl")
#Spec Dispatch Types
export Enthalpy, Entropy, InternalEnergy, 
Gibbs, Helmholtz, 
VolumeAmount, Pressure, Temperature, MaterialAmount,
MaterialCompounds, PhaseFractions, VaporFraction,
MolecularWeight,
PurePhase
#Spec types
export Spec,Specs

#spec functions
export specs,specs_grid

#spec utilities
export value,specification,values
include("specs_model.jl")
#FromSpecs Model
export FromSpecs
include("fromspecs_props.jl")
export pressure,temperature,mass,moles


end # module