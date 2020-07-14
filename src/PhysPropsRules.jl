module PhysPropsRules


using Unitful
import Unitful: @u_str

include("utilities.jl")
include("properties.jl")
include("units.jl")
include("spec.jl")
include("specs_model.jl")
include("specs_matter.jl")


export Spec,Specs,specs,specs_grid
value,specification,values,FromSpecs
export Enthalpy, Entropy, InternalEnergy, 
Gibbs, Helmholtz, 
VolumeAmount, Pressure, Temperature, MaterialAmount,
MaterialCompounds, PhaseFractions, VaporFraction,
MolecularWeight
PurePhase
export pressure,temperature,mass,moles


end # module