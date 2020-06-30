module PhysPropsRules


using Unitful
import Unitful: @u_str

include("utilities.jl")
include("properties.jl")
include("units.jl")
include("spec.jl")
export Spec
export Enthalpy, Entropy, InternalEnergy, 
Gibbs, Helmholtz, 
Volume, Pressure, Temperature, Mass, Moles,
MassNumbers, MolNumbers, MassFractions, MolFractions, 
PhaseFractions, VaporFraction,
MolecularWeight

end # module,