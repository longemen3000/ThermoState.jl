module PhysPropsRules


using Unitful
import Unitful: @u_str

include("utilities.jl")
include("docs.jl")
include("units.jl")
include("spec.jl")
include("show.jl")
include("specs_model.jl")
include("spec_compounds.jl")
include("fromspecs_props.jl")


#Spec Dispatch Types
#Spec types
export Spec,Specs

#spec functions
export spec,specs,specs_grid

#spec utilities
export value,specification,values
#FromSpecs Model
export FromSpecs
export pressure,temperature,mass,moles


end # module