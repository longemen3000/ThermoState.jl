module PhysPropsRules


using Unitful
import Unitful: @u_str

include("utilities.jl")
include("docs.jl")
include("types.jl")
include("spec.jl")
include("units.jl")
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
export value,specification,values,get_spec
#FromSpecs Model
export FromSpecs
export pressure,temperature,mass,moles
export @to_units


end # module