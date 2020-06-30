module PhysPropsRules


using Unitful
import Unitful: @u_str

include("utilities.jl")
include("properties.jl")
include("units.jl")
include("spec.jl")
export SpecReal

end # module
