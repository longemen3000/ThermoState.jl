using ThermoState
using Documenter

makedocs(;
    modules=[ThermoState],
    authors="longemen3000 <longemen3000@gmail.com> and contributors",
    repo="https://github.com/longemen3000/ThermoState.jl/blob/{commit}{path}#L{line}",
    sitename="ThermoState.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://longemen3000.github.io/ThermoState.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "The Specification Object (`Spec`)" => "spec.md",
        "Thermodynamic State Object (`ThermodynamicState`)" => "state.md",
        "Utilities" => "utils.md",
        "Examples" => [
            "new_model.md",
        ],
        "API"=> "api.md"
    ],
)

deploydocs(;
    repo="github.com/longemen3000/ThermoState.jl",
)
