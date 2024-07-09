# go to the docs directory
# ] pkg> activate .
# ] pkg> instantiate
# run include("make.jl")
# or from the root, run in the shell: julia --project=docs docs/make.jl

using Pkg 
push!(LOAD_PATH,"../../src/")
using Documenter
using SummationByParts

makedocs(;
    sitename = "SummationByParts.jl",
    modules = [SummationByParts],
    doctest = false,
    warnonly= true, 
    format = Documenter.HTML(
        size_threshold_ignore = ["reference.md"]
    ),
    pages = [
    "Home" => "index.md",
    "Reference" => "reference.md",
    ],
)
deploydocs(; repo = "https://github.com/OptimalDesignLab/SummationByParts.jl")
