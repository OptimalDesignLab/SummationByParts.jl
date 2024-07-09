# go to the docs directory
# ] pkg> activate .
# ] pkg> instantiate
# run include("make.jl")

using Pkg 
# pkg"activate .."
push!(LOAD_PATH,"../../src/")
using Documenter
using SummationByParts

# makedocs(;
#     sitename = "SummationByParts.jl",
#     modules = [SummationByParts],
#     format = Documenter.HTML(),
#     clean=true,
#     warnonly = true, 
#     doctest = false)

    
# DocMeta.setdocmeta!(SummationByParts, :DocTestSetup, :(using SummationByParts); recursive = true)

makedocs(;
    sitename = "SummationByParts.jl",
    modules = [SummationByParts],
    doctest = false,
    warnonly= true, 
    format = Documenter.HTML(),
    pages = [
    "Home" => "index.md",
    "Reference" => "reference.md",
    # "Build Face Operators" => "buildfaceoperators.md"
    ],
)

# makedocs(;
#     sitename = "SummationByParts.jl",
#     modules = [SummationByParts],
#     doctest = false,
#     warnonly= true, 
#     format = Documenter.HTML(),
#     pages = [
#     "Home" => "index.md",
#     "References" => Any[
#     "buildfaceoperators" => "buildfaceoperators.md",
#     "buildoperators" => "buildoperators.md",
#     "Cubature" => "cubature.md",
#     "derivecubature" =>"derivecubature.md",
#     "differentiate_rev" => "differentiate_rev.md",
#     "Optimizer" => "optimizer.md",
#     "OrthoPoly" => "orthopoly.md", 
#     "SymCubatures" => "symcubatures.md"],

#     # "Build Face Operators" => "buildfaceoperators.md"
#     ],

# deploydocs(; repo = "https://github.com/OptimalDesignLab/SummationByParts.jl", push_preview = true)

#"references.md",