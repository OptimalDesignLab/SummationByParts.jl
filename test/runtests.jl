using SummationByParts
using SummationByParts.OrthoPoly
using SummationByParts.Cubature
using SummationByParts.SymCubatures
using FactCheck

include("test_orthopoly.jl")
include("test_symcubatures.jl")
include("test_cubature.jl")
include("test_buildoperators.jl")

FactCheck.exitstatus()