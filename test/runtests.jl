using SummationByParts
#using Base.Test
using FactCheck

include("test_orthopoly.jl")
include("test_symcubatures.jl")
include("test_cubature.jl")
include("test_SummationByParts.jl")

FactCheck.exitstatus()