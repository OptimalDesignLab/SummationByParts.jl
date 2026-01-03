using SummationByParts
using SummationByParts.OrthoPoly
using SummationByParts.Cubature
using SummationByParts.SymCubatures
using LinearAlgebra
using Random
using Test
import SummationByParts.Add

include("test_weakdifferentiate.jl")
include("test_weakdifferentiate_rev.jl")
include("test_weakdifferentiate_jac.jl")
include("test_differentiate.jl")
include("test_differentiate_rev.jl")
include("test_volumeintegrate.jl")
include("test_volumeintegrate_rev.jl")
include("test_facenormal.jl")
include("test_facenormal_rev.jl")
include("test_mappingjacobian.jl")
include("test_mappingjacobian_rev.jl")
include("test_faceinterpolate.jl")
include("test_faceinterpolate_rev.jl")
include("test_faceintegrate.jl")
include("test_faceintegrate_rev.jl")
include("test_faceintegrate_jac.jl")
include("test_utils.jl")
include("test_orthopoly.jl")
include("test_symcubatures.jl") 
include("test_cubature.jl") 
include("test_buildfaceoperators.jl")
include("test_buildoperators.jl") # need to check last test
include("test_directionaldifferentiate.jl")
include("test_edgestabilize.jl")

println("All tests are completed.")