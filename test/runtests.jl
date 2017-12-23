using SummationByParts
using SummationByParts.OrthoPoly
using SummationByParts.Cubature
using SummationByParts.SymCubatures
using FactCheck
using ArrayViews
using ODLCommonTools
import ODLCommonTools.sview
import SummationByParts.Add

view = ArrayViews.view
#=
include("test_orthopoly.jl")
include("test_symcubatures.jl")
include("test_cubature.jl")
include("test_buildfaceoperators.jl")
include("test_buildoperators.jl")
include("test_weakdifferentiate.jl")
include("test_weakdifferentiate_rev.jl")
include("test_weakdifferentiate_jac.jl")
include("test_differentiate.jl")
include("test_differentiate_rev.jl")
include("test_directionaldifferentiate.jl")
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
include("test_edgestabilize.jl")
=#
include("test_utils.jl")

FactCheck.exitstatus()
