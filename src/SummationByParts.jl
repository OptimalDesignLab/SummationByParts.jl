module SummationByParts

using Optim
using LineSearches
#using ArrayViews
#using ODLCommonTools
#import ODLCommonTools.sview

# """
#     getComplexStep(T)

# returns an appropriate complex-step size for the given type
# """
# getComplexStep{T <: Float32}(::Type{T}) = 1f-20
# getComplexStep{T <: Float64}(::Type{T}) = 1e-60
# getComplexStep{T <: Complex64}(::Type{T}) = 1f-20
# getComplexStep{T <: Complex128}(::Type{T}) = 1e-60

# getComplexstep(::Type{T}) where {T<:Float32} = 1f-20
# getComplexstep(::Type{T}) where {T<:Float64} = 1e-60
# getComplexstep(::Type{T}) where {T<:Complex{Float64}} = 1e-60

include("compile_time.jl")
include("face_types.jl")
include("orthopoly.jl")
include("symcubatures.jl")
include("optimizer.jl")
include("cubature.jl")
include("../quadrature_data/write_quad_file.jl")
include("../visual/plot_nodes.jl")
include("../visual/plot_nodes_web.jl")

using .OrthoPoly
using .SymCubatures
using .Cubature
using .Optimizer

export AbstractSBP
export LineSegSBP
export TriSBP, SparseTriSBP
export TetSBP, SparseTetSBP
export getLineSegSBPLobbato, getLineSegSBPLegendre
export getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE, getTriCubatureForTetFaceDiagE
export getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE
export AbstractFace
export DenseFace, LineSegFace, TriFace, TetFace
export SparseFace, TriSparseFace, TetSparseFace
export getLineSegFace, getTriFaceForDiagE, getTetFaceForDiagE
export Boundary, Interface

export plotly_tri_nodes, plotly_tet_nodes, write_tri_data, write_tet_data, plot_tri_nodes, plot_tet_nodes
export pso, levenberg_marquardt, rosenbrock, rastrigin
export calcnodes, calcminnodedistance, getNumFaceNodes, truncErr, quadTruncErr, checkInteriorNodeLocaton
export weakdifferentiate!, weakDifferentiateElement!
export weakdifferentiate_rev!, weakDifferentiateElement_rev!
export weakDifferentiateElement_jac!
export differentiate!, differentiateElement!
export differentiate_rev!, differentiateElement_rev!
export directionalDifferentiateElement!
export volumeintegrate!, volumeIntegrateElement!
export volumeintegrate_rev!, volumeIntegrateElement_rev!
export calcMappingJacobian!, calcMappingJacobianElement!, mappingjacobian!
export calcMappingJacobian_rev!, mappingjacobian_rev!
export facenormal!, calcFaceNormals!
export facenormal_rev!, calcFaceNormals_rev!
export boundaryinterpolate!, boundaryFaceInterpolate!
export boundaryinterpolate_rev!, boundaryFaceInterpolate_rev!
export interiorfaceinterpolate!, interiorFaceInterpolate!
export interiorfaceinterpolate_rev!, interiorFaceInterpolate_rev!
export integratefunctional!, integrateBoundaryFunctional!
export integratefunctional_rev!, integrateBoundaryFunctional_rev!
export boundaryintegrate!, boundaryFaceIntegrate!
export boundaryintegrate_rev!, boundaryFaceIntegrate_rev!
export boundaryFaceIntegrate_jac!
export interiorfaceintegrate!, interiorFaceIntegrate!
export interiorfaceintegrate_rev!, interiorFaceIntegrate_rev!
export interiorFaceIntegrate_jac!
export edgestabilize!, permuteinterface!

export minTriCubatureDiagE, minTriCubatureOmega

"""
### SBP.UnaryFunctor

`UnaryFunctor` is abstract type that is used to specify different types of
residual updates in the useoperators.jl and usefaceoperators.jl files

"""
abstract type UnaryFunctor end

"""
### SBP.Add

`Add` is basically an identity operator.  It is needed as a default
UnaryFunctor.

"""
mutable struct Add <: UnaryFunctor
end
function (update::Add)(scalar)
  return scalar
end

"""
### SBP.Subtract

`Subtract` is used to negate values, which is useful in residual updates

"""
mutable struct Subtract <: UnaryFunctor
end
function (update::Subtract)(scalar)
  return -scalar
end

include("sbp_types.jl") #<--- holds all the abstract and concrete types
include("outerconstructors.jl") #<--- outer constructors and factories
include("buildfaceoperators.jl") #<--- functions related to building face ops
include("buildoperators.jl") #<--- functions related to building SBP operators
include("weakdifferentiate.jl") #<--- functions for weak differentiation
include("weakdifferentiate_rev.jl") #<--- reverse-diff of weak differentiation
include("weakdifferentiate_jac.jl") #<--- Jacobians for weak differentiation
include("differentiate.jl") #<--- functions for strong differentiation
include("differentiate_rev.jl") #<--- reverse-diff of strong differentiation
include("directionaldifferentiate.jl") #<--- directional differentiation
include("volumeintegrate.jl") #<--- volume integration against test functions
include("volumeintegrate_rev.jl") #<--- reverse-diff of volume integration
include("facenormal.jl") #<--- functions to compute scaled face normals
include("facenormal_rev.jl") #<--- reverse-diff of scaled face normals
include("mappingjacobian.jl") #<--- functions to compute the mapping jacobian
include("mappingjacobian_rev.jl") #<--- reverse-diff of mappingjacobian
include("faceinterpolate.jl") #<--- functions to interpolate data to faces
include("faceinterpolate_rev.jl") #<--- functions to interpolate data to faces
include("faceintegrate.jl") #<--- functions for integration over faces
include("faceintegrate_rev.jl") #<--- reverse mode of faceintegrate.jl
include("faceintegrate_jac.jl") #<--- Jacobians for combined face operations
include("edgestabilize.jl") #<--- functions related to edge stabilization
include("utils.jl") # <--- miscellaneous functions
include("derivecubature.jl") # <--- functions to derive cubature rules

end # module
