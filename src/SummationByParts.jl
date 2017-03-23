__precompile__(false)
module SummationByParts

using ArrayViews
using ODLCommonTools
import ODLCommonTools.sview

include("orthopoly.jl")
include("symcubatures.jl")
include("cubature.jl")

using .OrthoPoly
using .SymCubatures
using .Cubature

export AbstractSBP, TriSBP, TetSBP, SparseTriSBP, SparseTetSBP
export getTriSBPGamma, getTriSBPOmega, getTriSBPWithDiagE
export getTetSBPGamma, getTetSBPOmega
export AbstractFace, TriFace, TetFace
export getTriFaceForDiagE
export Boundary, Interface

export calcnodes, calcminnodedistance, getNumFaceNodes
export weakdifferentiate!, weakDifferentiateElement!
export weakdifferentiate_rev!, weakDifferentiateElement_rev!
export differentiate!, differentiateElement!
export differentiate_rev!, differentiateElement_rev!
export directionalDifferentiateElement!
export volumeintegrate!, volumeIntegrateElement!
export volumeintegrate_rev!, volumeIntegrateElement_rev!
export calcMappingJacobian!, calcMappingJacobianElement!, mappingjacobian!
export calcMappingJacobian_rev!
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
export interiorfaceintegrate!, interiorFaceIntegrate!
export interiorfaceintegrate_rev!, interiorFaceIntegrate_rev!
export edgestabilize!, permuteinterface!

@doc """
### SBP.AbstractSBP

`AbstractSBP` is a parametric abstract type that defines summation-by-parts
finite-difference operators.

"""->
abstract AbstractSBP{T<:Number}
#abstract AbstractSBP{T<:AbstractFloat}

@doc """
### SBP.TriSBP

Defines diagonal-norm SBP first-derivative operators on a right-triangle.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the triangle required for these operators
* `cub` : a symmetric cubature type for triangles
* `vtx` : vertices of the reference element in computational space
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""->
immutable TriSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  cub::TriSymCub{T}
  vtx::Array{T,2}
  w::Array{T,1}
  Q::Array{T,3}

  # inner constructor
  function TriSBP(degree::Int, cub::TriSymCub{T}, vtx::Array{T,2},
                  w::Array{T,1}, Q::Array{T,3})
    @assert( degree >= 1 && degree <= 4)
    numnodes = cub.numnodes
    @assert( size(Q,1) == size(Q,2) == size(w,1) == numnodes )
    @assert( size(Q,3) == 2 )
    new(degree, numnodes, cub, vtx, w, Q)
  end
end

@doc """
### SBP.SparseTriSBP

Defines diagonal-norm SBP first-derivative operators on a right-triangle using a
cubature rule that is greater than 2*p-1.  This provides additional flexiblity
in the SBP operator that is used to make a sparse S.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the triangle required for these operators
* `cub` : a symmetric cubature type for triangles
* `vtx` : vertices of the reference element in computational space
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""->
immutable SparseTriSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  cub::TriSymCub{T}
  vtx::Array{T,2}
  w::Array{T,1}
  Q::Array{T,3}

  function SparseTriSBP(;degree::Int=1, faceorder::Array{Int,1}=[1;2;3], 
                        internal=false, vertices=true,
                        cubdegree::Int=2*degree+1)
    @assert( degree >= 1 && degree <= 4 )
    cub, vtx = tricubature(cubdegree, T, internal=internal, vertices=vertices)
    numnodes = cub.numnodes
    Q = zeros(T, (numnodes, numnodes, 2))
    w, Q = SummationByParts.buildsparseoperators(cub, vtx, degree)
    new(degree, numnodes, cub, vtx, w, Q)
  end
end

@doc """
### SBP.TetSBP

Defines diagonal-norm SBP first-derivative operators on a right-tetrahedron.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the tetrahedron required for these operators
* `cub` : a symmetric cubature type for tetrahedra
* `vtx` : vertices of the reference element in computational space
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""->
immutable TetSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  cub::TetSymCub{T}
  vtx::Array{T,2}
  w::Array{T,1}
  Q::Array{T,3}
  
  # inner constructor
  function TetSBP(degree::Int, cub::TetSymCub{T}, vtx::Array{T,2},
                  w::Array{T,1}, Q::Array{T,3})
    @assert( degree >= 1 && degree <= 4)
    numnodes = cub.numnodes
    @assert( size(Q,1) == size(Q,2) == size(w,1) == numnodes )
    @assert( size(Q,3) == 3 )
    new(degree, numnodes, cub, vtx, w, Q)
  end
end

@doc """
### SBP.SparseTetSBP

Defines diagonal-norm SBP first-derivative operators on a right-tetrahedron
using a cubature rule that is greater than 2*p-1.  This provides additional
flexiblity in the SBP operator that is used to make a sparse S.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the tetrahedron required for these operators
* `cub` : a symmetric cubature type for tetrahedra
* `vtx` : vertices of the reference element in computational space
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""->
immutable SparseTetSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  cub::TetSymCub{T}
  vtx::Array{T,2}
  w::Array{T,1}
  Q::Array{T,3}

  function SparseTetSBP(;degree::Int=1, faceorder::Array{Int,1}=[1;2;3;4],
                        internal=false, cubdegree::Int=2*degree-1)
    @assert( degree >= 1 && degree <= 3 )
    cub, vtx = tetcubature(cubdegree, T, internal=internal)
    numnodes = cub.numnodes
    Q = zeros(T, (numnodes, numnodes, 3))
    w, Q = SummationByParts.buildsparseoperators(cub, vtx, degree)
    new(degree, numnodes, cub, vtx, w, Q)
  end
end

@doc """
### SBP.AbstractFace

`AbstractFace` is a parametric abstract type that defines face-based data and
operations (e.g. volume-to-face reconstruction, face integration, etc) for
summation-by-parts finite-difference operators.

"""->
abstract AbstractFace{T<:Number}

@doc """
### SBP.TriFace

Defines a face between two TriSBP operators with the same cubature nodes

**Fields**

* `degree` : face integration is exact for polys of degree 2*`degree`
* `numnodes` : number of cubature nodes
* `stencilsize` : number of nodes in the reconstruction stencil
* `dstencilsize` : number of nodes in the derivative operator stencils
* `cub` : a symmetric cubature type for triangle faces (i.e. edges)
* `vtx` : the vertices of the face in reference space, [-1,1]
* `wface` : mass matrix (quadrature) for the face
* `interp[:,:]` : volume-to-face-nodes reconstruction operator
* `perm[:,:]` : permutation for volume nodes so `interp` can be used on all sides
* `deriv[:,:]` : derivative operators for face-based coordinate system
* `dperm[:,:]` : permutation for volume nodes so `deriv` can be used on all sides
* `nbrperm[:,:]` : permutation for face nodes on neighbour element

"""->
immutable TriFace{T} <: AbstractFace{T}
  degree::Int
  numnodes::Int
  stencilsize::Int
  dstencilsize::Int
  cub::LineSymCub{T}
  vtx::Array{T,2}
  wface::Array{T,1}
  normal::Array{T,2}
  interp::Array{T,2}
  perm::Array{Int,2}
  deriv::Array{T,3}
  dperm::Array{Int,2}
  nbrperm::Array{Int,2}

  # inner constructor
  function TriFace(degree::Int, facecub::LineSymCub{T}, facevtx::Array{T,2},
                   interp::Array{T,2}, perm::Array{Int,2},
                   deriv::Array{T,3}, dperm::Array{Int,2})
    @assert( degree >= 1 && degree <= 5 )
    numnodes = facecub.numnodes
    @assert( size(interp,2) == size(deriv,2) == numnodes )
    normal = T[0 -1; 1 1; -1 0].'
    nbrperm = SymCubatures.getneighbourpermutation(facecub)
    wface = SymCubatures.calcweights(facecub)
    stencilsize = size(interp,1)
    dstencilsize = size(deriv,1)
    new(degree, facecub.numnodes, stencilsize, dstencilsize, facecub, facevtx, 
        wface, normal, interp, perm, deriv, dperm, nbrperm)
  end
end

@doc """
### SBP.TetFace

Defines a face between two TetSBP operators with the same cubature nodes

**Fields**

* `degree` : face integration is exact for polys of degree 2*`degree`
* `numnodes` : number of cubature nodes
* `stencilsize` : number of nodes in the reconstruction stencil
* `dstencilsize` : number of nodes in the derivative operator stencils
* `cub` : a symmetric cubature type for tetrahedral faces (i.e. triangles)
* `vtx` : the vertices of the face in the reference space of the face
* `wface` : mass matrix (quadrature) for the face
* `interp[:,:]` : volume-to-face-nodes reconstruction operator
* `perm[:,:]` : permutation for volume nodes so `interp` can be used on all sides
* `deriv[:,:]` : derivative operators for face-based coordinate system
* `dperm[:,:]` : permutation for volume nodes so `deriv` can be used on all sides
* `nbrperm[:,:]` : permutation for face nodes on neighbour element

"""->
immutable TetFace{T} <: AbstractFace{T}
  degree::Int
  numnodes::Int
  stencilsize::Int
  #dstencilsize::Int
  cub::TriSymCub{T}
  vtx::Array{T,2}
  wface::Array{T,1}
  normal::Array{T,2}
  interp::Array{T,2}
  perm::Array{Int,2}
  #deriv::Array{T,3}
  #dperm::Array{Int,2}
  nbrperm::Array{Int,2}

  # inner constructor
  function TetFace(degree::Int, facecub::TriSymCub{T}, facevtx::Array{T,2},
                   interp::Array{T,2}, perm::Array{Int,2})
    @assert( degree >= 1 && degree <= 4 )
    numnodes = facecub.numnodes
    @assert( size(interp,2) == numnodes )
    normal = T[0 0 -1; 0 -1 0; 1 1 1; -1 0 0].'
    nbrperm = SymCubatures.getneighbourpermutation(facecub)
    wface = SymCubatures.calcweights(facecub)
    stencilsize = size(interp,1)
    new(degree, numnodes, stencilsize, facecub, facevtx, wface, normal, interp,
        perm, nbrperm)
  end
end

@doc """
### SBP.UnaryFunctor

`UnaryFunctor` is abstract type that is used to specify different types of
residual updates in the useoperators.jl and usefaceoperators.jl files

"""->
abstract UnaryFunctor

@doc """
### SBP.Add

`Add` is basically an identity operator.  It is needed as a default
UnaryFunctor.

"""->
type Add <: UnaryFunctor
end
function call(update::Add, scalar)
  return scalar
end

@doc """
### SBP.Subtract

`Subtract` is used to negate values, which is useful in residual updates

"""->
type Subtract <: UnaryFunctor
end
function call(update::Subtract, scalar)
  return -scalar
end

include("outerconstructors.jl") #<--- outer constructors and factories
include("buildfaceoperators.jl") #<--- functions related to building face operators
include("buildoperators.jl") #<--- functions related to building SBP operators
include("weakdifferentiate.jl") #<--- functions for weak differentiation
include("weakdifferentiate_rev.jl") #<--- reverse-diff of weak differentiation
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
include("edgestabilize.jl") #<--- functions related to edge stabilization
include("utils.jl") # <--- miscillaneous functions

end # module
