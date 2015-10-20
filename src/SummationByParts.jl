module SummationByParts

include("orthopoly.jl")
include("symcubatures.jl")
include("cubature.jl")

using ArrayViews
using PDESolverCommon
using .OrthoPoly
using .SymCubatures
using .Cubature

export AbstractSBP, TriSBP, TetSBP, Boundary, Interface, calcnodes, calcminnodedistance,
  weakdifferentiate!, differentiate!, directionaldifferentiate!, 
  volumeintegrate!, mappingjacobian!, boundaryintegrate!,
  interiorfaceintegrate!

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
* `numbndry` : number of boundary nodes for the operators
* `numfacenodes` : number of nodes on an individual face of the domain
* `facenodes` : indices of the nodes that lie on each face
* `facenormal` : unit normal to the faces, in reference space
* `cub` : a symmetric cubature type for triangles
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `wface` : mass matrix (dense) for the faces
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""->
immutable TriSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  numbndry::Int
  numfacenodes::Int
  facenodes::Array{Int,2}
  facenormal::Array{T,2}
  cub::TriSymCub{T}
  w::Array{T,1}
  wface::Array{T,2}
  Q::Array{T,3}

  function TriSBP(;degree::Int=1, faceorder::Array{Int,1}=[1;2;3], 
                  bubble::Int=-1, reorder=true)
    @assert( degree >= 1 && degree <= 4 )
    cub, vtx = tricubature(2*degree-1, T)
    numnodes = cub.numnodes
    numbndry = SymCubatures.getnumboundarynodes(cub)
    numfacenodes = SymCubatures.getnumfacenodes(cub)
    wface, facenodes = boundarymassmatrix(cub, vtx, degree)
    facenormal = T[0 -1; 1 1; -1 0].'
    Q = zeros(T, (numnodes, numnodes, 2))
    if bubble < 0
      w, Q[:,:,1], Q[:,:,2] = SummationByParts.buildoperators(cub, vtx, degree)
    else
      w, Q[:,:,1], Q[:,:,2] = SummationByParts.buildoperators(cub, vtx, degree,
                                                              bubble)
    end
    if reorder
      # reorder the nodes
      perm, faceperm = SummationByParts.getnodepermutation(cub, degree)
      w = w[perm]
      Q[:,:,1] = Q[perm,perm,1]
      Q[:,:,2] = Q[perm,perm,2]
      perminv = invperm(perm)
      for k = 1:3
        facenodes[:,k] = perminv[facenodes[:,k]]
      end
      wface = wface[faceperm,faceperm]
      # reorder the faces
      facenodes[:,:] = facenodes[faceperm,faceorder]
      facenormal = facenormal[:,faceorder]
    end

    new(degree, numnodes, numbndry, numfacenodes, facenodes, facenormal, cub,
        w, wface, Q)
  end
end

@doc """
### SBP.TetSBP

Defines diagonal-norm SBP first-derivative operators on a right-tetrahedron.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the tetrahedron required for these operators
* `numbndry` : number of boundary nodes for the operators 
* `numfacenodes` : number of nodes on an individual face of the domain
* `facenodes` : indices of the nodes that lie on each face
* `facenormal` : unit normal to the faces, in reference space
* `cub` : a symmetric cubature type for tetrahedra
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `wface` : mass matrix (dense) for the faces
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""->
immutable TetSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  numbndry::Int
  numfacenodes::Int
  facenodes::Array{Int,2}
  facenormal::Array{T,2}
  cub::TetSymCub{T}
  w::Array{T,1}
  wface::Array{T,2}
  Q::Array{T,3}

  function TetSBP(;degree::Int=1, faceorder::Array{Int,1}=[1;2;3;4],
                  reorder=true)
    @assert( degree >= 1 && degree <= 4 )
    cub, vtx = tetcubature(2*degree-1, T)
    numnodes = cub.numnodes
    numbndry = SymCubatures.getnumboundarynodes(cub)
    numfacenodes = SymCubatures.getnumfacenodes(cub)
    wface, facenodes = boundarymassmatrix(cub, vtx, degree)
    facenormal = T[0 0 -1; 1 1 1; -1 0 0; 0 -1 0].'
    Q = zeros(T, (numnodes, numnodes, 3))
    w, Q[:,:,1], Q[:,:,2], Q[:,:,3] = SummationByParts.buildoperators(cub, vtx, degree)

    if reorder
      # reorder the nodes
      perm, faceperm = SummationByParts.getnodepermutation(cub, degree)
      w = w[perm]
      Q[:,:,1] = Q[perm,perm,1]
      Q[:,:,2] = Q[perm,perm,2]
      Q[:,:,3] = Q[perm,perm,3]
      perminv = invperm(perm)
      for k = 1:4
        facenodes[:,k] = perminv[facenodes[:,k]]
      end
      wface = wface[faceperm,faceperm]
      # reorder the faces
      facenodes[:,:] = facenodes[faceperm,faceorder]
      facenormal = facenormal[:,faceorder]
    end

    new(degree, numnodes, numbndry, numfacenodes, facenodes, facenormal, cub,
        w, wface, Q)
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

Defines a face between two TriSBP operators of the same order.

**Fields**

* `degree` : maximum polynomial degree for which face integration is exact
* `cub` : a symmetric cubature type for triangle faces (i.e. edges)
* `wface` : mass matrix (quadrature) for the face
* `R[:,:]` : volume-to-face-nodes reconstruction operator
* `perm[:,:]` : permutation for volume nodes so `R` can be used on all sides

"""->
immutable TriFace{T} <: AbstractFace{T}
  degree::Int
  cub::TetSymCub{T}
  wface::Array{T}
  R::Array{T,2}
  perm::Array{Int,2}
end


include("buildoperators.jl") #<--- functions related to building SBP operators
include("useoperators.jl") #<--- functions for applying SBP operators

end # module
