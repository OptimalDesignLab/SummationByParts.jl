module SummationByParts

include("orthopoly.jl")
include("symcubatures.jl")
include("cubature.jl")

using .OrthoPoly
using .SymCubatures
using .Cubature

export SBPOperator, TriSBP, TetSBP, applyQ!, applyD!, applyH!

@doc """
### SBP.SBPOperator

`SBPOperator` is a parametric abstract type that defines summation-by-parts
finite-difference operators.

"""->
abstract SBPOperator{T<:FloatingPoint}

@doc """
### SBP.TriSBP

Defines diagonal-norm SBP first-derivative operators on a right-triangle.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the triangle required for these operators
* `numbndry` : number of boundary nodes for the operators 
* `x` : node locations
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""->
type TriSBP{T} <: SBPOperator{T}
  degree::Int
  numnodes::Int
  numbndry::Int
  #x::Array{T}
  cub::TriSymCub{T}
  w::Array{T}
  Q::Array{T}

  function TriSBP(;degree::Int=1)
    @assert( degree >= 1 && degree <= 4 )
    cub, vtx = tricubature(2*degree-1, T)
    numnodes = cub.numnodes
    numbndry = SymCubatures.getnumboundarynodes(cub)
    Q = zeros(T, (numnodes, numnodes, 2))
    w, Q[:,:,1], Q[:,:,2] = SummationByParts.buildoperators(cub, vtx, degree)
    x = zeros(T, (2, numnodes))
    #x[1,:], x[2,:] = SymCubatures.calcnodes(cub, vtx)
    new(degree, numnodes, numbndry, cub, w, Q)
  end
end

@doc """
### SBP.TetSBP

Defines diagonal-norm SBP first-derivative operators on a right-tetrahedron.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the tetrahedron required for these operators
* `numbndry` : number of boundary nodes for the operators 
* `x` : node locations
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""->
type TetSBP{T} <: SBPOperator{T}
  degree::Int
  numnodes::Int
  numbndry::Int
  #x::Array{T}
  cub::TetSymCub{T}
  w::Array{T}
  Q::Array{T}

  function TetSBP(;degree::Int=1)
    @assert( degree >= 1 && degree <= 4 )
    cub, vtx = tetcubature(2*degree-1, T)
    numnodes = cub.numnodes
    numbndry = SymCubatures.getnumboundarynodes(cub)
    Q = zeros(T, (numnodes, numnodes, 3))
    w, Q[:,:,1], Q[:,:,2], Q[:,:,3] = SummationByParts.buildoperators(cub, vtx, degree)
    x = zeros(T, (3, numnodes))
    x[1,:], x[2,:], x[3,:] = SymCubatures.calcnodes(cub, vtx)
    new(degree, numnodes, numbndry, cub, w, Q)
  end
end

include("buildoperators.jl") #<--- functions related to building SBP operators
include("useoperators.jl") #<--- functions for applying SBP operators

end # module
