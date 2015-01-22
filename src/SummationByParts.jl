module SummationByParts

include("orthopoly.jl")
include("symcubatures.jl")
include("cubature.jl")

export SBPOperator, TriSBP, TetSBP

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
* `x` : node locations
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Qx` : discrete x-derivative stiffness matrix operator
* `Qy` : discrete y-derivative stiffness matrix operator

"""->
type TriSBP{T} <: SBPOperator{T}
  degree::Int
  numnodes::Int
  x::Array{T}
  w::Array{T}
  Qx::Array{T}
  Qy::Array{T}

  function TriSBP(;degree::Int=1)
    @assert( degree >= 1 && degree <= 4 )
    w, x = tricubature(2*degree-1, T)
    numnodes = length(w)
    new(degree, numnodes, x, w, zeros(numnodes), zeros(numnodes))
  end
end

@doc """
### SBP.TetSBP

Defines diagonal-norm SBP first-derivative operators on a right-tetrahedron.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the tetrahedron required for these operators
* `x` : node locations
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Qx` : discrete x-derivative stiffness matrix operator
* `Qy` : discrete y-derivative stiffness matrix operator
* `Qz` : discrete z-derivative stiffness matrix operator

"""->
type TetSBP{T} <: SBPOperator{T}
  degree::Int
  numnodes::Int
  x::Array{T}
  w::Array{T}
  Qx::Array{T}
  Qy::Array{T}
  Qz::Array{T}
end

end # module
