module SummationByParts

include("orthopoly.jl")
include("symcubatures.jl")
include("cubature.jl")

include("experiment.jl")

export SBPOperator, TriSBP, TetSBP, cubatureerror

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
* `Qx` : discrete x-derivative stiffness matrix operator
* `Qy` : discrete y-derivative stiffness matrix operator

"""->
type TriSBP{T} <: SBPOperator{T}
  degree::Int
  numnodes::Int
  numbndry::Int
  x::Array{T}
  w::Array{T}
  Qx::Array{T}
  Qy::Array{T}
  Ex::Array{T}
  Ey::Array{T}

  function TriSBP(;degree::Int=1)
    @assert( degree >= 1 && degree <= 4 )
    cub, vtx = tricubature(2*degree-1, T)
    numnodes = cub.numnodes
    w = SymCubatures.calcweights(cub)
    
    new(degree, numnodes, numbndry, x, w, zeros(numnodes), zeros(numnodes))
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

@doc """
### SummationByParts.bndrynodalexpansion

Computes the transformation matrix that maps the Proriol orthogonal polynomials
to polynomials that are nodal on the boundary nodes, i.e. if E is the
transformation matrix and P is the matrix of Proriol polys, with the polynomials
listed by row, then P*E = I, when restricted to the boundary nodes.

**Inputs**

* `cub`: symmetric cubature rule for a right triangle
* `vtx`: vertices of the right triangle
* `d`: maximum total degree for the Proriol polynomials

**Outputs**

* `E`: transformation matrix

"""->
function bndrynodalexpansion{T}(cub::TriSymCub{T}, vtx::Array{T,2}, d::Int)
  numbndry = SymCubatures.getnumboundarynodes(cub)
  N = convert(Int64, (d+1)*(d+2)/2)
  xaug = zeros(T, N)
  yaug = zeros(T, N)
  x, y = SymCubatures.calcnodes(cub, vtx)
  xaug[1:numbndry] = x[1:numbndry]
  yaug[1:numbndry] = y[1:numbndry]
  # set the augmented interior nodes that make a unisolvent set for polys; use
  # uniform points in the interior for now
  ptr = numbndry+1
  for j = 1:(d-2)
    eta = 2.*j/d-1
    for i = 1:(d-j-1)
      xi = 2.*i/d-1
      xaug[ptr] = xi
      yaug[ptr] = eta
      ptr += 1
    end
  end
  V = zeros(T, (N, N))
  ptr = 1
  for r = 0:d
    for j = 0:r
      i = r-j
      V[:,ptr] = OrthoPoly.proriolpoly(xaug, yaug, i, j)
      ptr += 1
    end
  end
  return inv(V)
end

function bndrynodalexpansion{T}(cub::TetSymCub{T}, vtx::Array{T,2}, d::Int)
  numbndry = SymCubatures.getnumboundarynodes(cub)
  N = convert(Int64, (d+1)*(d+2)*(d+3)/6)
  xaug = zeros(T, N)
  yaug = zeros(T, N)
  zaug = zeros(T, N)
  x, y, z = SymCubatures.calcnodes(cub, vtx)
  xaug[1:numbndry] = x[1:numbndry]
  yaug[1:numbndry] = y[1:numbndry]
  zaug[1:numbndry] = z[1:numbndry]
  # set the augmented interior nodes that make a unisolvent set for polys; use
  # uniform points in the interior for now
  ptr = numbndry+1
  for k = 1:(d-3)
    zeta = 2.*k/d-1
    for j = 1:(d-k-2)
      eta = 2.*j/d-1
      for i = 1:(d-j-k-1)
        xi = 2.*i/d-1
        xaug[ptr] = xi
        yaug[ptr] = eta
        zaug[ptr] = zeta
        ptr += 1
      end
    end
  end
  V = zeros(T, (N, N))
  ptr = 1
  for r = 0:d
    for k = 0:r
      for j = 0:r-k
        i = r-j-k
        V[:,ptr] = OrthoPoly.proriolpoly(xaug, yaug, zaug, i, j, k)
        ptr += 1
      end
    end
  end
  return inv(V)
end

@doc """
### SummationByParts.massmatrices

Finds and returns the element mass matrix, `w`, as well as the symmetric part of
the SBP operators, `Ex`, `Ey` (`Ez`).  The latter operators coorespond to
boundary integrals in the divergence theorem, hence they are mass matrices of
the boundary facets.

**Inputs**

* `cub`: symmetric cubature rule for a right triangle
* `vtx`: vertices of the right triangle
* `d`: maximum total degree for the Proriol polynomials

**Outputs**

* `w`: diagonal element mass matrix, stored as a 1D array
* `Ex`, `Ey` (`Ez`): symmetric parts of the SBP first derivative operators

"""->
function massmatrices{T}(cub::TriSymCub{T}, vtx::Array{T,2}, d::Int)
  numbndry = SymCubatures.getnumboundarynodes(cub)
  N = convert(Int, (d+1)*(d+2)/2 )
  # compute the derivatives of the ortho polys
  P = zeros(T, (cub.numnodes,N) )
  dPdx = zeros(P)
  dPdy = zeros(P)
  x, y = SymCubatures.calcnodes(cub, vtx)
  ptr = 1
  for r = 0:d
    for j = 0:r
      i = r-j
      P[:,ptr] = OrthoPoly.proriolpoly(x, y, i, j)
      dPdx[:,ptr], dPdy[:,ptr] = OrthoPoly.diffproriolpoly(x, y, i, j)
      ptr += 1
    end
  end
  E = SummationByParts.bndrynodalexpansion(cub, vtx, d)
  P *= E
  dPdx *= E
  dPdy *= E
  # compute and return the mass matrices
  w = SymCubatures.calcweights(cub)
  Ex = zeros(T, (cub.numnodes,cub.numnodes) )
  Ey = zeros(Ex)
  Q = P.'*diagm(w)*dPdx
  Ex[1:numbndry,1:numbndry] = Q[1:numbndry,1:numbndry] + Q[1:numbndry,1:numbndry].'
  Q = P.'*diagm(w)*dPdy
  Ey[1:numbndry,1:numbndry] = Q[1:numbndry,1:numbndry] + Q[1:numbndry,1:numbndry].'
  return w, Ex, Ey
end

function massmatrices{T}(cub::TetSymCub{T}, vtx::Array{T,2}, d::Int)
  numbndry = SymCubatures.getnumboundarynodes(cub)
  N = convert(Int, (d+1)*(d+2)*(d+3)/6 )
  # compute the derivatives of the ortho polys
  P = zeros(T, (cub.numnodes,N) )
  dPdx = zeros(P)
  dPdy = zeros(P)
  dPdz = zeros(P)
  x, y, z = SymCubatures.calcnodes(cub, vtx)
  ptr = 1
  for r = 0:d
    for k = 0:r
      for j = 0:r-k
        i = r-j-k
        P[:,ptr] = OrthoPoly.proriolpoly(x, y, z, i, j, k)
        dPdx[:,ptr], dPdy[:,ptr], dPdz[:,ptr] =
          OrthoPoly.diffproriolpoly(x, y, z, i, j, k)
        ptr += 1
      end
    end
  end
  E = SummationByParts.bndrynodalexpansion(cub, vtx, d)
  P *= E
  dPdx *= E
  dPdy *= E
  dPdz *= E
  # compute and return the mass matrices
  w = SymCubatures.calcweights(cub)
  Ex = zeros(T, (cub.numnodes,cub.numnodes) )
  Ey = zeros(Ex)
  Ez = zeros(Ex)
  Q = P.'*diagm(w)*dPdx
  Ex[1:numbndry,1:numbndry] = Q[1:numbndry,1:numbndry] + Q[1:numbndry,1:numbndry].'
  Q = P.'*diagm(w)*dPdy
  Ey[1:numbndry,1:numbndry] = Q[1:numbndry,1:numbndry] + Q[1:numbndry,1:numbndry].'
  Q = P.'*diagm(w)*dPdz
  Ez[1:numbndry,1:numbndry] = Q[1:numbndry,1:numbndry] + Q[1:numbndry,1:numbndry].'
  return w, Ex, Ey, Ez
end  

# function accuracyequations()
#   x, y = SymCubatures.calcnodes(cub, vtx)  
#   # the number of unknowns for both skew-symmetric matrices Qx and Qy
#   numQvars = cub.numnodes*(cub.numnodes-1)
#   # the number of accuracy equations
#   numeqns = cub.numnodes*(d+1)*(d+2)
#   A = zeros(T, (numeqns, numQvars))
#   bx = zeros(T, numeqns)
#   by = zeros(T, numeqns)

#   # loop over ortho polys up to degree d
#   ptr = 0
#   for r = 0:d
#     for j = 0:r
#       i = r-j
#       P = OrthoPoly.proriolpoly(x, y, i, j)
#       dPdx, dPdy = OrthoPoly.diffproriolpoly(x, y, i, j)
#       # loop over the lower part of skew-symmetric matrix Qx
#       for row = 2:cub.numnodes
#         offset = convert(Int, (row-1)*(row-2)/2)
#         for col = 1:row-1
#           A(ptr+row, offset+col) += P[col]
#           A(ptr+col, offset+col) -= P[row]
#         end
#       end
#       bx[ptr+1:ptr+cub.numnodes] = diagm(w)*dPdx - Ex*P
#       by[ptr+1:ptr+cub.numnodes] = diagm(w)*dPdy - Ey*P
#       ptr += cub.numnodes
#     end
#   end
# end


end # module
