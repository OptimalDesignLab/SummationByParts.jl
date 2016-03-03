# This file gathers together functions used to build the SBP operators

@doc """
### SummationByParts.bndrynodalexpansion

Computes the transformation matrix that maps the Proriol orthogonal polynomials
to polynomials that are nodal on the boundary nodes, i.e. if E is the
transformation matrix and P is the matrix of Proriol polys, with the polynomials
listed by row, then P*E = I, when restricted to the boundary nodes.

**Inputs**

* `cub`: symmetric cubature rule
* `vtx`: vertices of the right simplex
* `d`: maximum total degree for the Proriol polynomials

**Outputs**

* `E`: transformation matrix

"""->
function bndrynodalexpansion{T}(cub::TriSymCub{T}, vtx::Array{T,2}, d::Int)
  numbndry = SymCubatures.getnumboundarynodes(cub)
  N = convert(Int64, (d+1)*(d+2)/2)
  xaug = zeros(T, N)
  yaug = zeros(T, N)
  x = SymCubatures.calcnodes(cub, vtx)
  # get the unique boundary node indices
  bndryindices = SymCubatures.getbndrynodeindices(cub)
  xaug[1:numbndry] = x[1,bndryindices]
  yaug[1:numbndry] = x[2,bndryindices]
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
  x = SymCubatures.calcnodes(cub, vtx)
  # get the unique boundary node indices
  bndryindices = SymCubatures.getbndrynodeindices(cub)
  xaug[1:numbndry] = x[1,bndryindices]
  yaug[1:numbndry] = x[2,bndryindices]
  zaug[1:numbndry] = x[3,bndryindices]
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
### SummationByParts.nodalexpansion

Computes the transformation matrix that maps the Proriol orthogonal polynomials
to polynomials that are nodal on the cubature nodes, i.e. if C is the
transformation matrix and P is the matrix of Proriol polys, with the polynomials
listed by row, then P*C = I, when restricted to the nodes.  This function can be
used to construct the operators in the spectral-element method of Giraldo and
Tayler.

**Inputs**

* `cub`: symmetric cubature rule
* `vtx`: vertices of the right simplex
* `d`: maximum total degree on the edges
* `e`: maximum degree of bubble functions on the interior

**Outputs**

* `C`: transformation matrix

"""->
function nodalexpansion{T}(cub::TriSymCub{T}, vtx::Array{T,2}, d::Int, e::Int)
  #error("This method was broken by nodal reordering; it is currently unused")
  numbndry = SymCubatures.getnumboundarynodes(cub)
  numbub = cub.numnodes - numbndry
  N = convert(Int64, (d+1)*(d+2)/2)

  # Step 1: compute the degree d expansion for the edge nodes.  Note that
  # (xaug[1:N],yaug[1:N]) is composed of the boundary nodes + interior nodes that
  # make Vandermonde unisolvent.  The last numint nodes in (xaug,yaug) are the
  # actual interior nodes.
  xaug = zeros(T, N+numbub)
  yaug = zeros(T, N+numbub)
  x = SymCubatures.calcnodes(cub, vtx)
  # get the unique boundary node indices
  bndryindices = SymCubatures.getbndrynodeindices(cub)
  xaug[1:numbndry] = x[1,bndryindices]
  yaug[1:numbndry] = x[2,bndryindices]
  ptr = numbndry+1
  # set uniform nodes on interior to make Vandermonde unisolvent
  for j = 1:(d-2)
    eta = 2.*j/d-1
    for i = 1:(d-j-1)
      xi = 2.*i/d-1
      xaug[ptr] = xi
      yaug[ptr] = eta
      ptr += 1
    end
  end
  # these are the actual interior nodes used for the cubature
  indices = SymCubatures.getinteriornodeindices(cub)
  xaug[N+1:end] = x[1,indices]
  yaug[N+1:end] = x[2,indices]

  Vbndry = zeros(T, (N, N))
  Pbub = zeros(T, (numbub, N))
  ptr = 1
  for r = 0:d
    for j = 0:r
      i = r-j
      P = OrthoPoly.proriolpoly(xaug, yaug, i, j)
      Vbndry[:,ptr] = P[1:N]
      Pbub[:,ptr] = P[N+1:end]
      ptr += 1
    end
  end

  # Step 2: compute error in boundary basis at interior nodes
  invVbndry = inv(Vbndry)
  Perr = Pbub*invVbndry

  # Step 3: construct an augmented set of nodes that will be used to correct the
  # boundary basis expansion and create the interior basis expansion
  Naug = convert(Int, (e+1)*(e+2)/2)
  xaug = zeros(T, (Naug) )
  yaug = zeros(T, (Naug) )
  xaug[1:numbub] = x[1,indices] 
  yaug[1:numbub] = x[2,indices] 
  ptr = numbub+1
  for j = 0:e-1
    xi = 2.*j/e-1
    # bottom side
    xaug[ptr] = xi
    yaug[ptr] = -1.0
    ptr += 1
    # left side
    xaug[ptr] = -1.0
    yaug[ptr] = -xi
    ptr += 1
    # hypotenuse
    xaug[ptr] = -xi
    yaug[ptr] = xi
    ptr += 1
  end
  
  # Step 4: find the boundary corrections and interior basis expansion
  Vbub = zeros(T, (Naug, Naug) )
  ptr = 1
  for r = 0:e
    for j = 0:r
      i = r-j
      Vbub[:,ptr] = OrthoPoly.proriolpoly(xaug, yaug, i, j)
      ptr += 1
    end
  end
  invVbub = inv(Vbub)
  C = zeros(T, (Naug, cub.numnodes))
  C[1:N,bndryindices] = invVbndry[:,1:numbndry]
  C[:,indices] += invVbub[:,1:numbub]
  # corrections
  for i = 1:numbndry
    C[:,bndryindices[i]] -= invVbub[:,1:numbub]*Perr[:,i]
  end
  return C
end

@doc """
### SummationByParts.boundaryoperators

Finds the symmetric part of the SBP operators, `Ex`, `Ey` (`Ez`).  These
operators coorespond to boundary integrals in the divergence theorem, and are
related to the mass matrices of the boundary faces.

**Inputs**

* `cub`: symmetric cubature rule
* `vtx`: vertices of the right simplex
* `d`: maximum total degree for the Proriol polynomials

**Outputs**

* `Ex`, `Ey` (`Ez`): symmetric parts of the SBP first derivative operators

"""->
function boundaryoperators{T}(cub::TriSymCub{T}, vtx::Array{T,2}, d::Int)
  numbndry = SymCubatures.getnumboundarynodes(cub)
  N = convert(Int, (d+1)*(d+2)/2 )
  # compute the derivatives of the ortho polys
  P = zeros(T, (cub.numnodes,N) )
  dPdx = zeros(P)
  dPdy = zeros(P)
  x = SymCubatures.calcnodes(cub, vtx)
  ptr = 1
  for r = 0:d
    for j = 0:r
      i = r-j
      P[:,ptr] = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      dPdx[:,ptr], dPdy[:,ptr] = 
      OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      ptr += 1
    end
  end
  E = SummationByParts.bndrynodalexpansion(cub, vtx, d)
  P *= E
  dPdx *= E
  dPdy *= E
  # compute and return the boundary operators
  w = SymCubatures.calcweights(cub)
  indices = SymCubatures.getbndrynodeindices(cub)  
  Ex = zeros(T, (cub.numnodes,cub.numnodes) )
  Ey = zeros(Ex)
  Q = P.'*diagm(w)*dPdx
  Ex[indices,indices] = Q[1:numbndry,1:numbndry] + Q[1:numbndry,1:numbndry].'
  Q = P.'*diagm(w)*dPdy
  Ey[indices,indices] = Q[1:numbndry,1:numbndry] + Q[1:numbndry,1:numbndry].'
  return Ex, Ey
end

function boundaryoperators{T}(cub::TetSymCub{T}, vtx::Array{T,2}, d::Int)
  numbndry = SymCubatures.getnumboundarynodes(cub)
  N = convert(Int, (d+1)*(d+2)*(d+3)/6 )
  # compute the derivatives of the ortho polys
  P = zeros(T, (cub.numnodes,N) )
  dPdx = zeros(P)
  dPdy = zeros(P)
  dPdz = zeros(P)
  x = SymCubatures.calcnodes(cub, vtx)
  ptr = 1
  for r = 0:d
    for k = 0:r
      for j = 0:r-k
        i = r-j-k
        P[:,ptr] = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), vec(x[3,:]),
                                         i, j, k)
        dPdx[:,ptr], dPdy[:,ptr], dPdz[:,ptr] =
          OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), vec(x[3,:]),
                                    i, j, k)
        ptr += 1
      end
    end
  end
  E = SummationByParts.bndrynodalexpansion(cub, vtx, d)
  P *= E
  dPdx *= E
  dPdy *= E
  dPdz *= E
  # compute and return the boundary operators
  w = SymCubatures.calcweights(cub)
  indices = SymCubatures.getbndrynodeindices(cub)  
  Ex = zeros(T, (cub.numnodes,cub.numnodes) )
  Ey = zeros(Ex)
  Ez = zeros(Ex)
  Q = P.'*diagm(w)*dPdx
  Ex[indices,indices] = Q[1:numbndry,1:numbndry] + Q[1:numbndry,1:numbndry].'
  Q = P.'*diagm(w)*dPdy
  Ey[indices,indices] = Q[1:numbndry,1:numbndry] + Q[1:numbndry,1:numbndry].'
  Q = P.'*diagm(w)*dPdz
  Ez[indices,indices] = Q[1:numbndry,1:numbndry] + Q[1:numbndry,1:numbndry].'
  return Ex, Ey, Ez
end  

@doc """
The following will eventually supercede current version and is based on the face operators

Finds the symmetric part of the SBP operators, e.g. `Ex`.  These operators
coorespond to boundary integrals in the divergence theorem, and are related to
the mass matrices of the boundary faces.

**Inputs**

* `face`: SBP face operator for a particular element
* `di`: desired coordinate direction of operator

**InOuts**

* `E`: symmetric part of the SBP first derivative operator in direction di

"""->
function boundaryoperator!{T}(face::AbstractFace{T}, di::Int,
                              E::AbstractArray{T,2})
  @assert( size(face.interp,1) <= size(E,1) == size(E,2) )
  @assert( di >= 1 && di <= 3)
  fill!(E, zero(T))
  # loop over faces of the element
  for findex = 1:size(face.perm,2) 
    E[face.perm[:,findex],face.perm[:,findex]] += 
      face.interp*diagm(face.wface.*face.normal[di,findex])*face.interp.'
  end
  return E
end

@doc """
### SummationByParts.boundarymassmatrix

Returns the (dense) mass matrix for a set of nodes on a reference boundary.
This mass matrix can be used for boundary integration or to impose boundary
conditions weakly.  The array `bndryindices` is also returned, which is a list
of element-node indices for each boundary (see also
SymCubatures.getfacenodeindices).

**Inputs**

* `cub`: symmetric cubature rule
* `vtx`: vertices of the right simplex
* `d`: maximum total degree for the polynomials

**Outputs**

* `Hbndry`: reference boundary mass matrix
* `bndryindices`: list of nodes that lie on each boundary

"""->
function boundarymassmatrix{T}(cub::TriSymCub{T}, vtx::Array{T,2}, d::Int)
  bndryindices = SymCubatures.getfacenodeindices(cub)
  numbndrynodes = size(bndryindices,1)
  x = SymCubatures.calcnodes(cub, vtx)
  xbndry = vec(x[1,bndryindices[:,1]])
  P = zeros(T, (numbndrynodes,numbndrynodes))
  for j = 0:d
    P[:,j+1] = OrthoPoly.jacobipoly(xbndry, 0.0, 0.0, j)
  end
  A = inv(P)
  Hbndry = A.'*A
  return Hbndry, bndryindices
end

function boundarymassmatrix{T}(cub::TetSymCub{T}, vtx::Array{T,2}, d::Int)
  bndryindices = SymCubatures.getfacenodeindices(cub)
  numbndrynodes = size(bndryindices,1)
  x = SymCubatures.calcnodes(cub, vtx)
  xbndry = vec(x[1,bndryindices[:,1]])
  ybndry = vec(x[2,bndryindices[:,1]])
  P = zeros(T, (numbndrynodes,numbndrynodes))
  ptr = 1
  for r = 0:d
    for j = 0:r
      i = r-j
      P[:,ptr] = OrthoPoly.proriolpoly(xbndry, ybndry, i, j)
      ptr += 1
    end
  end
  A = inv(P)
  Hbndry = A.'*A
  return Hbndry, bndryindices
end

@doc """
### SummationByParts.accuracyconstraints

Returns the accuracy constraints on the asymmetric part of the SBP stiffness
matrices.  These constraints are linear, and for each coordinate-direction
operator (i.e. Qx, Qy,...) the system matrix `A` is the same; only the
right-hand side changes.

The columns in `A` are ordered assuming only the strictly lower triangular part
of the operators are the unknowns.  These unknowns are ordered by row and then
column.  For example, entry Q_21 = -Q_12 is the number 1 variable, and
Q_32 = -Q_23 is the number 3 variable.

**Inputs**

* `cub`: symmetric cubature rule
* `vtx`: vertices of the right simplex
* `d`: maximum total degree for the Proriol polynomials
* `E`: the symmetric part of the SBP stiffness matrices

**Outputs**

* `A`: the system matrix for the linear accuracy constraints
* `bx`,`by` (`bz`): the right-hand-sides of the accuracy constraints

"""->
function accuracyconstraints{T}(cub::TriSymCub{T}, vtx::Array{T,2}, d::Int,
                                E::Array{T,3})
  x = SymCubatures.calcnodes(cub, vtx)
#  Ex, Ey = SummationByParts.boundaryoperators(cub, vtx, d)
  w = SymCubatures.calcweights(cub)
#  Ex *= 0.5; Ey *= 0.5
  # the number of unknowns for in the skew-symmetric matrices
  numQvars = convert(Int, cub.numnodes*(cub.numnodes-1)/2)
  # the number of accuracy equations
  numeqns = convert(Int, cub.numnodes*(d+1)*(d+2)/2)
  A = zeros(T, (numeqns, numQvars))
  bx = zeros(T, numeqns)
  by = zeros(T, numeqns)
  # loop over ortho polys up to degree d
  ptr = 0
  for r = 0:d
    for j = 0:r
      i = r-j
      P = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      dPdx, dPdy = OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      # loop over the lower part of the skew-symmetric matrices
      for row = 2:cub.numnodes
        offset = convert(Int, (row-1)*(row-2)/2)
        for col = 1:row-1
          A[ptr+row, offset+col] += P[col]
          A[ptr+col, offset+col] -= P[row]
        end
      end
      bx[ptr+1:ptr+cub.numnodes] = diagm(w)*dPdx - E[:,:,1]*P
      by[ptr+1:ptr+cub.numnodes] = diagm(w)*dPdy - E[:,:,2]*P
      ptr += cub.numnodes
    end
  end
  return A, bx, by
end

function accuracyconstraints{T}(cub::TetSymCub{T}, vtx::Array{T,2}, d::Int)
  x = SymCubatures.calcnodes(cub, vtx) 
  Ex, Ey, Ez = SummationByParts.boundaryoperators(cub, vtx, d)
  w = SymCubatures.calcweights(cub)
  Ex *= 0.5; Ey *= 0.5; Ez *= 0.5
  # the number of unknowns for both skew-symmetric matrices Qx, Qy, and Qz
  numQvars = convert(Int, cub.numnodes*(cub.numnodes-1)/2)
  # the number of accuracy equations
  numeqns = convert(Int, cub.numnodes*(d+1)*(d+2)*(d+3)/6)
  A = zeros(T, (numeqns, numQvars))
  bx = zeros(T, numeqns)
  by = zeros(T, numeqns)
  bz = zeros(T, numeqns)
  # loop over ortho polys up to degree d
  ptr = 0
  for r = 0:d
    for k = 0:r
      for j = 0:r-k
        i = r-j-k
        P = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), vec(x[3,:]), i, j, k)
        dPdx, dPdy, dPdz = OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]),
                                                     vec(x[3,:]), i, j, k)
        # loop over the lower part of the skew-symmetric matrices
        for row = 2:cub.numnodes
          offset = convert(Int, (row-1)*(row-2)/2)
          for col = 1:row-1
            A[ptr+row, offset+col] += P[col]
            A[ptr+col, offset+col] -= P[row]
          end
        end
        bx[ptr+1:ptr+cub.numnodes] = diagm(w)*dPdx - Ex*P
        by[ptr+1:ptr+cub.numnodes] = diagm(w)*dPdy - Ey*P
        bz[ptr+1:ptr+cub.numnodes] = diagm(w)*dPdz - Ez*P
        ptr += cub.numnodes
      end
    end
  end
  return A, bx, by, bz
end

@doc """
### SummationByParts.commuteerror

Returns the commute-error objective value.  For 2D SBP operators, this is
defined as ||H*(Dx*Dy - Dy*Dx)||^2, where the norm is the Frobenius norm.  For
3D operators, the error is defined as ||H*(Dx*Dy - Dy*Dx)||^2 + ||H*(Dx*Dz -
Dz*Dx)||^2 + ||H*(Dy*Dz - Dz*Dx||^2.

**Inputs**

* `w`: cubature rule weights
* `Qxpart`,`Qypart` (`Qzpart`): Q operators that satisfy the accuracy conditions
* `Z`: basis for the null space of the accuracy constraints (may be empty)
* `reducedsol`: the weights for `Z`; the first [1:numnodes] elements are for Qx

**Outputs**

* `f`: commute-error objective value

"""->
function commuteerror{T,T2}(w::Array{T}, Qxpart::AbstractArray{T,2},
                            Qypart::AbstractArray{T,2},
                            Z::Array{T,2}, reducedsol::Array{T2})
  # build Qx and Qy
  Qx = convert(Array{T2,2}, Qxpart)
  Qy = convert(Array{T2,2}, Qypart)
  numnodes = length(w)
  Qxnull = Z*reducedsol[1:size(Z,2)]
  Qynull = Z*reducedsol[size(Z,2)+1:end]
  for row = 2:numnodes
    offset = convert(Int, (row-1)*(row-2)/2)
    for col = 1:row-1
      Qx[row,col] += Qxnull[offset+col]
      Qx[col,row] -= Qxnull[offset+col]
      Qy[row,col] += Qynull[offset+col]
      Qy[col,row] -= Qynull[offset+col]
    end
  end
  # compute f = Frobenius_norm(H*(Dx*Dy - Dy*Dx)), and its derivatives
  f = zero(T2)
  dfdQx = zeros(Qx)
  dfdQy = zeros(Qy)
  for row = 1:numnodes
    for col = 1:numnodes
      Aij = zero(T2)
      for k = 1:numnodes
        Aij += (Qx[row,k]*Qy[k,col] - Qy[row,k]*Qx[k,col])/w[k]
      end
      f += 0.5*Aij*Aij
      for k = 1:numnodes
        dfdQx[row,k] += Aij*Qy[k,col]/w[k]
        dfdQx[k,col] -= Aij*Qy[row,k]/w[k]
        dfdQy[k,col] += Aij*Qx[row,k]/w[k]
        dfdQy[row,k] -= Aij*Qx[k,col]/w[k]
      end
    end
  end
  dfdQxnull = zeros(Qxnull)
  dfdQynull = zeros(Qynull)
  for row = 2:numnodes
    offset = convert(Int, (row-1)*(row-2)/2)
    for col = 1:row-1
      dfdQxnull[offset+col] += dfdQx[row,col] - dfdQx[col,row]
      dfdQynull[offset+col] += dfdQy[row,col] - dfdQy[col,row]
    end
  end
  dfdreducedsol = [Z.'*dfdQxnull; Z.'*dfdQynull]
  return f, dfdreducedsol
end

@doc """
### SummationByParts.buildoperators

Construct and return the SBP matrix operators, specifically the diagonal norm
matrix and the stiffness matrices.

**Inputs**

* `cub`: symmetric cubature rule
* `vtx`: vertices of the right simplex
* `d`: maximum total degree for the Proriol polynomials

**Outputs**

* `w`: the diagonal norm stored as a 1D array
* `Q`: the stiffness matrices

"""->
function buildoperators{T}(cub::TriSymCub{T}, vtx::Array{T,2}, d::Int;
                           internal::Bool=false)
  w = SymCubatures.calcweights(cub)
  face = TriFace{T}(d, cub, vtx)
  Q = zeros(T, (cub.numnodes,cub.numnodes,2) )
  SummationByParts.boundaryoperator!(face, 1, slice(Q,:,:,1))
  SummationByParts.boundaryoperator!(face, 2, slice(Q,:,:,2))
  scale!(Q, 0.5)
  A, bx, by = SummationByParts.accuracyconstraints(cub, vtx, d, Q)
  # use the minimum norm least-squares solution
  #Afact = qrfact(A)
  #x = Afact\bx; y = Afact\by
  #x = A\bx; y = A\by
  # use the minimum norm least-squares solution
  Ainv = pinv(A)
  x = Ainv*bx; y = Ainv*by
  for row = 2:cub.numnodes
    offset = convert(Int, (row-1)*(row-2)/2)
    for col = 1:row-1
      Q[row,col,1] += x[offset+col]
      Q[col,row,1] -= x[offset+col]
      Q[row,col,2] += y[offset+col]
      Q[col,row,2] -= y[offset+col]
    end
  end
  return w, Q
end

function buildoperators{T}(cub::TetSymCub{T}, vtx::Array{T,2}, d::Int)
  w = SymCubatures.calcweights(cub)
  Qx, Qy, Qz = SummationByParts.boundaryoperators(cub, vtx, d)
  A, bx, by, bz = SummationByParts.accuracyconstraints(cub, vtx, d)
  # use the minimum norm least-squares solution
  Afact = qrfact(A)
  x = Afact\bx; y = Afact\by; z = Afact\bz
  #Ainv = pinv(A)
  #x = Ainv*bx; y = Ainv*by; z = Ainv*bz
  scale!(Qx, 0.5)
  scale!(Qy, 0.5)
  scale!(Qz, 0.5)
  for row = 2:cub.numnodes
    offset = convert(Int, (row-1)*(row-2)/2)
    for col = 1:row-1
      Qx[row,col] += x[offset+col]
      Qx[col,row] -= x[offset+col]
      Qy[row,col] += y[offset+col]
      Qy[col,row] -= y[offset+col]
      Qz[row,col] += z[offset+col]
      Qz[col,row] -= z[offset+col]
    end
  end
  return w, Qx, Qy, Qz
end

@doc """
  
**Note**: when a fifth input is included, it is interpreted as the degree of the
  interior bubble functions, and the Q operators returned are the
  spectral-element differentiation operators of Giraldo and Taylor

"""->
function buildoperators{T}(cub::TriSymCub{T}, vtx::Array{T,2}, d::Int, e::Int)
  w = SymCubatures.calcweights(cub)
  N = convert(Int, (e+1)*(e+2)/2 )
  # compute the derivatives of the ortho polys
  P = zeros(T, (cub.numnodes, N) )
  dPdx = zeros(P)
  dPdy = zeros(P)
  x = SymCubatures.calcnodes(cub, vtx)
  ptr = 1
  for r = 0:e
    for j = 0:r
      i = r-j
      P[:,ptr] = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      dPdx[:,ptr], dPdy[:,ptr] =
        OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      ptr += 1
    end
  end
  C = SummationByParts.nodalexpansion(cub, vtx, d, e)
  P *= C
  dPdx *= C
  dPdy *= C
  # compute and return the operators
  w = SymCubatures.calcweights(cub)
  Q = zeros(T, (cub.numnodes,cub.numnodes,2) )
  Q[:,:,1] = P.'*diagm(w)*dPdx
  Q[:,:,2] = P.'*diagm(w)*dPdy

  # # The following forms the skew-symmetric form of the SEM operators
  # face = TriFace{T}(degree=d, faceonly=true)
  # SummationByParts.boundaryoperator!(face, 1, slice(Q,:,:,1))
  # SummationByParts.boundaryoperator!(face, 2, slice(Q,:,:,2))
  # Q[:,:,1] += P.'*diagm(w)*dPdx - dPdx.'*diagm(w)*P
  # Q[:,:,2] += P.'*diagm(w)*dPdy - dPdy.'*diagm(w)*P
  # scale!(Q, 0.5)

  return w, Q  
end

@doc """
### SummationByParts.getnodepermutation

The node ordering produced by SymCubature is not convenient for mapping local to
global node indices in the global residual assembly procedure.  This function
returns a reordering that is more suited for local-to-global mapping.

*Note*: the edge parameters of `cub` are assumed to be less than 0.5.

*WARNING*: the `faceperm` array has not been thoroughly tested

**Inputs**

* `cub`: symmetric cubature rule
* `d`: maximum polynomial degree for which the SBP operator is exact

**Outputs**

* `perm`: a permutation vector of indices
* `faceperm`: a permutation vector for the face indices

"""->
function getnodepermutation{T}(cub::TriSymCub{T}, d::Int)
  perm = zeros(Int, (cub.numnodes))
  ptr = 0 # index pointer for TriSymCub nodes
  permptr = 0 # index pointer for reordered nodes
  paramptr = 0 # index pointer for free-node parameters
  permptrint = SymCubatures.getnumboundarynodes(cub)
  
  # first, deal with 3-symmetries
  if cub.vertices
    perm[permptr+1:permptr+3] = [1;2;3] # vertices are unchanged
    ptr += 3
    permptr += 3
  end
  if cub.midedges
    # if midedges present there must be an odd number of nodes along the edge
    perm[permptr + div(d-1,2) + 1] = ptr+1
    perm[permptr + (d-1) + div(d-1,2) + 1] = ptr+2
    perm[permptr + 2*(d-1) + div(d-1,2) + 1] = ptr+3
    ptr += 3
  end
  for i = 1:cub.numS21
    # set S21 orbit nodes
    perm[permptrint+1:permptrint+3] = [ptr+1;ptr+2;ptr+3]
    permptrint += 3
    ptr += 3
    paramptr += 1
  end
  # next, deal with 6-symmetries

  # edge nodes in sequence along an edge, which requires that we order the
  # respective edge parameters, accounting for symmetry about alpha = 1/2
  edgeparam = cub.params[paramptr+1:paramptr+cub.numedge]
  #for i = 1:cub.numedge
  #  edgeparam[i] > 0.5 ? edgeparam[i] = 1 - edgeparam[i] : nothing
  #end
  # smaller parameters are further from the nodes, so do a reverse sort
  edgeperm = sortperm(edgeparam, rev=true)
  for i = 1:cub.numedge
    perm[permptr + i] = ptr + (edgeperm[i]-1)*6 + 2
    perm[permptr + (d-1) - i + 1] = ptr + (edgeperm[i]-1)*6 + 1
    perm[permptr + (d-1) + i] = ptr + (edgeperm[i]-1)*6 + 4
    perm[permptr + 2*(d-1) - i + 1] = ptr + (edgeperm[i]-1)*6 + 3
    #perm[permptr + 2*(d-1) + i] = ptr + (edgeperm[i]-1)*6 + 5
    #perm[permptr + 3*(d-1) - i + 1] = ptr + (edgeperm[i]-1)*6 + 6
    perm[permptr + 3*(d-1) - i + 1] = ptr + (edgeperm[i]-1)*6 + 5
    perm[permptr + 2*(d-1) + i] = ptr + (edgeperm[i]-1)*6 + 6
  end
  ptr += 6*cub.numedge
  permptr += 3*(d-1)
  paramptr += cub.numedge
 
  for i = 1:cub.numS111
    perm[permptrint+1:permptrint+6] = [ptr+1:ptr+6;]
    ptr += 6
    permptrint + 6
    paramptr += 2
  end
  if cub.centroid
    perm[permptrint+1] = ptr+1
    permptrint += 1
    ptr += 1
  end
  @assert(ptr == permptrint == cub.numnodes)

  # Next, find the permutation for the face nodes
  numbndrynodes = SymCubatures.getnumfacenodes(cub)
  faceperm = zeros(Int, (numbndrynodes))
  faceperm[1] = 1
  faceperm[numbndrynodes] = 2
  if cub.midedges
    # if midedges present there must be an odd number of nodes along the edge
    midnode = div(d+2,2)
    faceperm[midnode] = 3
    faceperm[2:(midnode-1)] = edgeperm + size(edgeperm,1) + 3
    faceperm[midnode+1:end-1] = midnode + edgeperm
  else
    faceperm[2:size(edgeperm,1)+1] = edgeperm + size(edgeperm,1) + 2
    faceperm[size(edgeperm,1)+2:end-1] = edgeperm + size(edgeperm,1) + 1
  end

  return perm, faceperm
end

function getnodepermutation{T}(cub::TetSymCub{T}, d::Int)
  @assert(d >= 1 && d <= 4, "implemented for d in [1,4] only")
  perm = zeros(Int, (cub.numnodes))
  ptr = 0 # index pointer for TetSymCub nodes
  permptr = 0 # index pointer for reordered nodes
  paramptr = 0 # index pointer for free-node parameters
  permptrint = SymCubatures.getnumboundarynodes(cub)
  numface = div((d-1)*(d-2),2) # this would need to change for d > 4

  # first, deal with 4-symmetries
  if cub.vertices
    perm[permptr+1:permptr+4] = [ptr+1:ptr+4;] # vertices are unchanged
    ptr += 4
    permptr += 4
  end
  fc = 0
  if cub.facecentroid
    # first, shift permptr to account for the edge nodes
    permptr += 12*cub.numedge
    cub.midedges ? permptr += 6 : nothing
    # assume face centroids are first in the face ordering
    for face = 1:4
      perm[permptr + (face-1)*numface + 1] = ptr + face
    end
    ptr += 4
    permptr -= 12*cub.numedge
    cub.midedges ? permptr -= 6 : nothing
    fc = 1
  end
  for i = 1:cub.numS31
    # set S31 orbit nodes
    perm[permptrint+1:permptrint+4] = [ptr+1:ptr+4;]
    permptrint += 4
    ptr += 4
    paramptr += 1
  end
  # next, deal with 6-symmetries
  if cub.midedges
    # if midedges present there must be an odd number of nodes along the edge
    for edge = 1:6
      perm[permptr + (edge-1)*(d-1) + div(d-1,2) + 1] = ptr + edge
    end
    ptr += 6
  end
  for i = 1:cub.numS22
    # set S22 oribt nodes
    perm[permptrint+1:permptrint+6] = [ptr+1:ptr+6;]
    permptrint += 6
    ptr += 6
    paramptr += 1
  end
  # next, deal with 12-symmetries
  # permute edge nodes; there are d-1 nodes on an edge, if vertices are excluded.
  # Edge nodes in sequence along an edge, which requires that we order the
  # respective edge parameters, accounting for symmetry about alpha = 1/2
  edgeparam = cub.params[paramptr+1:paramptr+cub.numedge]
  #for i = 1:cub.numedge
  #  edgeparam[i] > 0.5 ? edgeparam[i] = 1.0 - edgeparam[i] : nothing
  #end
  # smaller parameters are further from the nodes, so do a reverse sort
  edgeperm = sortperm(edgeparam, rev=true)
  for i = 1:cub.numedge
    for edge = 1:6
      perm[permptr + (edge-1)*(d-1) + i] = ptr + (edgeperm[i]-1)*12 + 2*edge 
      perm[permptr + edge*(d-1) - i + 1] = ptr + (edgeperm[i]-1)*12 + 2*edge - 1
    end
  end
  ptr += 12*cub.numedge
  permptr += 6*(d-1)
  paramptr += cub.numedge
  # edge nodes are now completed

  for i = 1:cub.numfaceS21
    # set face-based S21 orbit nodes
    for face = 1:4
      for j = 1:3
        # the + fc after permptr is for the face centroid
        perm[permptr + fc + (face-1)*numface + 3*(i-1) + j] = ptr + j
      end
      ptr += 3
    end
    paramptr += 1
  end
  permptr += 4*numface

  # next, deal with 24-symmetries
  # ...

  # finally, deal with 1-symmetry (centroid)
  if cub.centroid
    perm[permptrint+1] = ptr+1
    permptrint += 1
    ptr += 1
  end

  # Next, find the permutation for the face nodes
  numbndrynodes = SymCubatures.getnumfacenodes(cub)
  faceperm = zeros(Int, (numbndrynodes))

  # TEMP: this must be corrected
  faceperm = [1:numbndrynodes;]

  return perm, faceperm
end