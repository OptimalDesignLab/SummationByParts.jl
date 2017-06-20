# This file gathers together a hodge-podge of functions that are not easily
# categorized

@doc """
### SummationByParts.getNumFaceNodes

Returns the number of SBP element nodes on a face.

**Inputs**

* `sbp`: an SBP operator

**Returns**

* `numfacenodes`: number of nodes on (one) face of the element

"""->
function getNumFaceNodes{T}(sbp::AbstractSBP{T})
  return SymCubatures.getnumfacenodes(sbp.cub)
end

@doc """
### SummationByParts.getnbrnodeindex

Returns the face-node index on `face.faceR` equivalent to index `i` on
`face.faceL`.

**Inputs**

* `sbp`: an SBP operator
* `face`: an element interface
* `i`: face-node index on `face.faceL`

**Returns**

* `j`: face-node index on `face.faceR`

"""->
function getnbrnodeindex{T}(sbp::TriSBP{T}, face::Interface, i::Int)
  return [2; 1; sbp.numfacenodes:-1:3][i]
end

@doc """
### SummationByParts.calcnodes

This function returns the node coordinates for an SBP operator.  It basically
calls calcnodes for the underlying SymCubature.  This function assumes the
element mapping is linear, i.e. edges are lines.

**Inputs**

* `sbp`: an SBP operator
* `vtx`: the vertices that define the element

**Outputs**

* `x`: the node coordinates; 1st dimension is the coordinate, the second the node

**Example**
```
  # define a third-order accurate SBP on triangles
  sbp = TriSBP{Float64}(degree=2)
  # build a simple 2-element grid on a square domain
  x = zeros(Float64, (2,sbp.numnodes,2))
  vtx = [0. 0.; 1. 0.; 0. 1.]
  x[:,:,1] = calcnodes(sbp, vtx)
  vtx = [1. 0.; 1. 1.; 0. 1.]
  x[:,:,2] = calcnodes(sbp, vtx)
```
"""->
function calcnodes{T}(sbp::AbstractSBP{T}, vtx::Array{T}=sbp.vtx) 
  return SymCubatures.calcnodes(sbp.cub, vtx)
end

@doc """
### SummationByParts.calcminnodedistance

Returns the minimum distance between distinct nodes on an element with straight sides

**Inputs**

* `sbp`: an SBP operator
* `vtx`: the vertices that define the element

**Returns**

* `mindist`: the minimum distance between distinct nodes

"""->
function calcminnodedistance{T}(sbp::AbstractSBP{T}, vtx::Array{T}=sbp.vtx)
  x = SymCubatures.calcnodes(sbp.cub, vtx)
  mindist = convert(T, Inf)
  for i = 1:size(x,2)
    for j = i+1:size(x,2)
      mindist = min(mindist, norm(x[:,i] - x[:,j]))
    end
  end
  return mindist
end

@doc """
### SummationByParts.buildinterpolation

Builds a matrix operator that can reconstruct a field located at the sbp nodes
to an auxlliary set of nodes.

**Inputs**

* `sbp`: an SBP operator
* `xinterp`: points to interpolate to in ref coords, size = [ndim,numpoints]
* `d=sbp.degree`: (optional) interpolation is exact for degree d polys

**Returns**

* `R`: the interpolation operator, size = [numpoints, sbp.numnodes]

"""->
function buildinterpolation{T}(sbp::TriSBP{T}, xinterp::AbstractArray{T,2};
                               d::Int=sbp.degree)
  # evaluate the basis at the SBP nodes and the interpolation points
  @assert size(xinterp, 1) == 2
  N = convert(Int, (d+1)*(d+2)/2 )
  Psbp = zeros(T, (sbp.numnodes,N) )
  Pinterp = zeros(T, (size(xinterp,2),N) )
  xsbp = calcnodes(sbp, sbp.vtx)
  ptr = 1
  for r = 0:d
    for j = 0:r
      i = r-j
      Psbp[:,ptr] = OrthoPoly.proriolpoly(vec(xsbp[1,:]), vec(xsbp[2,:]), i, j)
      Pinterp[:,ptr] = OrthoPoly.proriolpoly(vec(xinterp[1,:]),
                                             vec(xinterp[2,:]), i, j)
      ptr += 1
    end
  end
  R = (pinv(Psbp.')*Pinterp.').'
  return R
end

function buildinterpolation{T}(sbp::TetSBP{T}, xinterp::AbstractArray{T,2};
                               d::Int=sbp.degree)
  # evaluate the basis at the SBP nodes and the interpolation points
  @assert size(xinterp, 1) == 3
  N = convert(Int, (d+1)*(d+2)*(d+3)/6)
  Psbp = zeros(T, (sbp.numnodes,N) )
  Pinterp = zeros(T, (size(xinterp,2),N) )
  xsbp = calcnodes(sbp, sbp.vtx)
  ptr = 1
  for r = 0:d
    for k = 0:r
      for j = 0:r-k
        i = r-j-k
        Psbp[:,ptr] = OrthoPoly.proriolpoly(vec(xsbp[1,:]), vec(xsbp[2,:]),
                                            vec(xsbp[3,:]), i, j, k)
        Pinterp[:,ptr] = OrthoPoly.proriolpoly(vec(xinterp[1,:]),
                                               vec(xinterp[2,:]), 
                                               vec(xinterp[3,:]),
                                               i, j, k)
        ptr += 1
      end
    end
  end
  R = (pinv(Psbp.')*Pinterp.').'
  return R
end

@doc """
### SummationByParts.permuteinterface!

For a given array of values on a faces, permutes the node values (in place) to 
be in the orientation specified by the orient field of the corresponding 
Interface.  Methods are available for scalar and vector fields.

**Inputs**

* `sbp`: an SBPFace operator
* `ifaces`: an array of Interfaces

**In/Outs**

* `uface`: the array of face values, must have dimensions
           [n x sbpface.numnodes x length(ifaces] if 3D array (where n
           is arbitrary) or [sbpface.numnodes x length(ifaces] if 2D.  
           The permutation is applied to the second array dimension for the 
           3D case or the first dimension in the 2D case.
"""->
function permuteinterface!{Tsbp, Tsol}(sbpface::AbstractFace{Tsbp}, 
                                       ifaces::AbstractArray{Interface}, 
                                       uface::AbstractArray{Tsol, 3})
  @assert length(ifaces) == size(uface, 3)
  @assert sbpface.numnodes == size(uface, 2)

  dofpernode, numfacenodes, numfaces = size(uface)
  # temporary array needed during permutation
  workarr = Array(Tsol, dofpernode, numfacenodes)

  for iface =1:length(ifaces)
    orient = ifaces[iface].orient
    permvec = sview(sbpface.nbrperm, :, orient)
    facedata = sview(uface, :, :, iface)
    permuteface!(permvec, workarr, facedata)
  end

  return nothing
end

function permuteinterface!{Tsbp, Tsol}(sbpface::AbstractFace{Tsbp}, 
                                       ifaces::AbstractArray{Interface}, 
                                       uface::AbstractArray{Tsol, 2})
  @assert length(ifaces) == size(uface, 2)
  @assert sbpface.numnodes == size(uface, 1)

  numfacenodes, numfaces = size(uface)
  # temporary array needed during permutation
  workarr = Array(Tsol, numfacenodes)

  for iface =1:length(ifaces)
    orient = ifaces[iface].orient
    permvec = sview(sbpface.nbrperm, :, orient)
    facedata = sview(uface, :, iface)
    permuteface!(permvec, workarr, facedata)
  end

  return nothing
end

@doc """
### SummationByParts.permuteface!

This function applys a permutation to the data on a particular face.

**Inputs**

* `permvec`: vector specifying the permutation to apply
* `workarr`: a temporary array, same size as face_data, that is overwritten
             during the computation
**In/Outs**

* `face_data`: an N-D array containing the data to be pemuted, where the
              permutation is applied to the second dimension of the array
              for N=2 and the first dimension if N=1.  N > 2 is not
              currently supported
"""->
function permuteface!{Ti <: Integer, Tsol}(permvec::AbstractArray{Ti, 1},
                                           workarr::AbstractArray{Tsol, 2},
                                           facedata::AbstractArray{Tsol, 2})
  # copy to temporary array, applying permutation
  for i=1:size(facedata, 2)  # loop over nodes on the face
    idx = permvec[i]
    for j=1:size(facedata, 1)  # all dofs on the node
      workarr[j, idx] = facedata[j, i]
    end
  end

  # copy back, using linear indexing
  for i=1:length(facedata)
    facedata[i] = workarr[i]
  end

  return nothing
end

function permuteface!{Ti <: Integer, Tsol}(permvec::AbstractArray{Ti, 1},
                                           workarr::AbstractArray{Tsol, 1},
                                           facedata::AbstractArray{Tsol, 1})
  # copy to temporary array, applying permutation
  for i=1:size(facedata, 1)  # loop over nodes on the face
    idx = permvec[i]
    workarr[idx] = facedata[i]
  end

  # copy back, using linear indexing
  for i=1:length(facedata)
    facedata[i] = workarr[i]
  end

  return nothing
end

@doc """
### SummationByParts.basispursuit!

Finds an approximate solution to the underdetermined problem `Ax = b` that is
sparse using the alternating direction method of multipliers (ADMM).

**Inputs**

* `A`: matrix in the linear equation that must be satisfied
* `b`: vector in the linear equation that must be satisfied
* `rho` (optional) : augmented Lagrangian parameter
* `alpha` (optional) : over-relaxation parameter (usually between 1.0 and 1.8)
* `hist` (optional): output the convergence history to the terminal
* `abstol` (optional): absolute tolerance
* `reltol` (optional): relative tolerance

**In/Outs**

* `x`: sparse solution of the problem

**Notes**

This is a direct translation of Boyd et al.'s Matlab implementation of ADMM for
basis pursuit.  This method is not well suited to high-accuracy solutions, so,
for our purposes, it is best used as a means of identifying the sparsity
pattern.

"""->
function basispursuit!(A::AbstractArray{Float64,2}, b::AbstractVector{Float64},
                       x::AbstractVector{Float64}; rho::Float64=1.0,
                       alpha::Float64=1.0, hist::Bool=false,
                       abstol::Float64=1e-4, reltol::Float64=1e-2)
  @assert( size(A,1) == size(b,1) )
  @assert( size(A,2) == size(x,1) )
  @assert( rho > 0.0 )
  @assert( alpha >= 1.0 && alpha <= 1.8 )
  @assert( abstol > 0.0 && reltol > 0.0 )
  function objective(x)
    return norm(x,1)    
  end
  function shrinkage(a, kappa)
    return max(0, a-kappa) - max(0, -a-kappa)
  end

  maxiter = 10000
  (m, n) = size(A)
  fill!(x, 0.0)
  z = zeros(x)
  u = zeros(x)
  
  if hist
    @printf("%10s %10s %10s %10s %10s %10s\n", "iter", "r norm", "eps pri",
            "s norm", "eps dual", "objective")
  end

  # precompute some static variables for x-update (projection on to Ax=b)
  #AAt = A*A.'
  #P = eye(n,n) - A.'*(AAt\A) # eye(n) - V*V.' where V are the right sing. vecs
  #q = A.'*(AAt \ b)
  Ainv = pinv(A)
  P = eye(n,n) - Ainv*A
  q = Ainv*b

  for k = 1:maxiter
    # x-update
    x[:] = P*(z - u) + q

    # z-update
    zold = z
    x_hat = alpha*x + (1 - alpha)*zold
    z = shrinkage(x_hat + u, 1/rho)

    u += (x_hat - z)
    
    # diagnostics, reporting, termination checks
    objval = objective(x)
    r_norm = norm(x - z)
    s_norm = norm(-rho*(z - zold))
    
    eps_pri = sqrt(n)*abstol + reltol*max(norm(x), norm(-z))
    eps_dual = sqrt(n)*abstol + reltol*norm(rho*u)

    if hist
      @printf("%10d %10f %10f %10f %10f %10f\n", k, r_norm, eps_pri, s_norm, eps_dual,
              objval)
    end
    if r_norm < eps_pri && s_norm < eps_dual
      break
    end
  end
end

@doc """
### SummationByParts.calcSparseSolution!

Finds a solution to the underdetermined problem `Ax = b` that is sparse.  Uses
basispursit! to find an approximate solution, which is used to eliminate columns
from `A`, from which an accurate solution is found.

**Inputs**

* `A`: matrix in the linear equation that must be satisfied
* `b`: vector in the linear equation that must be satisfied

**In/Outs**

* `x`: sparse solution of the problem

"""->
function calcSparseSolution!(A::AbstractArray{Float64,2},
                             b::AbstractVector{Float64},
                             x::AbstractVector{Float64})
  @assert( size(A,1) == size(b,1) )
  @assert( size(A,2) == size(x,1) )
  # find an approximate solution
  basispursuit!(A, b, x, rho=1.5, alpha=1.0, hist=false, abstol=1e-6,
                reltol=1e-6)
  # use the approximate solution to identify the non-zero entries
  rankA = rank(A)
  P = zeros(size(A,2),rankA)
  idx = sortperm(abs(x), rev=true)
  for i = 1:rankA
    P[idx[i],i] = 1.0
  end
  # find the reduced (full-rank) matrix and invert
  AP = A*P
  #println("size(AP) = ",size(AP))
  #println("rank(AP) = ",rank(AP))
  #x[:] = P*(AP\b)
  x[:] = P*(pinv(AP)*b)
end

function compareEigs(x::Float64, y::Float64)
  if isless(x,y)
    return true
  else
    return false
  end
end

function compareEigs(x::Complex128, y::Complex128)
  if abs(abs(x) - abs(y)) < 10*eps(Float64)
    # moduli are close, sort by imaginary part
    if imag(x) < imag(y)
      return true
    else
      return false
    end
  elseif abs(x) < abs(y)
    return true
  else
    return false
  end
end
  
@doc """
### SummationByParts.calcMatrixEigs!

Finds the eigenvalues of the given matrix in order of increasing modulus; that
is, this works for symmetric and non-symmetric square matrices.  This method is
basically a front end for `eigfact`.

**Inputs**

* `A`: matrix whose eigenvalues are desired

**In/Outs**

* `λ`: eigenvalues of `A` sorted in increasing modulus

"""->
function calcMatrixEigs!{T}(A::AbstractArray{T,2},
                            λ::AbstractVector{T})
  @assert( size(A,1) == size(A,2) )
  @assert( length(λ) == size(A,1) )
  n = size(A,1)  
  fac = eigfact(A)
  idx = zeros(Int,n)
  sortperm!(idx, fac[:values], lt=compareEigs)
  for i = 1:n
    λ[i] = fac[:values][idx[i]]
  end
end

@doc """
### SummationByParts.calcMatrixEigs_rev!

The reverse-mode differentiated version of `calcMatrixEigs!`.  The math behind
this is from Mike Giles report "An extended collection of matrix derivative
results for forward and reverse mode algorithmic differentiation."

**Inputs**

* `A`: matrix whose eigenvalues are desired
* `λ_bar`: ∂f/∂λ vector that multiplies derivatives from left

**In/Outs**

* `λ`: eigenvalues of `A` sorted in increasing modulus
* `A_bar`: derivatives ∂f/∂A

"""->
function calcMatrixEigs_rev!{T}(A::AbstractArray{T,2},
                                λ::AbstractVector{T},
                                λ_bar::AbstractVector{T},
                                A_bar::AbstractArray{T,2})
  @assert( size(A,1) == size(A,2) == size(A_bar,1) == size(A_bar,2) )
  @assert( length(λ) == length(λ_bar) == size(A,1) )
  n = size(A,1)  
  fac = eigfact(A)
  idx = zeros(Int,n)
  sortperm!(idx, fac[:values], lt=compareEigs)
  for i = 1:n
    λ[i] = fac[:values][idx[i]]
  end
  A_bar[:,:] = (fac[:vectors][:,idx]')\(diagm(λ_bar)*fac[:vectors][:,idx]') 
end

@doc """
### SummationByParts.conditionObj

Let `x` = `Znull`*`xred` + `xperp` be the (unique) entries in a skew symmetric
matrix `S`, and let `E` be a symmetric matrix.  This routine computes an
approximation of the condition number of the matrix `A` = (`S` + |`E`|), where
|⋅| is the elementwise absolute value.  The matrix `A` corresponds to a
weak-form discretization of linear advection.  The condition number is
approximated using KS aggregation to ensure the objective is differentiable.

**Inputs**

* `xred`: a reduced-space for the entries in the skew-symmetric matrix
* `p`: defines the KS parameter; as p tends to infinity, we get obj = kappa(A)
* `xperp`: a particular solution that satisfies the SBP accuracy conditions
* `Znull`: matrix that defines the null-space of the SBP accuracy conditions
* `E`: symmetric matrix, usually the boundary operator for an SBP matrix

**Returns**

* `obj`: approximate condition number of `A` as defined above.

"""->
function conditionObj(xred::AbstractVector{Float64}, p,
                      xperp::AbstractVector{Float64},
                      Znull::AbstractArray{Float64,2},
                      E::AbstractArray{Float64,2})
  @assert( length(xperp) == size(Znull,1) )
  @assert( length(xred) == size(Znull,2) )
  @assert( size(E,1) == size(E,2) )
  @assert( p >= 1.0 )
  n = size(E,1)
  x = Znull*xred + xperp
  # insert into (S + 1/2 |E|)
  A = zeros(n,n)
  for i = 1:n
    for j = 1:n
      A[i,j] = abs(E[i,j])
    end
  end
  for i = 2:n
    offset = convert(Int, (i-1)*(i-2)/2)
    for j = 1:i-1
      A[i,j] += x[offset+j]
      A[j,i] -= x[offset+j]
    end
  end

  # compute the SVD of A
  U, S, V = svd(A)
  maxKS = 0.0
  minKS = 0.0
  for i = 1:n
    maxKS += exp(p*(S[i] - S[1]))
    minKS += exp(p*(1/S[i] - 1/S[end]))
  end
  return (S[1] + (1/p)*log(maxKS/n))*(1/S[end] + (1/p)*log(minKS/n))
end

@doc """
### SummationByParts.eigenvalueObjGrad!

Computes the gradient of the function `conditionObj` with respect to `xred`,
and returns it in the array `g`.

**Inputs**

* `xred`: a reduced-space for the entries in the skew-symmetric matrix
* `p`: defines the KS parameter; as p tends to infinity, we get obj = kappa(A)
* `xperp`: a particular solution that satisfies the SBP accuracy conditions
* `Znull`: matrix that defines the null-space of the SBP accuracy conditions
* `E`: symmetric matrix, usually the boundary operator for an SBP matrix

**In/Outs**

* `g`: gradient of the objective `conditionObj` with respect to `xred`

"""->
function conditionObjGrad!(xred::AbstractVector{Float64}, p,
                           xperp::AbstractVector{Float64},
                           Znull::AbstractArray{Float64,2},
                           E::AbstractArray{Float64,2},
                           g::AbstractVector{Float64})
  @assert( length(xperp) == size(Znull,1) )
  @assert( length(xred) == size(Znull,2) )
  @assert( size(E,1) == size(E,2) )
  @assert( length(g) == length(xred) )
  @assert( p >= 1.0 )
  n = size(E,1)
  x = zeros(xperp)
  x = Znull*xred + xperp
  # insert into (S + 1/2 |E|)
  A = zeros(n,n)
  for i = 1:n
    for j = 1:n
      A[i,j] = abs(E[i,j])
    end
  end
  for i = 2:n
    offset = convert(Int, (i-1)*(i-2)/2)
    for j = 1:i-1
      A[i,j] += x[offset+j]
      A[j,i] -= x[offset+j]
    end
  end

  # partial forward sweep; get the SVD of A and some intermediate vars
  U, S, V = svd(A)
  maxKS = 0.0
  minKS = 0.0
  for i = 1:n
    maxKS += exp(p*(S[i] - S[1]))
    minKS += exp(p*(1/S[i] - 1/S[end]))
  end

  # start reverse sweep
  S_bar = zeros(S)
  #return (S[1] + (1/p)*log(maxKS/n))*(1/S[end] + (1/p)*log(minKS/n))
  maxKS_bar = (1/S[end] + (1/p)*log(minKS/n))/(p*maxKS)
  minKS_bar = (S[1] + (1/p)*log(maxKS/n))/(p*minKS)
  S_bar[1] = (1/S[end] + (1/p)*log(minKS/n))
  S_bar[end] = -(S[1] + (1/p)*log(maxKS/n))/(S[end]^2)
  for i = 1:n
    # maxKS += exp(p*(S[i] - S[1]))
    S_bar[i] += maxKS_bar*exp(p*(S[i] - S[1]))*p
    S_bar[1] -= maxKS_bar*exp(p*(S[i] - S[1]))*p
    # minKS += exp(p*(1/S[i] - 1/S[end]))
    S_bar[i] -= minKS_bar*exp(p*(1/S[i] - 1/S[end]))*p/(S[i]^2)
    S_bar[end] += minKS_bar*exp(p*(1/S[i] - 1/S[end]))*p/(S[end]^2)
  end
  A_bar = U*diagm(S_bar)*V'

  x_bar = zeros(x)
  for i = 2:n
    offset = convert(Int, (i-1)*(i-2)/2)
    for j = 1:i-1
      # A[i,j] = +complex(x[offset+j], 0.0)
      # A[j,i] = -complex(x[offset+j], 0.0)
      x_bar[offset+j] += A_bar[i,j] - A_bar[j,i]
    end
  end
  # x = Znull*xred + xperp
  g[:] = Znull.'*x_bar
end

@doc """
### SummationByParts.eigenvalueObj

Let `x` = `Znull`*`xred` + `xperp` be the (unique) entries in a skew symmetric
matrix `S`, and let `E` be a symmetric matrix.  This routine returns the
spectral radius of the matrix `A` = diagm(1./`w`)*(`S` + |`E`|), where |⋅| is
the elementwise absolute value.  The matrix `A` corresponds to a strong-form
discretization of linear advection.

**Inputs**

* `xred`: a reduced-space for the entries in the skew-symmetric matrix
* `p`: not used presently
* `xperp`: a particular solution that satisfies the SBP accuracy conditions
* `Znull`: matrix that defines the null-space of the SBP accuracy conditions
* `w`: diagonal norm entries in an SBP operator
* `E`: symmetric matrix, usually the boundary operator for an SBP matrix

**Returns**

* `obj`: the 2p-norm of the moduli of the eigenvalues of `A` as defined above.

"""->
function eigenvalueObj(xred::AbstractVector{Float64}, p::Int,
                       xperp::AbstractVector{Float64},
                       Znull::AbstractArray{Float64,2},
                       w::AbstractVector{Float64},
                       E::AbstractArray{Float64,2})
  @assert( length(xperp) == size(Znull,1) )
  @assert( length(xred) == size(Znull,2) )
  @assert( size(E,1) == size(E,2) == length(w) )
  @assert( p >= 1 )
  n = size(E,1)
  x = Znull*xred + xperp
  # insert into H^-1(S + 1/2 |E|)
  A = zeros(Complex128, (n,n))
  for i = 1:n
    for j = 1:n
      A[i,j] = complex(abs(E[i,j]), 0.0)
    end
  end
  for i = 2:n
    offset = convert(Int, (i-1)*(i-2)/2)
    for j = 1:i-1
      A[i,j] += complex(x[offset+j], 0.0)
      A[j,i] -= complex(x[offset+j], 0.0)
    end
  end
  for i = 1:n
    fac = 1./w[i]
    for j = 1:n
      A[i,j] *= fac
    end
  end
  # compute the (sorted) eigenvalues of A, and the objective
  λ = zeros(Complex128, n)
  calcMatrixEigs!(A, λ)
  return abs(λ[end])
end

@doc """
### SummationByParts.eigenvalueObjGrad!

Computes the gradient of the function `eigenvalueObj` with respect to `xred`,
and returns it in the array `g`.

**Inputs**

* `xred`: a reduced-space for the entries in the skew-symmetric matrix
* `p`: not used at present
* `xperp`: a particular solution that satisfies the SBP accuracy conditions
* `Znull`: matrix that defines the null-space of the SBP accuracy conditions
* `w`: diagonal norm entries in an SBP operator
* `E`: symmetric matrix, usually the boundary operator for an SBP matrix

**In/Outs**

* `g`: gradient of the objective `eigenvalueObj` with respect to `xred`

"""->
function eigenvalueObjGrad!(xred::AbstractVector{Float64}, p::Int,
                            xperp::AbstractVector{Float64},
                            Znull::AbstractArray{Float64,2},
                            w::AbstractVector{Float64},
                            E::AbstractArray{Float64,2},
                            g::AbstractVector{Float64})
  @assert( length(xperp) == size(Znull,1) )
  @assert( length(xred) == size(Znull,2) )
  @assert( size(E,1) == size(E,2) == length(w) )
  @assert( length(g) == length(xred) )
  @assert( p >= 1 )
  n = size(E,1)
  x = zeros(xperp)
  x = Znull*xred + xperp
  # insert into H^-1(S + 1/2 |E|)
  A = zeros(Complex128, (n,n))
  for i = 1:n
    for j = 1:n
      A[i,j] = complex(abs(E[i,j]), 0.0)
    end
  end
  for i = 2:n
    offset = convert(Int, (i-1)*(i-2)/2)
    for j = 1:i-1
      A[i,j] += complex(x[offset+j], 0.0)
      A[j,i] -= complex(x[offset+j], 0.0)
    end
  end
  for i = 1:n
    fac = 1./w[i]
    for j = 1:n
      A[i,j] *= fac
    end
  end
  # compute the (sorted) eigenvalues of A and the objective
  λ = zeros(Complex128, n)
  calcMatrixEigs!(A, λ)

  λ_bar = zeros(λ)
  λ_bar[end] = λ[end]/abs(λ[end])  
  A_bar = zeros(A)
  calcMatrixEigs_rev!(A, λ, λ_bar, A_bar)
  for i = 1:n
    fac = 1./w[i]
    for j = 1:n
      A_bar[i,j] *= fac
    end
  end
  # place A_bar into x_bar
  x_bar = zeros(x)
  for i = 2:n
    offset = convert(Int, (i-1)*(i-2)/2)
    for j = 1:i-1
      # A[i,j] = +complex(x[offset+j], 0.0)
      # A[j,i] = -complex(x[offset+j], 0.0)
      x_bar[offset+j] += real(A_bar[i,j]) - real(A_bar[j,i])
    end
  end
  # x = Znull*xred + xperp
  g[:] = Znull.'*x_bar
end
