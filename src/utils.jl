# This file gathers together a hodge-podge of functions that are not easily
# categorized

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
  if sbp.reorder
    perm, faceperm = SummationByParts.getnodepermutation(sbp.cub, sbp.degree)
    x = SymCubatures.calcnodes(sbp.cub, vtx)
    return x[:,perm]
  else
    return SymCubatures.calcnodes(sbp.cub, vtx)
  end
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

**Returns**

* `R`: the interpolation operator, size = [numpoints, sbp.numnodes]

"""->
function buildinterpolation{T}(sbp::TriSBP{T}, xinterp::AbstractArray{T,2})
  # evaluate the basis at the SBP nodes and the interpolation points
  @assert size(xinterp, 1) == 2
  d = sbp.degree
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

function buildinterpolation{T}(sbp::TetSBP{T}, xinterp::AbstractArray{T,2})
  # evaluate the basis at the SBP nodes and the interpolation points
  @assert size(xinterp, 1) == 3

  d = sbp.degree
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

Finds a solution to the underdetermined problem `Ax = b` that is sparse using
the alternating direction method of multipliers (ADMM).

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

  maxiter = 1000
  (m, n) = size(A)
  fill!(x, 0.0)
  z = zeros(x)
  u = zeros(x)
  
  if hist
    @printf("%10s %10s %10s %10s %10s %10s\n", "iter", "r norm", "eps pri",
            "s norm", "eps dual", "objective")
  end

  # precompute some static variables for x-update (projection on to Ax=b)
  AAt = A*A.'
  P = eye(n,n) - A.'*(AAt\A)
  q = A.'*(AAt \ b)
  #Ainv = pinv(A)
  #P = eye(n,n) - Ainv*Ainv.'
  #q = Ainv*b

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

