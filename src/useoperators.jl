# This file gathers together functions related to using the SBP operators

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
### SummationByParts.weakdifferentiate!

Applies the SBP stiffness matrix (or its transpose) to data in `flux` and
**adds** the result to `res`.  Different methods are available depending on the
rank of `flux`:

* For *scalar* fields, it is assumed that `flux` is a rank-2 array, with the
first dimension for the local-node index, and the second dimension for the
element index.
* For *vector* fields, `flux` is a rank-3 array, with the first dimension for
the index of the vector field, the second dimension for the local-node index,
and the third dimension for the element index.

Naturally, the number of entries in the dimension of `flux` (and `res`)
corresponding to the nodes must be equal to the number of nodes in the SBP
operator sbp.

**Inputs**

* `sbp`: an SBP operator type
* `di`: direction index of the operator that is desired (di=1 for Qx, etc)
* `flux`: the array that the operator is applied to
* `±` : PlusFunctor to add to res, MinusFunctor to subract
* `trans` (optional): if true, the transpose operation is applied

**In/Outs**

* `res`: where the result of applying Q[:,:,di] to u is stored

"""->
function weakdifferentiate!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int, 
                                            flux::AbstractArray{Tflx,2},
                                            res::AbstractArray{Tres,2},
                                            (±)::UnaryFunctor=Add();
                                            trans::Bool=false)
  @assert( sbp.numnodes == size(flux,1) && sbp.numnodes == size(res,1) )
  @assert( length(flux) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  if trans # apply transposed Q
    for elem = 1:size(flux,2)
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          res[i,elem] += ±(sbp.Q[j,i,di]*flux[j,elem])
        end
      end
    end
  else # apply Q
    for elem = 1:size(flux,2)
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          res[i,elem] += ±(sbp.Q[i,j,di]*flux[j,elem])
        end
      end
    end
  end
end

function weakdifferentiate!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                                            flux::AbstractArray{Tflx,3},
                                            res::AbstractArray{Tres,3},
                                            (±)::UnaryFunctor=Add();
                                            trans::Bool=false)
  @assert( sbp.numnodes == size(flux,2) && sbp.numnodes == size(res,2) )
  @assert( length(flux) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  if trans # apply transposed Q
    for elem = 1:size(flux,3)
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          for field = 1:size(flux,1)
            res[field,i,elem] += ±(sbp.Q[j,i,di]*flux[field,j,elem])
          end
        end
      end
    end
  else # apply Q
    for elem = 1:size(flux,3)
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          for field = 1:size(flux,1)
            res[field,i,elem] += ±(sbp.Q[i,j,di]*flux[field,j,elem])
          end
        end
      end
    end
  end
end

@doc """
### SummationByParts.differentiate!

Applies the SBP differentiation matrix operator, D, to data in `flux` and
**adds** the result to `res`.  Different methods are available depending on the
rank of `flux`:

* For *scalar* fields, it is assumed that `flux` is a rank-2 array, with the
first dimension for the local-node index, and the second dimension for the
element index.
* For *vector* fields, `flux` is a rank-3 array, with the first dimension for
the index of the vector field, the second dimension for the local-node index,
and the third dimension for the element index.

Naturally, the number of entries in the dimension of `flux` (and `res`)
corresponding to the nodes must be equal to the number of nodes in the SBP
operator sbp.

**Inputs**

* `sbp`: an SBP operator type
* `di`: direction index of the operator that is desired (di=1 for Dx, etc)
* `flux`: the array that the operator is applied to
* `±`: PlusFunctor to add to res, MinusFunctor to subract

**In/Outs**

* `res`: where the result of applying inv(H)*Q[:,:,di] to u is stored

"""->
function differentiate!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                                        flux::AbstractArray{Tflx,2},
                                        res::AbstractArray{Tres,2},
                                        (±)::UnaryFunctor=Add())
  @assert( sbp.numnodes == size(flux,1) && sbp.numnodes == size(res,1) )
  @assert( length(flux) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  Hinv = 1./sbp.w
  for elem = 1:size(flux,2)
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        res[i,elem] += ±(sbp.Q[i,j,di]*flux[j,elem])
      end
      res[i,elem] *= Hinv[i]
    end
  end
end

function differentiate!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                                        flux::AbstractArray{Tflx,3},
                                        res::AbstractArray{Tres,3},
                                        (±)::UnaryFunctor=Add())
  @assert( sbp.numnodes == size(flux,2) && sbp.numnodes == size(res,2) )
  @assert( length(flux) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  Hinv = 1./sbp.w
  for elem = 1:size(flux,3)
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        for field = 1:size(flux,1)
          res[field,i,elem] += ±(sbp.Q[i,j,di]*flux[field,j,elem])
        end
      end
      for field = 1:size(flux,1)
        res[field,i,elem] *= Hinv[i]
      end
    end
  end
end

@doc """
### SummationByParts.directionaldifferentiate!

Performs a directional derivative (in reference space) at a given node.  In
constrast with other differentiate routines, the input field `u` is for **a
single element**, not a collection of elements.

**WARNING**: In the case of a vector field u, the directional derivative is
  added to the output Ddir; the user must zero this before.

**Inputs**

* `sbp`: an SBP operator type
* `dir`: a direction vector for the directional derivative
* `u`: the field that is being differentiated (either a scalar or vector)
* `i`: index of the node at which the derivative is desired

**Returns or In/Outs**

* `Ddir`: derivative of `u` in direction `dir`

"""->
function directionaldifferentiate!{Tsbp,Tmsh,Tsol}(sbp::AbstractSBP{Tsbp},
                                                   dir::Array{Tmsh,1},
                                                   u::AbstractArray{Tsol,1},
                                                   i::Int)
  @assert( size(sbp.Q, 3) == size(dir,1) )
  Ddir = zero(Tmsh)*zero(Tsol)
  for di = 1:size(sbp.Q, 3)
    tmp = zero(Tmsh)*zero(Tsol)
    for j = 1:sbp.numnodes
      tmp += sbp.Q[i,j,di]*u[j]
    end
    Ddir += dir[di]*tmp/sbp.w[i]
  end
  return Ddir
end

function directionaldifferentiate!{Tsbp,Tmsh,Tsol,Tres}(sbp::AbstractSBP{Tsbp},
                                                        dir::Array{Tmsh,1}, 
                                                        u::AbstractArray{Tsol,2},
                                                        i::Int,
                                                        Ddir::Array{Tres,1})
  @assert( size(sbp.Q, 3) == size(dir,1) )
  for di = 1:size(sbp.Q, 3)
    tmp = zeros(Ddir)
    for j = 1:sbp.numnodes
      for field = 1:size(u,1)
        tmp[field] += sbp.Q[i,j,di]*u[field,j]
      end
    end
    for field = 1:size(u,1)
      Ddir[field] += dir[di]*tmp[field]/sbp.w[i]
    end
  end
end

@doc """
### SummationByParts.volumeintegrate!

Applies the SBP mass matrix operator, H, to data in `u` and stores
the result in `res`.  Different methods are available depending on the rank of
`u`:

* For *scalar* fields, it is assumed that `u` is a rank-2 array, with the first
dimension for the local-node index, and the second dimension for the element
index.
* For *vector* fields, `u` is a rank-3 array, with the first dimension for the
index of the vector field, the second dimension for the local-node index, and
the third dimension for the element index.

Naturally, the number of entries in the dimension of `u` (and `res`)
corresponding to the nodes must be equal to the number of nodes in the SBP
operator sbp.

**Inputs**

* `sbp`: an SBP operator type
* `u`: the array that the operator is applied to
* `±`: PlusFunctor to add to res, MinusFunctor to subract

**In/Outs**

* `res`: where the result of applying H to u is stored

"""->
function volumeintegrate!{Tsbp,Tsol,Tres}(sbp::AbstractSBP{Tsbp},
                                          u::AbstractArray{Tsol,2},
                                          res::AbstractArray{Tres,2},
                                          (±)::UnaryFunctor=Add())
  @assert( sbp.numnodes == size(u,1) && sbp.numnodes == size(res,1) )
  @assert( length(u) == length(res) )
  H = sbp.w
  for elem = 1:size(u,2)
    for i = 1:sbp.numnodes
      res[i,elem] += ±(H[i]*u[i,elem])
    end
  end
end

function volumeintegrate!{Tsbp,Tsol,Tres}(sbp::AbstractSBP{Tsbp},
                                          u::AbstractArray{Tsol,3},
                                          res::AbstractArray{Tres,3},
                                          (±)::UnaryFunctor=Add())
  @assert( sbp.numnodes == size(u,2) && sbp.numnodes == size(res,2) )
  @assert( length(u) == length(res) )
  H = sbp.w
  for elem = 1:size(u,3)
    for i = 1:sbp.numnodes
      for field = 1:size(u,1)
        res[field,i,elem] += ±(H[i]*u[field,i,elem])
      end
    end
  end
end

@doc """
### SummationByParts.mappingjacobian!

Evaluates the (scaled) Jacobian of the mapping from reference coordinates to
physical coordinates, as well as the determinant of the Jacobian.  The values
returned in dξdx are scaled by the determinant, so they have the same units as
the boundary measure (i.e. length in 2D, or length^2 in 3D).  This scaling is
adopted, because conservation laws written in conservative form in the reference
frame use the scaled Jacobian.

**Inputs**

* `sbp`: an SBP operator type
* `x`: the physical coordinates; 1st dim = coord, 2nd dim = node, 3rd dim = elem

**In/Outs**

* `dξdx`: the scaled Jacobian of the mapping; 1st dim = ref coord, 2nd dim =
  phys coord, 3rd dim = node, 3rd dim = elem
* `jac`: the determinant of the Jacobian; 1st dim = node, 2nd dim = elem

"""->
function mappingjacobian!{Tsbp,Tmsh}(sbp::TriSBP{Tsbp},
                                     x::AbstractArray{Tmsh,3},
                                     dξdx::AbstractArray{Tmsh,4},
                                     jac::AbstractArray{Tmsh,2})
  @assert( sbp.numnodes == size(x,2) && sbp.numnodes == size(dξdx,3) )
  @assert( size(x,3) == size(dξdx,4) )
  @assert( size(x,1) == 2 && size(dξdx,1) == 2 && size(dξdx,2) == 2 )
  fill!(dξdx, zero(Tmsh))
  dxdξ = zeros(Tmsh, (2,sbp.numnodes,size(x,3)))
  # compute d(x,y)/dxi and set deta/dx and deta/dy
  differentiate!(sbp, 1, x, dxdξ)
  for elem = 1:size(x,3)
    for i = 1:size(x,2)
      dξdx[2,1,i,elem] = -dxdξ[2,i,elem]
      dξdx[2,2,i,elem] = dxdξ[1,i,elem]
    end
  end
  # compute d(x,y)/deta and set dxi/dx and dxi/dy
  fill!(dxdξ, zero(Tmsh))
  differentiate!(sbp, 2, x, dxdξ)
  for elem = 1:size(x,3)
    for i = 1:size(x,2)
      dξdx[1,2,i,elem] = -dxdξ[1,i,elem]
      dξdx[1,1,i,elem] = dxdξ[2,i,elem]
    end
  end
  # compute the determinant of the Jacobian
  for elem = 1:size(x,3)
    for i = 1:sbp.numnodes
      jac[i,elem] = one(Tmsh)/(dξdx[1,1,i,elem]*dξdx[2,2,i,elem] - 
                               dξdx[1,2,i,elem]*dξdx[2,1,i,elem])
    end
  end
  # check for negative jac here?
end

function mappingjacobian!{Tsbp,Tmsh}(sbp::SparseTriSBP{Tsbp},
                                     x::AbstractArray{Tmsh,3},
                                     dξdx::AbstractArray{Tmsh,4},
                                     jac::AbstractArray{Tmsh,2})
  @assert( sbp.numnodes == size(x,2) && sbp.numnodes == size(dξdx,3) )
  @assert( size(x,3) == size(dξdx,4) )
  @assert( size(x,1) == 2 && size(dξdx,1) == 2 && size(dξdx,2) == 2 )
  fill!(dξdx, zero(Tmsh))
  dxdξ = zeros(Tmsh, (2,sbp.numnodes,size(x,3)))
  # compute d(x,y)/dxi and set deta/dx and deta/dy
  differentiate!(sbp, 1, x, dxdξ)
  for elem = 1:size(x,3)
    for i = 1:size(x,2)
      dξdx[2,1,i,elem] = -dxdξ[2,i,elem]
      dξdx[2,2,i,elem] = dxdξ[1,i,elem]
    end
  end
  # compute d(x,y)/deta and set dxi/dx and dxi/dy
  fill!(dxdξ, zero(Tmsh))
  differentiate!(sbp, 2, x, dxdξ)
  for elem = 1:size(x,3)
    for i = 1:size(x,2)
      dξdx[1,2,i,elem] = -dxdξ[1,i,elem]
      dξdx[1,1,i,elem] = dxdξ[2,i,elem]
    end
  end
  # compute the determinant of the Jacobian
  for elem = 1:size(x,3)
    for i = 1:sbp.numnodes
      jac[i,elem] = one(Tmsh)/(dξdx[1,1,i,elem]*dξdx[2,2,i,elem] - 
                               dξdx[1,2,i,elem]*dξdx[2,1,i,elem])
    end
  end
  # check for negative jac here?
end

function mappingjacobian!{Tsbp,Tmsh}(sbp::TetSBP{Tsbp},
                                     x::AbstractArray{Tmsh,3},
                                     dξdx::AbstractArray{Tmsh,4},
                                     jac::AbstractArray{Tmsh,2})
  @assert( sbp.numnodes == size(x,2) && sbp.numnodes == size(dξdx,3) )
  @assert( size(x,3) == size(dξdx,4) )
  @assert( size(x,1) == 3 && size(dξdx,1) == 3 && size(dξdx,2) == 3 )
  numelem = size(x,3)
  dxdξ = zeros(Tmsh, (3,sbp.numnodes,numelem,3))
  # calculate the derivative of the coordinates with respect to (xi,eta,zeta)
  # using the SBP operator
  for di = 1:3
    differentiate!(sbp, di, x, sub(dxdξ,:,:,:,di)) 
  end
  fill!(dξdx, zero(Tmsh))
  # calculate the metrics: the outer loop calculates the derivatives of
  # coordinates-times-coordinate-derivatives, for example, d/deta( z * d(y)/d
  # zeta), and this contribution is added to the relevant metric.
  for di = 1:3
    it1 = mod(di,3)+1
    it2 = mod(di+1,3)+1
    for elem = 1:numelem
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          # the 0.5 factor in the SBP entry accounts for the averaging over the
          # two possible analytical formulations
          coeff = 0.5*sbp.Q[i,j,di]
          # this second set of indices (di2, it21, it22) denote the x-, y-,
          # z-coordinate indices, while the first set (di, it1, it2) denotes the xi-,
          # eta-, zeta-indices
          for di2 = 1:3
            it21 = mod(di2,3)+1
            it22 = mod(di2+1,3)+1
            dξdx[it1,di2,i,elem] += 
            coeff*(x[it22,j,elem]*dxdξ[it21,j,elem,it2] -
                   x[it21,j,elem]*dxdξ[it22,j,elem,it2])            
            dξdx[it2,di2,i,elem] += 
            coeff*(x[it21,j,elem]*dxdξ[it22,j,elem,it1] -
                   x[it22,j,elem]*dxdξ[it21,j,elem,it1])
          end
        end
      end
    end
  end
  # scale metrics by norm and calculate the determinant of the mapping
  for elem = 1:numelem
    for i = 1:sbp.numnodes
      dξdx[:,:,i,elem] /= sbp.w[i]
      jac[i,elem] = one(Tmsh)/(dxdξ[1,i,elem,1]*dxdξ[2,i,elem,2]*dxdξ[3,i,elem,3] +
                               dxdξ[1,i,elem,2]*dxdξ[2,i,elem,3]*dxdξ[3,i,elem,1] +
                               dxdξ[1,i,elem,3]*dxdξ[2,i,elem,1]*dxdξ[3,i,elem,2] -
                               dxdξ[1,i,elem,1]*dxdξ[2,i,elem,3]*dxdξ[3,i,elem,2] -
                               dxdξ[1,i,elem,2]*dxdξ[2,i,elem,1]*dxdξ[3,i,elem,3] -
                               dxdξ[1,i,elem,3]*dxdξ[2,i,elem,2]*dxdξ[3,i,elem,1])
    end
  end
  # check for negative jac here?
end

@doc """
### SummationByParts.calcmappingjacobian!

Uses a given set of Lagrangian element nodes to determine an analytical
(polynomial) mapping, and then uses this mapping to determine the Jacobian of
the mapping.  The approach varies depending on the dimension of the problem:

* For *2-dimensional problems* the exact Jacobian of the mapping is used;
* For *3-dimensional problems* an optimization problem is solved

**Inputs**

* `sbp`: an SBP operator type
* `mapdegree`: the polynomial degree of the mapping
* `xlag`: Lagrangian nodes in physical space; [coord, Lagrangian node]
* `xref`: Lagrangian nodes in reference space; [coord, Lagrangian node]
* `Eone`: Ex*one, Ey*one (Ez*one); [sbp node, coord] see notes below

**In/Outs**

* `xsbp`: location of the SBP nodes in physical space; [coord, sbp node]
* `dξdx`: the scaled Jacobian of the mapping; [ref coord, phys coord, sbp node]
* `jac`: the determinant of the Jacobian; [sbp node]

**Notes**

The array `Eone` is the product of the boundary operators, in physical space,
with the vector of ones (see Crean et al., Entropy-Conservative,
Multidimensional Summation-By-Parts Discretization of the Euler Equations, as
well as the test in `test/test_useoperators.jl`).  These products are used to
define the metric invariants.  *They are only needed by the 2-dimensional code,
and so this array can be passed empty in that case*.

"""->
function calcmappingjacobian!{Tsbp,Tmsh}(sbp::TriSBP{Tsbp}, mapdegree::Int,
                                         xlag::AbstractArray{Tmsh,2},
                                         xref::AbstractArray{Tmsh,2},
                                         xsbp::AbstractArray{Tmsh,2},
                                         dξdx::AbstractArray{Tmsh,3},
                                         jac::AbstractArray{Tmsh},
                                         Eone::AbstractArray{Tmsh,2}=
                                           Array(Tmsh,0,0))
  @assert( sbp.numnodes == size(xsbp,2) == size(dξdx,3) == size(jac,1) )
  @assert( size(xlag,1) == size(xref,1) == size(xsbp,1) == size(dξdx,1) 
           == size(dξdx,2) == 2 )
  numdof = binomial(mapdegree+2,2)
  @assert( size(xlag,2) == size(xref,2) == numdof )
  # Step 1: find the polynomial mapping using xlag
  V = zeros(Tmsh, (numdof,numdof) )
  ptr = 1
  for r = 0:mapdegree
    for j = 0:r
      i = r-j
      V[:,ptr] = OrthoPoly.proriolpoly(vec(xref[1,:]), vec(xref[2,:]), i, j)
      ptr += 1
    end
  end
  coeff = zeros(Tmsh, (numdof,2))
  coeff = V\(xlag.')
  # Step 2: compute the mapped SBP nodes and the analytical Jacobian at sbp nodes
  x = calcnodes(sbp) # <-- SBP nodes in reference space
  dxdξ = zeros(Tmsh, (2,2,sbp.numnodes))
  fill!(xsbp, zero(Tmsh))
  ptr = 1
  for r = 0:mapdegree
    for j = 0:r
      i = r-j
      P = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      dPdξ, dPdη = OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      for di = 1:2
        for nd = 1:sbp.numnodes
          xsbp[di,nd] += coeff[ptr,di]*P[nd]
          dxdξ[di,1,nd] += coeff[ptr,di]*dPdξ[nd]
          dxdξ[di,2,nd] += coeff[ptr,di]*dPdη[nd]
        end
      end
      ptr += 1
    end
  end
  # Step 3: form the scaled metrics and determinant
  for i = 1:sbp.numnodes
    dξdx[1,1,i] = dxdξ[2,2,i]
    dξdx[1,2,i] = -dxdξ[1,2,i]
    dξdx[2,1,i] = -dxdξ[2,1,i]
    dξdx[2,2,i] = dxdξ[1,1,i]
    jac[i] = one(Tmsh)/(dξdx[1,1,i]*dξdx[2,2,i] - dξdx[1,2,i]*dξdx[2,1,i])
  end
end

function calcmappingjacobian!{Tsbp,Tmsh}(sbp::TetSBP{Tsbp},
                                         mapdegree::Int,
                                         xlag::AbstractArray{Tmsh,2},
                                         xref::AbstractArray{Tmsh,2},
                                         xsbp::AbstractArray{Tmsh,2},
                                         dξdx::AbstractArray{Tmsh,3},
                                         jac::AbstractArray{Tmsh},
                                         Eone::AbstractArray{Tmsh,2})
  @assert( sbp.numnodes == size(xsbp,2) == size(dξdx,3) == size(jac,1) 
           == size(Eone,1) )
  @assert( size(xlag,1) == size(xref,1) == size(xsbp,1) == size(dξdx,1) 
           == size(dξdx,2) == size(Eone,2) == 3 )
  numdof = binomial(mapdegree+3,3)
  @assert( size(xlag,2) == size(xref,2) == numdof )
  # Step 1: find the polynomial mapping using xlag
  V = zeros(Tmsh, (numdof,numdof) )
  ptr = 1
  for r = 0:mapdegree
    for k = 0:r
      for j = 0:r-k
        i = r-j-k
        V[:,ptr] = OrthoPoly.proriolpoly(vec(xref[1,:]), vec(xref[2,:]),
                                         vec(xref[3,:]), i, j, k)
        ptr += 1
      end
    end
  end
  coeff = zeros(Tmsh, (numdof,3))
  coeff = V\(xlag.')
  # Step 2: compute the mapped SBP nodes and the analytical Jacobian at sbp nodes
  x = calcnodes(sbp) # <-- SBP nodes in reference space
  dxdξ = zeros(Tmsh, (3,3,sbp.numnodes))
  fill!(xsbp, zero(Tmsh))
  ptr = 1
  for r = 0:mapdegree
    for k = 0:r
      for j = 0:r-k
        i = r-j-k
        P = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), vec(x[3,:]), i, j, k)
        dPdξ, dPdη, dPdζ = OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]),
                                                     vec(x[3,:]), i, j, k)
        for di = 1:3
          for nd = 1:sbp.numnodes
            xsbp[di,nd] += coeff[ptr,di]*P[nd]
            dxdξ[di,1,nd] += coeff[ptr,di]*dPdξ[nd]
            dxdξ[di,2,nd] += coeff[ptr,di]*dPdη[nd]
            dxdξ[di,3,nd] += coeff[ptr,di]*dPdζ[nd]
          end
        end
        ptr += 1
      end
    end
  end
  # Step 3: find the minimum-norm solution that satisfies the metric invariants
  dξdx_targ = zeros(dxdξ)  
  for i = 1:sbp.numnodes
    for di = 1:3
      it1 = mod(di,3)+1
      it2 = mod(di+1,3)+1
      for di2 = 1:3
        it21 = mod(di2,3)+1
        it22 = mod(di2+1,3)+1
        dξdx_targ[di2,di,i] = (dxdξ[it1,it21,i]*dxdξ[it2,it22,i] - 
                               dxdξ[it1,it22,i]*dxdξ[it2,it21,i])
      end
    end
    jac[i] = one(Tmsh)/(dxdξ[1,1,i]*dxdξ[2,2,i]*dxdξ[3,3,i] +
                        dxdξ[1,2,i]*dxdξ[2,3,i]*dxdξ[3,1,i] +
                        dxdξ[1,3,i]*dxdξ[2,1,i]*dxdξ[3,2,i] -
                        dxdξ[1,1,i]*dxdξ[2,3,i]*dxdξ[3,2,i] -
                        dxdξ[1,2,i]*dxdξ[2,1,i]*dxdξ[3,3,i] -
                        dxdξ[1,3,i]*dxdξ[2,2,i]*dxdξ[3,1,i])
  end
  Qt = zeros(Tsbp, (sbp.numnodes, 3*sbp.numnodes) )
  Qt = [sbp.Q[:,:,1].' sbp.Q[:,:,2].' sbp.Q[:,:,3].']
  Qtinv = pinv(Qt)
  targ = zeros(Tmsh, (3*sbp.numnodes))
  sol = zeros(targ)
  for di = 1:3
    for di2 = 1:3
      for i = 1:sbp.numnodes      
        targ[i + (di2-1)*sbp.numnodes] = dξdx_targ[di2,di,i]
      end
    end
    b = Qt*targ - Eone[:,di]
    sol = Qtinv*b
    for di2 = 1:3
      for i = 1:sbp.numnodes
        dξdx[di2,di,i] = dξdx_targ[di2,di,i] - sol[i + (di2-1)*sbp.numnodes]
      end
    end  
  end
end


# @doc """
# ### SummationByParts.calcmappingjacobianreverse!

# Forms the reverse-mode of algorithmic differentiation product for the method
# `calcmappingjacobian`.  Specifically, it computes `xlag_r` = `xsbp_r`^T *
# ∂`xsbp`/∂`xlag` + `dξdx_r`^T * ∂`dξdx`/∂`xlag` + `jac_r`^T * ∂`jac`/∂`xlag`
# where the various quantities are defined below (the `_r` denotes the reverse
# mode) .  Note that the input and output order of the arguments follows that used
# in `calcmappingjacobian`.

# **Inputs**

# * `sbp`: an SBP operator type
# * `mapdegree`: the polynomial degree of the mapping
# * `xsbp_r`: gradient of location of the SBP nodes; [coord, sbp node]
# * `dξdx_r`: gradient of the scaled Jacobian; [ref coord, phys coord, sbp node]
# * `jac_r`: gradient of the determinant of the Jacobian; [sbp node]
# * `xref`: Lagrangian nodes in reference space; [coord, Lagrangian node]

# **In/Outs**

# * `xlag_r`: gradient of Lagrangian nodes; [coord, Lagrangian node]
# * `Eone_r`: gradient of Ex*one, Ey*one (Ez*one); [sbp node, coord]

# **Notes**

# See `calcmappingjacobian!` for an explanation of `Eone`; it is only needed in
# the 3D case.

# """->
# function calcmappingjacobianreverse!{Tsbp,Tmsh}(sbp::TriSBP{Tsbp}, 
#                                                 mapdegree::Int,
#                                                 xlag_r::AbstractArray{Tmsh,2},
#                                                 xref::AbstractArray{Tmsh,2},
#                                                 xsbp_r::AbstractArray{Tmsh,2},
#                                                 dξdx_r::AbstractArray{Tmsh,3},
#                                                 jac_r::AbstractArray{Tmsh},
#                                                 Eone_r::AbstractArray{Tmsh,2})
#   @assert( sbp.numnodes == size(xsbp_r,2) == size(dξdx_r,3) == size(jac_r,1) )
#   @assert( size(xlag_r,1) == size(xref,1) == size(xsbp_r,1) == size(dξdx_r,1) 
#            == size(dξdx_R,2) == 2 )
#   numdof = binomial(mapdegree+2,2)
#   @assert( size(xlag_r,2) == size(xref,2) == numdof )
#   # Step 3: compute dxdξ_r = dξdx_r^T*∂(dξdx)/∂(dxdξ) + jac_r^T*∂(jac)/∂(dxdξ)
#   dxdξ_r = zeros(Tmsh, (2,2,sbp.numnodes))
#   for i = 1:sbp.numnodes
#     # first, we need to account for dependence of determinant on dξdx
#     fac = jac_r[i]*jac[i]*jac[i]
#     dξdx_r[1,1,i] -= fac*dξdx[2,2,i]
#     dξdx_r[2,2,i] -= fac*dξdx[1,1,i]
#     dξdx_r[1,2,i] += fac*dξdx[2,1,i]
#     dξdx_r[2,1,i] += fac*dξdx[1,2,i]
#     # now dxdξ_r = (dξdx_r^T + jac_r^T*∂(jac)/∂(dξdx))*∂(dξdx)/∂(dxdξ)
#     dxdξ_r[2,2,i] += dξdx_r[1,1,i] 
#     dxdξ_r[1,2,i] -= dξdx_r[1,2,i] 
#     dxdξ_r[2,1,i] -= dξdx_r[2,1,i]
#     dxdξ_r[1,1,i] += dξdx_r[2,2,i]
#   end
#   # Step 2: coeff_r = dxdξ_r^T*∂(dxdξ)/∂(coeff) + xsbp_r^T*∂(xsbp)/∂(coeff)
#   x = calcnodes(sbp) # <-- SBP nodes in reference space
#   coeff_r = zeros(Tmsh, (numdof,2))
#   ptr = 1
#   for r = 0:mapdegree
#     for j = 0:r
#       i = r-j
#       P = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
#       dPdξ, dPdη = OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
#       for di = 1:2
#         for nd = 1:sbp.numnodes
#           coeff_r[ptr,di] += (xsbp_r[di,nd]*P[nd] + dxdξ_r[di,1,nd]*dPdξ[nd] +
#                               dxdξ_r[di,2,nd]*dPdη[nd])
#         end
#       end
#       ptr += 1
#     end
#   end
#   # Step 1: xlag_r = coeff_r^T*∂(coeff)/∂(xlag)
#   V = zeros(Tmsh, (numdof,numdof) )
#   ptr = 1
#   for r = 0:mapdegree
#     for j = 0:r
#       i = r-j
#       V[:,ptr] = OrthoPoly.proriolpoly(vec(xref[1,:]), vec(xref[2,:]), i, j)
#       ptr += 1
#     end
#   end
#   xlag_r = ((V.')\coeff_r).'
#   Eone_r = Array(Tmsh,0,0)
# end