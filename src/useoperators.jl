# This file gathers together functions related to using the SBP operators

@doc """
### SummationByParts.applyQ!

Applies the SBP Q matrix operator to data in `u` and stores the result in `res`.
Different methods are available depending on the rank of `u`:

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
* `di`: direction index of the operator that is desired (di=1 for Qx, etc)
* `u`: the array that the operator is applied to

**In/Outs**

* `res`: where the result of applying Q[:,:,di] to u is stored

"""->
function applyQ!{T}(sbp::SBPOperator{T}, di::Int, u::Array{T,2}, res::Array{T,2})
  @assert( sbp.numnodes == size(u,1) && sbp.numnodes == size(res,1) )
  @assert( length(u) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  for elem = 1:size(u,2)
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        res[i,elem] += sbp.Q[i,j,di]*u[j,elem]
      end
    end
  end
end

function applyQ!{T}(sbp::SBPOperator{T}, di::Int, u::Array{T,3}, res::Array{T,3})
  @assert( sbp.numnodes == size(u,2) && sbp.numnodes == size(res,2) )
  @assert( length(u) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  for elem = 1:size(u,3)
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        for field = 1:size(u,1)
          res[field,i,elem] += sbp.Q[i,j,di]*u[field,j,elem]
        end
      end
    end
  end
end

@doc """
### SummationByParts.applyD!

Applies the SBP differentiation matrix operator, D, to data in `u` and stores
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
* `di`: direction index of the operator that is desired (di=1 for Dx, etc)
* `u`: the array that the operator is applied to

**In/Outs**

* `res`: where the result of applying inv(H)*Q[:,:,di] to u is stored

"""->
function applyD!{T}(sbp::SBPOperator{T}, di::Int, u::Array{T,2}, res::Array{T,2})
  @assert( sbp.numnodes == size(u,1) && sbp.numnodes == size(res,1) )
  @assert( length(u) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  Hinv = 1./sbp.w
  for elem = 1:size(u,2)
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        res[i,elem] += sbp.Q[i,j,di]*u[j,elem]
      end
      res[i,elem] *= Hinv[i]
    end
  end
end

function applyD!{T}(sbp::SBPOperator{T}, di::Int, u::Array{T,3}, res::Array{T,3})
  @assert( sbp.numnodes == size(u,2) && sbp.numnodes == size(res,2) )
  @assert( length(u) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  Hinv = 1./sbp.w
  for elem = 1:size(u,3)
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        for field = 1:size(u,1)
          res[field,i,elem] += sbp.Q[i,j,di]*u[field,j,elem]
        end
      end
      for field = 1:size(u,1)
        res[field,i,elem] *= Hinv[i]
      end
    end
  end
end

@doc """
### SummationByParts.applyH!

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

**In/Outs**

* `res`: where the result of applying H to u is stored

"""->
function applyH!{T}(sbp::SBPOperator{T}, u::Array{T,2}, res::Array{T,2})
  @assert( sbp.numnodes == size(u,1) && sbp.numnodes == size(res,1) )
  @assert( length(u) == length(res) )
  for elem = 1:size(u,2)
    for i = 1:sbp.numnodes
      res[i,elem] += sbp.w[i]*u[i,elem]
    end
  end
end

function applyH!{T}(sbp::SBPOperator{T}, u::Array{T,3}, res::Array{T,3})
  @assert( sbp.numnodes == size(u,2) && sbp.numnodes == size(res,2) )
  @assert( length(u) == length(res) )
  for elem = 1:size(u,3)
    for i = 1:sbp.numnodes
      for field = 1:size(u,1)
        res[field,i,elem] += sbp.w[i]*u[field,i,elem]
      end
    end
  end
end

@doc """
### SummationByParts.mappingjacobian!

Evaluates the (scaled) Jacobian of the mapping from reference coordinates to
physical coordinates, as well as the determinant of the Jacobian.  The values
returned in dxidx are scaled by the determinant, so they have the same units as
the boundary measure (i.e. length in 2D, or length^2 in 3D).  This scaling is
adopted, because conservation laws written in conservative form in the reference
frame use the scaled Jacobian.

**Inputs**

* `sbp`: an SBP operator type
* `x`: the physical coordinates; 1st dim = coord, 2nd dim = node, 3rd dim = elem

**In/Outs**

* `dxidx`: the scaled Jacobian of the mapping; 1st dim = ref coord, 2nd dim =
  phys coord, 3rd dim = node, 3rd dim = elem.
* `jac`: the determinant of the Jacobian

"""->
function mappingjacobian!{T}(sbp::TriSBP{T}, x::Array{T,3}, dxidx::Array{T,4},
                             jac::Array{T,2})
  @assert( sbp.numnodes == size(x,2) && sbp.numnodes == size(dxidx,3) )
  @assert( size(x,3) == size(dxidx,4) )
  @assert( size(x,1) == 2 && size(dxidx,1) == 2 && size(dxidx,2) == 2 )
  fill!(dxidx, zero(T))
  dxdxi = zeros(T, (2,sbp.numnodes,size(x,3)))
  # compute d(x,y)/dxi and set deta/dx and deta/dy
  applyD!(sbp, 1, x, dxdxi)
  dxidx[2,1,:,:] = -dxdxi[2,:,:]
  dxidx[2,2,:,:] = dxdxi[1,:,:]
  # compute d(x,y)/deta and set dxi/dx and dxi/dy
  fill!(dxdxi, zero(T))
  applyD!(sbp, 2, x, dxdxi)
  dxidx[1,2,:,:] = -dxdxi[1,:,:]
  dxidx[1,1,:,:] = dxdxi[2,:,:]
  # compute the determinant of the Jacobian
  for elem = 1:size(x,3)
    for i = 1:sbp.numnodes
      jac[i,elem] = dxidx[1,1,i,elem]*dxidx[2,2,i,elem] - 
      dxidx[1,2,i,elem]*dxidx[2,1,i,elem]
    end
  end
  # check for negative jac here?
end

function mappingjacobian!{T}(sbp::TetSBP{T}, x::Array{T,3}, dxidx::Array{T,4},
                             jac::Array{T,2})
  @assert( sbp.numnodes == size(x,2) && sbp.numnodes == size(dxidx,3) )
  @assert( size(x,3) == size(dxidx,4) )
  @assert( size(x,1) == 3 && size(dxidx,1) == 3 && size(dxidx,2) == 3 )
  error("not implemented yet")
end