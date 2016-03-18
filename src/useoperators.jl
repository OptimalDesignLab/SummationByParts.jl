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

This function returns the node coordinates for an SBP operator.  The nodes are
ordered first by vertex, then by edge (with nodes ordered in sequence along the
directed edge), then by face (if appropriate), and then finally by volumn nodes.
This function assumes the element mapping is linear, i.e. edges are lines.

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
function calcnodes{T}(sbp::TriSBP{T}, vtx::Array{T}=sbp.vtx)
  perm, faceperm = SummationByParts.getnodepermutation(sbp.cub, sbp.degree)
  x = zeros(T, (2, sbp.numnodes))
  x = SymCubatures.calcnodes(sbp.cub, vtx)
  return x[:,perm]
end

function calcnodes{T}(sbp::TetSBP{T}, vtx::Array{T}=sbp.vtx)
  perm, faceperm = SummationByParts.getnodepermutation(sbp.cub, sbp.degree)
  x = zeros(T, (3, sbp.numnodes))
  x = SymCubatures.calcnodes(sbp.cub, vtx)
  return x[:,perm]
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
function calcminnodedistance{T}(sbp::AbstractSBP{T}, vtx::Array{T})
  x = calcnodes(sbp, vtx)
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
* `trans` (optional): if true, the transpose operation is applied

**In/Outs**

* `res`: where the result of applying Q[:,:,di] to u is stored

"""->
function weakdifferentiate!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int, 
                                            flux::AbstractArray{Tflx,2},
                                            res::AbstractArray{Tres,2};
                                            trans::Bool=false)
  @assert( sbp.numnodes == size(flux,1) && sbp.numnodes == size(res,1) )
  @assert( length(flux) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  if trans # apply transposed Q
    @inbounds begin
      for elem = 1:size(flux,2)
        for i = 1:sbp.numnodes
          for j = 1:sbp.numnodes
            res[i,elem] += sbp.Q[j,i,di]*flux[j,elem] 
          end
        end
      end
    end
  else # apply Q
    @inbounds begin
      for elem = 1:size(flux,2)
        for i = 1:sbp.numnodes
          for j = 1:sbp.numnodes
            res[i,elem] += sbp.Q[i,j,di]*flux[j,elem] 
          end
        end
      end
    end
  end
end

function weakdifferentiate!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                                            flux::AbstractArray{Tflx,3},
                                            res::AbstractArray{Tres,3};
                                            trans::Bool=false)
  @assert( sbp.numnodes == size(flux,2) && sbp.numnodes == size(res,2) )
  @assert( length(flux) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  if trans # apply transposed Q
    @inbounds begin
      for elem = 1:size(flux,3)
        for i = 1:sbp.numnodes
          for j = 1:sbp.numnodes
            for field = 1:size(flux,1)
              res[field,i,elem] += sbp.Q[j,i,di]*flux[field,j,elem]
            end
          end
        end
      end
    end
  else # apply Q
    @inbounds begin
      for elem = 1:size(flux,3)
        for i = 1:sbp.numnodes
          for j = 1:sbp.numnodes
            for field = 1:size(flux,1)
              res[field,i,elem] += sbp.Q[i,j,di]*flux[field,j,elem]
            end
          end
        end
      end
    end
  end
end

# macro update(x,op,y)
#   return :($x=$op($x,$y))
# end
# function weakdifferentiate!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int, 
#                                             flux::AbstractArray{Tflx,2},
#                                             res::AbstractArray{Tres,2},
#                                             op::Function; trans::Bool=false)
#   @assert( sbp.numnodes == size(flux,1) && sbp.numnodes == size(res,1) )
#   @assert( length(flux) == length(res) )
#   @assert( di > 0 && di <= size(sbp.Q,3) )
#   if trans # apply transposed Q
#     @inbounds begin
#       for elem = 1:size(flux,2)
#         for i = 1:sbp.numnodes
#           for j = 1:sbp.numnodes
#             # for op = + the following is equivalent to
#             # res[i,elem] += sbp.Q[j,i,di]*flux[j,elem] 
#             @update(res[i,elem], op, sbp.Q[j,i,di]*flux[j,elem])
#           end
#         end
#       end
#     end
#   else # apply Q
#     @inbounds begin
#       for elem = 1:size(flux,2)
#         for i = 1:sbp.numnodes
#           for j = 1:sbp.numnodes
#             # for op = + the following is equivalent to
#             # res[i,elem] += sbp.Q[i,j,di]*flux[j,elem] 
#             @update(res[i,elem], op, sbp.Q[i,j,di]*flux[j,elem])
#           end
#         end
#       end
#     end
#   end
# end

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

**In/Outs**

* `res`: where the result of applying inv(H)*Q[:,:,di] to u is stored

"""->
function differentiate!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                                        flux::AbstractArray{Tflx,2},
                                        res::AbstractArray{Tres,2})
  @assert( sbp.numnodes == size(flux,1) && sbp.numnodes == size(res,1) )
  @assert( length(flux) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  Hinv = 1./sbp.w
  @inbounds begin
    for elem = 1:size(flux,2)
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          res[i,elem] += sbp.Q[i,j,di]*flux[j,elem]
        end
        res[i,elem] *= Hinv[i]
      end
    end
  end
end

function differentiate!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                                        flux::AbstractArray{Tflx,3},
                                        res::AbstractArray{Tres,3})
  @assert( sbp.numnodes == size(flux,2) && sbp.numnodes == size(res,2) )
  @assert( length(flux) == length(res) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  Hinv = 1./sbp.w
  @inbounds begin
    for elem = 1:size(flux,3)
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          for field = 1:size(flux,1)
            res[field,i,elem] += sbp.Q[i,j,di]*flux[field,j,elem]
          end
        end
        for field = 1:size(flux,1)
          res[field,i,elem] *= Hinv[i]
        end
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
  @inbounds begin
    for di = 1:size(sbp.Q, 3)
      tmp = zero(Tmsh)*zero(Tsol)
      for j = 1:sbp.numnodes
        tmp += sbp.Q[i,j,di]*u[j]
      end
      Ddir += dir[di]*tmp/sbp.w[i]
    end
  end
  return Ddir
end

function directionaldifferentiate!{Tsbp,Tmsh,Tsol,Tres}(sbp::AbstractSBP{Tsbp},
                                                        dir::Array{Tmsh,1}, 
                                                        u::AbstractArray{Tsol,2},
                                                        i::Int,
                                                        Ddir::Array{Tres,1})
  @assert( size(sbp.Q, 3) == size(dir,1) )
  @inbounds begin
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

**In/Outs**

* `res`: where the result of applying H to u is stored

"""->
function volumeintegrate!{Tsbp,Tsol,Tres}(sbp::AbstractSBP{Tsbp},
                                          u::AbstractArray{Tsol,2},
                                          res::AbstractArray{Tres,2})
  @assert( sbp.numnodes == size(u,1) && sbp.numnodes == size(res,1) )
  @assert( length(u) == length(res) )
  H = sbp.w
  @inbounds begin
    for elem = 1:size(u,2)
      for i = 1:sbp.numnodes
        res[i,elem] += H[i]*u[i,elem]
      end
    end
  end
end

function volumeintegrate!{Tsbp,Tsol,Tres}(sbp::AbstractSBP{Tsbp},
                                          u::AbstractArray{Tsol,3},
                                          res::AbstractArray{Tres,3})
  @assert( sbp.numnodes == size(u,2) && sbp.numnodes == size(res,2) )
  @assert( length(u) == length(res) )
  H = sbp.w
  @inbounds begin
    for elem = 1:size(u,3)
      for i = 1:sbp.numnodes
        for field = 1:size(u,1)
          res[field,i,elem] += H[i]*u[field,i,elem]
        end
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

