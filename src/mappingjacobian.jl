# This file gathers together functions related to the calculation of the
# coordinate mapping Jacobian on elements and faces

@doc """
### SummationByParts.calcMappingJacobian!

Uses a given set of Lagrangian element nodes to determine an analytical
(polynomial) mapping, and then uses this mapping to determine the Jacobian of
the mapping.  The approach varies depending on the dimension of the problem:

* For *2-dimensional problems* the exact Jacobian of the mapping is used;
* For *3-dimensional problems* an optimization problem is solved

**Inputs**

* `sbp`: an SBP operator type
* `mapdegree`: the polynomial degree of the mapping
* `xref`: Lagrangian nodes in reference space; [coord, Lagrangian node]
* `xlag`: Lagrangian nodes in physical space; [coord, Lagrangian node, element]
* `Eone`: Ex*one, Ey*one (Ez*one); [sbp node, coord, element] see notes below

**In/Outs**

* `xsbp`: location of the SBP nodes in physical space; [coord, sbp node, element]
* `dξdx`: scaled Jacobian of mapping; [ref coord, phys coord, sbp node, element]
* `jac`: the determinant of the Jacobian; [sbp node, element]

**Notes**

The array `Eone` is the product of the boundary operators, in physical space,
with the vector of ones (see Crean et al., Entropy-Conservative,
Multidimensional Summation-By-Parts Discretization of the Euler Equations, as
well as the test in `test/test_mappingjacobian.jl`).  These products are used to
define the metric invariants.  *They are not needed by the 2-dimensional code*,
and so this array can be passed empty in that case.

"""->
function calcMappingJacobian!{Tsbp,Tmsh}(sbp::TriSBP{Tsbp}, mapdegree::Int,
                                         xref::AbstractArray{Tmsh,2},
                                         xlag::AbstractArray{Tmsh,3},
                                         xsbp::AbstractArray{Tmsh,3},
                                         dξdx::AbstractArray{Tmsh,4},
                                         jac::AbstractArray{Tmsh,2},
                                         Eone::AbstractArray{Tmsh,3}=
                                           Array(Tmsh,0,0,0))
  @assert( sbp.numnodes == size(xsbp,2) == size(dξdx,3) == size(jac,1) )
  @assert( size(xlag,1) == size(xref,1) == size(xsbp,1) == size(dξdx,1) 
           == size(dξdx,2) == 2 )
  @assert( size(xsbp,3) == size(xlag,3) == size(dξdx,4) == size(jac,2) )
  numdof = binomial(mapdegree+2,2)
  @assert( size(xlag,2) == size(xref,2) == numdof )
  # find the inverse of the Vandermonde matrix
  V = zeros(Tmsh, (numdof,numdof) )
  ptr = 1
  for r = 0:mapdegree
    for j = 0:r
      i = r-j
      V[:,ptr] = OrthoPoly.proriolpoly(vec(xref[1,:]), vec(xref[2,:]), i, j)
      ptr += 1
    end
  end
  Vinv = inv(V)
  # get the SBP nodes in reference space in order to find the orthogonal
  # polynomials and their derivatives at these nodes
  x = calcnodes(sbp)
  P = zeros(Tmsh, (sbp.numnodes, numdof))
  dPdξ = zeros(Tmsh, (sbp.numnodes, numdof))
  dPdη = zeros(Tmsh, (sbp.numnodes, numdof))
  ptr = 1
  for r = 0:mapdegree
    for j = 0:r
      i = r-j
      P[:,ptr] = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      dPdξ[:,ptr], dPdη[:,ptr] =
        OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      ptr += 1
    end
  end
  # loop over each element...
  fill!(xsbp, zero(Tmsh))
  coeff = zeros(Tmsh, (numdof,2))
  dxdξ = zeros(Tmsh, (2,2,sbp.numnodes))
  for e = 1:size(xlag,3)  
    # find the coefficents of the polynomial mapping using xlag and Vinv
    for di = 1:2
      for i = 1:numdof
        coeff[i,di] = zero(Tmsh)
        for j = 1:numdof
          coeff[i,di] += Vinv[i,j]*xlag[di,j,e]
        end
      end
    end
    # compute the mapped SBP nodes and the coordinate derivatives at these nodes
    fill!(dxdξ, zero(Tmsh))
    for i = 1:numdof
      for di = 1:2
        for nd = 1:sbp.numnodes
          xsbp[di,nd,e] += coeff[i,di]*P[nd,i]
          dxdξ[di,1,nd] += coeff[i,di]*dPdξ[nd,i]
          dxdξ[di,2,nd] += coeff[i,di]*dPdη[nd,i]
        end
      end
    end
    # form the scaled metrics and determinant
    for i = 1:sbp.numnodes
      dξdx[1,1,i,e] = dxdξ[2,2,i]
      dξdx[1,2,i,e] = -dxdξ[1,2,i]
      dξdx[2,1,i,e] = -dxdξ[2,1,i]
      dξdx[2,2,i,e] = dxdξ[1,1,i]
      jac[i,e] = one(Tmsh)/(dξdx[1,1,i,e]*dξdx[2,2,i,e] -
                            dξdx[1,2,i,e]*dξdx[2,1,i,e])
    end
  end
end

function calcMappingJacobian!{Tsbp,Tmsh}(sbp::TetSBP{Tsbp}, mapdegree::Int,
                                         xref::AbstractArray{Tmsh,2},
                                         xlag::AbstractArray{Tmsh,3},
                                         xsbp::AbstractArray{Tmsh,3},
                                         dξdx::AbstractArray{Tmsh,4},
                                         jac::AbstractArray{Tmsh,2},
                                         Eone::AbstractArray{Tmsh,3})
  @assert( sbp.numnodes == size(xsbp,2) == size(dξdx,3) == size(jac,1) 
           == size(Eone,1) )
  @assert( size(xlag,1) == size(xref,1) == size(xsbp,1) == size(dξdx,1) 
           == size(dξdx,2) == size(Eone,2) == 3 )
  @assert( size(xsbp,3) == size(xlag,3) == size(dξdx,4) == size(jac,2) ==
           size(Eone,3) )
  numdof = binomial(mapdegree+3,3)
  @assert( size(xlag,2) == size(xref,2) == numdof )
  # find the inverse of the Vandermonde matrix
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
  Vinv = inv(V)
  # get the SBP nodes in reference space in order to find the orthogonal
  # polynomials and their derivatives at these nodes
  x = calcnodes(sbp)
  P = zeros(Tmsh, (sbp.numnodes, numdof))
  dPdξ = zeros(Tmsh, (sbp.numnodes, numdof))
  dPdη = zeros(Tmsh, (sbp.numnodes, numdof))
  dPdζ = zeros(Tmsh, (sbp.numnodes, numdof))
  ptr = 1
  for r = 0:mapdegree
    for k = 0:r
      for j = 0:r-k
        i = r-j-k
        P[:,ptr] = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), vec(x[3,:]),
                                         i, j, k)
        dPdξ[:,ptr], dPdη[:,ptr] , dPdζ[:,ptr] =
          OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), vec(x[3,:]),
                                    i, j, k)
        ptr += 1
      end
    end
  end

  fill!(xsbp, zero(Tmsh))
  coeff = zeros(Tmsh, (numdof,3))
  dxdξ = zeros(Tmsh, (3,3,sbp.numnodes))
  dξdx_targ = zeros(dxdξ)
  Qt = zeros(Tsbp, (sbp.numnodes, 3*sbp.numnodes) )
  Qt = [sbp.Q[:,:,1].' sbp.Q[:,:,2].' sbp.Q[:,:,3].']
  Qtinv = pinv(Qt)
  targ = zeros(Tmsh, (3*sbp.numnodes))
  sol = zeros(targ)
  # loop over each element...  
  for e = 1:size(xlag,3)
    # find the coefficents of the polynomial mapping using xlag and Vinv
    for di = 1:3
      for i = 1:numdof
        coeff[i,di] = zero(Tmsh)
        for j = 1:numdof
          coeff[i,di] += Vinv[i,j]*xlag[di,j,e]
        end
      end
    end
    # compute the mapped SBP nodes and the coordinate derivatives at these nodes
    fill!(dxdξ, zero(Tmsh))
    for i = 1:numdof
      for di = 1:3
        for nd = 1:sbp.numnodes
          xsbp[di,nd,e] += coeff[i,di]*P[nd,i]
          dxdξ[di,1,nd] += coeff[i,di]*dPdξ[nd,i]
          dxdξ[di,2,nd] += coeff[i,di]*dPdη[nd,i]
          dxdξ[di,3,nd] += coeff[i,di]*dPdζ[nd,i]
        end
      end
    end
    # find the target values for the mapping Jacobian; also, find the Jacobian
    # determinant
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
      jac[i,e] = one(Tmsh)/(dxdξ[1,1,i]*dxdξ[2,2,i]*dxdξ[3,3,i] +
                            dxdξ[1,2,i]*dxdξ[2,3,i]*dxdξ[3,1,i] +
                            dxdξ[1,3,i]*dxdξ[2,1,i]*dxdξ[3,2,i] -
                            dxdξ[1,1,i]*dxdξ[2,3,i]*dxdξ[3,2,i] -
                            dxdξ[1,2,i]*dxdξ[2,1,i]*dxdξ[3,3,i] -
                            dxdξ[1,3,i]*dxdξ[2,2,i]*dxdξ[3,1,i])
    end
    # find the minimum-norm solution that satisfies the metric invariants    
    for di = 1:3
      # check that Eone sums to zero
      @assert( abs(sum(Eone[:,di,e])) < 1e-14 )
      for di2 = 1:3
        for i = 1:sbp.numnodes      
          targ[i + (di2-1)*sbp.numnodes] = dξdx_targ[di2,di,i]
        end
      end
      b = Qt*targ - Eone[:,di,e]
      sol = Qtinv*b
      for di2 = 1:3
        for i = 1:sbp.numnodes
          dξdx[di2,di,i,e] = dξdx_targ[di2,di,i] - sol[i + (di2-1)*sbp.numnodes]
        end
      end  
    end
  end
end

@doc """
### SummationByParts.calcMappingJacobianElement!

Single element variant of calcMappingJacobian!.  Uses a given set of Lagrangian
element nodes to determine an analytical (polynomial) mapping, and then uses
this mapping to determine the Jacobian of the mapping.  The approach varies
depending on the dimension of the problem:

* For *2-dimensional problems* the exact Jacobian of the mapping is used;
* For *3-dimensional problems* an optimization problem is solved

**Inputs**

* `sbp`: an SBP operator type
* `mapdegree`: the polynomial degree of the mapping
* `xref`: Lagrangian nodes in reference space; [coord, Lagrangian node]
* `xlag`: Lagrangian nodes in physical space; [coord, Lagrangian node]
* `Eone`: Ex*one, Ey*one (Ez*one); [sbp node, coord] see notes below

**In/Outs**

* `xsbp`: location of the SBP nodes in physical space; [coord, sbp node]
* `dξdx`: the scaled Jacobian of the mapping; [ref coord, phys coord, sbp node]
* `jac`: the determinant of the Jacobian; [sbp node]

**Notes**

The array `Eone` is the product of the boundary operators, in physical space,
with the vector of ones (see Crean et al., Entropy-Conservative,
Multidimensional Summation-By-Parts Discretization of the Euler Equations, as
well as the test in `test/test_mappingjacobian.jl`).  These products are used to
define the metric invariants.  *They are not needed by the 2-dimensional code,
and so this array can be passed empty in that case*.

"""->
function calcMappingJacobianElement!{
  Tsbp,Tmsh}(sbp::TriSBP{Tsbp}, mapdegree::Int, xref::AbstractArray{Tmsh,2},
             xlag::AbstractArray{Tmsh,2}, xsbp::AbstractArray{Tmsh,2},
             dξdx::AbstractArray{Tmsh,3}, jac::AbstractArray{Tmsh},
             Eone::AbstractArray{Tmsh,2}=Array(Tmsh,0,0))
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

function calcMappingJacobianElement!{
  Tsbp,Tmsh}(sbp::TetSBP{Tsbp}, mapdegree::Int, xref::AbstractArray{Tmsh,2}, 
             xlag::AbstractArray{Tmsh,2}, xsbp::AbstractArray{Tmsh,2},
             dξdx::AbstractArray{Tmsh,3}, jac::AbstractArray{Tmsh},
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

@doc """
### SummationByParts.mappingjacobian!

Evaluates the Jacobian of the mapping from face-reference coordinates to
physical coordinates, as well as the determinant of the Jacobian.

**Inputs**

* `sbpface`: an SBP face operator type
* `ifaces`: list of element interfaces stored as an array of `Interface`s
* `x`: the physical coordinates in [coord, node, elem] format

**In/Outs**

* `dξdx`: the Jacobian in [ξ coord, x coord, face node, L/R, face] format
* `jac`: the determinant in [face node, L/R, face] format

"""->
function mappingjacobian!{Tsbp,Tmsh}(sbpface::AbstractFace{Tsbp},
                                     ifaces::Array{Interface},
                                     x::AbstractArray{Tmsh,3},
                                     dξdx::AbstractArray{Tmsh,5},
                                     jac::AbstractArray{Tmsh,3})
  @assert( size(ifaces,1) == size(dξdx,5) == size(jac,3) )
  @assert( size(dξdx,4) == size(jac,2) == 2 )
  @assert( size(dξdx,3) == size(jac,1) == sbpface.numnodes )
  @assert( size(x,1) == 2 && size(dξdx,1) == 2 && size(dξdx,2) == 2 )
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]
      # compute the derivatives dxdξ
      dxdξL = zeros(Tmsh, (2,2))
      dxdξR = zeros(Tmsh, (2,2))
      for di = 1:2
        for di2 = 1:2
          for j = 1:sbpface.dstencilsize
            dxdξL[di,di2] +=
              sbpface.deriv[j,i,di2]*x[di,sbpface.dperm[j,face.faceL],
                                       face.elementL]
            dxdξR[di,di2] +=
              sbpface.deriv[j,iR,di2]*x[di,sbpface.dperm[j,face.faceR],
                                        face.elementR]
          end
        end
      end
      # compute the Jacobian determinants
      jac[i,1,findex] = one(Tmsh)/(dxdξL[1,1]*dxdξL[2,2] - 
                                   dxdξL[1,2]*dxdξL[2,1])
      jac[iR,2,findex] = one(Tmsh)/(dxdξR[1,1]*dxdξR[2,2] - 
                                    dxdξR[1,2]*dxdξR[2,1])

      # compute the derivatives dξdx
      dξdx[1,1,i,1,findex] = dxdξL[2,2]*jac[i,1,findex]
      dξdx[1,2,i,1,findex] = -dxdξL[1,2]*jac[i,1,findex]
      dξdx[2,2,i,1,findex] = dxdξL[1,1]*jac[i,1,findex]
      dξdx[2,1,i,1,findex] = -dxdξL[2,1]*jac[i,1,findex]
      
      dξdx[1,1,iR,2,findex] = dxdξR[2,2]*jac[iR,2,findex]
      dξdx[1,2,iR,2,findex] = -dxdξR[1,2]*jac[iR,2,findex]
      dξdx[2,2,iR,2,findex] = dxdξR[1,1]*jac[iR,2,findex]
      dξdx[2,1,iR,2,findex] = -dxdξR[2,1]*jac[iR,2,findex]
    end
  end
end

@doc """
### SummationByParts.mappingjacobian!

**Deprecated**:

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
  # calculate the metrics
  for di = 1:3
    it1 = mod(di,3)+1
    it2 = mod(di+1,3)+1
    for elem = 1:numelem
      for i = 1:sbp.numnodes
        for di2 = 1:3
          it21 = mod(di2,3)+1
          it22 = mod(di2+1,3)+1
          dξdx[it1,di2,i,elem] = (dxdξ[it22,i,elem,di]*dxdξ[it21,i,elem,it2] -
                                  dxdξ[it21,i,elem,di]*dxdξ[it22,i,elem,it2])            
        end
      end
    end
  end    
  # calculate the determinant of the mapping
  for elem = 1:numelem
    for i = 1:sbp.numnodes
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
