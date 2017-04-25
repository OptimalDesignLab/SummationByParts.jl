# This file gathers together functions related to the reverse-mode differentiation of mapping Jacobian

@doc """
### SummationByParts.calcMappingJacobianRevDiff!

Forms the reverse-mode of algorithmic differentiation product for the method
`calcMappingJacobian!`.  Specifically, it computes `xlag_r` = `xsbp_r`^T *
∂`xsbp`/∂`xlag` + `dξdx_r`^T * ∂`dξdx`/∂`xlag` + `jac_r`^T * ∂`jac`/∂`xlag`
and `Eone_r` = `dξdx_r`^T * ∂`dξdx`/∂`xlag` where the various quantities are defined below (the `_r` denotes the reverse mode).  Note that the input and output order of the arguments follows that used in `calcMappingJacobian`.

**Inputs**

* `sbp`: an SBP operator type
* `mapdegree`: the polynomial degree of the mapping
* `xref`: Lagrangian nodes in reference space; [coord, Lagrangian node]
* `xsbp_r`: gradient w.r.t. SBP nodes; [coord, sbp node]
* `dξdx`: scaled Jacobian of mapping (as output from calcMappingJacobian!)
* `dξdx_r`: gradient w.r.t. Jacobian; [ref coord, phys coord, sbp node, element]
* `jac`: the determinant of the Jacobian (as output from calcMappingJacobian!)
* `jac_r`: gradient w.r.t. determinant of the Jacobian; [sbp node, element]

**In/Outs**

* `xlag_r`: gradient w.r.t. Lagrangian nodes; [coord, Lagrangian node, element]
* `Eone_r`: gradient w.r.t. Ex*one, Ey*one (Ez*one); [sbp node, coord, element]

**Notes**

See `calcMappingJacobian!` for an explanation of `Eone`; it is only needed in
the 3D case, but Eone_r needs to be supplied in both 2D and 3D.

"""->
function calcMappingJacobianRevDiff!{
  Tsbp,Tmsh, T2}(sbp::TriSBP{Tsbp}, mapdegree::Int, xref::AbstractArray{T2,2},
             xlag_r::AbstractArray{Tmsh,3}, xsbp_r::AbstractArray{Tmsh,3},
             dξdx::AbstractArray{Tmsh,4}, dξdx_r::AbstractArray{Tmsh,4},
             jac::AbstractArray{Tmsh,2}, jac_r::AbstractArray{Tmsh,2},
             Eone_r::AbstractArray{Tmsh,3})

  @assert( sbp.numnodes == size(xsbp_r,2) == size(dξdx,3) ==size(dξdx_r,3) ==
           size(jac_r,1) == size(Eone_r,1) )
  @assert( size(xlag_r,1) == size(xref,1) == size(xsbp_r,1) == size(dξdx,1) ==
           size(dξdx_r,1) == size(dξdx,2) == size(dξdx_r,2) == size(Eone_r,2) ==
           2 )
  @assert( size(xsbp_r,3) == size(xlag_r,3) == size(dξdx,4) == size(dξdx_r,4) ==
           size(jac,2) == size(jac_r,2) == size(Eone_r,3) )
  numdof = binomial(mapdegree+2,2)
  @assert( size(xlag_r,2) == size(xref,2) == numdof )
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
  fill!(xlag_r, zero(Tmsh))
  fill!(Eone_r, zero(Tmsh))
  coeff_r = zeros(Tmsh, (numdof,2))
  dxdξ_r = zeros(Tmsh, (2,2,sbp.numnodes))
  for e = 1:size(xlag_r,3)
    # compute dxdξ_r = dξdx_r^T*∂(dξdx)/∂(dxdξ) + jac_r^T*∂(jac)/∂(dxdξ)
    fill!(dxdξ_r,zero(Tmsh))
    for i = 1:sbp.numnodes
      # first, we need to account for dependence of determinant on dξdx
      fac = jac_r[i,e]*jac[i,e]*jac[i,e]
      dξdx_r[1,1,i,e] -= fac*dξdx[2,2,i,e]
      dξdx_r[2,2,i,e] -= fac*dξdx[1,1,i,e]
      dξdx_r[1,2,i,e] += fac*dξdx[2,1,i,e]
      dξdx_r[2,1,i,e] += fac*dξdx[1,2,i,e]
      # now dxdξ_r = (dξdx_r^T + jac_r^T*∂(jac)/∂(dξdx))*∂(dξdx)/∂(dxdξ)
      dxdξ_r[2,2,i] += dξdx_r[1,1,i,e]
      dxdξ_r[1,2,i] -= dξdx_r[1,2,i,e]
      dxdξ_r[2,1,i] -= dξdx_r[2,1,i,e]
      dxdξ_r[1,1,i] += dξdx_r[2,2,i,e]
    end
    # coeff_r = dxdξ_r^T*∂(dxdξ)/∂(coeff) + xsbp_r^T*∂(xsbp)/∂(coeff)
    fill!(coeff_r, zero(Tmsh))
    for i = 1:numdof
      for di = 1:2
        for nd = 1:sbp.numnodes
          coeff_r[i,di] += (xsbp_r[di,nd,e]*P[nd,i] + dxdξ_r[di,1,nd]*dPdξ[nd,i] + dxdξ_r[di,2,nd]*dPdη[nd,i])
        end
      end
    end
    # xlag_r = coeff_r^T*∂(coeff)/∂(xlag)
    for di = 1:2
      for i = 1:numdof
        for j = 1:numdof
          # coeff[i,di] += Vinv[i,j]*xlag[di,j,e]
          xlag_r[di,j,e] += coeff_r[i,di]*Vinv[i,j]
        end
      end
    end
  end # loop over elements
end

# function calcMappingJacobian!{Tsbp,Tmsh}(sbp::TriSBP{Tsbp},
#                                          mapdegree::Int,
#                                          xref::AbstractArray{Tmsh,2},
#                                          xlag_r::AbstractArray{Tmsh,3},
#                                          xsbp_r::AbstractArray{Tmsh,3},
#                                          dξdx::AbstractArray{Tmsh,4},
#                                          dξdx_r::AbstractArray{Tmsh,4},
#                                          jac::AbstractArray{Tmsh,2},
#                                          jac_r::AbstractArray{Tmsh,2},
#                                          Eone_r::AbstractArray{Tmsh,3})
#   @assert( sbp.numnodes == size(xsbp_r,2) == size(dξdx,3) == size(dξdx_r,3) == size(jac,1) == size(jac_r,1) == size(Eone,1) )
#   @assert( size(xlag_r,1) == size(xref,1) == size(xsbp_r,1) == size(dξdx,1)
#            == size(dξdx_r,1) == size(dξdx,2) == size(dξdx_r,2) == size(Eone_r,2) == 3 )
#   @assert( size(xsbp_r,3) == size(xlag_r,3) == size(dξdx,4) == size(dξdx,4) == size(jac,2) == size(jac_r,2) == size(Eone,3) )
#   numdof = binomial(mapdegree+3,3)
#   @assert( size(xlag_r,2) == size(xref,2) == numdof )
#   # find the inverse of the Vandermonde matrix
#   V = zeros(Tmsh, (numdof,numdof) )
#   ptr = 1
#   for r = 0:mapdegree
#     for k = 0:r
#       for j = 0:r-k
#         i = r-j-k
#         V[:,ptr] = OrthoPoly.proriolpoly(vec(xref[1,:]), vec(xref[2,:]),
#                                          vec(xref[3,:]), i, j, k)
#         ptr += 1
#       end
#     end
#   end
#   Vinv = inv(V)
#   # get the SBP nodes in reference space in order to find the orthogonal
#   # polynomials and their derivatives at these nodes
#   x = calcnodes(sbp)
#   P = zeros(Tmsh, (sbp.numnodes, numdof))
#   dPdξ = zeros(Tmsh, (sbp.numnodes, numdof))
#   dPdη = zeros(Tmsh, (sbp.numnodes, numdof))
#   dPdζ = zeros(Tmsh, (sbp.numnodes, numdof))
#   ptr = 1
#   for r = 0:mapdegree
#     for k = 0:r
#       for j = 0:r-k
#         i = r-j-k
#         P[:,ptr] = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), vec(x[3,:]),
#                                          i, j, k)
#         dPdξ[:,ptr], dPdη[:,ptr] , dPdζ[:,ptr] =
#           OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), vec(x[3,:]),
#                                     i, j, k)
#         ptr += 1
#       end
#     end
#   end

#   fill!(xsbp, zero(Tmsh))
#   coeff_r = zeros(Tmsh, (numdof,3))
#   dxdξ_r = zeros(Tmsh, (3,3,sbp.numnodes))
#   dξdx_targ_r = zeros(dxdξ_r)
#   Qt = zeros(Tsbp, (sbp.numnodes, 3*sbp.numnodes) )
#   Qt = [sbp.Q[:,:,1].' sbp.Q[:,:,2].' sbp.Q[:,:,3].']
#   Qtinv = pinv(Qt)
#   targ_r = zeros(Tmsh, (3*sbp.numnodes))
#   sol_r = zeros(targ_r)
#   # loop over each element...
#   for e = 1:size(xlag,3)

#     # compute dξdx_targ_r = dξdx_r^T ∂(dξdx)/∂(dξdx_r) and Eone_r
#     fill!(targ_r, zero(Tmsh))
#     fill!(sol_r, zero(Tmsh))
#     for di = 1:3
#       fill!(sol_r, zero(Tmsh))
#       for di2 = 1:3
#         for i = 1:sbp.numnodes
#           #dξdx[di2,di,i] = dξdx_targ[di2,di,i] - sol[i + (di2-1)*sbp.numnodes]
#           dξdx_targ_r[di2,di,i] += dξdx_r[di2,di,i]
#           sol_r[i + (di2-1)*sbp.numnodes] -= dξdx_r[di2,di,i]
#         end
#       end
#       # b = Qt*targ - Eone[:,di]
#       # sol = Qtinv*b
#       b_r = Qtinv'*sol_r
#       targ_r -= Qt'*b_r
#       Eone_r[:,di] -=  b_r
#       for di2 = 1:3
#         for i = 1:sbp.numnodes
#           #targ[i + (di2-1)*sbp.numnodes] = dξdx_targ[di2,di,i]
#           dξdx_targ_r[di2,di,i] += targ_r[i + (di2-1)*sbp.numnodes]
#         end
#       end



#       # check that Eone sums to zero
#       @assert( abs(sum(Eone[:,di,e])) < 1e-14 )
#       for di2 = 1:3
#         for i = 1:sbp.numnodes
#           targ[i + (di2-1)*sbp.numnodes] = dξdx_targ[di2,di,i]
#         end
#       end
#       b = Qt*targ - Eone[:,di,e]
#       sol = Qtinv*b
#       for di2 = 1:3
#         for i = 1:sbp.numnodes
#           dξdx[di2,di,i,e] = dξdx_targ[di2,di,i] - sol[i + (di2-1)*sbp.numnodes]
#         end
#       end
#     end
#     end





#     # find the coefficents of the polynomial mapping using xlag and Vinv
#     for di = 1:3
#       for i = 1:numdof
#         coeff[i,di] = zero(Tmsh)
#         for j = 1:numdof
#           coeff[i,di] += Vinv[i,j]*xlag[di,j,e]
#         end
#       end
#     end
#     # compute the mapped SBP nodes and the coordinate derivatives at these nodes
#     fill!(dxdξ, zero(Tmsh))
#     for i = 1:numdof
#       for di = 1:3
#         for nd = 1:sbp.numnodes
#           xsbp[di,nd,e] += coeff[i,di]*P[nd,i]
#           dxdξ[di,1,nd] += coeff[i,di]*dPdξ[nd,i]
#           dxdξ[di,2,nd] += coeff[i,di]*dPdη[nd,i]
#           dxdξ[di,3,nd] += coeff[i,di]*dPdζ[nd,i]
#         end
#       end
#     end
#     # find the target values for the mapping Jacobian; also, find the Jacobian
#     # determinant
#     for i = 1:sbp.numnodes
#       for di = 1:3
#         it1 = mod(di,3)+1
#         it2 = mod(di+1,3)+1
#         for di2 = 1:3
#           it21 = mod(di2,3)+1
#           it22 = mod(di2+1,3)+1
#           dξdx_targ[di2,di,i] = (dxdξ[it1,it21,i]*dxdξ[it2,it22,i] -
#                                  dxdξ[it1,it22,i]*dxdξ[it2,it21,i])
#         end
#       end
#       jac[i] = one(Tmsh)/(dxdξ[1,1,i]*dxdξ[2,2,i]*dxdξ[3,3,i] +
#                           dxdξ[1,2,i]*dxdξ[2,3,i]*dxdξ[3,1,i] +
#                           dxdξ[1,3,i]*dxdξ[2,1,i]*dxdξ[3,2,i] -
#                           dxdξ[1,1,i]*dxdξ[2,3,i]*dxdξ[3,2,i] -
#                           dxdξ[1,2,i]*dxdξ[2,1,i]*dxdξ[3,3,i] -
#                           dxdξ[1,3,i]*dxdξ[2,2,i]*dxdξ[3,1,i])
#     end
#     # find the minimum-norm solution that satisfies the metric invariants
#     for di = 1:3
#       # check that Eone sums to zero
#       @assert( abs(sum(Eone[:,di])) < 1e-14 )
#       for di2 = 1:3
#         for i = 1:sbp.numnodes
#           targ[i + (di2-1)*sbp.numnodes] = dξdx_targ[di2,di,i]
#         end
#       end
#       b = Qt*targ - Eone[:,di]
#       sol = Qtinv*b
#       for di2 = 1:3
#         for i = 1:sbp.numnodes
#           dξdx[di2,di,i] = dξdx_targ[di2,di,i] - sol[i + (di2-1)*sbp.numnodes]
#         end
#       end
#     end
#   end
# end

