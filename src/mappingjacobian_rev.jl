# This file gathers together functions related to the reverse-mode differentiation of mapping Jacobian

@doc """
### SummationByParts.calcMappingJacobian_rev!

Forms the reverse-mode of algorithmic differentiation product for the method
`calcMappingJacobian!`.  Specifically, it computes `xlag_bar` = `xsbp_bar`^T *
∂`xsbp`/∂`xlag` + `dξdx_bar`^T * ∂`dξdx`/∂`xlag` + `jac_bar`^T * ∂`jac`/∂`xlag` and
`Eone_bar` = `dξdx_bar`^T * ∂`dξdx`/∂`xlag` where the various quantities are defined
below (the `_bar` denotes the reverse mode).  Note that the input and output order
of the arguments follows that used in `calcMappingJacobian!`.

**Inputs**

* `sbp`: an SBP operator type
* `mapdegree`: the polynomial degree of the mapping
* `xref`: Lagrangian nodes in reference space; [coord, Lagrangian node]
* `xsbp_bar`: gradient w.r.t. SBP nodes; [coord, sbp node]
* `dξdx`: scaled Jacobian of mapping (as output from calcMappingJacobian!)
* `dξdx_bar`: gradient w.r.t. Jacobian; [ref coord, phys coord, sbp node, element]
* `jac`: the determinant of the Jacobian (as output from calcMappingJacobian!)
* `jac_bar`: gradient w.r.t. determinant of the Jacobian; [sbp node, element]

**In/Outs**

* `xlag_bar`: gradient w.r.t. Lagrangian nodes; [coord, Lagrangian node, element]
* `Eone_bar`: gradient w.r.t. Ex*one, Ey*one (Ez*one); [sbp node, coord, element]

**Notes**

See `calcMappingJacobian!` for an explanation of `Eone`; it is only needed in
the 3D case, but Eone_bar needs to be supplied in both 2D and 3D.

"""->
function calcMappingJacobian_rev!{
  Tsbp,Tmsh}(sbp::TriSBP{Tsbp}, mapdegree::Int, xref::AbstractArray{Tmsh,2},
             xlag::AbstractArray{Tmsh,3}, xlag_bar::AbstractArray{Tmsh,3},
             xsbp_bar::AbstractArray{Tmsh,3}, dξdx_bar::AbstractArray{Tmsh,4},
             jac_bar::AbstractArray{Tmsh,2}, Eone_bar::AbstractArray{Tmsh,3})
  @assert( sbp.numnodes == size(xsbp_bar,2) == size(dξdx_bar,3) ==
           size(jac_bar,1) == size(Eone_bar,1) )
  @assert( size(xlag,1) == size(xlag_bar,1) == size(xref,1) == size(xsbp_bar,1)
           == size(dξdx_bar,1) == size(dξdx_bar,2) == size(Eone_bar,2) == 2 )
  @assert( size(xsbp_bar,3) == size(xlag,3) == size(xlag_bar,3) ==
           size(dξdx_bar,4) == size(jac_bar,2) == size(Eone_bar,3) )
  numdof = binomial(mapdegree+2,2)
  @assert( size(xlag,2) == size(xlag_bar,2) == size(xref,2) == numdof )
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
  fill!(xlag_bar, zero(Tmsh))
  fill!(Eone_bar, zero(Tmsh))
  coeff = zeros(Tmsh, (numdof,2))
  coeff_bar = zeros(Tmsh, (numdof,2))
  dxdξ = zeros(Tmsh, (2,2,sbp.numnodes))
  dxdξ_bar = zeros(Tmsh, (2,2,sbp.numnodes))
  dξdx = zeros(Tmsh, (2,2))
  for e = 1:size(xlag_bar,3)
    # The mapping determinant depends nonlinearly on dxdξ, so we need to
    # recompute these terms
    # find the coefficents of the polynomial mapping using xlag and Vinv
    for di = 1:2
      for i = 1:numdof
        coeff[i,di] = zero(Tmsh)
        for j = 1:numdof
          coeff[i,di] += Vinv[i,j]*xlag[di,j,e]
        end
      end
    end
    # compute the coordinate derivatives at the SBP nodes
    fill!(dxdξ, zero(Tmsh))
    for i = 1:numdof
      for di = 1:2
        for nd = 1:sbp.numnodes
          dxdξ[di,1,nd] += coeff[i,di]*dPdξ[nd,i]
          dxdξ[di,2,nd] += coeff[i,di]*dPdη[nd,i]
        end
      end
    end
    
    # start reverse sweep
    
    # compute dxdξ_bar = dξdx_bar^T*∂(dξdx)/∂(dxdξ) + jac_bar^T*∂(jac)/∂(dxdξ)
    fill!(dxdξ_bar,zero(Tmsh))
    for i = 1:sbp.numnodes
      dξdx[1,1] = dxdξ[2,2,i]
      dξdx[1,2] = -dxdξ[1,2,i]
      dξdx[2,1] = -dxdξ[2,1,i]
      dξdx[2,2] = dxdξ[1,1,i]
      jac = one(Tmsh)/(dξdx[1,1]*dξdx[2,2] - dξdx[1,2]*dξdx[2,1])
      # first, we need to account for dependence of determinant on dξdx
      fac = jac_bar[i,e]*jac*jac
      dξdx_bar[1,1,i,e] -= fac*dξdx[2,2]
      dξdx_bar[2,2,i,e] -= fac*dξdx[1,1]
      dξdx_bar[1,2,i,e] += fac*dξdx[2,1]
      dξdx_bar[2,1,i,e] += fac*dξdx[1,2]
      # now dxdξ_bar = (dξdx_bar^T + jac_bar^T*∂(jac)/∂(dξdx))*∂(dξdx)/∂(dxdξ)
      dxdξ_bar[2,2,i] += dξdx_bar[1,1,i,e]
      dxdξ_bar[1,2,i] -= dξdx_bar[1,2,i,e]
      dxdξ_bar[2,1,i] -= dξdx_bar[2,1,i,e]
      dxdξ_bar[1,1,i] += dξdx_bar[2,2,i,e]
    end
    # coeff_bar = dxdξ_bar^T*∂(dxdξ)/∂(coeff) + xsbp_bar^T*∂(xsbp)/∂(coeff)
    fill!(coeff_bar, zero(Tmsh))
    for i = 1:numdof
      for di = 1:2
        for nd = 1:sbp.numnodes
          coeff_bar[i,di] += (xsbp_bar[di,nd,e]*P[nd,i] +
                              dxdξ_bar[di,1,nd]*dPdξ[nd,i] +
                              dxdξ_bar[di,2,nd]*dPdη[nd,i])
        end
      end
    end
    # xlag_bar = coeff_bar^T*∂(coeff)/∂(xlag)
    for di = 1:2
      for i = 1:numdof
        for j = 1:numdof
          # coeff[i,di] += Vinv[i,j]*xlag[di,j,e]
          xlag_bar[di,j,e] += coeff_bar[i,di]*Vinv[i,j]
        end
      end
    end
  end # loop over elements
end

function calcMappingJacobian_rev!{
  Tsbp,Tmsh}(sbp::TetSBP{Tsbp}, mapdegree::Int, xref::AbstractArray{Tmsh,2},
             xlag::AbstractArray{Tmsh,3}, xlag_bar::AbstractArray{Tmsh,3},
             xsbp_bar::AbstractArray{Tmsh,3}, dξdx_bar::AbstractArray{Tmsh,4},
             jac_bar::AbstractArray{Tmsh,2}, Eone_bar::AbstractArray{Tmsh,3})
  @assert( sbp.numnodes == size(xsbp_bar,2) == size(dξdx_bar,3) ==
           size(jac_bar,1) == size(Eone_bar,1) )
  @assert( size(xlag,1) == size(xlag_bar,1) == size(xref,1) == size(xsbp_bar,1)
           == size(dξdx_bar,1) == size(dξdx_bar,2) == size(Eone_bar,2) == 3 )
  @assert( size(xsbp_bar,3) == size(xlag,3) == size(xlag_bar,3) ==
           size(dξdx_bar,4) == size(jac_bar,2) == size(Eone_bar,3) )
  numdof = binomial(mapdegree+3,3)
  @assert( size(xlag,2) == size(xlag_bar,2) == size(xref,2) == numdof )
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
  coeff = zeros(Tmsh, (numdof,3))
  coeff_bar = zeros(Tmsh, (numdof,3))
  dxdξ = zeros(Tmsh, (3,3,sbp.numnodes))
  dxdξ_bar = zeros(Tmsh, (3,3,sbp.numnodes))
  dξdx_targ_bar = zeros(Tmsh, (3,3,sbp.numnodes))
  Qt = zeros(Tsbp, (sbp.numnodes, 3*sbp.numnodes) )
  Qt = [sbp.Q[:,:,1].' sbp.Q[:,:,2].' sbp.Q[:,:,3].']
  Qtinv = pinv(Qt)
  targ_bar = zeros(Tmsh, (3*sbp.numnodes))
  sol_bar = zeros(targ_bar)
  # loop over each element...
  for e = 1:size(xlag,3)
    # unlike the 2D case, the jacobian terms dξdx depend nonlinearly on dxdξ, so
    # we need to recompute these quantities
    # find the coefficents of the polynomial mapping using xlag and Vinv
    for di = 1:3
      for i = 1:numdof
        coeff[i,di] = zero(Tmsh)
        for j = 1:numdof
          coeff[i,di] += Vinv[i,j]*xlag[di,j,e]
        end
      end
    end
    # compute the coordinate derivatives at the SBP nodes
    fill!(dxdξ, zero(Tmsh))
    for i = 1:numdof
      for di = 1:3
        for nd = 1:sbp.numnodes
          dxdξ[di,1,nd] += coeff[i,di]*dPdξ[nd,i]
          dxdξ[di,2,nd] += coeff[i,di]*dPdη[nd,i]
          dxdξ[di,3,nd] += coeff[i,di]*dPdζ[nd,i]
        end
      end
    end

    # Start the reverse sweep
    
    # compute dξdx_targ_bar = dξdx_bar^T ∂(dξdx)/∂(dξdx_bar) and Eone_bar
    fill!(dξdx_targ_bar, zero(Tmsh))
    for di = 1:3
      fill!(targ_bar, zero(Tmsh))
      fill!(sol_bar, zero(Tmsh))
      for di2 = 1:3
        for i = 1:sbp.numnodes
          #dξdx[di2,di,i,e] = dξdx_targ[di2,di,i] - sol[i + (di2-1)*sbp.numnodes]
          dξdx_targ_bar[di2,di,i] += dξdx_bar[di2,di,i,e]
          sol_bar[i + (di2-1)*sbp.numnodes] -= dξdx_bar[di2,di,i,e]
        end
      end
      # sol = Qtinv*b
      b_bar = Qtinv'*sol_bar
      # b = Qt*targ - Eone[:,di]
      targ_bar += Qt'*b_bar
      Eone_bar[:,di] -=  b_bar
      for di2 = 1:3
        for i = 1:sbp.numnodes
          #targ[i + (di2-1)*sbp.numnodes] = dξdx_targ[di2,di,i]
          dξdx_targ_bar[di2,di,i] += targ_bar[i + (di2-1)*sbp.numnodes]
        end
      end
    end 
    # compute dxdξ_bar = dξdx_targ_bar*∂(dξdx_targ)/∂(dxdξ) +
    # jac_bar*∂(jac)/∂(dxdξ)
    fill!(dxdξ_bar, zero(Tmsh))
    for i = 1:sbp.numnodes
      jac = one(Tmsh)/(dxdξ[1,1,i]*dxdξ[2,2,i]*dxdξ[3,3,i] +
                       dxdξ[1,2,i]*dxdξ[2,3,i]*dxdξ[3,1,i] +
                       dxdξ[1,3,i]*dxdξ[2,1,i]*dxdξ[3,2,i] -
                       dxdξ[1,1,i]*dxdξ[2,3,i]*dxdξ[3,2,i] -
                       dxdξ[1,2,i]*dxdξ[2,1,i]*dxdξ[3,3,i] -
                       dxdξ[1,3,i]*dxdξ[2,2,i]*dxdξ[3,1,i])
      jac2 = jac_bar[i,e]*jac*jac
      dxdξ_bar[1,1,i] -= jac2*(dxdξ[2,2,i]*dxdξ[3,3,i] -
                               dxdξ[2,3,i]*dxdξ[3,2,i])
      dxdξ_bar[2,2,i] -= jac2*(dxdξ[1,1,i]*dxdξ[3,3,i] -
                               dxdξ[1,3,i]*dxdξ[3,1,i])
      dxdξ_bar[3,3,i] -= jac2*(dxdξ[1,1,i]*dxdξ[2,2,i] -
                               dxdξ[1,2,i]*dxdξ[2,1,i])
      dxdξ_bar[1,2,i] -= jac2*(dxdξ[2,3,i]*dxdξ[3,1,i] -
                               dxdξ[2,1,i]*dxdξ[3,3,i])
      dxdξ_bar[2,3,i] -= jac2*(dxdξ[1,2,i]*dxdξ[3,1,i] -
                               dxdξ[1,1,i]*dxdξ[3,2,i])
      dxdξ_bar[3,1,i] -= jac2*(dxdξ[1,2,i]*dxdξ[2,3,i] -
                               dxdξ[1,3,i]*dxdξ[2,2,i])
      dxdξ_bar[1,3,i] -= jac2*(dxdξ[2,1,i]*dxdξ[3,2,i] -
                               dxdξ[2,2,i]*dxdξ[3,1,i])
      dxdξ_bar[2,1,i] -= jac2*(dxdξ[1,3,i]*dxdξ[3,2,i] -
                               dxdξ[1,2,i]*dxdξ[3,3,i])
      dxdξ_bar[3,2,i] -= jac2*(dxdξ[1,3,i]*dxdξ[2,1,i] -
                               dxdξ[1,1,i]*dxdξ[2,3,i])
      for di = 1:3
        it1 = mod(di,3)+1
        it2 = mod(di+1,3)+1
        for di2 = 1:3
          it21 = mod(di2,3)+1
          it22 = mod(di2+1,3)+1
          # dξdx_targ[di2,di,i] = (dxdξ[it1,it21,i]*dxdξ[it2,it22,i] -
          #                        dxdξ[it1,it22,i]*dxdξ[it2,it21,i])
          dxdξ_bar[it1,it21,i] += dξdx_targ_bar[di2,di,i]*dxdξ[it2,it22,i]
          dxdξ_bar[it2,it22,i] += dξdx_targ_bar[di2,di,i]*dxdξ[it1,it21,i]
          dxdξ_bar[it1,it22,i] -= dξdx_targ_bar[di2,di,i]*dxdξ[it2,it21,i]
          dxdξ_bar[it2,it21,i] -= dξdx_targ_bar[di2,di,i]*dxdξ[it1,it22,i]
        end
      end
    end
    # compute coeff_bar = dxdξ_bar*∂(dxdξ)/∂(coeff) + xsbp_bar*∂(xsbp)/∂(coeff)
    fill!(coeff_bar, zero(Tmsh))
    for i = 1:numdof
      for di = 1:3
        for nd = 1:sbp.numnodes
          # xsbp[di,nd,e] += coeff[i,di]*P[nd,i]
          coeff_bar[i,di] += xsbp_bar[di,nd,e]*P[nd,i]
          # dxdξ[di,1,nd] += coeff[i,di]*dPdξ[nd,i]
          coeff_bar[i,di] += dxdξ_bar[di,1,nd]*dPdξ[nd,i]
          # dxdξ[di,2,nd] += coeff[i,di]*dPdη[nd,i]
          coeff_bar[i,di] += dxdξ_bar[di,2,nd]*dPdη[nd,i]
          # dxdξ[di,3,nd] += coeff[i,di]*dPdζ[nd,i]
          coeff_bar[i,di] += dxdξ_bar[di,3,nd]*dPdζ[nd,i]
        end
      end
    end
    # xlag_bar = coeff_bar^T*∂(coeff)/∂(xlag)
    for di = 1:3
      for i = 1:numdof
        for j = 1:numdof
          # coeff[i,di] += Vinv[i,j]*xlag[di,j,e]
          xlag_bar[di,j,e] += Vinv[i,j]*coeff_bar[i,di]
        end
      end
    end
  end
end

@doc """
### SummationByParts.mappingjacobian_rev!

**Deprecated**:

Forms the reverse-mode of algorithmic differentiation product for the method
`mappingjacobian!`.  Specifically, it computes `x_bar` = `dξdx_bar`^T *
∂`dξdx`/∂`x` + `jac_bar`^T * ∂`jac`/∂`x` where the various quantities are
defined below (the `_bar` denotes the reverse mode).  Note that the input and
output order of the arguments follows that used in `mappingjacobian!`.

**Inputs**

* `sbp`: an SBP operator type
* `x`: the physical coordinates; 1st dim = coord, 2nd dim = node, 3rd dim = elem
* `dξdx_bar`: gradient w.r.t. Jacobian; [ref coord, phys coord, sbp node, element]
* `jac_bar`: gradient w.r.t. determinant of the Jacobian; [sbp node, element]

**In/Outs**

* `x_bar`: gradient w.r.t. SBP nodes; [coord, Lagrangian node, element]

  """->
function mappingjacobian_rev!{Tsbp,Tmsh}(sbp::TriSBP{Tsbp},
                                         x::AbstractArray{Tmsh,3},
                                         x_bar::AbstractArray{Tmsh,3},
                                         dξdx_bar::AbstractArray{Tmsh,4},
                                         jac_bar::AbstractArray{Tmsh,2})
  @assert( sbp.numnodes == size(x,2) == size(x_bar,2) == size(dξdx_bar,3)
           == size(jac_bar,1) )
  @assert( size(x,3) == size(x_bar,3) == size(dξdx_bar,4) == size(jac_bar,2) )
  @assert( size(x,1) == size(x_bar,1) == size(dξdx_bar,1) == size(dξdx_bar,2)
           == 2 )
  work = zeros(Tmsh, (2,sbp.numnodes,2))
  dξdx = zeros(Tmsh, (2,2,sbp.numnodes))
  for elem = 1:size(x,3)
    # compute the coordinate derivatives
    fill!(work, zero(Tmsh))
    for di = 1:2
      differentiateElement!(sbp, di, sub(x,:,:,elem), sub(work,:,:,di)) 
    end
    # compute the scaled metrics: could also pass these in to avoid recomputing
    # them...
    for i = 1:sbp.numnodes
      dξdx[1,1,i] = work[2,i,2]
      dξdx[1,2,i] = -work[1,i,2]
      dξdx[2,1,i] = -work[2,i,1]
      dξdx[2,2,i] = work[1,i,1]
    end    
    # start the reverse sweep
    for i = 1:sbp.numnodes
      jac = one(Tmsh)/(dξdx[1,1,i]*dξdx[2,2,i] - dξdx[1,2,i]*dξdx[2,1,i])
      jac2 = jac_bar[i,elem]*jac*jac
      dξdx_bar[1,1,i,elem] -= jac2*dξdx[2,2,i]
      dξdx_bar[1,2,i,elem] += jac2*dξdx[2,1,i]
      dξdx_bar[2,1,i,elem] += jac2*dξdx[1,2,i]
      dξdx_bar[2,2,i,elem] -= jac2*dξdx[1,1,i]
    end
    fill!(work, zero(Tmsh))
    for i = 1:sbp.numnodes
      work[2,i,2] += dξdx_bar[1,1,i,elem]
      work[1,i,2] -= dξdx_bar[1,2,i,elem]
      work[2,i,1] -= dξdx_bar[2,1,i,elem]
      work[1,i,1] += dξdx_bar[2,2,i,elem]
    end
    for di = 1:2
      differentiateElement_rev!(sbp, di, sub(x_bar,:,:,elem), sub(work,:,:,di)) 
    end
  end
end


function mappingjacobian_rev!{Tsbp,Tmsh}(sbp::TetSBP{Tsbp},
                                         x::AbstractArray{Tmsh,3},
                                         x_bar::AbstractArray{Tmsh,3},
                                         dξdx_bar::AbstractArray{Tmsh,4},
                                         jac_bar::AbstractArray{Tmsh,2})
  @assert( sbp.numnodes == size(x,2) == size(x_bar,2) == size(dξdx_bar,3)
           == size(jac_bar,1) )
  @assert( size(x,3) == size(x_bar,3) == size(dξdx_bar,4) == size(jac_bar,2) )
  @assert( size(x,1) == size(x_bar,1) == size(dξdx_bar,1) == size(dξdx_bar,2)
           == 3 )
  work = zeros(Tmsh, (3,sbp.numnodes,3))
  dxdξ = zeros(Tmsh, (3,3,sbp.numnodes))
  dxdξ_bar = zeros(dxdξ)
  for elem = 1:size(x,3)
    # compute the coordinate derivatives
    fill!(work, zero(Tmsh))
    for di = 1:3
      differentiateElement!(sbp, di, sub(x,:,:,elem), sub(work,:,:,di)) 
    end
    permutedims!(dxdξ, work, [1,3,2]) # probably slow

    # Start the reverse sweep
    fill!(dxdξ_bar, zero(Tmsh))
    
    # Get the contribution of the Jacobian product to dxdξ
    for i = 1:sbp.numnodes
      jac = one(Tmsh)/(dxdξ[1,1,i]*dxdξ[2,2,i]*dxdξ[3,3,i] +
                       dxdξ[1,2,i]*dxdξ[2,3,i]*dxdξ[3,1,i] +
                       dxdξ[1,3,i]*dxdξ[2,1,i]*dxdξ[3,2,i] -
                       dxdξ[1,1,i]*dxdξ[2,3,i]*dxdξ[3,2,i] -
                       dxdξ[1,2,i]*dxdξ[2,1,i]*dxdξ[3,3,i] -
                       dxdξ[1,3,i]*dxdξ[2,2,i]*dxdξ[3,1,i])
      jac2 = jac_bar[i,elem]*jac*jac
      dxdξ_bar[1,1,i] -= jac2*(dxdξ[2,2,i]*dxdξ[3,3,i] -
                               dxdξ[2,3,i]*dxdξ[3,2,i])
      dxdξ_bar[2,2,i] -= jac2*(dxdξ[1,1,i]*dxdξ[3,3,i] -
                               dxdξ[1,3,i]*dxdξ[3,1,i])
      dxdξ_bar[3,3,i] -= jac2*(dxdξ[1,1,i]*dxdξ[2,2,i] -
                               dxdξ[1,2,i]*dxdξ[2,1,i])
      dxdξ_bar[1,2,i] -= jac2*(dxdξ[2,3,i]*dxdξ[3,1,i] -
                               dxdξ[2,1,i]*dxdξ[3,3,i])
      dxdξ_bar[2,3,i] -= jac2*(dxdξ[1,2,i]*dxdξ[3,1,i] -
                               dxdξ[1,1,i]*dxdξ[3,2,i])
      dxdξ_bar[3,1,i] -= jac2*(dxdξ[1,2,i]*dxdξ[2,3,i] -
                               dxdξ[1,3,i]*dxdξ[2,2,i])
      dxdξ_bar[1,3,i] -= jac2*(dxdξ[2,1,i]*dxdξ[3,2,i] -
                               dxdξ[2,2,i]*dxdξ[3,1,i])
      dxdξ_bar[2,1,i] -= jac2*(dxdξ[1,3,i]*dxdξ[3,2,i] -
                               dxdξ[1,2,i]*dxdξ[3,3,i])
      dxdξ_bar[3,2,i] -= jac2*(dxdξ[1,3,i]*dxdξ[2,1,i] -
                               dxdξ[1,1,i]*dxdξ[2,3,i])
    end
    
    # contribution of metrics to dxdξ
    for di = 1:3
      it1 = mod(di,3)+1
      it2 = mod(di+1,3)+1
      for i = 1:sbp.numnodes
        for di2 = 1:3
          it21 = mod(di2,3)+1
          it22 = mod(di2+1,3)+1
          #dξdx[it1,di2,i,elem] += (dxdξ[it22,i,elem,di]*dxdξ[it21,i,elem,it2] -
          #                         dxdξ[it21,i,elem,di]*dxdξ[it22,i,elem,it2])
          dxdξ_bar[it22,di,i] += dxdξ[it21,it2,i]*dξdx_bar[it1,di2,i,elem]
          dxdξ_bar[it21,it2,i] += dxdξ[it22,di,i]*dξdx_bar[it1,di2,i,elem]
          dxdξ_bar[it21,di,i] -= dxdξ[it22,it2,i]*dξdx_bar[it1,di2,i,elem]
          dxdξ_bar[it22,it2,i] -= dxdξ[it21,di,i]*dξdx_bar[it1,di2,i,elem]
        end
      end # i loop
    end # di loop

    permutedims!(work, dxdξ_bar, [1,3,2]) # probably slow    
    for di = 1:3
      #differentiate!(sbp, di, x, sub(dxdξ,:,:,:,di))
      differentiateElement_rev!(sbp, di, sub(x_bar,:,:,elem),
                                sub(work,:,:,di))
    end
  end # loop over elements
end
