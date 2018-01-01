# This file gathers together methods related to computing the Jacobian of
# face-based integral terms; it combines operations from both faceinterpolate
# and faceintegrate.

@doc """
### SummationByParts.interiorFaceIntegrate_jac!

Given the face-node flux Jacobians, this method computes the Jacobian of the
chain interiorFaceInterpolate! --> face-flux Jacobian --> interiorFaceIntegrate!
and adds the contributions to the adjacent elements' Jacobians matrices.
Different methods are available depending on the rank of `dfluxduL` and
`dfluxduR`:

* For *scalar* fields, it is assumed that `dfluxduL` and `dfluxduR` are rank-1
arrays, with the only dimension for the face local-node index.  `dresLduL` is the Jacobian of the left-element residual with respect to the left-element solution, and similarly for `dresLduR`, `dresRduL`, `dresRduR`.
* For *vector* fields, `dfluxduL` and `dfluxduR` are rank-3 arrays, with the
first and second dimensions for indices of the vector field (at a particular
face node) and the second dimension for the face local-node index.  `dresLduL`
is a rank-4 array that is the Jacobian of the left-element residual with respect
to the left-element solution; the first 2 dimensions are of size nvar =
size(dfluxduL,1), while the third and fourth dimensions are of size
`sbp.numnodes`.  Similarly for `dresLduR`, `dresRduL`, `dresRduR`.

**Inputs**

* `sbpface`: an SBP AbstractFace type
* `ifaces`: list of element interfaces stored as an array of `Interface`s
* `dfluxduL`: array of the derivative of the flux w.r.t. the left trace
* `dfluxduR`: array of the derivative of the flux w.r.t. the right trace
* `±`: PlusFunctor to add to res, MinusFunctor to subract

**In/Outs**

* `dresLduL`: Jacobian of left residual w.r.t. left state
* `dresLduR`: Jacobian of left residual w.r.t. right state
* `dresRduL`: Jacobian of right residual w.r.t. left state
* `dresRduR`: Jacobian of right residual w.r.t. right state

"""->
function interiorFaceIntegrate_jac!{
  Tsbp,Tflx,Tres}(sbpface::DenseFace{Tsbp}, iface::Interface,
                  dfluxduL::AbstractArray{Tflx,1},
                  dfluxduR::AbstractArray{Tflx,1},
                  dresLduL::AbstractArray{Tres,2},
                  dresLduR::AbstractArray{Tres,2},
                  dresRduL::AbstractArray{Tres,2},
                  dresRduR::AbstractArray{Tres,2},
                  (±)::UnaryFunctor=Add())
  @asserts_enabled begin
    @assert( size(sbpface.interp,1) <= size(dresLduL,1) )
    @assert( size(sbpface.interp,2) == size(dfluxduL,1) == size(dfluxduR,1) )
    @assert( size(dresLduL,1) == size(dresLduL,2) ==
             size(dresLduR,1) == size(dresLduR,2) ==
             size(dresRduL,1) == size(dresRduL,2) ==
             size(dresRduR,1) == size(dresRduR,2) )
  end
  
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    for j1 = 1:sbpface.stencilsize
      # This is the loop over the columns
      colL = sbpface.perm[j1,iface.faceL]
      colR = sbpface.perm[j1,iface.faceR]      
      for j2 = 1:sbpface.stencilsize
        # This is the loop over the rows
        rowL = sbpface.perm[j2,iface.faceL]
        rowR = sbpface.perm[j2,iface.faceR]

        dresLduL[rowL,colL] += ±(sbpface.interp[j2,i]*sbpface.wface[i]
                                 *dfluxduL[i]*sbpface.interp[j1,i])
        dresLduR[rowL,colR] += ±(sbpface.interp[j2,i]*sbpface.wface[i]
                                 *dfluxduR[i]*sbpface.interp[j1,iR])
        dresRduR[rowR,colR] -= ±(sbpface.interp[j2,iR]*sbpface.wface[iR]
                                 *dfluxduR[i]*sbpface.interp[j1,iR])
        dresRduL[rowR,colL] -= ±(sbpface.interp[j2,iR]*sbpface.wface[iR]
                                 *dfluxduL[i]*sbpface.interp[j1,i])
      end
    end
  end
end

function interiorFaceIntegrate_jac!{
  Tsbp,Tflx,Tres}(sbpface::DenseFace{Tsbp}, iface::Interface,
                  dfluxduL::AbstractArray{Tflx,3},
                  dfluxduR::AbstractArray{Tflx,3},
                  dresLduL::AbstractArray{Tres,4},
                  dresLduR::AbstractArray{Tres,4},
                  dresRduL::AbstractArray{Tres,4},
                  dresRduR::AbstractArray{Tres,4},
                  (±)::UnaryFunctor=Add())
  @asserts_enabled begin
    @assert( size(sbpface.interp,1) <= size(dresLduL,3) )
    @assert( size(sbpface.interp,2) == size(dfluxduL,3) == size(dfluxduR,3) )
    @assert( size(dresLduL,3) == size(dresLduL,4) ==
             size(dresLduR,3) == size(dresLduR,4) ==
             size(dresRduL,3) == size(dresRduL,4) ==
             size(dresRduR,3) == size(dresRduR,4) )
    @assert( size(dfluxduL,1) == size(dfluxduL,2) ==
             size(dfluxduR,1) == size(dfluxduR,2) ==
             size(dresLduL,1) == size(dresLduL,2) ==
             size(dresLduR,1) == size(dresLduR,2) ==
             size(dresRduL,1) == size(dresRduL,2) ==
             size(dresRduR,1) == size(dresRduR,2) )
  end
  
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    for j1 = 1:sbpface.stencilsize
      # This is the loop over the columns
      colL = sbpface.perm[j1,iface.faceL]
      colR = sbpface.perm[j1,iface.faceR]      
      for j2 = 1:sbpface.stencilsize
        # This is the loop over the rows
        rowL = sbpface.perm[j2,iface.faceL]
        rowR = sbpface.perm[j2,iface.faceR]
        
        for q = 1:size(dfluxduL,2)
          for p = 1:size(dfluxduL,1)
            dresLduL[p,q,rowL,colL] += ±(sbpface.interp[j2,i]*sbpface.wface[i]
                                         *dfluxduL[p,q,i]*sbpface.interp[j1,i])
            dresLduR[p,q,rowL,colR] += ±(sbpface.interp[j2,i]*sbpface.wface[i]
                                         *dfluxduR[p,q,i]*sbpface.interp[j1,iR])
            dresRduR[p,q,rowR,colR] -= ±(sbpface.interp[j2,iR]*sbpface.wface[iR]
                                         *dfluxduR[p,q,i]*sbpface.interp[j1,iR])
            dresRduL[p,q,rowR,colL] -= ±(sbpface.interp[j2,iR]*sbpface.wface[iR]
                                         *dfluxduL[p,q,i]*sbpface.interp[j1,i])
          end
        end
      end
    end
  end
end

function interiorFaceIntegrate_jac!{
  Tsbp,Tflx,Tres}(sbpface::SparseFace{Tsbp}, iface::Interface,
                  dfluxduL::AbstractArray{Tflx,1},
                  dfluxduR::AbstractArray{Tflx,1},
                  dresLduL::AbstractArray{Tres,2},
                  dresLduR::AbstractArray{Tres,2},
                  dresRduL::AbstractArray{Tres,2},
                  dresRduR::AbstractArray{Tres,2},
                  (±)::UnaryFunctor=Add())
  @asserts_enabled begin
    @assert( size(dresLduL,1) == size(dresLduL,2) ==
             size(dresLduR,1) == size(dresLduR,2) ==
             size(dresRduL,1) == size(dresRduL,2) ==
             size(dresRduR,1) == size(dresRduR,2) )
  end
  
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    colL = sbpface.perm[i,iface.faceL]
    colR = sbpface.perm[iR,iface.faceR]
    rowL = colL
    rowR = colR
    dresLduL[rowL,colL] += ±(sbpface.wface[i]*dfluxduL[i])
    dresLduR[rowL,colR] += ±(sbpface.wface[i]*dfluxduR[i])
    dresRduR[rowR,colR] -= ±(sbpface.wface[iR]*dfluxduR[i])
    dresRduL[rowR,colL] -= ±(sbpface.wface[iR]*dfluxduL[i])
  end
end

function interiorFaceIntegrate_jac!{
  Tsbp,Tflx,Tres}(sbpface::SparseFace{Tsbp}, iface::Interface,
                  dfluxduL::AbstractArray{Tflx,3},
                  dfluxduR::AbstractArray{Tflx,3},
                  dresLduL::AbstractArray{Tres,4},
                  dresLduR::AbstractArray{Tres,4},
                  dresRduL::AbstractArray{Tres,4},
                  dresRduR::AbstractArray{Tres,4},
                  (±)::UnaryFunctor=Add())
  @asserts_enabled begin
    @assert( size(dresLduL,3) == size(dresLduL,4) ==
             size(dresLduR,3) == size(dresLduR,4) ==
             size(dresRduL,3) == size(dresRduL,4) ==
             size(dresRduR,3) == size(dresRduR,4) )
    @assert( size(dfluxduL,1) == size(dfluxduL,2) ==
             size(dfluxduR,1) == size(dfluxduR,2) ==
             size(dresLduL,1) == size(dresLduL,2) ==
             size(dresLduR,1) == size(dresLduR,2) ==
             size(dresRduL,1) == size(dresRduL,2) ==
             size(dresRduR,1) == size(dresRduR,2) )
  end
  
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    colL = sbpface.perm[i,iface.faceL]
    colR = sbpface.perm[iR,iface.faceR]      
    rowL = colL
    rowR = colR   
    for q = 1:size(dfluxduL,2)
      for p = 1:size(dfluxduL,1)
        dresLduL[p,q,rowL,colL] += ±(sbpface.wface[i]*dfluxduL[p,q,i])
        dresLduR[p,q,rowL,colR] += ±(sbpface.wface[i]*dfluxduR[p,q,i])
        dresRduR[p,q,rowR,colR] -= ±(sbpface.wface[iR]*dfluxduR[p,q,i])
        dresRduL[p,q,rowR,colL] -= ±(sbpface.wface[iR]*dfluxduL[p,q,i])
      end
    end
  end
end

