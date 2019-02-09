# This file gathers together functions related to forming Jacobians of "weak"
# differentiation using the SBP operators

"""
### SummationByParts.weakDifferentiateElement_jac!

This function applies the SBP stiffness matrix (or its transpose) to a diagonal,
or block diagonal, matrix of flux jacobians `dfluxdu`.  The resulting matrix is
added to or subtracted from `dresdu`.  Different methods are available depending
on the rank of `dfluxdu`:

* For *scalar* fields, it is assumed that `dfluxdu` is a rank-1 array, with the
only dimension for the local-node index.  `dresdu` is a rank-2 array, with size
`sbp.numnodes` x `sbp.numnodes`.
* For *vector* fields, `dfluxdu` is a rank-3 array, with the first and second
dimensions for indices of the vector field (at a particular node) and the second
dimension for the local-node index.  `dresdu` is a rank-4 array; the first 2
dimensions are of size nvar = size(dfluxdu,1), while the third and fourth
dimensions are of size `sbp.numnodes`.

Naturally, the number of entries in the dimension of `dfluxdu` (and `dresdu`)
corresponding to the nodes must be equal to the number of nodes in the SBP
operator sbp.

**Inputs**

* `sbp`: an SBP operator type
* `di`: direction index of the operator that is desired (di=1 for Qx, etc)
* `dfluxdu`: array of the derivative of the flux w.r.t. the state
* `±` : PlusFunctor to add to res, MinusFunctor to subract
* `trans` (optional): if true, the transpose operation is applied

**In/Outs**

* `dresdu`: stores Q[:,:,di]*diag(dfluxdu) or Q[:,:,di]'*diag(dfluxdu)

"""
function weakDifferentiateElement_jac!{
  Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                  dfluxdu::AbstractArray{Tflx,1},
                  dresdu::AbstractArray{Tres,2},
                  (±)::UnaryFunctor=Add(),
                  trans::Bool=false)
  @asserts_enabled begin
    @assert( sbp.numnodes == size(dfluxdu,1) == size(dresdu,1) ==
             size(dresdu,2) )
    @assert( di > 0 && di <= size(sbp.Q,3) )
  end
  if trans # apply transposed Q
    for j = 1:sbp.numnodes
      for i = 1:sbp.numnodes                
        # res[i] += ±(sbp.Q[j,i,di]*flux[j])
        dresdu[i,j] += ±(sbp.Q[j,i,di]*dfluxdu[j])
      end
    end
  else # apply Q
    for j = 1:sbp.numnodes    
      for i = 1:sbp.numnodes
        # res[i] += ±(sbp.Q[i,j,di]*flux[j])
        dresdu[i,j] += ±(sbp.Q[i,j,di]*dfluxdu[j])
      end
    end
  end
end

function weakDifferentiateElement_jac!{
  Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                  dfluxdu::AbstractArray{Tflx,3},
                  dresdu::AbstractArray{Tres,4},
                  (±)::UnaryFunctor=Add(),
                  trans::Bool=false)
  @asserts_enabled begin
    @assert( size(dfluxdu,1) == size(dfluxdu,2) == size(dresdu,1) ==
             size(dresdu,2) )
    @assert( sbp.numnodes == size(dfluxdu,3) == size(dresdu,3) ==
             size(dresdu,4) )
    @assert( di > 0 && di <= size(sbp.Q,3) )
  end
  if trans # apply transposed Q
    for j = 1:sbp.numnodes
      for i = 1:sbp.numnodes
        for q = 1:size(dfluxdu,2)
          for p = 1:size(dfluxdu,1)
            # res[field,i] += ±(sbp.Q[j,i,di]*flux[field,j])
            dresdu[p,q,i,j] += ±(sbp.Q[j,i,di]*dfluxdu[p,q,j])
          end
        end
      end
    end
  else # apply Q
    for j = 1:sbp.numnodes
      for i = 1:sbp.numnodes
        for q = 1:size(dfluxdu,2)
          for p = 1:size(dfluxdu,1)
            # res[field,i] += ±(sbp.Q[i,j,di]*flux[field,j])
            dresdu[p,q,i,j] += ±(sbp.Q[i,j,di]*dfluxdu[p,q,j])
          end
        end
      end
    end
  end
end
