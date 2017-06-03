# This file gathers together functions related to "weak" differentiation using
# the SBP operators

@doc """
### SummationByParts.weakdifferentiate!

Applies the SBP stiffness matrix (or its transpose) to data in `flux` and adds
to or subtracts from `res`.  Different methods are available depending on the
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
### SummationByParts.weakDifferentiateElement!

This is the single-element variant of weakdifferentiate!.  Applies the SBP
stiffness matrix (or its transpose) to data in `flux` and adds to or subtracts
from `res`.  Different methods are available depending on the rank of `flux`:

* For *scalar* fields, it is assumed that `flux` is a rank-1 array, with the only
dimension for the local-node index.
* For *vector* fields, `flux` is a rank-2 array, with the first dimension for
the index of the vector field, and the second dimension for the local-node index.

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
function weakDifferentiateElement!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp},
                                                   di::Int,
                                                   flux::AbstractArray{Tflx,1},
                                                   res::AbstractArray{Tres,1},
                                                   (±)::UnaryFunctor=Add();
                                                   trans::Bool=false)
  @assert( sbp.numnodes == size(flux,1) == size(res,1) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  if trans # apply transposed Q
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        res[i] += ±(sbp.Q[j,i,di]*flux[j])
      end
    end
  else # apply Q
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        res[i] += ±(sbp.Q[i,j,di]*flux[j])
      end
    end
  end
end

function weakDifferentiateElement!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp},
                                                   di::Int,
                                                   flux::AbstractArray{Tflx,2},
                                                   res::AbstractArray{Tres,2},
                                                   (±)::UnaryFunctor=Add(),
                                                   trans::Bool=false)
#  @assert( sbp.numnodes == size(flux,2) == size(res,2) )
#  @assert( length(flux) == length(res) )
#  @assert( di > 0 && di <= size(sbp.Q,3) )
  if trans # apply transposed Q
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        for field = 1:size(flux,1)
          res[field,i] += ±(sbp.Q[j,i,di]*flux[field,j])
        end
      end
    end
  else # apply Q
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        for field = 1:size(flux,1)
          res[field,i] += ±(sbp.Q[i,j,di]*flux[field,j])
        end
      end
    end
  end
end
