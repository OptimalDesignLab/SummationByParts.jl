# This file gathers together functions related to volumne integration over a
# test function using SBP operators

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
### SummationByParts.volumeIntegrateElement!

This is the single-element variant of volumeIntegrate!.  Applies the SBP mass
matrix operator, H, to data in `u` and stores the result in `res`.  Different
methods are available depending on the rank of `u`:

* For *scalar* fields, it is assumed that `u` is a rank-1 array, with the first
and only dimension for the local-node index.
* For *vector* fields, `u` is a rank-2 array, with the first dimension for the
index of the vector field, and the second dimension for the local-node index.

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
function volumeIntegrateElement!{Tsbp,Tsol,Tres}(sbp::AbstractSBP{Tsbp},
                                                 u::AbstractArray{Tsol,1},
                                                 res::AbstractArray{Tres,1},
                                                 (±)::UnaryFunctor=Add())
  @asserts_enabled begin
    @assert( sbp.numnodes == size(u,1) == size(res,1) )
  end

  for i = 1:sbp.numnodes
    res[i] += ±(sbp.w[i]*u[i])
  end
end

function volumeIntegrateElement!{Tsbp,Tsol,Tres}(sbp::AbstractSBP{Tsbp},
                                                 u::AbstractArray{Tsol,2},
                                                 res::AbstractArray{Tres,2},
                                                 (±)::UnaryFunctor=Add())
  @asserts_enabled begin
    @assert( sbp.numnodes == size(u,2) == size(res,2) )
    @assert( length(u) == length(res) )
  end

  for i = 1:sbp.numnodes
    for field = 1:size(u,1)
      res[field,i] += ±(sbp.w[i]*u[field,i])
    end
  end
end
