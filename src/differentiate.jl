# This file gathers together functions related to strong differentiation using
# the SBP operators

"""
### SummationByParts.differentiate!

Applies the SBP differentiation matrix operator, D, to data in `flux` and stores
the result in `res`.  Different methods are available depending on the rank of
`flux`:

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

"""
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

"""
### SummationByParts.differentiateElement!

This is the single-element variant of differentiate!  Applies the SBP
differentiation matrix operator, D, to data in `flux` and stores the result in
`res`.  Different methods are available depending on the rank of `flux`:

* For *scalar* fields, it is assumed that `flux` is a rank-1 array, with first
and only dimension for the local-node index.
* For *vector* fields, `flux` is a rank-2 array, with the first dimension for
the index of the vector field, and the second dimension for the local-node index.

Naturally, the number of entries in the dimension of `flux` (and `res`)
corresponding to the nodes must be equal to the number of nodes in the SBP
operator sbp.

**Inputs**

* `sbp`: an SBP operator type
* `di`: direction index of the operator that is desired (di=1 for Dx, etc)
* `flux`: the array that the operator is applied to
* `±`: PlusFunctor to add to res, MinusFunctor to subract
* `trans` (optional): if true, the transpose operation is applied

**In/Outs**

* `res`: where the result of applying inv(H)*Q[:,:,di] to u is stored

"""
function differentiateElement!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                                               flux::AbstractArray{Tflx,1},
                                               res::AbstractArray{Tres,1},
                                               (±)::UnaryFunctor=Add(),
                                               trans::Bool=false)
  @asserts_enabled begin
    @assert( sbp.numnodes == size(flux,1) == size(res,1) )
    @assert( di > 0 && di <= size(sbp.Q,3) )
  end
  
  if trans # apply transposed D
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        res[i] += ±(sbp.Q[j,i,di]*flux[j]/sbp.w[j])
      end
    end    
  else # apply D
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        res[i] += ±(sbp.Q[i,j,di]*flux[j])
      end
      res[i] /= sbp.w[i]
    end
  end
end

function differentiateElement!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                                               flux::AbstractArray{Tflx,2},
                                               res::AbstractArray{Tres,2},
                                               (±)::UnaryFunctor=Add(),
                                               trans::Bool=false)
  @asserts_enabled begin
    @assert( sbp.numnodes == size(flux,2) == size(res,2) )
    @assert( length(flux) == length(res) )
    @assert( di > 0 && di <= size(sbp.Q,3) )
  end

  if trans # apply transposed D
#    Hinvflux = zeros(Tres, (size(flux,1)))
    for i = 1:sbp.numnodes
      Hinv = 1./sbp.w[i]
#      for field = 1:size(flux,1)
#        Hinvflux[field] = Hinv*flux[field,i]
#      end
      for j = 1:sbp.numnodes
        for field = 1:size(flux,1)
#          res[field,j] += ±(sbp.Q[i,j,di]*Hinvflux[field])
          res[field,j] += ±(sbp.Q[i,j,di]*Hinv*flux[field, i])
        end
      end
    end
  else # apply D
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        for field = 1:size(flux,1)
          res[field,i] += ±(sbp.Q[i,j,di]*flux[field,j])
        end
      end
      Hinv = 1./sbp.w[i]
      for field = 1:size(flux,1)
        res[field,i] *= Hinv
      end
    end
  end
end
