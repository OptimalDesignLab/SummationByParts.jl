# This file contains the reverse-mode version of the methods in differentiate.jl

@doc """
### SummationByParts.differentiate_rev!

This is the reverse differentiated version of differentiate!.  See
differentiate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `flux` variable.

**Inputs**

* `sbp`: an SBP operator type
* `di`: direction index of the operator that is desired (di=1 for Dx, etc)
* `res_bar`: vector applied to the left of the D operator
* `±`: PlusFunctor to add to res_bar, MinusFunctor to subract

**In/Outs**

* `flux_bar`: the result of the vector matrix product between D and `res_bar`

"""->
function differentiate_rev!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                                            flux_bar::AbstractArray{Tflx,2},
                                            res_bar::AbstractArray{Tres,2},
                                            (±)::UnaryFunctor=Add())
  @assert( sbp.numnodes == size(flux_bar,1) && sbp.numnodes == size(res_bar,1) )
  @assert( length(flux_bar) == length(res_bar) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  tmp_bar = zero(Tres)
  for elem = 1:size(flux_bar,2)
    for i = 1:sbp.numnodes
      #res[i,elem] *= Hinv[i] 
      tmp_bar = res_bar[i,elem]./sbp.w[i] # this is to avoid overwritting res_bar
      for j = 1:sbp.numnodes
        #res[i,elem] += ±(sbp.Q[i,j,di]*flux[j,elem])
        flux_bar[j,elem] += ±(sbp.Q[i,j,di]*tmp_bar)
      end
    end
  end
end

function differentiate_rev!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                                        flux_bar::AbstractArray{Tflx,3},
                                        res_bar::AbstractArray{Tres,3},
                                        (±)::UnaryFunctor=Add())
  @assert( sbp.numnodes == size(flux_bar,2) && sbp.numnodes == size(res_bar,2) )
  @assert( length(flux_bar) == length(res_bar) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  Hinv = zero(Tsbp)
  tmp_bar = zeros(Tres, size(flux_bar,1))
  for elem = 1:size(flux_bar,3)
    for i = 1:sbp.numnodes
      Hinv = 1./sbp.w[i]
      for field = 1:size(flux_bar,1)
        tmp_bar[field] = Hinv*res_bar[field,i,elem]
      end      
      for j = 1:sbp.numnodes
        for field = 1:size(flux_bar,1)
          #res[field,i,elem] += ±(sbp.Q[i,j,di]*flux[field,j,elem])
          flux_bar[field,j,elem] += ±(sbp.Q[i,j,di]*tmp_bar[field])
        end
      end
    end
  end
end

@doc """
### SummationByParts.differentiateElement_rev!

This is the reverse differentiated version of differentiateElement!.  See
differentiate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `flux` variable.

**Inputs**

* `sbp`: an SBP operator type
* `di`: direction index of the operator that is desired (di=1 for Dx, etc)
* `res_bar`: vector applied to the left of the D operator
* `±`: PlusFunctor to add to res_bar, MinusFunctor to subract

**In/Outs**

* `flux_bar`: the result of the vector matrix product between D and `res_bar`

"""->
function differentiateElement_rev!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                                               flux_bar::AbstractArray{Tflx,1},
                                               res_bar::AbstractArray{Tres,1},
                                               (±)::UnaryFunctor=Add())
  @assert( sbp.numnodes == size(flux_bar,1) == size(res_bar,1) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  tmp_bar = zero(Tres)
  for i = 1:sbp.numnodes
    # res[i] /= sbp.w[i]    
    tmp_bar = res_bar[i]./sbp.w[i]    
    for j = 1:sbp.numnodes
      #res[i] += ±(sbp.Q[i,j,di]*flux[j])
      flux_bar[j] += ±(sbp.Q[i,j,di]*tmp_bar)
    end
  end
end

function differentiateElement_rev!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                                               flux_bar::AbstractArray{Tflx,2},
                                               res_bar::AbstractArray{Tres,2},
                                               (±)::UnaryFunctor=Add())
  @assert( sbp.numnodes == size(flux_bar,2) == size(res_bar,2) )
  @assert( length(flux_bar) == length(res_bar) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  Hinv = zero(Tsbp)
  tmp_bar = zeros(Tres, size(flux_bar,1))
  for i = 1:sbp.numnodes
    Hinv = 1./sbp.w[i]
    for field = 1:size(flux_bar,1)
      #res[field,i] *= Hinv
      tmp_bar[field] = Hinv*res_bar[field,i]
    end    
    for j = 1:sbp.numnodes
      for field = 1:size(flux_bar,1)
        #res[field,i] += ±(sbp.Q[i,j,di]*flux[field,j])
        flux_bar[field,j] += ±(sbp.Q[i,j,di]*tmp_bar[field])
      end
    end
  end
end
