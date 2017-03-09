# This file contains the reverse-mode version of the methods in
# weakdifferentiate.jl

@doc """
### SummationByParts.weakdifferentiate_rev!

This is the reverse differentiated version of weakdifferentiate!.  See
weakdifferentiate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `flux` variable.

**Inputs**

* `sbp`: an SBP operator type
* `di`: direction index of the operator that is desired (di=1 for Qx, etc)
* `res_bar`: vector applied to the left of the Q operator
* `±` : PlusFunctor to add to res_bar, MinusFunctor to subract
* `trans` (optional): if true, the transpose operation is applied

**In/Outs**

* `flux_bar`: the result of the vector matrix product between Q and `res_bar`

"""->
function weakdifferentiate_rev!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int, 
                                                flux_bar::AbstractArray{Tflx,2},
                                                res_bar::AbstractArray{Tres,2},
                                                (±)::UnaryFunctor=Add();
                                                trans::Bool=false)
  @assert( sbp.numnodes == size(flux_bar,1) && sbp.numnodes == size(res_bar,1) )
  @assert( length(flux_bar) == length(res_bar) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  if trans # apply transposed Q
    for elem = 1:size(flux_bar,2)
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          #res[i,elem] += ±(sbp.Q[j,i,di]*flux[j,elem])
          flux_bar[j,elem] += ±(sbp.Q[j,i,di]*res_bar[i,elem])
        end
      end
    end
  else # apply Q
    for elem = 1:size(flux_bar,2)
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          #res[i,elem] += ±(sbp.Q[i,j,di]*flux[j,elem])
          flux_bar[j,elem] += ±(sbp.Q[i,j,di]*res_bar[i,elem])
        end
      end
    end
  end
end

function weakdifferentiate_rev!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                                                flux_bar::AbstractArray{Tflx,3},
                                                res_bar::AbstractArray{Tres,3},
                                                (±)::UnaryFunctor=Add();
                                                trans::Bool=false)
  @assert( sbp.numnodes == size(flux_bar,2) && sbp.numnodes == size(res_bar,2) )
  @assert( length(flux_bar) == length(res_bar) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  if trans # apply transposed Q
    for elem = 1:size(flux_bar,3)
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          for field = 1:size(flux_bar,1)
            #res[field,i,elem] += ±(sbp.Q[j,i,di]*flux[field,j,elem])
            flux_bar[field,j,elem] += ±(sbp.Q[j,i,di]*res_bar[field,i,elem])
          end
        end
      end
    end
  else # apply Q
    for elem = 1:size(flux_bar,3)
      for i = 1:sbp.numnodes
        for j = 1:sbp.numnodes
          for field = 1:size(flux_bar,1)
            #res[field,i,elem] += ±(sbp.Q[i,j,di]*flux[field,j,elem])
            flux_bar[field,j,elem] += ±(sbp.Q[i,j,di]*res_bar[field,i,elem])
          end
        end
      end
    end
  end
end

@doc """
### SummationByParts.weakDifferentiateElement_rev!

This is the reverse differentiated version of weakDifferentiateElement!.  See
weakdifferentiate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `flux` variable.

**Inputs**

* `sbp`: an SBP operator type
* `di`: direction index of the operator that is desired (di=1 for Qx, etc)
* `res_bar`: vector applied to the left of the Q operator
* `±` : PlusFunctor to add to res_bar, MinusFunctor to subract
* `trans` (optional): if true, the transpose operation is applied

**In/Outs**

* `flux_bar`: the result of the vector matrix product between Q and `res_bar`

"""->
function weakDifferentiateElement_rev!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                                                       flux_bar::AbstractArray{Tflx,1},
                                                       res_bar::AbstractArray{Tres,1},
                                                       (±)::UnaryFunctor=Add();
                                                       trans::Bool=false)
  @assert( sbp.numnodes == size(flux_bar,1) == size(res_bar,1) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  if trans # apply transposed Q
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        #res[i] += ±(sbp.Q[j,i,di]*flux[j])
        flux_bar[j] += ±(sbp.Q[j,i,di]*res_bar[i])
      end
    end
  else # apply Q
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        #res[i] += ±(sbp.Q[i,j,di]*flux[j])
        flux_bar[j] += ±(sbp.Q[i,j,di]*res_bar[i])
      end
    end
  end
end

function weakDifferentiateElement_rev!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp}, di::Int,
                                                       flux_bar::AbstractArray{Tflx,2},
                                                       res_bar::AbstractArray{Tres,2},
                                                       (±)::UnaryFunctor=Add();
                                                       trans::Bool=false)
  @assert( sbp.numnodes == size(flux_bar,2) == size(res_bar,2) )
  @assert( length(flux_bar) == length(res_bar) )
  @assert( di > 0 && di <= size(sbp.Q,3) )
  if trans # apply transposed Q
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        for field = 1:size(flux_bar,1)
          #res[field,i] += ±(sbp.Q[j,i,di]*flux[field,j])
          flux_bar[field,j] += ±(sbp.Q[j,i,di]*res_bar[field,i])
        end
      end
    end
  else # apply Q
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        for field = 1:size(flux_bar,1)
          #res[field,i] += ±(sbp.Q[i,j,di]*flux[field,j])
          flux_bar[field,j] += ±(sbp.Q[i,j,di]*res_bar[field,i])
        end
      end
    end
  end
end
