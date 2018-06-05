# This file contains the reverse-mode version of the methods in
# volumeintegrate.jl

"""
### SummationByParts.volumeintegrate_rev!

This is the reverse differentiated version of volumeintegrate!.  See
volumeintegrate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `u` variable.

**Inputs**

* `sbp`: an SBP operator type
* `res_bar`: vector applied to the left of the H operator
* `±`: PlusFunctor to add to res, MinusFunctor to subract

**In/Outs**

* `u_bar`: the result of the vector matrix product between H and `res_bar`

"""
function volumeintegrate_rev!{Tsbp,Tsol,Tres}(sbp::AbstractSBP{Tsbp},
                                              u_bar::AbstractArray{Tsol,2},
                                              res_bar::AbstractArray{Tres,2},
                                              (±)::UnaryFunctor=Add())
  @assert( sbp.numnodes == size(u_bar,1) && sbp.numnodes == size(res_bar,1) )
  @assert( length(u_bar) == length(res_bar) )
  H = sbp.w
  for elem = 1:size(u_bar,2)
    for i = 1:sbp.numnodes
      # res[i,elem] += ±(H[i]*u[i,elem])
      u_bar[i,elem] += ±(H[i]*res_bar[i,elem])
    end
  end
end

function volumeintegrate_rev!{Tsbp,Tsol,Tres}(sbp::AbstractSBP{Tsbp},
                                              u_bar::AbstractArray{Tsol,3},
                                              res_bar::AbstractArray{Tres,3},
                                              (±)::UnaryFunctor=Add())
  @assert( sbp.numnodes == size(u_bar,2) && sbp.numnodes == size(res_bar,2) )
  @assert( length(u_bar) == length(res_bar) )
  H = sbp.w
  for elem = 1:size(u_bar,3)
    for i = 1:sbp.numnodes
      for field = 1:size(u_bar,1)
        # res[field,i,elem] += ±(H[i]*u[field,i,elem])
        u_bar[field,i,elem] += ±(H[i]*res_bar[field,i,elem])
      end
    end
  end
end

"""
### SummationByParts.volumeIntegrateElement_rev!

This is the reverse differentiated version of volumeIntegrateElement!.  See
volumeintegrate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `u` variable.

**Inputs**

* `sbp`: an SBP operator type
* `res_bar`: vector applied to the left of the H operator
* `±`: PlusFunctor to add to res, MinusFunctor to subract

**In/Outs**

* `u_bar`: the result of the vector matrix product between H and `res_bar`

"""
function volumeIntegrateElement_rev!{Tsbp,Tsol,Tres}(sbp::AbstractSBP{Tsbp},
                                                     u_bar::AbstractArray{Tsol,1},
                                                     res_bar::AbstractArray{Tres,1},
                                                     (±)::UnaryFunctor=Add())
  @asserts_enabled begin
    @assert( sbp.numnodes == size(u_bar,1) == size(res_bar,1) )
  end

  for i = 1:sbp.numnodes
    # res[i] += ±(sbp.w[i]*u[i])
    u_bar[i] += ±(sbp.w[i]*res_bar[i])
  end
end

function volumeIntegrateElement_rev!{Tsbp,Tsol,Tres}(sbp::AbstractSBP{Tsbp},
                                                     u_bar::AbstractArray{Tsol,2},
                                                     res_bar::AbstractArray{Tres,2},
                                                     (±)::UnaryFunctor=Add())

  @asserts_enabled begin
    @assert( sbp.numnodes == size(u_bar,2) == size(res_bar,2) )
    @assert( length(u_bar) == length(res_bar) )
  end

  for i = 1:sbp.numnodes
    for field = 1:size(u_bar,1)
      # res[field,i] += ±(sbp.w[i]*u[field,i])
      u_bar[field,i] += ±(sbp.w[i]*res_bar[field,i])
    end
  end
end
