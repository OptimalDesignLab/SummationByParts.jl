# This file contains the reverse-mode version of the methods in faceintegrate.jl

@doc """
### SummationByParts.integratefunctional!

This is the reverse differentiated version of integratefunctional!.  See
faceintegrate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `flux` variable.

**Inputs**

* `sbpface`: an SBP face operator type
* `bndryfaces`: list of boundary faces stored as an array of `Boundary`s
* `fun_bar`: incoming quantity that multiplies the functional from the left
* `±`: PlusFunctor to add to fun, MinusFunctor to subract

**In/Outs**

* `flux_bar`: result of the vector-Jacobian product.

**Returns**

* `fun`: in the case of the scalar version, the functional value is returned

"""->
function integratefunctional_rev!{Tsbp,Tflx,Tfun}(sbpface::AbstractFace{Tsbp},
                                                  bndryfaces::Array{Boundary},
                                                  flux_bar::AbstractArray{Tflx,3},
                                                  fun_bar::AbstractArray{Tfun,1},
                                                  (±)::UnaryFunctor=Add())
  @assert( size(sbpface.interp,2) == size(flux_bar,2) )
  @assert( size(flux_bar,1) == size(fun_bar,1) )
  @assert( size(bndryfaces,1) == size(flux_bar,3) )
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      for field = 1:size(fun_bar,1)
        #fun[field] += ±(sbpface.wface[i]*flux[field,i,bindex])
        flux_bar[field,i,bindex] += ±(sbpface.wface[i])*fun_bar[field]
      end
    end
  end
end

