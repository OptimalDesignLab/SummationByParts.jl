# This file contains the reverse-mode version of the methods in faceintegrate.jl

@doc """
### SummationByParts.integratefunctional_rev!

This is the reverse differentiated version of integratefunctional!.  See
faceintegrate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `flux` variable.

**Inputs**

* `sbpface`: an SBP face operator type
* `bndryfaces`: list of boundary faces stored as an array of `Boundary`s
* `fun_bar`: incoming quantity that multiplies the functional from the left
* `±`: PlusFunctor to add to `flux_bar`, MinusFunctor to subract

**In/Outs**

* `flux_bar`: result of the vector-Jacobian product.

"""->
function integratefunctional_rev!{Tsbp,Tflx,Tfun}(sbpface::DenseFace{Tsbp},
                                                  bndryfaces::AbstractArray{Boundary},
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

function integratefunctional_rev!{Tsbp,Tflx,Tfun}(sbpface::SparseFace{Tsbp},
                                                  bndryfaces::AbstractArray{Boundary},
                                                  flux_bar::AbstractArray{Tflx,3},
                                                  fun_bar::AbstractArray{Tfun,1},
                                                  (±)::UnaryFunctor=Add())
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

@doc """
### SummationByParts.integrateBoundaryFunctional_rev!

This is the reverse differentiated version of integrateBoundaryFunctional!.  See
faceintegrate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `flux` variable.

**Inputs**

* `sbpface`: an SBP face operator type
* `face`: the face of the element to integrate
* `fun_bar`: incoming quantity that multiplies the functional from the left
* `±`: PlusFunctor to add to `flux_bar`, MinusFunctor to subract

**In/Outs**

* `flux_bar`: result of the vector-Jacobian product.

"""->
function integrateBoundaryFunctional_rev!{
  Tsbp,Tflx,Tfun}(sbpface::DenseFace{Tsbp}, face::Integer,
                  flux_bar::AbstractArray{Tflx,2}, fun_bar::AbstractArray{Tfun,1},
                  (±)::UnaryFunctor=Add())
  @assert( size(sbpface.interp,2) == size(flux_bar,2) )
  @assert( size(flux_bar,1) == size(fun_bar,1) )
  for i = 1:sbpface.numnodes
    for field = 1:size(fun_bar,1)
      # fun[field] += ±(sbpface.wface[i]*flux[field,i])
      flux_bar[field,i] += ±(sbpface.wface[i]*fun_bar[field])
    end
  end
end

function integrateBoundaryFunctional_rev!{
  Tsbp,Tflx,Tfun}(sbpface::SparseFace{Tsbp}, face::Integer,
                  flux_bar::AbstractArray{Tflx,2}, fun_bar::AbstractArray{Tfun,1},
                  (±)::UnaryFunctor=Add())
  @assert( size(flux_bar,1) == size(fun_bar,1) )
  for i = 1:sbpface.numnodes
    for field = 1:size(fun_bar,1)
      # fun[field] += ±(sbpface.wface[i]*flux[field,i])
      flux_bar[field,i] += ±(sbpface.wface[i]*fun_bar[field])
    end
  end
end

@doc """
### SummationByParts.boundaryintegrate_rev!

This is the reverse differentiated version of boundaryintegrate!.  See
faceintegrate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `flux` variable.

**Inputs**

* `sbpface`: an SBP AbstractFace type
* `bndryfaces`: list of boundary faces stored as an array of `Boundary`s
* `res_bar`: vector applied to the left of the (R^T*B) operator
* `±`: PlusFunctor to add to `flux_bar`, MinusFunctor to subract

**In/Outs**

* `flux_bar`: result of the vector matrix product between (R^T*B) and `res_bar`

"""->
function boundaryintegrate_rev!{Tsbp,Tflx,Tres}(sbpface::DenseFace{Tsbp},
                                                bndryfaces::AbstractArray{Boundary},
                                                flux_bar::AbstractArray{Tflx,2},
                                                res_bar::AbstractArray{Tres,2},
                                                (±)::UnaryFunctor=Add())
  @assert( size(sbpface.interp,1) <= size(res_bar,1) )
  @assert( size(sbpface.interp,2) == size(flux_bar,1) )
  @assert( size(bndryfaces,1) == size(flux_bar,2) )
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      wflux_bar = zero(Tres)
      for j = 1:sbpface.stencilsize
        # res[sbpface.perm[j,bndry.face],bndry.element] +=
        #   ±(sbpface.interp[j,i]*wflux)
        wflux_bar += ±(sbpface.interp[j,i]*
                       res_bar[sbpface.perm[j,bndry.face],bndry.element])           
      end
      # wflux = sbpface.wface[i]*flux[i,bindex]
      flux_bar[i,bindex] += sbpface.wface[i]*wflux_bar
    end
  end
end

function boundaryintegrate_rev!{Tsbp,Tflx,Tres}(sbpface::DenseFace{Tsbp},
                                            bndryfaces::AbstractArray{Boundary},
                                            flux_bar::AbstractArray{Tflx,3},
                                            res_bar::AbstractArray{Tres,3},
                                            (±)::UnaryFunctor=Add())
  @assert( size(sbpface.interp,1) <= size(res_bar,2) )
  @assert( size(sbpface.interp,2) == size(flux_bar,2) )
  @assert( size(bndryfaces,1) == size(flux_bar,3) )
  wflux_bar = zeros(Tres, size(res_bar,1))
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      fill!(wflux_bar, zero(Tres))
      for j = 1:sbpface.stencilsize
        for field = 1:size(res_bar,1)
          # res[field,sbpface.perm[j,bndry.face],bndry.element] +=
          #   ±(sbpface.interp[j,i]*wflux[field])
          wflux_bar[field] +=
            ±(sbpface.interp[j,i]*
              res_bar[field,sbpface.perm[j,bndry.face],bndry.element])
        end
      end
      for field = 1:size(res_bar,1)
        # wflux[field] = sbpface.wface[i]*flux[field,i,bindex]
        flux_bar[field,i,bindex] += sbpface.wface[i]*wflux_bar[field]
      end
    end
  end
end

function boundaryintegrate_rev!{Tsbp,Tflx,Tres}(sbpface::SparseFace{Tsbp},
                                                bndryfaces::AbstractArray{Boundary},
                                                flux_bar::AbstractArray{Tflx,2},
                                                res_bar::AbstractArray{Tres,2},
                                                (±)::UnaryFunctor=Add())
  @assert( size(bndryfaces,1) == size(flux_bar,2) )
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      # res[sbpface.perm[i,bndry.face],bndry.element] +=
      #   ±(sbpface.wface[i]*flux[i,bindex])
      flux_bar[i,bindex] +=
        ±(sbpface.wface[i]*res_bar[sbpface.perm[i,bndry.face],bndry.element])
    end
  end
end

function boundaryintegrate_rev!{Tsbp,Tflx,Tres}(sbpface::SparseFace{Tsbp},
                                            bndryfaces::AbstractArray{Boundary},
                                            flux_bar::AbstractArray{Tflx,3},
                                            res_bar::AbstractArray{Tres,3},
                                            (±)::UnaryFunctor=Add())
  @assert( size(bndryfaces,1) == size(flux_bar,3) )
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      for field = 1:size(res_bar,1)
        # res[field,sbpface.perm[i,bndry.face],bndry.element] +=
        #   ±(sbpface.wface[i]*flux[field,i,bindex])
        flux_bar[field,i,bindex] +=
          ±(sbpface.wface[i]*
            res_bar[field,sbpface.perm[i,bndry.face],bndry.element])
      end
    end
  end
end

@doc """
### SummationByParts.boundaryFaceIntegrate_rev!

This is the reverse differentiated version of boundaryFaceIntegrate!.  See
faceintegrate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `flux` variable.

**Inputs**

* `sbpface`: an SBP AbstractFace type
* `face`: the face of the element to integrate and project back to the element
* `res_bar`: vector applied to the left of the (R^T*B) operator
* `±`: PlusFunctor to add to `flux_bar`, MinusFunctor to subract

**In/Outs**

* `flux_bar`: result of the vector matrix product between (R^T*B) and `res_bar`

"""->
function boundaryFaceIntegrate_rev!{Tsbp,Tflx,Tres}(sbpface::AbstractFace{Tsbp},
                                                    face::Integer,
                                                    flux_bar::AbstractArray{Tflx,1},
                                                    res_bar::AbstractArray{Tres,1},
                                                    (±)::UnaryFunctor=Add())
  @assert( size(sbpface.interp,1) <= size(res_bar,1) )
  @assert( size(sbpface.interp,2) == size(flux_bar,1) )
  for i = 1:sbpface.numnodes
    wflux_bar = zero(Tres)
    for j = 1:sbpface.stencilsize
      # res[sbpface.perm[j,face]] += ±(sbpface.interp[j,i]*wflux)
      wflux_bar += ±(sbpface.interp[j,i]*res_bar[sbpface.perm[j,face]])
    end
    # wflux = sbpface.wface[i]*flux[i]
    flux_bar[i] += sbpface.wface[i]*wflux_bar
  end
end

function boundaryFaceIntegrate_rev!{Tsbp,Tflx,Tres}(sbpface::DenseFace{Tsbp},
                                                    face::Integer,
                                                    flux_bar::AbstractArray{Tflx,2},
                                                    res_bar::AbstractArray{Tres,2},
                                                    (±)::UnaryFunctor=Add())
  @assert( size(sbpface.interp,1) <= size(res_bar,2) )
  @assert( size(sbpface.interp,2) == size(flux_bar,2) )
  @assert( size(flux_bar,1) == size(res_bar,1) )
  wflux_bar = zeros(Tres, size(res_bar,1))
  for i = 1:sbpface.numnodes
    fill!(wflux_bar, zero(Tres))
    for j = 1:sbpface.stencilsize
      for field = 1:size(res_bar,1)
        # res[field,sbpface.perm[j,face]] += ±(sbpface.interp[j,i]*wflux[field])
        wflux_bar[field] += ±(sbpface.interp[j,i]*
                              res_bar[field,sbpface.perm[j,face]])
      end
    end
    for field = 1:size(res_bar,1)
      # wflux[field] = sbpface.wface[i]*flux[field,i]
      flux_bar[field,i] += sbpface.wface[i]*wflux_bar[field]
    end
  end
end

function boundaryFaceIntegrate_rev!{Tsbp,Tflx,Tres}(sbpface::SparseFace{Tsbp},
                                                    face::Integer,
                                                    flux_bar::AbstractArray{Tflx,1},
                                                    res_bar::AbstractArray{Tres,1},
                                                    (±)::UnaryFunctor=Add())
  for i = 1:sbpface.numnodes
    # res[sbpface.perm[i,face]] += ±(sbpface.wface[i]*flux[i])
    flux_bar[i] += ±(sbpface.wface[i]*res_bar[sbpface.perm[i,face]])
  end
end

function boundaryFaceIntegrate_rev!{Tsbp,Tflx,Tres}(sbpface::SparseFace{Tsbp},
                                                    face::Integer,
                                                    flux_bar::AbstractArray{Tflx,2},
                                                    res_bar::AbstractArray{Tres,2},
                                                    (±)::UnaryFunctor=Add())
  @assert( size(flux_bar,1) == size(res_bar,1) )
  for i = 1:sbpface.numnodes
    for field = 1:size(res_bar,1)
      # res[field,sbpface.perm[i,face]] += ±(sbpface.wface[i]*flux[field,i])
      flux_bar[field,i] += ±(sbpface.wface[i]*
                             res_bar[field,sbpface.perm[i,face]])
    end
  end
end

@doc """
### SummationByParts.interiorfaceintegrate_rev!

This is the reverse differentiated version of interiorfaceintegrate!.  See
faceintegrate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `flux` variable.

**Inputs**

* `sbpface`: an SBP AbstractFace type
* `ifaces`: list of element interfaces stored as an array of `Interface`s
* `res_bar`: vector applied to the left of the (R^T*B) operator
* `±`: PlusFunctor to add to `flux_bar`, MinusFunctor to subract

**In/Outs**

* `flux_bar`: result of the vector matrix product between (R^T*B) and `res_bar`

"""->
function interiorfaceintegrate_rev!{Tsbp,Tflx,Tres}(sbpface::DenseFace{Tsbp},
                                                    ifaces::AbstractArray{Interface},
                                                    flux_bar::AbstractArray{Tflx,2},
                                                    res_bar::AbstractArray{Tres,2},
                                                    (±)::UnaryFunctor=Add())
  @assert( size(sbpface.interp,1) <= size(res_bar,1) )
  @assert( size(sbpface.interp,2) == size(flux_bar,1) )
  @assert( size(ifaces,1) == size(flux_bar,2) )
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]
      for j = 1:sbpface.stencilsize
        # res[sbpface.perm[j,face.faceL],face.elementL] += 
        #  ±(sbpface.interp[j,i]*sbpface.wface[i]*flux[i,findex])
        flux_bar[i,findex] += ±(sbpface.interp[j,i]*sbpface.wface[i]*
                                res_bar[sbpface.perm[j,face.faceL],face.elementL])
        # res[sbpface.perm[j,face.faceR],face.elementR] -=
        #  ±(sbpface.interp[j,iR]*sbpface.wface[iR]*flux[i,findex])
        flux_bar[i,findex] -= ±(sbpface.interp[j,iR]*sbpface.wface[iR]*
                                res_bar[sbpface.perm[j,face.faceR],face.elementR])
      end
    end
  end
end

function interiorfaceintegrate_rev!{Tsbp,Tflx,Tres}(sbpface::DenseFace{Tsbp},
                                                ifaces::AbstractArray{Interface},
                                                flux_bar::AbstractArray{Tflx,3},
                                                res_bar::AbstractArray{Tres,3},
                                                (±)::UnaryFunctor=Add())
  @assert( size(res_bar,1) == size(flux_bar,1) )  
  @assert( size(sbpface.interp,1) <= size(res_bar,2) )
  @assert( size(sbpface.interp,2) == size(flux_bar,2) )
  @assert( size(ifaces,1) == size(flux_bar,3) )
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]
      for j = 1:sbpface.stencilsize
        for field = 1:size(res_bar,1)
          # res[field,sbpface.perm[j,face.faceL],face.elementL] += 
          #   ±(sbpface.interp[j,i]*sbpface.wface[i]*flux[field,i,findex])
          flux_bar[field,i,findex] +=
            ±(sbpface.interp[j,i]*sbpface.wface[i]*
              res_bar[field,sbpface.perm[j,face.faceL],face.elementL])
          # res[field,sbpface.perm[j,face.faceR],face.elementR] -=
          #   ±(sbpface.interp[j,iR]*sbpface.wface[iR]*flux[field,i,findex])
          flux_bar[field,i,findex] -=
            ±(sbpface.interp[j,iR]*sbpface.wface[iR]*
              res_bar[field,sbpface.perm[j,face.faceR],face.elementR])
        end
      end
    end
  end
end

function interiorfaceintegrate_rev!{Tsbp,Tflx,Tres}(sbpface::SparseFace{Tsbp},
                                                    ifaces::AbstractArray{Interface},
                                                    flux_bar::AbstractArray{Tflx,2},
                                                    res_bar::AbstractArray{Tres,2},
                                                    (±)::UnaryFunctor=Add())
  @assert( size(ifaces,1) == size(flux_bar,2) )
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]
      # res[sbpface.perm[i,face.faceL],face.elementL] +=
      #   ±(sbpface.wface[i]*flux[i,findex])
      flux_bar[i,findex] += ±(sbpface.wface[i]*
                              res_bar[sbpface.perm[i,face.faceL],face.elementL])
      # res[sbpface.perm[iR,face.faceR],face.elementR] -=
      #   ±(sbpface.wface[iR]*flux[i,findex])
      flux_bar[i,findex] -= ±(sbpface.wface[i]*
                              res_bar[sbpface.perm[iR,face.faceR],face.elementR])
    end
  end
end

function interiorfaceintegrate_rev!{Tsbp,Tflx,Tres}(sbpface::SparseFace{Tsbp},
                                                    ifaces::AbstractArray{Interface},
                                                    flux_bar::AbstractArray{Tflx,3},
                                                    res_bar::AbstractArray{Tres,3},
                                                    (±)::UnaryFunctor=Add())
  @assert( size(res_bar,1) == size(flux_bar,1) )  
  @assert( size(ifaces,1) == size(flux_bar,3) )
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]
      for field = 1:size(res_bar,1)
        # res[field,sbpface.perm[i,face.faceL],face.elementL] +=
        #   ±(sbpface.wface[i]*flux[field,i,findex])
        flux_bar[field,i,findex] +=
          ±(sbpface.wface[i]*
            res_bar[field,sbpface.perm[i,face.faceL],face.elementL])
        # res[field,sbpface.perm[iR,face.faceR],face.elementR] -=
        #   ±(sbpface.wface[iR]*flux[field,i,findex])
        flux_bar[field,i,findex] -=
          ±(sbpface.wface[i]*
            res_bar[field,sbpface.perm[iR,face.faceR],face.elementR])
      end
    end
  end
end

@doc """
### SummationByParts.interiorFaceIntegrate_rev!

This is the reverse differentiated version of interiorFaceIntegrate!.  See
faceintegrate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `flux` variable.

**Inputs**

* `sbpface`: an SBP AbstractFace type
* `ifaces`: list of element interfaces stored as an array of `Interface`s
* `resL_bar`: vector (left element) applied to the left of the (R^T*B) operator 
* `resR_bar`: vector (right element) applied to the left of the (R^T*B) operator 
* `±`: PlusFunctor to add to `flux_bar`, MinusFunctor to subract

**In/Outs**

* `flux_bar`: result of the vector-matrix product

"""->
function interiorFaceIntegrate_rev!{Tsbp,Tflx,Tres}(sbpface::DenseFace{Tsbp},
                                                    iface::Interface,
                                                    flux_bar::AbstractArray{Tflx,1},
                                                    resL_bar::AbstractArray{Tres,1},
                                                    resR_bar::AbstractArray{Tres,1},
                                                    (±)::UnaryFunctor=Add())
  @assert( size(sbpface.interp,1) <= size(resL_bar,1) )
  @assert( size(sbpface.interp,1) <= size(resR_bar,1) )
  @assert( size(sbpface.interp,2) == size(flux_bar,1) )
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    for j = 1:sbpface.stencilsize
      # resL[sbpface.perm[j,iface.faceL]] += 
      #   ±(sbpface.interp[j,i]*sbpface.wface[i]*flux[i])
      flux_bar[i] += ±(sbpface.interp[j,i]*sbpface.wface[i]*
                       resL_bar[sbpface.perm[j,iface.faceL]])
      # resR[sbpface.perm[j,iface.faceR]] -=
      #   ±(sbpface.interp[j,iR]*sbpface.wface[iR]*flux[i])
      flux_bar[i] -= ±(sbpface.interp[j,iR]*sbpface.wface[iR]*
                       resR_bar[sbpface.perm[j,iface.faceR]])
    end
  end
end

function interiorFaceIntegrate_rev!{Tsbp,Tflx,Tres}(sbpface::DenseFace{Tsbp},
                                                    iface::Interface,
                                                    flux_bar::AbstractArray{Tflx,2},
                                                    resL_bar::AbstractArray{Tres,2},
                                                    resR_bar::AbstractArray{Tres,2},
                                                    (±)::UnaryFunctor=Add())
  @assert( size(resL_bar,1) == size(resR_bar,1) == size(flux_bar,1) )  
  @assert( size(sbpface.interp,1) <= size(resL_bar,2) )
  @assert( size(sbpface.interp,1) <= size(resR_bar,2) )
  @assert( size(sbpface.interp,2) == size(flux_bar,2) )
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    for j = 1:sbpface.stencilsize
      for field = 1:size(flux_bar,1)
        # resL[field,sbpface.perm[j,iface.faceL]] += 
        #   ±(sbpface.interp[j,i]*sbpface.wface[i]*flux[field,i])
        flux_bar[field,i] += ±(sbpface.interp[j,i]*sbpface.wface[i]*
                               resL_bar[field,sbpface.perm[j,iface.faceL]])
        # resR[field,sbpface.perm[j,iface.faceR]] -=
        #   ±(sbpface.interp[j,iR]*sbpface.wface[iR]*flux[field,i])
        flux_bar[field,i] -= ±(sbpface.interp[j,iR]*sbpface.wface[iR]*
                               resR_bar[field,sbpface.perm[j,iface.faceR]])
      end
    end
  end
end

function interiorFaceIntegrate_rev!{Tsbp,Tflx,Tres}(sbpface::SparseFace{Tsbp},
                                                    iface::Interface,
                                                    flux_bar::AbstractArray{Tflx,1},
                                                    resL_bar::AbstractArray{Tres,1},
                                                    resR_bar::AbstractArray{Tres,1},
                                                    (±)::UnaryFunctor=Add())
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    # resL[sbpface.perm[i,iface.faceL]] += ±(sbpface.wface[i]*flux[i])
    flux_bar[i] += ±(sbpface.wface[i]*resL_bar[sbpface.perm[i,iface.faceL]])
    # resR[sbpface.perm[iR,iface.faceR]] -= ±(sbpface.wface[i]*flux[i])
    flux_bar[i] -= ±(sbpface.wface[i]*resR_bar[sbpface.perm[iR,iface.faceR]])
  end
end

function interiorFaceIntegrate_rev!{Tsbp,Tflx,Tres}(sbpface::SparseFace{Tsbp},
                                                    iface::Interface,
                                                    flux_bar::AbstractArray{Tflx,2},
                                                    resL_bar::AbstractArray{Tres,2},
                                                    resR_bar::AbstractArray{Tres,2},
                                                    (±)::UnaryFunctor=Add())
  @assert( size(resL_bar,1) == size(resR_bar,1) == size(flux_bar,1) )  
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    for field = 1:size(flux_bar,1)
      # resL[field,sbpface.perm[i,iface.faceL]] +=
      #   ±(sbpface.wface[i]*flux[field,i])
      flux_bar[field,i] += ±(sbpface.wface[i]*
                             resL_bar[field,sbpface.perm[i,iface.faceL]])
      # resR[field,sbpface.perm[iR,iface.faceR]] -=
      #   ±(sbpface.wface[i]*flux[field,i])
      flux_bar[field,i] -= ±(sbpface.wface[i]*
                             resR_bar[field,sbpface.perm[iR,iface.faceR]])
    end
  end
end
