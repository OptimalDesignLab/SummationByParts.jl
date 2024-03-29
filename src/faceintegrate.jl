# This file gathers together methods related to integration over faces, both
# against test functions and for integral functionals.

"""
### SummationByParts.integratefunctional!

Integrates a given scalar (or vector) field over the boundary faces.

* For *scalar* fields, the dimensions of `uface` correspond to [face-node index,
  boundary index] and the scalar functional is a return value
* For *vector* fields, the dimensions of `uface` correspond to [field index,
  face-node index, boundary index] and the dimensions of `fun` correspond to
  [field index].

**Inputs**

* `sbpface`: an SBP face operator type
* `bndryfaces`: list of boundary faces stored as an array of `Boundary`s
* `flux`:  array of field data that is being integrated
* `±`: PlusFunctor to add to fun, MinusFunctor to subract

**In/Outs**

* `fun`: functional value (or vector) being *contributed* to by the integration

**Returns**

* `fun`: in the case of the scalar version, the functional value is returned

"""
function integratefunctional!(sbpface::DenseFace{Tsbp},
                                         bndryfaces::AbstractArray{Boundary},
                                         flux::AbstractArray{Tflx,2},
                                         (±)::UnaryFunctor=Add()) where {Tsbp,Tflx}
  @assert( size(sbpface.interp,2) == size(flux,1) )
  @assert( size(bndryfaces,1) == size(flux,2) )
  fun = zero(Tflx)
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      fun += ±(sbpface.wface[i]*flux[i,bindex])
    end
  end
  return fun
end

function integratefunctional!(sbpface::DenseFace{Tsbp},
                                              bndryfaces::AbstractArray{Boundary},
                                              flux::AbstractArray{Tflx,3},
                                              fun::AbstractArray{Tfun,1},
                                              (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tfun}
  @assert( size(sbpface.interp,2) == size(flux,2) )
  @assert( size(flux,1) == size(fun,1) )
  @assert( size(bndryfaces,1) == size(flux,3) )
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      for field = 1:size(fun,1)
        fun[field] += ±(sbpface.wface[i]*flux[field,i,bindex])
      end
    end
  end
end

function integratefunctional!(sbpface::SparseFace{Tsbp},
                                         bndryfaces::AbstractArray{Boundary},
                                         flux::AbstractArray{Tflx,2},
                                         (±)::UnaryFunctor=Add()) where {Tsbp,Tflx}
  @assert( size(bndryfaces,1) == size(flux,2) )
  fun = zero(Tflx)
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      fun += ±(sbpface.wface[i]*flux[i,bindex])
    end
  end
  return fun
end

function integratefunctional!(sbpface::SparseFace{Tsbp},
                                              bndryfaces::AbstractArray{Boundary},
                                              flux::AbstractArray{Tflx,3},
                                              fun::AbstractArray{Tfun,1},
                                              (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tfun}
  @assert( size(flux,1) == size(fun,1) )
  @assert( size(bndryfaces,1) == size(flux,3) )
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      for field = 1:size(fun,1)
        fun[field] += ±(sbpface.wface[i]*flux[field,i,bindex])
      end
    end
  end
end

"""
### SummationByParts.integrateBoundaryFunctional!

This is the single-face variant of integratefunctional!.  Integrates a given
scalar (or vector) field over the boundary faces.

* For *scalar* fields, the dimensions of `uface` correspond to [face-node index]
  and the scalar functional is a return value

* For *vector* fields, the dimensions of `uface` correspond to [field index,
  face-node index] and the dimensions of `fun` correspond to [field index].

**Inputs**

* `sbpface`: an SBP face operator type
* `face`: the face of the element to integrate
* `flux`:  array of field data that is being integrated
* `±`: PlusFunctor to add to fun, MinusFunctor to subract

**In/Outs**

* `fun`: functional value (or vector) being *contributed* to by the integration

**Returns**

* `fun`: in the case of the scalar version, the functional value is returned

"""
function integrateBoundaryFunctional!(sbpface::DenseFace{Tsbp},
                                                 face::Integer, 
                                                 flux::AbstractArray{Tflx,1},
                                                 (±)::UnaryFunctor=Add()) where {Tsbp,Tflx}
  @assert( size(sbpface.interp,2) == size(flux,1) )
  fun = zero(Tflx)
  for i = 1:sbpface.numnodes
      fun += ±(sbpface.wface[i]*flux[i])
  end
  return fun
end

function integrateBoundaryFunctional!(sbpface::DenseFace{Tsbp}, face::Integer,
                  flux::AbstractArray{Tflx,2}, fun::AbstractArray{Tfun,1},
                  (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tfun}
  @assert( size(sbpface.interp,2) == size(flux,2) )
  @assert( size(flux,1) == size(fun,1) )
  for i = 1:sbpface.numnodes
    for field = 1:size(fun,1)
      fun[field] += ±(sbpface.wface[i]*flux[field,i])
    end
  end
end

function integrateBoundaryFunctional!(sbpface::SparseFace{Tsbp},
                                                 face::Integer, 
                                                 flux::AbstractArray{Tflx,1},
                                                 (±)::UnaryFunctor=Add()) where {Tsbp,Tflx}
  fun = zero(Tflx)
  for i = 1:sbpface.numnodes
      fun += ±(sbpface.wface[i]*flux[i])
  end
  return fun
end

function integrateBoundaryFunctional!(sbpface::SparseFace{Tsbp}, face::Integer,
                  flux::AbstractArray{Tflx,2}, fun::AbstractArray{Tfun,1},
                  (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tfun}
  @assert( size(flux,1) == size(fun,1) )
  for i = 1:sbpface.numnodes
    for field = 1:size(fun,1)
      fun[field] += ±(sbpface.wface[i]*flux[field,i])
    end
  end
end

"""
### SummationByParts.boundaryintegrate!

Scales flux values at boundary cubature points by cubature weights, and then
performs transposed interpolation/extrapolation back to volume nodes. Different
methods are available depending on the rank of `flux`:

* For *scalar* fields, it is assumed that `flux` is a rank-2 array, with the
first dimension for the face-node index, and the second dimension for the
boundary index.
* For *vector* fields, `flux` is a rank-3 array, with the first dimension for
the index of the vector field, the second dimension for the face-node index, and
the third dimension for the boundary index.

The dimensions of `res` are still based on elements; the last dimension is for
the element index and the second-last dimension is for the element-local node
index.

**Inputs**

* `sbpface`: an SBP AbstractFace type
* `bndryfaces`: list of boundary faces stored as an array of `Boundary`s
* `flux`: array of flux data that is being integrated
* `±`: PlusFunctor to add to res, MinusFunctor to subract

**In/Outs**

* `res`: where the result of the integration is stored

**WARNING**: the order of the boundaries in `bndryfaces` and `flux` must be
  consistent.

"""
function boundaryintegrate!(sbpface::DenseFace{Tsbp},
                                            bndryfaces::AbstractArray{Boundary},
                                            flux::AbstractArray{Tflx,2},
                                            res::AbstractArray{Tres,2},
                                            (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  @assert( size(sbpface.interp,1) <= size(res,1) )
  @assert( size(sbpface.interp,2) == size(flux,1) )
  @assert( size(bndryfaces,1) == size(flux,2) )
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      wflux = sbpface.wface[i]*flux[i,bindex]
      for j = 1:sbpface.stencilsize
        res[sbpface.perm[j,bndry.face],bndry.element] +=
          ±(sbpface.interp[j,i]*wflux)
      end
    end
  end
end

function boundaryintegrate!(sbpface::DenseFace{Tsbp},
                                            bndryfaces::AbstractArray{Boundary},
                                            flux::AbstractArray{Tflx,3},
                                            res::AbstractArray{Tres,3},
                                            (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  @assert( size(sbpface.interp,1) <= size(res,2) )
  @assert( size(sbpface.interp,2) == size(flux,2) )
  @assert( size(bndryfaces,1) == size(flux,3) )
  wflux = zeros(Tflx, size(res,1))
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      for field = 1:size(res,1)
        wflux[field] = sbpface.wface[i]*flux[field,i,bindex]
      end
      for j = 1:sbpface.stencilsize
        for field = 1:size(res,1)
          res[field,sbpface.perm[j,bndry.face],bndry.element] +=
            ±(sbpface.interp[j,i]*wflux[field])
        end
      end
    end
  end
end

function boundaryintegrate!(sbpface::SparseFace{Tsbp},
                                            bndryfaces::AbstractArray{Boundary},
                                            flux::AbstractArray{Tflx,2},
                                            res::AbstractArray{Tres,2},
                                            (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  @assert( size(bndryfaces,1) == size(flux,2) )
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      res[sbpface.perm[i,bndry.face],bndry.element] +=
        ±(sbpface.wface[i]*flux[i,bindex])
    end
  end
end

function boundaryintegrate!(sbpface::SparseFace{Tsbp},
                                            bndryfaces::AbstractArray{Boundary},
                                            flux::AbstractArray{Tflx,3},
                                            res::AbstractArray{Tres,3},
                                            (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  @assert( size(bndryfaces,1) == size(flux,3) )
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      for field = 1:size(res,1)
        res[field,sbpface.perm[i,bndry.face],bndry.element] +=
          ±(sbpface.wface[i]*flux[field,i,bindex])
      end
    end
  end
end

"""
### SummationByParts.boundaryFaceIntegrate!

This is the single-face variant of boundaryintegrate!. Scales flux values at
boundary cubature points by cubature weights, and then performs transposed
interpolation/extrapolation back to volume nodes. Different methods are
available depending on the rank of `flux`:

* For *scalar* fields, it is assumed that `flux` is a rank-1 array, with the
first and only dimension for the face-node index.
* For *vector* fields, `flux` is a rank-2 array, with the first dimension for
the index of the vector field, and the second dimension for the face-node index.

The dimensions of `res` are still based on elements; the last dimension (in the
scalar case, the only dimension) is for the element-local node index.

**Inputs**

* `sbpface`: an SBP AbstractFace type
* `face`: the face of the element to integrate and project back to the element
* `flux`: array of flux data that is being integrated
* `±`: PlusFunctor to add to res, MinusFunctor to subract

**In/Outs**

* `res`: where the result of the integration is stored

"""
function boundaryFaceIntegrate!(sbpface::DenseFace{Tsbp},
                                                face::Integer,
                                                flux::AbstractArray{Tflx,1},
                                                res::AbstractArray{Tres,1},
                                                (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  @assert( size(sbpface.interp,1) <= size(res,1) )
  @assert( size(sbpface.interp,2) == size(flux,1) )
  for i = 1:sbpface.numnodes
    wflux = sbpface.wface[i]*flux[i]
    for j = 1:sbpface.stencilsize
      res[sbpface.perm[j,face]] += ±(sbpface.interp[j,i]*wflux)
    end
  end
end

function boundaryFaceIntegrate!(sbpface::DenseFace{Tsbp},
                                                face::Integer,
                                                flux::AbstractArray{Tflx,2},
                                                res::AbstractArray{Tres,2},
                                                (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  @assert( size(sbpface.interp,1) <= size(res,2) )
  @assert( size(sbpface.interp,2) == size(flux,2) )
  @assert( size(flux,1) == size(res,1) )
  wflux = zeros(Tflx, size(res,1))
  for i = 1:sbpface.numnodes
    for field = 1:size(res,1)
      wflux[field] = sbpface.wface[i]*flux[field,i]
    end
    for j = 1:sbpface.stencilsize
      for field = 1:size(res,1)
        res[field,sbpface.perm[j,face]] += ±(sbpface.interp[j,i]*wflux[field])
      end
    end
  end
end

function boundaryFaceIntegrate!(sbpface::SparseFace{Tsbp},
                                                face::Integer,
                                                flux::AbstractArray{Tflx,1},
                                                res::AbstractArray{Tres,1},
                                                (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  for i = 1:sbpface.numnodes
    res[sbpface.perm[i,face]] += ±(sbpface.wface[i]*flux[i])
  end
end

function boundaryFaceIntegrate!(sbpface::SparseFace{Tsbp},
                                                face::Integer,
                                                flux::AbstractArray{Tflx,2},
                                                res::AbstractArray{Tres,2},
                                                (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  @assert( size(flux,1) == size(res,1) )
  for i = 1:sbpface.numnodes
    for field = 1:size(res,1)
      res[field,sbpface.perm[i,face]] += ±(sbpface.wface[i]*flux[field,i])
    end
  end
end

"""
### SummationByParts.interiorfaceintegrate!

Scales flux values at element-interface cubature points by cubature weights, and
then performs transposed interpolation/extrapolation back to volume nodes.
Different methods are available depending on the rank of `flux`:

* For *scalar* fields, it is assumed that `flux` is a rank-2 array, with the
first dimension for the face-node index, and the second dimension for the
interface index.
* For *vector* fields, `flux` is a rank-3 array, with the first dimension for
the index of the vector field, the second dimension for the face-node index, and
the third dimension for the interface index.

The dimensions of `res` are still based on elements; the last dimension is for
the element index and the second-last dimension is for the element-local node
index.

**Inputs**

* `sbpface`: an SBP AbstractFace type
* `ifaces`: list of element interfaces stored as an array of `Interface`s
* `flux`: array of flux data that is being integrated
* `±`: PlusFunctor to add to res, MinusFunctor to subract

**In/Outs**

* `res`: where the result of the integration is stored

**WARNING**: the order of the interfaces in `ifaces` and `flux` must be
  consistent.

"""
function interiorfaceintegrate!(sbpface::DenseFace{Tsbp},
                                                ifaces::AbstractArray{Interface},
                                                flux::AbstractArray{Tflx,2},
                                                res::AbstractArray{Tres,2},
                                                (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  @assert( size(sbpface.interp,1) <= size(res,1) )
  @assert( size(sbpface.interp,2) == size(flux,1) )
  @assert( size(ifaces,1) == size(flux,2) )
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]
      for j = 1:sbpface.stencilsize
        res[sbpface.perm[j,face.faceL],face.elementL] += 
        ±(sbpface.interp[j,i]*sbpface.wface[i]*flux[i,findex])
        res[sbpface.perm[j,face.faceR],face.elementR] -=
          ±(sbpface.interp[j,iR]*sbpface.wface[iR]*flux[i,findex])
      end
    end
  end
end

function interiorfaceintegrate!(sbpface::DenseFace{Tsbp},
                                                ifaces::AbstractArray{Interface},
                                                flux::AbstractArray{Tflx,3},
                                                res::AbstractArray{Tres,3},
                                                (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  @assert( size(res,1) == size(flux,1) )  
  @assert( size(sbpface.interp,1) <= size(res,2) )
  @assert( size(sbpface.interp,2) == size(flux,2) )
  @assert( size(ifaces,1) == size(flux,3) )
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]
      for j = 1:sbpface.stencilsize
        for field = 1:size(res,1)
          res[field,sbpface.perm[j,face.faceL],face.elementL] += 
          ±(sbpface.interp[j,i]*sbpface.wface[i]*flux[field,i,findex])
          res[field,sbpface.perm[j,face.faceR],face.elementR] -=
            ±(sbpface.interp[j,iR]*sbpface.wface[iR]*flux[field,i,findex])
        end
      end
    end
  end
end

function interiorfaceintegrate!(sbpface::SparseFace{Tsbp},
                                                ifaces::AbstractArray{Interface},
                                                flux::AbstractArray{Tflx,2},
                                                res::AbstractArray{Tres,2},
                                                (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  @assert( size(ifaces,1) == size(flux,2) )
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]
      res[sbpface.perm[i,face.faceL],face.elementL] +=
        ±(sbpface.wface[i]*flux[i,findex])
      res[sbpface.perm[iR,face.faceR],face.elementR] -=
        ±(sbpface.wface[iR]*flux[i,findex])
    end
  end
end

function interiorfaceintegrate!(sbpface::SparseFace{Tsbp},
                                                ifaces::AbstractArray{Interface},
                                                flux::AbstractArray{Tflx,3},
                                                res::AbstractArray{Tres,3},
                                                (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  @assert( size(res,1) == size(flux,1) )  
  @assert( size(ifaces,1) == size(flux,3) )
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]
      for field = 1:size(res,1)
        res[field,sbpface.perm[i,face.faceL],face.elementL] +=
          ±(sbpface.wface[i]*flux[field,i,findex])
        res[field,sbpface.perm[iR,face.faceR],face.elementR] -=
          ±(sbpface.wface[iR]*flux[field,i,findex])
      end
    end
  end
end

"""
### SummationByParts.interiorFaceIntegrate!

This is the single-face variant of interiorfaceintegrate!.  Scales flux values
at element-interface cubature points by cubature weights, and then performs
transposed interpolation/extrapolation back to volume nodes.  Different methods
are available depending on the rank of `flux`:

* For *scalar* fields, it is assumed that `flux` is a rank-1 array, with the
first and only dimension for the face-node index
* For *vector* fields, `flux` is a rank-2 array, with the first dimension for
the index of the vector field, and the second dimension for the face-node index

The dimensions of `resL` and `resR` are still based on elements; the last
dimension (in the scalar case, the only dimension) is for the element-local node
index.

**Inputs**

* `sbpface`: an SBP AbstractFace type
* `ifaces`: list of element interfaces stored as an array of `Interface`s
* `flux`: array of flux data that is being integrated
* `±`: PlusFunctor to add to res, MinusFunctor to subract

**In/Outs**

* `resL`: where the result of the integration is stored for the *left* element
* `resR`: where the result of the integration is stored for the *right* element

"""
function interiorFaceIntegrate!(sbpface::DenseFace{Tsbp},
                                                iface::Interface,
                                                flux::AbstractArray{Tflx,1},
                                                resL::AbstractArray{Tres,1},
                                                resR::AbstractArray{Tres,1},
                                                (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  @assert( size(sbpface.interp,1) <= size(resL,1) )
  @assert( size(sbpface.interp,1) <= size(resR,1) )
  @assert( size(sbpface.interp,2) == size(flux,1) )
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    for j = 1:sbpface.stencilsize
      resL[sbpface.perm[j,iface.faceL]] += 
      ±(sbpface.interp[j,i]*sbpface.wface[i]*flux[i])
      resR[sbpface.perm[j,iface.faceR]] -=
        ±(sbpface.interp[j,iR]*sbpface.wface[iR]*flux[i])
    end
  end
end

function interiorFaceIntegrate!(sbpface::DenseFace{Tsbp},
                                                iface::Interface,
                                                flux::AbstractArray{Tflx,2},
                                                resL::AbstractArray{Tres,2},
                                                resR::AbstractArray{Tres,2},
                                                (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  @assert( size(resL,1) == size(resR,1) == size(flux,1) )  
  @assert( size(sbpface.interp,1) <= size(resL,2) )
  @assert( size(sbpface.interp,1) <= size(resR,2) )
  @assert( size(sbpface.interp,2) == size(flux,2) )
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    for j = 1:sbpface.stencilsize
      for field = 1:size(flux,1)
        resL[field,sbpface.perm[j,iface.faceL]] += 
        ±(sbpface.interp[j,i]*sbpface.wface[i]*flux[field,i])
        resR[field,sbpface.perm[j,iface.faceR]] -=
          ±(sbpface.interp[j,iR]*sbpface.wface[iR]*flux[field,i])
      end
    end
  end
end

function interiorFaceIntegrate!(sbpface::SparseFace{Tsbp},
                                                iface::Interface,
                                                flux::AbstractArray{Tflx,1},
                                                resL::AbstractArray{Tres,1},
                                                resR::AbstractArray{Tres,1},
                                                (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    resL[sbpface.perm[i,iface.faceL]] += ±(sbpface.wface[i]*flux[i])
    resR[sbpface.perm[iR,iface.faceR]] -= ±(sbpface.wface[i]*flux[i])
  end
end

function interiorFaceIntegrate!(sbpface::SparseFace{Tsbp},
                                                iface::Interface,
                                                flux::AbstractArray{Tflx,2},
                                                resL::AbstractArray{Tres,2},
                                                resR::AbstractArray{Tres,2},
                                                (±)::UnaryFunctor=Add()) where {Tsbp,Tflx,Tres}
  @assert( size(resL,1) == size(resR,1) == size(flux,1) )  
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    for field = 1:size(flux,1)
      resL[field,sbpface.perm[i,iface.faceL]] +=
        ±(sbpface.wface[i]*flux[field,i])
      resR[field,sbpface.perm[iR,iface.faceR]] -=
        ±(sbpface.wface[i]*flux[field,i])
    end
  end
end
