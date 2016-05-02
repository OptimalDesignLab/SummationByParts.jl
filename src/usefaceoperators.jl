# This file gathers together functions related to element face operations

@doc """
### SummationByParts.boundaryinterpolate!

Interpolates element-node data to the element-face cubature nodes for the given
set of faces.  Different methods are available depending on the rank of `uvol`:

* For *scalar* fields, it is assumed that `uvol` (`uface`) is a rank-2 array,
with the first dimension for the node index, and the second dimension for the
element index (boundary index).
* For *vector* fields, `uvol` (`uface`) is a rank-3 array, with the first dimension for
the index of the vector field, the second dimension for the node index, and
the third dimension for the element index (boundary index).

**Inputs**

* `sbpface`: an SBP face operator type
* `bndryfaces`: list of boundary faces stored as an array of `Boundary`s
* `uvol`: array of field data that is being interpolated

**In/Outs**

* `uface`: field data interpolated to the faces

"""->
function boundaryinterpolate!{Tsbp,Tsol}(sbpface::AbstractFace{Tsbp},
                                         bndryfaces::Array{Boundary},
                                         uvol::AbstractArray{Tsol,2},
                                         uface::AbstractArray{Tsol,2})
  @assert( size(sbpface.interp,1) <= size(uvol,1) )
  @assert( size(sbpface.interp,2) == size(uface,1) )
  @inbounds begin
    for (bindex, bndry) in enumerate(bndryfaces)
      for i = 1:sbpface.numnodes
        uface[i,bindex] = zero(Tsol)
        for j = 1:sbpface.stencilsize
          uface[i,bindex] += sbpface.interp[j,i]*uvol[sbpface.perm[j,bndry.face],
                                                      bndry.element]
        end
      end
    end
  end
end

function boundaryinterpolate!{Tsbp,Tsol}(sbpface::AbstractFace{Tsbp},
                                         bndryfaces::Array{Boundary},
                                         uvol::AbstractArray{Tsol,3},
                                         uface::AbstractArray{Tsol,3})
  @assert( size(uvol,1) == size(uface,1) )
  @assert( size(sbpface.interp,1) <= size(uvol,2) )
  @assert( size(sbpface.interp,2) == size(uface,2) )
  @inbounds begin
    for (bindex, bndry) in enumerate(bndryfaces)
      for i = 1:sbpface.numnodes
        for field=1:size(uvol, 1)
          uface[field,i,bindex] = zero(Tsol)
        end
        for j = 1:sbpface.stencilsize
           for field = 1:size(uvol,1)
             uface[field,i,bindex] += sbpface.interp[j,i]*
             uvol[field,sbpface.perm[j,bndry.face],bndry.element]
           end
        end
      end
    end
  end
end

@doc """
### SummationByParts.boundaryinterpolate!

Interpolates vector field values at the nodes of a given element to a 
specified face of the element.

**Inputs**

* `sbpface`: an SBP face operator
* `face`: the face of the element to interpolate to
* `uvol`: the values at the nodes of the elements, dimensions numcomp x numnodes
           where numcomp is the number of components in the vector field and
           numnodes is the number of nodes in the element

**In/Outs**
* `uface`: the result of the interpolation, numcomp x sbpface.numnodes

"""->
function boundaryinterpolate!{Tsbp,Tsol}(sbpface::AbstractFace{Tsbp},
                                         face::Integer,
                                         uvol::AbstractArray{Tsol,2},
                                         uface::AbstractArray{Tsol,2})
  for i = 1:sbpface.numnodes
    for field=1:size(uvol, 1)
      uface[field,i] = zero(Tsol)
    end
    for j = 1:sbpface.stencilsize
       for field = 1:size(uvol,1)
         uface[field,i] += sbpface.interp[j,i]*uvol[field,sbpface.perm[j,face]]
       end
    end
  end
end

@doc """
### SummationByParts.boundaryintegrate!

Integrates a given flux over a boundary using appropriate mass matrices defined
on the element faces.  Different methods are available depending on the rank of
`flux`:

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

* `sbp`: an SBP operator type
* `bndryfaces`: list of boundary faces stored as an array of `Boundary`s
* `flux`: array of flux data that is being integrated

**In/Outs**

* `res`: where the result of the integration is stored

**WARNING**: the order of the boundaries in `bndryfaces` and `flux` must be
  consistent.

"""->
function boundaryintegrate!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp},
                                            bndryfaces::Array{Boundary},
                                            flux::AbstractArray{Tflx,2},
                                            res::AbstractArray{Tres,2})
  @assert( sbp.numnodes == size(res,1) )
  @assert( sbp.numfacenodes == size(flux,1) )
  @assert( size(bndryfaces,1) == size(flux,2) )
  @inbounds begin
    for (bindex, bndry) in enumerate(bndryfaces)
      for i = 1:sbp.numfacenodes        
        for j = 1:sbp.numfacenodes
          jB = sbp.facenodes[j, bndry.face]::Int # element index for jth node on face
          res[jB,bndry.element] += sbp.wface[j,i]*flux[i,bindex]
        end
      end
    end
  end
end

function boundaryintegrate!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp},
                                            bndryfaces::Array{Boundary},
                                            flux::AbstractArray{Tflx,3},
                                            res::AbstractArray{Tres,3})
  @assert( size(res,1) == size(flux,1) )
  @assert( sbp.numnodes == size(res,2) )
  @assert( sbp.numfacenodes == size(flux,2) )
  @assert( size(bndryfaces,1) == size(flux,3) )
  @inbounds begin
    for (bindex, bndry) in enumerate(bndryfaces)
      for i = 1:sbp.numfacenodes        
        for j = 1:sbp.numfacenodes
          jB = sbp.facenodes[j, bndry.face]::Int # element index for jth node on face
          for field = 1:size(res,1)
            res[field,jB,bndry.element] += sbp.wface[j,i]*flux[field,i,bindex]
          end
        end
      end
    end
  end
end

function boundaryintegrate!{Tsbp,Tflx,Tres}(sbpface::AbstractFace{Tsbp},
                                            bndryfaces::Array{Boundary},
                                            flux::AbstractArray{Tflx,2},
                                            res::AbstractArray{Tres,2})
  @assert( size(sbpface.interp,1) <= size(res,1) )
  @assert( size(sbpface.interp,2) == size(flux,1) )
  @assert( size(bndryfaces,1) == size(flux,2) )
  @inbounds begin
    for (bindex, bndry) in enumerate(bndryfaces)
      for i = 1:sbpface.numnodes
        wflux = sbpface.wface[i]*flux[i,bindex]
        for j = 1:sbpface.stencilsize
          res[sbpface.perm[j,bndry.face],bndry.element] +=
            sbpface.interp[j,i]*wflux
        end
      end
    end
  end
end

function boundaryintegrate!{Tsbp,Tflx,Tres}(sbpface::AbstractFace{Tsbp},
                                            bndryfaces::Array{Boundary},
                                            flux::AbstractArray{Tflx,3},
                                            res::AbstractArray{Tres,3})
  @assert( size(sbpface.interp,1) <= size(res,2) )
  @assert( size(sbpface.interp,2) == size(flux,2) )
  @assert( size(bndryfaces,1) == size(flux,3) )
  wflux = zeros(Tflx, size(res,1))
  @inbounds begin
    for (bindex, bndry) in enumerate(bndryfaces)
      for i = 1:sbpface.numnodes
        for field = 1:size(res,1)
          wflux[field] = sbpface.wface[i]*flux[field,i,bindex]
        end
        for j = 1:sbpface.stencilsize
          for field = 1:size(res,1)
            res[field,sbpface.perm[j,bndry.face],bndry.element] +=
              sbpface.interp[j,i]*wflux[field]
          end
        end
      end
    end
  end
end

@doc """
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

**In/Outs**

* `fun`: functional value (or vector) being *contributed* to by the integration

**Returns**

* `fun`: in the case of the scalar version, the functional value is returned

"""->

function integratefunctional!{Tsbp,Tflx}(sbpface::AbstractFace{Tsbp},
                                         bndryfaces::Array{Boundary},
                                         flux::AbstractArray{Tflx,2})
  @assert( size(sbpface.interp,2) == size(flux,1) )
  @assert( size(bndryfaces,1) == size(flux,2) )
  fun = zero(Tflx)
  @inbounds begin
    for (bindex, bndry) in enumerate(bndryfaces)
      for i = 1:sbpface.numnodes
        fun += sbpface.wface[i]*flux[i,bindex]
      end
    end
  end
  return fun
end

function integratefunctional!{Tsbp,Tflx,Tfun}(sbpface::AbstractFace{Tsbp},
                                              bndryfaces::Array{Boundary},
                                              flux::AbstractArray{Tflx,3},
                                              fun::AbstractArray{Tfun,1})
  @assert( size(sbpface.interp,2) == size(flux,2) )
  @assert( size(flux,1) == size(fun,1) )
  @assert( size(bndryfaces,1) == size(flux,3) )
  @inbounds begin
    for (bindex, bndry) in enumerate(bndryfaces)
      for i = 1:sbpface.numnodes
        for field = 1:size(fun,1)
          fun[field] += sbpface.wface[i]*flux[field,i,bindex]
        end
      end
    end
  end
end

# function boundaryintegrate!{Tsbp,Tsol,Tres}(sbpface::AbstractFace{Tsbp},
#                                             bndryfaces::Array{Boundary},
#                                             u::AbstractArray{Tsol,2},
#                                             fluxfunc::FluxFunction{Tsol},
#                                             res::AbstractArray{Tres,2})
#   @assert( size(sol,1) == size(res,1) )
#   @assert( size(sbpface.interp,2) <= size(sol,1) )
#   @inbounds begin
#     for (bindex, bndry) in enumerate(bndryfaces)
#       for i = 1:sbpface.numnodes
#         # interpolate to face
#         uface = 0.0
#         for j = 1:sbpface.stencilsize
#           uface += sbpface.interp[i,j]*u[sbpface.perm[j,bndry.face],
#                                          bndry.element]
#         end
#         # compute flux at this face node
#         flux = 0.0
#         fluxfunc(uface, sbpface.normal, flux)
#         # return to volume nodes
#         for j = 1:sbpface.numvolnodes
#           res[sbpface.perm[j,bndry.face],bndry.element] +=
#             sbpface.wface[i]*sbpface.interp[i,j]*flux
#         end
#       end
#     end
#   end
# end


@doc """
### SummationByParts.interiorfaceinterpolate!

Interpolates element-node data to the element-face cubature nodes for the given
set of faces.  Different methods are available depending on the rank of `uvol`:

* For *scalar* fields, the dimensions of `uvol` correspond to [node index,
  element index] and the dimensions of `uface` correspond to [-/+ side,
  face-node index, interface index]
* For *vector* fields, the dimensions of `uvol` correspond to [field index,
  node index, element index] and the dimensions of `uface` correspond to
  [field index, -/+ side, face-node index, interface index]

**Inputs**

* `sbpface`: an SBP face operator type
* `ifaces`: list of interior faces stored as an array of `Interface`s
* `uvol`: array of field data that is being interpolated

**In/Outs**

* `uface`: field data interpolated to the faces

"""->
function interiorfaceinterpolate!{Tsbp,Tsol}(sbpface::AbstractFace{Tsbp},
                                             ifaces::Array{Interface},
                                             uvol::AbstractArray{Tsol,2},
                                             uface::AbstractArray{Tsol,3})
  @assert( size(sbpface.interp,1) <= size(uvol,1) )
  @assert( size(sbpface.interp,2) == size(uface,2) )
  @inbounds begin
    for (findex, face) in enumerate(ifaces)
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,face.orient]
        uface[1,i,findex] = zero(Tsol)
        uface[2,i,findex] = zero(Tsol)
        for j = 1:sbpface.stencilsize
          uface[1,i,findex] += sbpface.interp[j,i]*uvol[sbpface.perm[j,face.faceL],
                                                        face.elementL]
          uface[2,i,findex] += sbpface.interp[j,iR]*uvol[sbpface.perm[j,face.faceR],
                                                         face.elementR]
        end
      end
    end
  end
end

function interiorfaceinterpolate!{Tsbp,Tsol}(sbpface::AbstractFace{Tsbp},
                                             ifaces::Array{Interface},
                                             uvol::AbstractArray{Tsol,3},
                                             uface::AbstractArray{Tsol,4})
  @assert( size(uvol,1) == size(uface,1) )
  @assert( size(sbpface.interp,1) <= size(uvol,2) )
  @assert( size(sbpface.interp,2) == size(uface,3) )
  @inbounds begin
    for (findex, face) in enumerate(ifaces)
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,face.orient]
        for field=1:size(uvol, 1)
          uface[field,1,i,findex] = zero(Tsol)
          uface[field,2,i,findex] = zero(Tsol)
        end
        for j = 1:sbpface.stencilsize
          for field = 1:size(uvol,1)
            uface[field,1,i,findex] += sbpface.interp[j,i]*
            uvol[field,sbpface.perm[j,face.faceL],face.elementL]
            uface[field,2,i,findex] += sbpface.interp[j,iR]*
            uvol[field,sbpface.perm[j,face.faceR],face.elementR]
          end
        end
      end
    end
  end
end

@doc """
### SummationByParts.interiorfaceintegrate!

Integrates a given flux over interior element interfaces using appropriate mass
matrices defined on the element faces.  Different methods are available
depending on the rank of `flux`:

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

* `sbp`: an SBP operator type
* `ifaces`: list of element interfaces stored as an array of `Interface`s
* `flux`: array of flux data that is being integrated

**In/Outs**

* `res`: where the result of the integration is stored

**WARNING**: the order of the interfaces in `ifaces` and `flux` must be
  consistent.

"""->
function interiorfaceintegrate!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp},
                                                ifaces::Array{Interface},
                                                flux::AbstractArray{Tflx,2},
                                                res::AbstractArray{Tres,2})
  @assert( sbp.numnodes == size(res,1) )
  @assert( sbp.numfacenodes == size(flux,1) )
  @assert( size(ifaces,1) == size(flux,2) )
  # JEH: temporary, until nbrnodeindex is part of sbp type
  nbrnodeindex = [sbp.numfacenodes:-1:1;]
  @inbounds begin
    for (findex, face) in enumerate(ifaces)
      for i = 1:sbp.numfacenodes    
        for j = 1:sbp.numfacenodes
          jL = sbp.facenodes[j, face.faceL]::Int
          jR = sbp.facenodes[nbrnodeindex[j], face.faceR]::Int
          res[jL,face.elementL] += sbp.wface[j,i]*flux[i,findex]
          res[jR,face.elementR] -= sbp.wface[j,i]*flux[i,findex]
        end
      end
    end
  end
end

function interiorfaceintegrate!{Tsbp,Tflx,Tres}(sbp::AbstractSBP{Tsbp},
                                                ifaces::Array{Interface},
                                                flux::AbstractArray{Tflx,3},
                                                res::AbstractArray{Tres,3})
  @assert( size(res,1) == size(flux,1) )
  @assert( sbp.numnodes == size(res,2) )
  @assert( sbp.numfacenodes == size(flux,2) )
  @assert( size(ifaces,1) == size(flux,3) )
  # JEH: temporary, until nbrnodeindex is part of sbp type
  nbrnodeindex = [sbp.numfacenodes:-1:1;]
  @inbounds begin
    for (findex, face) in enumerate(ifaces)
      for i = 1:sbp.numfacenodes    
        for j = 1:sbp.numfacenodes
          jL = sbp.facenodes[j, face.faceL]::Int
          jR = sbp.facenodes[nbrnodeindex[j], face.faceR]::Int
          for field = 1:size(res,1)
            res[field,jL,face.elementL] += sbp.wface[j,i]*flux[field,i,findex]
            res[field,jR,face.elementR] -= sbp.wface[j,i]*flux[field,i,findex]
          end
        end
      end
    end
  end
end

function interiorfaceintegrate!{Tsbp,Tflx,Tres}(sbpface::AbstractFace{Tsbp},
                                                ifaces::Array{Interface},
                                                flux::AbstractArray{Tflx,2},
                                                res::AbstractArray{Tres,2})
  @assert( size(sbpface.interp,1) <= size(res,1) )
  @assert( size(sbpface.interp,2) == size(flux,1) )
  @assert( size(ifaces,1) == size(flux,2) )
  @inbounds begin
    for (findex, face) in enumerate(ifaces)
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,face.orient]
        for j = 1:sbpface.stencilsize
          res[sbpface.perm[j,face.faceL],face.elementL] += 
            sbpface.interp[j,i]*sbpface.wface[i]*flux[i,findex]
          res[sbpface.perm[j,face.faceR],face.elementR] -=
            sbpface.interp[j,iR]*sbpface.wface[iR]*flux[i,findex]
        end
      end
    end
  end
end

function interiorfaceintegrate!{Tsbp,Tflx,Tres}(sbpface::AbstractFace{Tsbp},
                                                ifaces::Array{Interface},
                                                flux::AbstractArray{Tflx,3},
                                                res::AbstractArray{Tres,3})
  @assert( size(res,1) == size(flux,1) )  
  @assert( size(sbpface.interp,1) <= size(res,2) )
  @assert( size(sbpface.interp,2) == size(flux,2) )
  @assert( size(ifaces,1) == size(flux,3) )
  @inbounds begin
    for (findex, face) in enumerate(ifaces)
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,face.orient]
        for j = 1:sbpface.stencilsize
          for field = 1:size(res,1)
            res[field,sbpface.perm[j,face.faceL],face.elementL] += 
              sbpface.interp[j,i]*sbpface.wface[i]*flux[field,i,findex]
            res[field,sbpface.perm[j,face.faceR],face.elementR] -=
              sbpface.interp[j,iR]*sbpface.wface[iR]*flux[field,i,findex]
          end
        end
      end
    end
  end
end

@doc """
### SummationByParts.mappingjacobian!

Evaluates the Jacobian of the mapping from face-reference coordinates to
physical coordinates, as well as the determinant of the Jacobian.

**Inputs**

* `sbpface`: an SBP face operator type
* `ifaces`: list of element interfaces stored as an array of `Interface`s
* `x`: the physical coordinates in [coord, node, elem] format

**In/Outs**

* `dξdx`: the Jacobian in [ξ coord, x coord, face node, L/R, face] format
* `jac`: the determinant in [face node, L/R, face] format

"""->
function mappingjacobian!{Tsbp,Tmsh}(sbpface::AbstractFace{Tsbp},
                                     ifaces::Array{Interface},
                                     x::AbstractArray{Tmsh,3},
                                     dξdx::AbstractArray{Tmsh,5},
                                     jac::AbstractArray{Tmsh,3})
  @assert( size(ifaces,1) == size(dξdx,5) == size(jac,3) )
  @assert( size(dξdx,4) == size(jac,2) == 2 )
  @assert( size(dξdx,3) == size(jac,1) == sbpface.numnodes )
  @assert( size(x,1) == 2 && size(dξdx,1) == 2 && size(dξdx,2) == 2 )
  @inbounds begin
    for (findex, face) in enumerate(ifaces)
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,face.orient]
        # compute the derivatives dxdξ
        dxdξL = zeros(Tmsh, (2,2))
        dxdξR = zeros(Tmsh, (2,2))
        for di = 1:2
          for di2 = 1:2
            for j = 1:sbpface.dstencilsize
              dxdξL[di,di2] +=
                sbpface.deriv[j,i,di2]*x[di,sbpface.dperm[j,face.faceL],
                                         face.elementL]
              dxdξR[di,di2] +=
                sbpface.deriv[j,iR,di2]*x[di,sbpface.dperm[j,face.faceR],
                                          face.elementR]
            end
          end
        end
        # compute the Jacobian determinants
        jac[i,1,findex] = one(Tmsh)/(dxdξL[1,1]*dxdξL[2,2] - 
                                     dxdξL[1,2]*dxdξL[2,1])
        jac[iR,2,findex] = one(Tmsh)/(dxdξR[1,1]*dxdξR[2,2] - 
                                      dxdξR[1,2]*dxdξR[2,1])

        # compute the derivatives dξdx
        dξdx[1,1,i,1,findex] = dxdξL[2,2]*jac[i,1,findex]
        dξdx[1,2,i,1,findex] = -dxdξL[1,2]*jac[i,1,findex]
        dξdx[2,2,i,1,findex] = dxdξL[1,1]*jac[i,1,findex]
        dξdx[2,1,i,1,findex] = -dxdξL[2,1]*jac[i,1,findex]

        dξdx[1,1,iR,2,findex] = dxdξR[2,2]*jac[iR,2,findex]
        dξdx[1,2,iR,2,findex] = -dxdξR[1,2]*jac[iR,2,findex]
        dξdx[2,2,iR,2,findex] = dxdξR[1,1]*jac[iR,2,findex]
        dξdx[2,1,iR,2,findex] = -dxdξR[2,1]*jac[iR,2,findex]
      end
    end
  end
end

@doc """
### SummationByParts.edgestabilize!

Applies edge stabilization to a given field, differentiating in the direction
specified by `dirvec`, and scaling by the `tau` field.

**Inputs**

* `sbpface`: an SBP face operator type
* `ifaces`: list of element interfaces stored as an array of `Interface`s
* `dirvec`: direction to differentiate in [xi coord, face node, L/R, face] format
* `tau`: scaling term in [face node, face] format
* `u`: field being stablized in [vol node, element] format

**In/Outs**

* `res`: where the result is stored in [vol node, element] format

"""->
function edgestabilize!{Tsbp,Tmsh,Tsol}(sbpface::AbstractFace{Tsbp},
                                        ifaces::Array{Interface},
                                        dirvec::Array{Tmsh,4},
                                        tau::Array{Tmsh,2},
                                        u::Array{Tsol,2},
                                        res::Array{Tsol,2})
  @assert( size(u) == size(res) )
  @assert( size(ifaces,1) == size(dirvec,4) == size(tau,2) )
  @assert( size(dirvec,2) == size(tau,1) == sbpface.numnodes )
  @assert( size(dirvec,3) == 2 )
  dudξ = zeros(Tsol, (2))
  const left = 1
  const right = 1
  @inbounds begin
    for (findex, face) in enumerate(ifaces)
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,face.orient]
        Du = zero(Tsol)
        for di = 1:2
          fill!(dudξ, zero(Tsol))
          # compute the derivatives in the ξ[di] direction
          for j = 1:sbpface.dstencilsize
            dudξ[left] += sbpface.deriv[j,i,di]*u[sbpface.dperm[j,face.faceL],
                                                  face.elementL]
            dudξ[right] += sbpface.deriv[j,iR,di]*u[sbpface.dperm[j,face.faceR],
                                                    face.elementR]
          end
          # contract with direction vector
          Du += dudξ[left]*dirvec[di,i,left,findex] + 
          dudξ[right]*dirvec[di,iR,right,findex]
        end
        # scale by tau and face cubature
        Du *= sbpface.wface[i]*tau[i,findex]
        # now apply transposed jump derivative
        for di = 1:2
          dudξ[left] = Du*dirvec[di,i,left,findex]
          dudξ[right] = Du*dirvec[di,iR,right,findex]
          for j = 1:sbpface.dstencilsize
            res[sbpface.dperm[j,face.faceL],face.elementL] +=
              sbpface.deriv[j,i,di]*dudξ[left]
            res[sbpface.dperm[j,face.faceR],face.elementR] += 
              sbpface.deriv[j,iR,di]*dudξ[right]
          end
        end
      end
    end
  end
end
