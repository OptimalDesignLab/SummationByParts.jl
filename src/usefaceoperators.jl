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

**WARNING**: the order of the boundaries in `bndryfaces` and `uface` must be
  consistent.

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
        uface[:,i,bindex] = zeros(Tsol, size(uvol,1))
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
        for j = 1:sbpface.stencilsize
          res[sbpface.perm[j,bndry.face],bndry.element] +=
            sbpface.interp[j,i]*sbpface.wface[i]*flux[j,bindex]
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
  @inbounds begin
    for (bindex, bndry) in enumerate(bndryfaces)
      for i = 1:sbpface.numnodes
        for j = 1:sbpface.stencilsize
          for field = 1:size(res,1)
            res[field,sbpface.perm[j,bndry.face],bndry.element] +=
              sbpface.interp[j,i]*sbpface.wface[i]*flux[field,j,bindex]
          end
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

* For *scalar* fields, it is assumed that `uvol` (`uface`) is a rank-2 array
(rank-3 array), with the first dimension for the node index, and the second
(third) dimension for the element index (interface index).
* For *vector* fields, `uvol` (`uface`) is a rank-3 array (rank-4 array), with
the first dimension for the index of the vector field, the second dimension for
the node index, and the third dimension (fourth dimension) for the element index
(interface index).

**Inputs**

* `sbpface`: an SBP face operator type
* `ifaces`: list of interior faces stored as an array of `Interface`s
* `uvol`: array of field data that is being interpolated

**In/Outs**

* `uface`: field data interpolated to the faces

**WARNING**: the order of the interfaces in `ifaces` and `uface` must be
  consistent.

"""->
function interiorfaceinterpolate!{Tsbp,Tsol}(sbpface::AbstractFace{Tsbp},
                                             ifaces::Array{Interface},
                                             uvol::AbstractArray{Tsol,2},
                                             uface::AbstractArray{Tsol,3})
  @assert( size(sbpface.interp,1) <= size(uvol,1) )
  @assert( size(sbpface.interp,2) == size(uface,1) )
  @inbounds begin
    for (findex, face) in enumerate(ifaces)
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,face.orient]
        uface[i,1,findex] = zero(Tsol)
        uface[iR,2,findex] = zero(Tsol)
        for j = 1:sbpface.stencilsize
          # note that uface[i,1,findex] and uface[i,2,findex] do not necessarily
          # correspond to the same node in this case.
          uface[i,1,findex] += sbpface.interp[j,i]*uvol[sbpface.perm[j,face.faceL],
                                                        face.elementL]
          uface[iR,2,findex] += sbpface.interp[j,iR]*uvol[sbpface.perm[j,face.faceR],
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
  @assert( size(sbpface.interp,2) == size(uface,2) )
  @inbounds begin
    for (findex, face) in enumerate(ifaces)
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,face.orient]
        uface[:,i,1,findex] = zeros(Tsol, size(uvol,1))
        uface[:,iR,2,findex] = zeros(Tsol, size(uvol,1))
        for j = 1:sbpface.stencilsize
          for field = 1:size(uvol,1)
            # note that uface[i,1,findex] and uface[i,2,findex] do not
            # necessarily correspond to the same node in this case.
            uface[field,i,1,findex] += sbpface.interp[j,i]*
            uvol[field,sbpface.perm[j,face.faceL],face.elementL]
            uface[field,iR,2,findex] += sbpface.interp[j,iR]*
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
              sbpface.interp[j,i]*sbpface.wface[i]*flux[i,findex]
            res[field,sbpface.perm[j,face.faceR],face.elementR] -=
              sbpface.interp[j,iR]*sbpface.wface[iR]*flux[i,findex]
          end
        end
      end
    end
  end
end