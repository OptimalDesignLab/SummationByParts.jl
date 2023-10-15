# This file gathers together functions related to interpolation of quantities to
# the element faces

"""
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

* `sbpface`: an SBP AbstractFace operator type
* `bndryfaces`: list of boundary faces stored as an array of `Boundary`s
* `uvol`: array of field data that is being interpolated

**In/Outs**

* `uface`: field data interpolated to the faces

"""
function boundaryinterpolate!(sbpface::DenseFace{Tsbp},
                                         bndryfaces::AbstractArray{Boundary},
                                         uvol::AbstractArray{Tsol,2},
                                         uface::AbstractArray{Tsol,2}) where {Tsbp,Tsol}
  @assert( size(sbpface.interp,1) <= size(uvol,1) )
  @assert( size(sbpface.interp,2) == size(uface,1) )
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

function boundaryinterpolate!(sbpface::DenseFace{Tsbp},
                                         bndryfaces::AbstractArray{Boundary},
                                         uvol::AbstractArray{Tsol,3},
                                         uface::AbstractArray{Tsol,3}) where {Tsbp,Tsol}
  @assert( size(uvol,1) == size(uface,1) )
  @assert( size(sbpface.interp,1) <= size(uvol,2) )
  @assert( size(sbpface.interp,2) == size(uface,2) )
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

function boundaryinterpolate!(sbpface::SparseFace{Tsbp},
                                         bndryfaces::AbstractArray{Boundary},
                                         uvol::AbstractArray{Tsol,2},
                                         uface::AbstractArray{Tsol,2}) where {Tsbp,Tsol}
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      uface[i,bindex] = uvol[sbpface.perm[i,bndry.face],bndry.element]
    end
  end
end

function boundaryinterpolate!(sbpface::SparseFace{Tsbp},
                                         bndryfaces::AbstractArray{Boundary},
                                         uvol::AbstractArray{Tsol,3},
                                         uface::AbstractArray{Tsol,3}) where {Tsbp,Tsol}
  @assert( size(uvol,1) == size(uface,1) )
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      for field=1:size(uvol, 1)
        uface[field,i,bindex] = uvol[field,sbpface.perm[i,bndry.face],
                                     bndry.element]
      end
    end
  end
end

"""
### SummationByParts.boundaryFaceInterpolate!

This is the single-face variant of boundaryinterpolate!.  Interpolates vector
field values at the nodes of a given element to a specified face of the
element. Different methods are available depending on the rank of `uvol`:

* For *scalar* fields, it is assumed that `uvol` (`uface`) is a rank-1 array,
with the first and only dimension for the node index.
* For *vector* fields, `uvol` (`uface`) is a rank-2 array, with the first
dimension for the index of the vector field, and the second dimension for the
node index.

**Inputs**

* `sbpface`: an SBP AbstractFace operator
* `face`: the face of the element to interpolate to
* `uvol`: array of field data that is being interpolated

**In/Outs**
  * `uface`: field data interpolated to the faces
  
"""
function boundaryFaceInterpolate!(sbpface::DenseFace{Tsbp},
                                             face::Integer,
                                             uvol::AbstractArray{Tsol,1},
                                             uface::AbstractArray{Tsol,1}) where {Tsbp,Tsol}
  @assert( size(sbpface.interp,1) <= size(uvol,1) )
  @assert( size(sbpface.interp,2) == size(uface,1) )
  for i = 1:sbpface.numnodes
    uface[i] = zero(Tsol)
    for j = 1:sbpface.stencilsize
      uface[i] += sbpface.interp[j,i]*uvol[sbpface.perm[j,face]]
    end
  end
end

function boundaryFaceInterpolate!(sbpface::DenseFace{Tsbp},
                                             face::Integer,
                                             uvol::AbstractArray{Tsol,2},
                                             uface::AbstractArray{Tsol,2}) where {Tsbp,Tsol}
  @assert( size(uvol,1) == size(uface,1) )
  @assert( size(sbpface.interp,1) <= size(uvol,2) )
  @assert( size(sbpface.interp,2) == size(uface,2) )
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

function boundaryFaceInterpolate!(sbpface::SparseFace{Tsbp},
                                             face::Integer,
                                             uvol::AbstractArray{Tsol,1},
                                             uface::AbstractArray{Tsol,1}) where {Tsbp,Tsol}
  for i = 1:sbpface.numnodes
    uface[i] = uvol[sbpface.perm[i,face]]
  end
end

function boundaryFaceInterpolate!(sbpface::SparseFace{Tsbp},
                                             face::Integer,
                                             uvol::AbstractArray{Tsol,2},
                                             uface::AbstractArray{Tsol,2}) where {Tsbp,Tsol}
  @assert( size(uvol,1) == size(uface,1) )
  for i = 1:sbpface.numnodes
    for field=1:size(uvol, 1)
      uface[field,i] = uvol[field,sbpface.perm[i,face]]
    end
  end
end

"""
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

"""
function interiorfaceinterpolate!(sbpface::DenseFace{Tsbp},
                                             ifaces::AbstractArray{Interface},
                                             uvol::AbstractArray{Tsol,2},
                                             uface::AbstractArray{Tsol,3}) where {Tsbp,Tsol}
  @assert( size(sbpface.interp,1) <= size(uvol,1) )
  @assert( size(sbpface.interp,2) == size(uface,2) )
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

function interiorfaceinterpolate!(sbpface::DenseFace{Tsbp},
                                             ifaces::AbstractArray{Interface},
                                             uvol::AbstractArray{Tsol,3},
                                             uface::AbstractArray{Tsol,4}) where {Tsbp,Tsol}
  @assert( size(uvol,1) == size(uface,1) )
  @assert( size(sbpface.interp,1) <= size(uvol,2) )
  @assert( size(sbpface.interp,2) == size(uface,3) )
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

function interiorfaceinterpolate!(sbpface::SparseFace{Tsbp},
                                             ifaces::AbstractArray{Interface},
                                             uvol::AbstractArray{Tsol,2},
                                             uface::AbstractArray{Tsol,3}) where {Tsbp,Tsol}
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]     
      uface[1,i,findex] = uvol[sbpface.perm[i,face.faceL],face.elementL]
      uface[2,i,findex] = uvol[sbpface.perm[iR,face.faceR],face.elementR]
    end
  end
end

function interiorfaceinterpolate!(sbpface::SparseFace{Tsbp},
                                             ifaces::AbstractArray{Interface},
                                             uvol::AbstractArray{Tsol,3},
                                             uface::AbstractArray{Tsol,4}) where {Tsbp,Tsol}
  @assert( size(uvol,1) == size(uface,1) )
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]
      for field=1:size(uvol, 1)
        uface[field,1,i,findex] =
          uvol[field,sbpface.perm[i,face.faceL],face.elementL]
        uface[field,2,i,findex] =
          uvol[field,sbpface.perm[iR,face.faceR],face.elementR]
      end
    end
  end
end

"""
### SummationByParts.interiorFaceInterpolate!

This is the single face variant of interiorFaceInterpolate!  Interpolates
element-node data to the element-face cubature nodes for the face; the nominal
*right* element will have its face data ordered consistently with that of the
left.  Different methods are available depending on the rank of `uvol`:

* For *scalar* fields, the only dimension of `uL` and `uR` corresponds to [node
  index] and the dimension of `ufaceL` and `ufaceR` corresponds to [face-node
  index]
* For *vector* fields, the dimensions of `uL` and `uR` correspond to [field
  index, node index] and the dimensions of `ufaceL` and `ufaceR` correspond to
  [field index, face-node index]

**Inputs**

* `sbpface`: an SBP face operator type
* `iface`: Interface type that specifies face indices and relative orientation
* `uL`: field data that is being interpolated from the nominal *left* element
* `uR`: field data that is being interpolated from the nominal *right* element

**In/Outs**

* `ufaceL`: field data from `uL` interpolated to the face
* `ufaceR`: field data from `uR` interpolated to the face

"""
function interiorFaceInterpolate!(sbpface::DenseFace{Tsbp},
                                             iface::Interface,
                                             uL::AbstractArray{Tsol,1},
                                             uR::AbstractArray{Tsol,1},
                                             ufaceL::AbstractArray{Tsol,1},
                                             ufaceR::AbstractArray{Tsol,1}) where {Tsbp,Tsol}
  @assert( size(sbpface.interp,1) <= size(uL,1) )
  @assert( size(sbpface.interp,1) <= size(uR,1) )
  @assert( size(sbpface.interp,2) == size(ufaceL,1) == size(ufaceR,1) )
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    ufaceL[i] = zero(Tsol)
    ufaceR[i] = zero(Tsol)
    for j = 1:sbpface.stencilsize
      ufaceL[i] += sbpface.interp[j,i]*uL[sbpface.perm[j,iface.faceL]]
      ufaceR[i] += sbpface.interp[j,iR]*uR[sbpface.perm[j,iface.faceR]]
    end
  end
end

function interiorFaceInterpolate!(sbpface::DenseFace{Tsbp},
                                             iface::Interface,
                                             uL::AbstractArray{Tsol,2},
                                             uR::AbstractArray{Tsol,2},
                                             ufaceL::AbstractArray{Tsol,2},
                                             ufaceR::AbstractArray{Tsol,2}) where {Tsbp,Tsol}
  @assert( size(uL,1) == size(ufaceL,1) == size(uR,1) == size(ufaceR,1) )
  @assert( size(sbpface.interp,1) <= size(uL,2) )
  @assert( size(sbpface.interp,1) <= size(uR,2) )
  @assert( size(sbpface.interp,2) == size(ufaceL,2) == size(ufaceR,2) )
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    @simd for field=1:size(uL,1)
      ufaceL[field,i] = zero(Tsol)
      ufaceR[field,i] = zero(Tsol)
    end
    for j = 1:sbpface.stencilsize
      for field = 1:size(uL,1)
        ufaceL[field,i] += sbpface.interp[j,i]*
        uL[field,sbpface.perm[j,iface.faceL]]
        ufaceR[field,i] += sbpface.interp[j,iR]*
        uR[field,sbpface.perm[j,iface.faceR]]
      end
    end
  end
end

function interiorFaceInterpolate!(sbpface::SparseFace{Tsbp},
                                             iface::Interface,
                                             uL::AbstractArray{Tsol,1},
                                             uR::AbstractArray{Tsol,1},
                                             ufaceL::AbstractArray{Tsol,1},
                                             ufaceR::AbstractArray{Tsol,1}) where {Tsbp,Tsol}
  @assert( size(ufaceL,1) == size(ufaceR,1) )
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    ufaceL[i] = uL[sbpface.perm[i,iface.faceL]]
    ufaceR[i] = uR[sbpface.perm[iR,iface.faceR]]
  end
end

function interiorFaceInterpolate!(sbpface::SparseFace{Tsbp},
                                             iface::Interface,
                                             uL::AbstractArray{Tsol,2},
                                             uR::AbstractArray{Tsol,2},
                                             ufaceL::AbstractArray{Tsol,2},
                                             ufaceR::AbstractArray{Tsol,2}) where {Tsbp,Tsol}
  @assert( size(uL,1) == size(ufaceL,1) == size(uR,1) == size(ufaceR,1) )
  @assert( size(ufaceL,2) == size(ufaceR,2) )
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    @simd for field=1:size(uL,1)
      ufaceL[field,i] = uL[field,sbpface.perm[i,iface.faceL]]
      ufaceR[field,i] = uR[field,sbpface.perm[iR,iface.faceR]]
    end
  end
end
