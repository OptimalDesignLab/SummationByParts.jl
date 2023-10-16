# This file contains the reverse-mode version of the methods in
# faceinterpolate.jl

"""
### SummationByParts.boundaryinterpolate_rev!

This is the reverse differentiated version of boundaryinterpolate!.  See
faceinterpolate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `uvol` variable.

**Inputs**

* `sbpface`: an SBP AbstractFace operator type
* `bndryfaces`: list of boundary faces stored as an array of `Boundary`s
* `uface_bar`: vector applied to the left of the interpolation operator (R)

**In/Outs**

* `uvol_bar`: the result of the vector matrix product between R and `uface_bar`

"""
function boundaryinterpolate_rev!(sbpface::DenseFace{Tsbp},
                                             bndryfaces::AbstractArray{Boundary},
                                             uvol_bar::AbstractArray{Tsol,2},
                                             uface_bar::AbstractArray{Tsol,2}) where {Tsbp,Tsol}
  @assert( size(sbpface.interp,1) <= size(uvol_bar,1) )
  @assert( size(sbpface.interp,2) == size(uface_bar,1) )
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      for j = 1:sbpface.stencilsize
        # uface[i,bindex] += sbpface.interp[j,i]*uvol[sbpface.perm[j,bndry.face],
        #                                             bndry.element]
        uvol_bar[sbpface.perm[j,bndry.face],bndry.element] +=
          sbpface.interp[j,i]*uface_bar[i,bindex]
      end
    end
  end
end

function boundaryinterpolate_rev!(sbpface::DenseFace{Tsbp},
                                             bndryfaces::AbstractArray{Boundary},
                                             uvol_bar::AbstractArray{Tsol,3},
                                             uface_bar::AbstractArray{Tsol,3}) where {Tsbp,Tsol}
  @assert( size(uvol_bar,1) == size(uface_bar,1) )
  @assert( size(sbpface.interp,1) <= size(uvol_bar,2) )
  @assert( size(sbpface.interp,2) == size(uface_bar,2) )
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      for j = 1:sbpface.stencilsize
         for field = 1:size(uvol_bar,1)
           # uface[field,i,bindex] += sbpface.interp[j,i]*
           # uvol[field,sbpface.perm[j,bndry.face],bndry.element]
           uvol_bar[field,sbpface.perm[j,bndry.face],bndry.element] +=
             sbpface.interp[j,i]*uface_bar[field,i,bindex]
         end
      end
    end
  end
end

function boundaryinterpolate_rev!(sbpface::SparseFace{Tsbp},
                                             bndryfaces::AbstractArray{Boundary},
                                             uvol_bar::AbstractArray{Tsol,2},
                                             uface_bar::AbstractArray{Tsol,2}) where {Tsbp,Tsol}
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      # uface[i,bindex] = uvol[sbpface.perm[i,bndry.face],bndry.element]
      uvol_bar[sbpface.perm[i,bndry.face],bndry.element] += uface_bar[i,bindex]
    end
  end
end

function boundaryinterpolate_rev!(sbpface::SparseFace{Tsbp},
                                             bndryfaces::AbstractArray{Boundary},
                                             uvol_bar::AbstractArray{Tsol,3},
                                             uface_bar::AbstractArray{Tsol,3}) where {Tsbp,Tsol}
  @assert( size(uvol_bar,1) == size(uface_bar,1) )
  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      for field=1:size(uvol_bar, 1)
        # uface[field,i,bindex] = uvol[field,sbpface.perm[i,bndry.face],
        #                              bndry.element]
        uvol_bar[field,sbpface.perm[i,bndry.face], bndry.element] +=
          uface_bar[field,i,bindex]
      end
    end
  end
end

"""
### SummationByParts.boundaryFaceInterpolate_rev!

This is the reverse differentiated version of boundaryFaceInterpolate!.  See
faceinterpolate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `uvol` variable.

**Inputs**

* `sbpface`: an SBP AbstractFace operator
* `face`: the face of the element to interpolate to
* `uface_bar`: vector applied to the left of the interpolation operator (R)

**In/Outs**

* `uvol_bar`: the result of the vector-matrix product between R and `uface_bar`
  
"""
function boundaryFaceInterpolate_rev!(sbpface::DenseFace{Tsbp},
                                                 face::Integer,
                                                 uvol_bar::AbstractArray{Tsol,1},
                                                 uface_bar::AbstractArray{Tsol,1}) where {Tsbp,Tsol}
  @assert( size(sbpface.interp,1) <= size(uvol_bar,1) )
  @assert( size(sbpface.interp,2) == size(uface_bar,1) )
  for i = 1:sbpface.numnodes
    for j = 1:sbpface.stencilsize
      # uface[i] += sbpface.interp[j,i]*uvol[sbpface.perm[j,face]]
      uvol_bar[sbpface.perm[j,face]] += sbpface.interp[j,i]*uface_bar[i]
    end
  end
end

function boundaryFaceInterpolate_rev!(sbpface::DenseFace{Tsbp},
                                                 face::Integer,
                                                 uvol_bar::AbstractArray{Tsol,2},
                                                 uface_bar::AbstractArray{Tsol,2}) where {Tsbp,Tsol}
  @assert( size(uvol_bar,1) == size(uface_bar,1) )
  @assert( size(sbpface.interp,1) <= size(uvol_bar,2) )
  @assert( size(sbpface.interp,2) == size(uface_bar,2) )
  for i = 1:sbpface.numnodes
    for j = 1:sbpface.stencilsize
       for field = 1:size(uvol_bar,1)
         # uface[field,i] += sbpface.interp[j,i]*uvol[field,sbpface.perm[j,face]]
         uvol_bar[field,sbpface.perm[j,face]] +=
           sbpface.interp[j,i]*uface_bar[field,i]
       end
    end
  end
end

function boundaryFaceInterpolate_rev!(sbpface::SparseFace{Tsbp},
                                                 face::Integer,
                                                 uvol_bar::AbstractArray{Tsol,1},
                                                 uface_bar::AbstractArray{Tsol,1}) where {Tsbp,Tsol}
  for i = 1:sbpface.numnodes
    # uface[i] = uvol[sbpface.perm[i,face]]
    uvol_bar[sbpface.perm[i,face]] += uface_bar[i] 
  end
end

function boundaryFaceInterpolate_rev!(sbpface::SparseFace{Tsbp},
                                                 face::Integer,
                                                 uvol_bar::AbstractArray{Tsol,2},
                                                 uface_bar::AbstractArray{Tsol,2}) where {Tsbp,Tsol}
  @assert( size(uvol_bar,1) == size(uface_bar,1) )
  for i = 1:sbpface.numnodes
    for field=1:size(uvol_bar, 1)
      # uface[field,i] = uvol[field,sbpface.perm[i,face]]
      uvol_bar[field,sbpface.perm[i,face]] += uface_bar[field,i]
    end
  end
end

"""
### SummationByParts.interiorfaceinterpolate_rev!

This is the reverse differentiated version of interiorfaceinterpolate!.  See
faceinterpolate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `uvol` variable.

**Inputs**

* `sbpface`: an SBP face operator type
* `ifaces`: list of interior faces stored as an array of `Interface`s
* `uface_bar`: vector applied to the left of the interpolation operator (R)

**In/Outs**

* `uvol_bar`: the result of the vector-matrix product between R and `uface_bar`

"""
function interiorfaceinterpolate_rev!(sbpface::DenseFace{Tsbp},
                                             ifaces::AbstractArray{Interface},
                                             uvol_bar::AbstractArray{Tsol,2},
                                             uface_bar::AbstractArray{Tsol,3}) where {Tsbp,Tsol}
  @assert( size(sbpface.interp,1) <= size(uvol_bar,1) )
  @assert( size(sbpface.interp,2) == size(uface_bar,2) )
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]
      for j = 1:sbpface.stencilsize
        # uface[1,i,findex] += sbpface.interp[j,i]*uvol[sbpface.perm[j,face.faceL],
        #                                               face.elementL]
        uvol_bar[sbpface.perm[j,face.faceL],face.elementL] +=
          sbpface.interp[j,i]*uface_bar[1,i,findex]
        # uface[2,i,findex] += sbpface.interp[j,iR]*uvol[sbpface.perm[j,face.faceR],
        #                                                face.elementR]
        uvol_bar[sbpface.perm[j,face.faceR],face.elementR] +=
          sbpface.interp[j,iR]*uface_bar[2,i,findex]
      end
    end
  end
end

function interiorfaceinterpolate_rev!(sbpface::DenseFace{Tsbp},
                                             ifaces::AbstractArray{Interface},
                                             uvol_bar::AbstractArray{Tsol,3},
                                             uface_bar::AbstractArray{Tsol,4}) where {Tsbp,Tsol}
  @assert( size(uvol_bar,1) == size(uface_bar,1) )
  @assert( size(sbpface.interp,1) <= size(uvol_bar,2) )
  @assert( size(sbpface.interp,2) == size(uface_bar,3) )
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]
      for j = 1:sbpface.stencilsize
        for field = 1:size(uvol_bar,1)
          # uface[field,1,i,findex] += sbpface.interp[j,i]*
          # uvol[field,sbpface.perm[j,face.faceL],face.elementL]
          uvol_bar[field,sbpface.perm[j,face.faceL],face.elementL] +=
            sbpface.interp[j,i]*uface_bar[field,1,i,findex]
          # uface[field,2,i,findex] += sbpface.interp[j,iR]*
          # uvol[field,sbpface.perm[j,face.faceR],face.elementR]
          uvol_bar[field,sbpface.perm[j,face.faceR],face.elementR] +=
            sbpface.interp[j,iR]*uface_bar[field,2,i,findex]
        end
      end
    end
  end
end

function interiorfaceinterpolate_rev!(sbpface::SparseFace{Tsbp},
                                                 ifaces::AbstractArray{Interface},
                                                 uvol_bar::AbstractArray{Tsol,2},
                                                 uface_bar::AbstractArray{Tsol,3}) where {Tsbp,Tsol}
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]     
      # uface[1,i,findex] = uvol[sbpface.perm[i,face.faceL],face.elementL]
      uvol_bar[sbpface.perm[i,face.faceL],face.elementL] += uface_bar[1,i,findex]
      # uface[2,i,findex] = uvol[sbpface.perm[iR,face.faceR],face.elementR]
      uvol_bar[sbpface.perm[iR,face.faceR],face.elementR] += uface_bar[2,i,findex]
    end
  end
end

function interiorfaceinterpolate_rev!(sbpface::SparseFace{Tsbp},
                                                 ifaces::AbstractArray{Interface},
                                                 uvol_bar::AbstractArray{Tsol,3},
                                                 uface_bar::AbstractArray{Tsol,4}) where {Tsbp,Tsol}
  @assert( size(uvol_bar,1) == size(uface_bar,1) )
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]
      for field=1:size(uvol_bar, 1)
        # uface[field,1,i,findex] =
        #   uvol[field,sbpface.perm[i,face.faceL],face.elementL]
        uvol_bar[field,sbpface.perm[i,face.faceL],face.elementL] +=
          uface_bar[field,1,i,findex]
        # uface[field,2,i,findex] =
        #   uvol[field,sbpface.perm[iR,face.faceR],face.elementR]
        uvol_bar[field,sbpface.perm[iR,face.faceR],face.elementR] +=
          uface_bar[field,2,i,findex]
      end
    end
  end
end

"""
### SummationByParts.interiorFaceInterpolate_rev!

This is the reverse differentiated version of interiorFaceInterpolate!.  See
faceinterpolate.jl for further details of the primal method.  This function is
differentiated with respect to the primal version's `uL` and `uR` variables.

**Inputs**

* `sbpface`: an SBP face operator type
* `iface`: Interface type that specifies face indices and relative orientation
* `ufaceL_bar`: vector applied to the left of the interpolation operator (R)
* `ufaceR_bar`: vector applied to the left of the interpolation operator (R)

**In/Outs**

* `uL_bar`: the result of the vector-matrix product between R and `ufaceL_bar`
* `uR_bar`: the result of the vector-matrix product between R and `ufaceR_bar`

"""
function interiorFaceInterpolate_rev!(sbpface::DenseFace{Tsbp},
                                             iface::Interface,
                                             uL_bar::AbstractArray{Tsol,1},
                                             uR_bar::AbstractArray{Tsol,1},
                                             ufaceL_bar::AbstractArray{Tsol,1},
                                             ufaceR_bar::AbstractArray{Tsol,1}) where {Tsbp,Tsol}
  @assert( size(sbpface.interp,1) <= size(uL_bar,1) )
  @assert( size(sbpface.interp,1) <= size(uR_bar,1) )
  @assert( size(sbpface.interp,2) == size(ufaceL_bar,1) == size(ufaceR_bar,1) )
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    for j = 1:sbpface.stencilsize
      # ufaceL[i] += sbpface.interp[j,i]*uL[sbpface.perm[j,iface.faceL]]
      uL_bar[sbpface.perm[j,iface.faceL]] += sbpface.interp[j,i]*ufaceL_bar[i]
      # ufaceR[i] += sbpface.interp[j,iR]*uR[sbpface.perm[j,iface.faceR]]
      uR_bar[sbpface.perm[j,iface.faceR]] += sbpface.interp[j,iR]*ufaceR_bar[i]
    end
  end
end

function interiorFaceInterpolate_rev!(sbpface::DenseFace{Tsbp},
                                             iface::Interface,
                                             uL_bar::AbstractArray{Tsol,2},
                                             uR_bar::AbstractArray{Tsol,2},
                                             ufaceL_bar::AbstractArray{Tsol,2},
                                             ufaceR_bar::AbstractArray{Tsol,2}) where {Tsbp,Tsol}
  @assert( size(uL_bar,1) == size(ufaceL_bar,1) == size(uR_bar,1) == size(ufaceR_bar,1) )
  @assert( size(sbpface.interp,1) <= size(uL_bar,2) )
  @assert( size(sbpface.interp,1) <= size(uR_bar,2) )
  @assert( size(sbpface.interp,2) == size(ufaceL_bar,2) == size(ufaceR_bar,2) )
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    for j = 1:sbpface.stencilsize
      for field = 1:size(uL_bar,1)
        # ufaceL[field,i] += sbpface.interp[j,i]*
        # uL[field,sbpface.perm[j,iface.faceL]]
        uL_bar[field,sbpface.perm[j,iface.faceL]] +=
          sbpface.interp[j,i]*ufaceL_bar[field,i]
        # ufaceR_bar[field,i] += sbpface.interp[j,iR]*
        # uR_bar[field,sbpface.perm[j,iface.faceR]]
        uR_bar[field,sbpface.perm[j,iface.faceR]] +=
          sbpface.interp[j,iR]*ufaceR_bar[field,i]        
      end
    end
  end
end

function interiorFaceInterpolate_rev!(sbpface::SparseFace{Tsbp},
                                                 iface::Interface,
                                                 uL_bar::AbstractArray{Tsol,1},
                                                 uR_bar::AbstractArray{Tsol,1},
                                                 ufaceL_bar::AbstractArray{Tsol,1},
                                                 ufaceR_bar::AbstractArray{Tsol,1}) where {Tsbp,Tsol}
  @assert( size(ufaceL_bar,1) == size(ufaceR_bar,1) )
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    # ufaceL[i] = uL[sbpface.perm[i,iface.faceL]]
    uL_bar[sbpface.perm[i,iface.faceL]] += ufaceL_bar[i]
    # ufaceR[i] = uR[sbpface.perm[iR,iface.faceR]]
    uR_bar[sbpface.perm[iR,iface.faceR]] += ufaceR_bar[i]
  end
end

function interiorFaceInterpolate_rev!(sbpface::SparseFace{Tsbp},
                                                 iface::Interface,
                                                 uL_bar::AbstractArray{Tsol,2},
                                                 uR_bar::AbstractArray{Tsol,2},
                                                 ufaceL_bar::AbstractArray{Tsol,2},
                                                 ufaceR_bar::AbstractArray{Tsol,2}) where {Tsbp,Tsol}
  @assert( size(uL_bar,1) == size(ufaceL_bar,1) == size(uR_bar,1) ==
           size(ufaceR_bar,1) )
  @assert( size(ufaceL_bar,2) == size(ufaceR_bar,2) )
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    @simd for field=1:size(uL_bar,1)
      # ufaceL[field,i] = uL[field,sbpface.perm[i,iface.faceL]]
      uL_bar[field,sbpface.perm[i,iface.faceL]] += ufaceL_bar[field,i]
      # ufaceR[field,i] = uR[field,sbpface.perm[iR,iface.faceR]]
      uR_bar[field,sbpface.perm[iR,iface.faceR]] += ufaceR_bar[field,i]
    end
  end
end
