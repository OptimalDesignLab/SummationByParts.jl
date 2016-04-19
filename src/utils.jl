# This file gathers together a hodge-podge of functions that are not easily
# categorized

# make sview point to either safe or unsafe views
global const use_safe_views = true
if use_safe_views
  global const sview = view
else
  global const sview = unsafe_view
end

@doc """
### SummationByParts.buildinterpolation

Builds a matrix operator that can reconstruct a field located at the sbp nodes
to an auxlliary set of nodes.

**Inputs**

* `sbp`: an SBP operator
* `xinterp`: points to interpolate to in ref coords, size = [ndim,numpoints]

**Returns**

* `R`: the interpolation operator, size = [numpoints, sbp.numnodes]

"""->
function buildinterpolation{T}(sbp::TriSBP{T}, xinterp::AbstractArray{T,2})
  # evaluate the basis at the SBP nodes and the interpolation points
  d = sbp.degree
  N = convert(Int, (d+1)*(d+2)/2 )
  Psbp = zeros(T, (sbp.numnodes,N) )  
  Pinterp = zeros(T, (size(xinterp,2),N) ) 
  xsbp = SymCubatures.calcnodes(sbp.cub, sbp.vtx)
  ptr = 1
  for r = 0:d
    for j = 0:r
      i = r-j
      Psbp[:,ptr] = OrthoPoly.proriolpoly(vec(xsbp[1,:]), vec(xsbp[2,:]), i, j)
      Pinterp[:,ptr] = OrthoPoly.proriolpoly(vec(xinterp[1,:]), 
                                             vec(xinterp[2,:]), i, j)
      ptr += 1
    end
  end
  R = (pinv(Psbp.')*Pinterp.').'
  return R
end

@doc """
### SummationByParts.permuteinterface!

For a given array of values on a faces, permutes the node values (in place) to 
be in the orientation specified by the orient field of the corresponding 
Interface.  Methods are available for scalar and vector fields.


**Inputs**

* `sbp`: an SBPFace operator
* `ifaces`: an array of Interfaces

**In/Outs**

* `uface`: the array of face values, must have dimensions
           [n x sbpface.numnodes x length(ifaces] if 3D array (where n
           is arbitrary) or [sbpface.numnodes x length(ifaces] if 2D.  
           The permutation is applied to the second array dimension for the 
           3D case or the first dimension in the 2D case.
"""->
function permuteinterface!{Tsbp, Tsol}(sbpface::AbstractFace{Tsbp}, 
                                       ifaces::AbstractArray{Interface}, 
                                       uface::AbstractArray{Tsol, 3})
  @assert length(ifaces) == size(uface, 3)
  @assert sbpface.numnodes == size(uface, 2)

  dofpernode, numfacenodes, numfaces = size(uface)
  # temporary array needed during permutation
  workarr = Array(Tsol, dofpernode, numfacenodes)

  for iface =1:length(ifaces)
    orient = ifaces[iface].orient
    permvec = sview(sbpface.nbrperm, :, orient)
    facedata = sview(uface, :, :, iface)
    permuteface!(permvec, workarr, facedata)
  end

  return nothing
end

function permuteinterface!{Tsbp, Tsol}(sbpface::AbstractFace{Tsbp}, 
                                       ifaces::AbstractArray{Interface}, 
                                       uface::AbstractArray{Tsol, 2})
  @assert length(ifaces) == size(uface, 2)
  @assert sbpface.numnodes == size(uface, 1)

  numfacenodes, numfaces = size(uface)
  # temporary array needed during permutation
  workarr = Array(Tsol, numfacenodes)

  for iface =1:length(ifaces)
    orient = ifaces[iface].orient
    permvec = sview(sbpface.nbrperm, :, orient)
    facedata = sview(uface, :, iface)
    permuteface!(permvec, workarr, facedata)
  end

  return nothing
end



@doc """
### SummationByParts.permuteface!

This function applys a permutation to the data on a particular face.

**Inputs**

* `permvec`: vector specifying the permutation to apply
* `workarr`: a temporary array, same size as face_data, that is overwritten
             during the computation
**In/Outs**

* `face_data`: an N-D array containing the data to be pemuted, where the
              permutation is applied to the second dimension of the array
              for N=2 and the first dimension if N=1.  N > 2 is not
              currently supported
"""->
function permuteface!{Ti <: Integer, Tsol}(permvec::AbstractArray{Ti, 1},
                                           workarr::AbstractArray{Tsol, 2},
                                           facedata::AbstractArray{Tsol, 2})

  # copy to temporary array, applying permutation
  for i=1:size(facedata, 2)  # loop over nodes on the face
    idx = permvec[i]
    for j=1:size(facedata, 1)  # all dofs on the node
      workarr[j, idx] = facedata[j, i]
    end
  end

  # copy back, using linear indexing
  for i=1:length(facedata)
    facedata[i] = workarr[i]
  end

  return nothing
end

function permuteface!{Ti <: Integer, Tsol}(permvec::AbstractArray{Ti, 1},
                                           workarr::AbstractArray{Tsol, 1},
                                           facedata::AbstractArray{Tsol, 1})

  # copy to temporary array, applying permutation
  for i=1:size(facedata, 1)  # loop over nodes on the face
    idx = permvec[i]
    workarr[idx] = facedata[i]
  end

  # copy back, using linear indexing
  for i=1:length(facedata)
    facedata[i] = workarr[i]
  end

  return nothing
end
