# This file gathers together functions related to the calculation of scaled face
# normals for general (curvilinear) faces

@doc """
### SummationByParts.calcFaceNormals!

Uses a given set of Lagrangian face nodes to determine an analytical
(polynomial) mapping, and then uses this mapping to determine the scaled
face-normal vector.

**Inputs**

* `sbpface`: an SBP face operator type
* `mapdegree`: the polynomial degree of the mapping
* `xref`: Lagrangian nodes in reference space; [coord, Lag node]
* `xlag`: Lagrangian nodes in physical space; [coord, Lag node, face]

**In/Outs**

* `xsbp`: SBP-face nodes in physical space; [coord, sbp node, face]
* `nrm`: scaled face-normal at the sbpface nodes; [component, sbp node, face]

"""->
function calcFaceNormals!{Tsbp,Tmsh}(sbpface::TriFace{Tsbp},
                                     mapdegree::Int,
                                     xref::AbstractArray{Tmsh,2},
                                     xlag::AbstractArray{Tmsh,3},
                                     xsbp::AbstractArray{Tmsh,3},
                                     nrm::AbstractArray{Tmsh,3})
  @assert( size(xlag,1) == size(xsbp,1) == size(nrm,1) == 2 )
  @assert( size(xref,1) == 1 )
  @assert( size(xsbp,2) == size(nrm,2) )
  numdof = (mapdegree+1)
  @assert( size(xlag,2) == size(xref,2) == numdof )
  @assert( size(xlag,3) == size(xsbp,3) == size(nrm,3) )
  # find the inverse of the Vandermonde matrix
  V = zeros(Tmsh, (numdof,numdof) )
  for i = 0:mapdegree
    V[:,i+1] = OrthoPoly.jacobipoly(vec(xref[1,:]), 0.0, 0.0, i)
  end
  Vinv = inv(V)
  # get the SBP nodes in reference space in order to find the orthogonal
  # polynomials and their derivatives at these nodes
  x = SymCubatures.calcnodes(sbpface.cub, sbpface.vtx)
  P = zeros(Tmsh, (sbpface.numnodes, numdof))
  dPdξ = zeros(Tmsh, (sbpface.numnodes, numdof))
  for i = 0:mapdegree
    P[:,i+1] = OrthoPoly.jacobipoly(vec(x[1,:]), 0.0, 0.0, i)
    dPdξ[:,i+1] = OrthoPoly.diffjacobipoly(vec(x[1,:]), 0.0, 0.0, i)
  end
  # loop over each face...
  fill!(nrm, zero(Tmsh))
  fill!(xsbp, zero(Tmsh))
  coeff = zeros(Tmsh, (numdof,2))
  for f = 1:size(xlag,3)
    # find the coefficents of the polynomial mapping using xlag and Vinv
    for i = 1:numdof
      coeff[i,1] = zero(Tmsh)
      coeff[i,2] = zero(Tmsh)
      for j = 1:numdof
        coeff[i,1] += Vinv[i,j]*xlag[1,j,f]
        coeff[i,2] += Vinv[i,j]*xlag[2,j,f]
      end
    end
    # compute the mapped SBP nodes and the analytical normal at these nodes
    for i = 1:numdof
      for nd = 1:sbpface.numnodes
        xsbp[1,nd,f] += coeff[i,1]*P[nd,i]
        xsbp[2,nd,f] += coeff[i,2]*P[nd,i]
        nrm[1,nd,f] += coeff[i,2]*dPdξ[nd,i]
        nrm[2,nd,f] -= coeff[i,1]*dPdξ[nd,i]
      end
    end
  end
end

function calcFaceNormals!{Tsbp,Tmsh}(sbpface::TetFace{Tsbp},
                                     mapdegree::Int,
                                     xref::AbstractArray{Tmsh,2},
                                     xlag::AbstractArray{Tmsh,3},
                                     xsbp::AbstractArray{Tmsh,3},
                                     nrm::AbstractArray{Tmsh,3})
  @assert( size(xlag,1) == size(xsbp,1) == size(nrm,1) == 3 )
  @assert( size(xref,1) == 2 )
  @assert( size(xsbp,2) == size(nrm,2) )
  numdof = binomial(mapdegree+2,2)
  @assert( size(xlag,2) == size(xref,2) == numdof )
  @assert( size(xlag,3) == size(xsbp,3) == size(nrm,3) )
  # find the inverse of the Vandermonde matrix
  V = zeros(Tmsh, (numdof,numdof) )
  ptr = 1
  for r = 0:mapdegree
    for j = 0:r
      i = r-j
      V[:,ptr] = OrthoPoly.proriolpoly(vec(xref[1,:]), vec(xref[2,:]), i, j)
      ptr += 1
    end
  end
  Vinv = inv(V)
  # get the SBP nodes in reference space in order to find the orthogonal
  # polynomials and their derivatives at these nodes
  x = SymCubatures.calcnodes(sbpface.cub, sbpface.vtx)
  P = zeros(Tmsh, (sbpface.numnodes, numdof))
  dPdξ = zeros(Tmsh, (sbpface.numnodes, numdof))
  dPdη = zeros(Tmsh, (sbpface.numnodes, numdof))
  ptr = 1
  for r = 0:mapdegree
    for j = 0:r
      i = r-j
      P[:,ptr] = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      dPdξ[:,ptr], dPdη[:,ptr] =
        OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      ptr += 1
    end
  end
  # loop over each face...
  fill!(nrm, zero(Tmsh))
  fill!(xsbp, zero(Tmsh))
  coeff = zeros(Tmsh, (numdof,3))
  dxdξ = zeros(Tmsh, (3,2,sbpface.numnodes))
  for f = 1:size(xlag,3)
    # find the coefficents of the polynomial mapping using xlag and Vinv
    for di = 1:3
      for i = 1:numdof
        coeff[i,di] = zero(Tmsh)
        for j = 1:numdof
          coeff[i,di] += Vinv[i,j]*xlag[di,j,f]
        end
      end
    end
    # compute the mapped SBP nodes and the tangent vectors at these nodes
    fill!(dxdξ, zero(Tmsh))
    for i = 1:numdof
      for di = 1:3
        for nd = 1:sbpface.numnodes
          xsbp[di,nd,f] += coeff[i,di]*P[nd,i]
          dxdξ[di,1,nd] += coeff[i,di]*dPdξ[nd,i]
          dxdξ[di,2,nd] += coeff[i,di]*dPdη[nd,i]
        end
      end
    end
    # compute the face-normal vector using the tangent vectors
    for di = 1:3
      it1 = mod(di,3)+1
      it2 = mod(di+1,3)+1
      for i = 1:sbpface.numnodes
        nrm[di,i,f] = dxdξ[it1,1,i]*dxdξ[it2,2,i] - dxdξ[it2,1,i]*dxdξ[it1,2,i]
      end
    end
  end
end

@doc """
### SummationByParts.facenormal!

This is the single-face variant of calcFaceNormals!.  Uses a given set of
Lagrangian face nodes to determine an analytical (polynomial) mapping, and then
uses this mapping to determine the scaled face-normal vector.

**Inputs**

* `sbpface`: an SBP face operator type
* `mapdegree`: the polynomial degree of the mapping
* `xref`: Lagrangian nodes in reference space; [coord, Lagrangian node]
* `xlag`: Lagrangian nodes in physical space; [coord, Lagrangian node]

**In/Outs**

* `xsbp`: location of the SBP-face nodes in physical space; [coord, sbp node]
* `nrm`: scaled face-normal at the sbpface nodes

"""->
function facenormal!{Tsbp,Tmsh}(sbpface::TriFace{Tsbp},
                                mapdegree::Int,
                                xref::AbstractArray{Tmsh,2},
                                xlag::AbstractArray{Tmsh,2},
                                xsbp::AbstractArray{Tmsh,2},
                                nrm::AbstractArray{Tmsh,2})
  @assert( size(xlag,1) == size(xsbp,1) == size(nrm,1) == 2 )
  @assert( size(xref,1) == 1 )
  @assert( size(xsbp,2) == size(nrm,2) )
  numdof = (mapdegree+1)
  @assert( size(xlag,2) == size(xref,2) == numdof )
  # Step 1: find the polynomial mapping using xlag
  V = zeros(Tmsh, (numdof,numdof) )
  for i = 0:mapdegree
    V[:,i+1] = OrthoPoly.jacobipoly(vec(xref[1,:]), 0.0, 0.0, i)
  end
  coeff = zeros(Tmsh, (numdof,2))
  coeff = V\(xlag.')
  # Step 2: compute the mapped SBP nodes and the analytical normal at sbp nodes
  x = SymCubatures.calcnodes(sbpface.cub, sbpface.vtx) # <-- SBP nodes, ref. spc
  fill!(nrm, zero(Tmsh))
  fill!(xsbp, zero(Tmsh))
  for i = 0:mapdegree
    P = OrthoPoly.jacobipoly(vec(x[1,:]), 0.0, 0.0, i)
    dPdξ = OrthoPoly.diffjacobipoly(vec(x[1,:]), 0.0, 0.0, i)
    for nd = 1:sbpface.numnodes
      xsbp[1,nd] += coeff[i+1,1]*P[nd]
      xsbp[2,nd] += coeff[i+1,2]*P[nd]
      nrm[1,nd] += coeff[i+1,2]*dPdξ[nd]
      nrm[2,nd] -= coeff[i+1,1]*dPdξ[nd]
    end
  end
end

function facenormal!{Tsbp,Tmsh}(sbpface::TetFace{Tsbp},
                                mapdegree::Int,
                                xref::AbstractArray{Tmsh,2},
                                xlag::AbstractArray{Tmsh,2},
                                xsbp::AbstractArray{Tmsh,2},
                                nrm::AbstractArray{Tmsh,2})
  @assert( size(xlag,1) == size(xsbp,1) == size(nrm,1) == 3 )
  @assert( size(xref,1) == 2 )
  @assert( size(xsbp,2) == size(nrm,2) )
  numdof = binomial(mapdegree+2,2)
  @assert( size(xlag,2) == size(xref,2) == numdof )
  # Step 1: find the polynomial mapping using xlag
  V = zeros(Tmsh, (numdof,numdof) )
  ptr = 1
  for r = 0:mapdegree
    for j = 0:r
      i = r-j
      V[:,ptr] = OrthoPoly.proriolpoly(vec(xref[1,:]), vec(xref[2,:]), i, j)
      ptr += 1
    end
  end
  coeff = zeros(Tmsh, (numdof,3))
  coeff = V\(xlag.')
  # Step 2: compute the mapped SBP nodes and the tangent vectors at sbp nodes
  x = SymCubatures.calcnodes(sbpface.cub, sbpface.vtx) # <-- SBP nodes, ref. spc
  fill!(xsbp, zero(Tmsh))
  dxdξ = zeros(Tmsh, (3,2,sbpface.numnodes))
  ptr = 1
  for r = 0:mapdegree
    for j = 0:r
      i = r-j
      P = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      dPdξ, dPdη = OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      for di = 1:3
        for nd = 1:sbpface.numnodes
          xsbp[di,nd] += coeff[ptr,di]*P[nd]
          dxdξ[di,1,nd] += coeff[ptr,di]*dPdξ[nd]
          dxdξ[di,2,nd] += coeff[ptr,di]*dPdη[nd]
        end
      end
      ptr += 1
    end
  end
  # Step 3: compute the face-normal vector using the tangent vectors
  for di = 1:3
    it1 = mod(di,3)+1
    it2 = mod(di+1,3)+1
    for i = 1:sbpface.numnodes
      nrm[di,i] = dxdξ[it1,1,i]*dxdξ[it2,2,i] - dxdξ[it2,1,i]*dxdξ[it1,2,i]
    end
  end
end
