# This file gathers together outer constructors for the SBP operators

"""
### SBP.getLineSegSBPLobbato

Returns Gauss-Lobbato type elements, that have nodes on the element boundary

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators

**Returns**

* `sbp`: a Gauss-Lobbato operator of the appropriate degree

"""
function getLineSegSBPLobbato(;degree::Int=1, Tsbp::Type=Float64)
  cub, vtx = quadrature(2*degree-1, Tsbp, internal=false)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 2))
  w, Q = SummationByParts.buildoperators(cub, vtx, degree)
  return LineSegSBP{Tsbp}(degree, cub, vtx, w, Q)
end

"""
### SBP.getLineSegSBPLegendre

Returns Gauss-Legendre type elements, that do not have nodes on the element
boundary

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators

**Returns**

* `sbp`: a Legendre-Gauss operator of the appropriate degree

"""
function getLineSegSBPLegendre(;degree::Int=1, Tsbp::Type=Float64)
  cub, vtx = quadrature(2*degree, Tsbp, internal=true)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 2))
  w, Q = SummationByParts.buildoperators(cub, vtx, degree)
  return LineSegSBP{Tsbp}(degree, cub, vtx, w, Q)
end

"""
### SBP.getTriSBPGamma

Returns SBP-Gamma type elements, that have nodes on the element boundary

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators

**Returns**

* `sbp`: an SBP-Gamma operator of the appropriate degree

"""
function getTriSBPGamma(;degree::Int=1, Tsbp::Type=Float64)
  cub, vtx = getTriCubatureGamma(2*degree-1, Tsbp)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 2))
  w, Q = SummationByParts.buildoperators(cub, vtx, degree)
  return TriSBP{Tsbp}(degree, cub, vtx, w, Q)
end

function getTriSBPOmega0(;degree::Int=1, Tsbp::Type=Float64)

  if degree <= 2
    cubdegree = 2*degree
  else
    cubdegree = 2*degree - 1
  end

  cub, vtx = getTriCubatureOmega(cubdegree, Tsbp)
  w, Q = SummationByParts.buildoperators(cub, vtx, degree)
  return TriSBP{Tsbp}(degree, cub, vtx, w, Q)
end


"""
### SBP.getTriSBPOmega

Returns SBP-Omega type elements, that have no nodes on the element boundary

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators

**Returns**

* `sbp`: an SBP-Omega operator of the appropriate degree

"""
function getTriSBPOmega(;degree::Int=1, Tsbp::Type=Float64)
  cub, vtx = getTriCubatureOmega(2*degree, Tsbp)
  w, Q = SummationByParts.buildoperators(cub, vtx, degree)
  return TriSBP{Tsbp}(degree, cub, vtx, w, Q)
end

"""
### SBP.getTriSBPOmega2

Like getTRISBPomega, but ensures the operator has a degree 2p cubature
rule

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators

**Returns**

* `sbp`: an SBP-Omega type operator of the appropriate degree

"""
function getTriSBPOmega2(;degree::Int=1, Tsbp::Type=Float64)

  cub, vtx = getTriCubatureOmega(2*degree, Tsbp)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 2))
#  w, Q = SummationByParts.buildoperators(cub, vtx, degree)
  w, Q = buildMinConditionOperators(cub, vtx, degree, vertices=false)

  return TriSBP{Tsbp}(degree, cub, vtx, w, Q)
end


"""
### SBP.getTriSBPDiagE

Returns SBP-DiagE type elements, whose boundary nodes are positioned at cubature
points

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators
* `vertices`: include vertices in the operator
* `mincond`: use the min-condition number operators

**Returns**

* `sbp`: an SBP-DiagE operator of the appropriate degree

"""
function getTriSBPDiagE(;degree::Int=1, Tsbp::Type=Float64,
                        vertices::Bool=true, mincond::Bool=true)
  if degree == 0
    cub, vtx = getTriCubatureGamma(1, Tsbp)
    Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 2))
    w, Q = SummationByParts.buildoperators(cub, vtx, 1)
    return TriSBP{Tsbp}(degree, cub, vtx, w, Q)
  end
  cub, vtx = getTriCubatureDiagE(2*degree, Tsbp, vertices=vertices)
 
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 2))
  if mincond
    # use the files written to disk
    w = SymCubatures.calcweights(cub)
    vstr = vertices ? "vert" : "novert"
    if degree >= 1 && degree <= 4
      Q = reshape( readdlm(joinpath(dirname(@__FILE__), "tri_diage_p$(degree)_$vstr.dat")), cub.numnodes, cub.numnodes, 2)

    else  # perform optimization
      w, Q = SummationByParts.buildMinConditionOperators(cub, vtx, degree,
                                                         vertices=vertices)
    end
  else
    w, Q = SummationByParts.buildMinFrobeniusOperators(cub, vtx, degree,
                                                       vertices=vertices)
  end

  return TriSBP{Tsbp}(degree, cub, vtx, w, Q)
end

"""
### SBP.getTetSBPGamma

Returns SBP-Gamma type elements, that have nodes on the element boundary

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators

**Returns**

* `sbp`: an SBP-Gamma operator of the appropriate degree

"""
function getTetSBPGamma(;degree::Int=1, Tsbp::Type=Float64)
  
  cub, vtx = getTetCubatureGamma(2*degree-1, Tsbp)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 3))
  w, Q = SummationByParts.buildoperators(cub, vtx, degree)
  return TetSBP{Tsbp}(degree, cub, vtx, w, Q)
end

"""
### SBP.getTetSBPOmega

Returns SBP-Omega type elements, that have no nodes on the element boundary

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators

**Returns**

* `sbp`: an SBP-Omega operator of the appropriate degree

"""
function getTetSBPOmega(;degree::Int=1, Tsbp::Type=Float64)
  cub, vtx = getTetCubatureOmega(2*degree-1, Tsbp)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 3))
  w, Q = SummationByParts.buildoperators(cub, vtx, degree)
  return TetSBP{Tsbp}(degree, cub, vtx, w, Q)  
end

"""
### SBP.getTetSBPDiagE

Returns SBP-DiagE type elements, whose boundary nodes are positioned at cubature
points

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators

**Returns**

* `sbp`: an SBP-DiagE operator of the appropriate degree

"""
function getTetSBPDiagE(;degree::Int=1, Tsbp::Type=Float64,
                        edges::Bool=false, vertices::Bool=false)
  @assert( degree >= 1 && degree <= 4 )
  # the operators saved to disk only support the default values for hte
  # node locations -> error if the user askes for something else
  @assert !edges  
  @assert !vertices
  cub, vtx = getTetCubatureDiagE(2*degree, Tsbp, vertices=vertices)
  w = zeros(Tsbp, (cub.numnodes))
  w = SymCubatures.calcweights(cub)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 3))
  if degree >= 1 && degree <= 4
    Q = reshape(readdlm(dirname(@__FILE__)"/tet_diage_p$(degree).dat", '\t', Tsbp),
                size(Q))
  else
    error("degree $degree Tet SBP diagE not supported")
  end
  return TetSBP{Tsbp}(degree, cub, vtx, w, Q)
end

"""
### SBP.getLineSegFace

Returns a trival face for line-segment elements.

**Inputs**

* `degree`: face integration is exact for polys of degree 2*`degree`
* `volcub`: cubature rule for the associated "volume"
* `vtx`: vertices of the line-segment

**Returns**

* `sbpface`: an SBP face type for line-segment elements

"""
function getLineSegFace{T}(degree::Int, volcub::LineSymCub{T}, vtx::Array{T,2})
  facecub, facevtx = pointCubature()
  R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx,
                                                     degree)
  D, Dperm = SummationByParts.buildfacederivatives(facecub, volcub, vtx, degree)
  wface = SymCubatures.calcweights(facecub)
  stencilsize = size(R,2)
  dstencilsize = size(D,1)
  return LineSegFace{T}(degree, facecub, facevtx, R.', perm, D, Dperm)
end

"""
### SBP.TriFace

Outer constructor for backward compatibility

**Inputs**

* `degree`: face integration is exact for polys of degree 2*`degree`
* `volcub`: cubature associated with the volume domain this face is part of
* `vtx`: vertices of the volume domain
* `vertices`: if true, element vertices are included in the nodes

**Returns**

* `sbpface`: an SBP face type for triangle elements

"""
function (::Type{TriFace{T}}){T}(degree::Int, volcub::TriSymCub{T},
                                 vtx::Array{T,2}; vertices::Bool=false)
  @assert( degree >= 1 && degree <= 5 )
  local R::Array{T,2}
  local perm::Array{Int,2}
  if vertices
    facecub, facevtx = quadrature(2*degree, T, internal=false)
    R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx,
                                                       degree+1)
  else
    facecub, facevtx = quadrature(2*degree, T, internal=true)
    R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx,
                                                       degree)    
  end
  D, Dperm = SummationByParts.buildfacederivatives(facecub, volcub, vtx,
                                                   degree)
  nbrperm = SymCubatures.getneighbourpermutation(facecub)
  wface = SymCubatures.calcweights(facecub)
  stencilsize = size(R,2)
  dstencilsize = size(D,1)
  TriFace{T}(degree, facecub, facevtx, R.', perm, D, Dperm)
end

"""
### SBP.getTriFaceForDiagE

Returns a quadrature that can be used to construct SBP operators with diagonal
boundary operators, E.

**Inputs**

* `degree`: face integration is exact for polys of degree 2*`degree`
* `volcub`: cubature rule for the associated volume
* `vtx`: vertices of the triangle
* `vertices`: (optional) if true, the face cubature includes the vertices

**Returns**

* `sbpface`: an SBP face type for triangle elements

"""
function getTriFaceForDiagE{T}(degree::Int, volcub::TriSymCub{T},
                               vtx::Array{T,2}; vertices::Bool=true)
  #@assert( degree >= 1 && degree <= 4 )
  if vertices
    facecub, facevtx = quadrature(2*degree, T, internal=false)
    R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx,
                                                       degree+1)
  else
    facecub, facevtx = quadrature(2*degree, T, internal=true)
    R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx,
                                                       degree)
  end
  perm_red = zeros(Int, (facecub.numnodes,3))
  for i = 1:facecub.numnodes
    node = indmax(R[i,:])
    perm_red[i,:] = perm[node,:]
  end  
  D, Dperm = SummationByParts.buildfacederivatives(facecub, volcub, vtx,
                                                   degree)
  nbrperm = SymCubatures.getneighbourpermutation(facecub)
  wface = SymCubatures.calcweights(facecub)
  return TriSparseFace{T}(degree, facecub, facevtx, perm_red, D, Dperm)
end

"""
### SBP.TetFace

Outer constructor for backward compatibility

**Inputs**

* `degree`: face integration is exact for polys of degree 2*`degree`
* `volcub`: cubature associated with the volume domain this face is part of
* `vtx`: vertices of the volume domain

**Returns**

* `sbpface`: an SBP face type for tetrahedral elements

"""
function (::Type{TetFace{T}}){T}(degree::Int, volcub::TetSymCub{T},
                                 vtx::Array{T,2})
  @assert( degree >= 1 && degree <= 4 )
  facecub, facevtx = getTriCubatureOmega(2*degree, T)
  normal = T[0 0 -1; 0 -1 0; 1 1 1; -1 0 0].'
  R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx,
                                                     degree)
  #D, Dperm = SummationByParts.buildfacederivatives(facecub, volcub, vtx,
  #                                                 degree)
  nbrperm = SymCubatures.getneighbourpermutation(facecub)
  wface = SymCubatures.calcweights(facecub)
  stencilsize = size(R,2)
  #dstencilsize = size(D,1)
  #new(degree, facecub.numnodes, stencilsize, dstencilsize, facecub, wface,
  #    normal, R.', perm, D, Dperm, nbrperm)
  TetFace{T}(degree, facecub, facevtx, R.', perm)
end

"""
### SBP.getTetFaceForDiagE

Returns a quadrature that can be used to construct SBP operators with diagonal
boundary operators, E.

**Inputs**

* `degree`: face integration is exact for polys of degree 2*`degree`+1
* `volcub`: cubature rule for the associated volume
* `vtx`: vertices of the triangle

**Returns**

* `sbpface`: an SBP face type for triangle elements

"""
function getTetFaceForDiagE{T}(degree::Int, volcub::TetSymCub{T},
                               vtx::Array{T,2})
  @assert( degree >= 1 && degree <= 4 )
  
  facecub, facevtx = getTriCubatureOmega(2*degree, T)
  
  #facecub = SymCubatures.TriSymCub{T}(vertices=true, numS21=1)
  #SymCubatures.setparams!(facecub, T[0.9055050463303657])
  #SymCubatures.setweights!(facecub, T[(1+sqrt(21))/60;(39+sqrt(21))/60])
  #facevtx = T[-1 -1; 1 -1; -1 1]
  R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx,
                                                     degree+1)
  #println("TEMP change in getTriFaceForDiagE!!!!")
  #facecub, facevtx = quadrature(2*degree, T, internal=true)
  #R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx,
  #                                                   degree)

  perm_red = zeros(Int, (facecub.numnodes,4))
  for i = 1:facecub.numnodes
    node = indmax(R[i,:])
    perm_red[i,:] = perm[node,:]
  end    
  nbrperm = SymCubatures.getneighbourpermutation(facecub)
  wface = SymCubatures.calcweights(facecub)
  return TetSparseFace{T}(degree, facecub, facevtx, perm_red)
end
