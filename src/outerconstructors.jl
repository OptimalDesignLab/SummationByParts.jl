# This file gathers together outer constructors for the SBP operators

@doc """
### SBP.getTriSBPGamma

Returns SBP-Gamma type elements, that have nodes on the element boundary

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators

**Returns**

* `sbp`: an SBP-Gamma operator of the appropriate degree

"""->
function getTriSBPGamma(;degree::Int=1, Tsbp::Type=Float64)
  cub, vtx = getTriCubatureGamma(2*degree-1, Tsbp)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 2))
  w, Q = SummationByParts.buildoperators(cub, vtx, degree)
  return TriSBP{Tsbp}(degree, cub, vtx, w, Q)
end

@doc """
### SBP.getTriSBPOmega

Returns SBP-Omega type elements, that have no nodes on the element boundary

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators

**Returns**

* `sbp`: an SBP-Omega operator of the appropriate degree

"""->
function getTriSBPOmega(;degree::Int=1, Tsbp::Type=Float64)
  cub, vtx = getTriCubatureOmega(2*degree, Tsbp)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 2))
  w, Q = SummationByParts.buildoperators(cub, vtx, degree)
  return TriSBP{Tsbp}(degree, cub, vtx, w, Q)
end

@doc """
### SBP.getTriSBPDiagE

Returns SBP-DiagE type elements, whose boundary nodes are positioned at cubature
points

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators

**Returns**

* `sbp`: an SBP-DiagE operator of the appropriate degree

"""->
function getTriSBPDiagE(;degree::Int=1, Tsbp::Type=Float64,
                        vertices::Bool=true)
  cub, vtx = getTriCubatureDiagE(2*degree, Tsbp, vertices=vertices)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 2))
  w, Q = SummationByParts.buildMinConditionOperators(cub, vtx, degree,
                                                     vertices=vertices)
  return TriSBP{Tsbp}(degree, cub, vtx, w, Q)
end

@doc """
### SBP.getTetSBPGamma

Returns SBP-Gamma type elements, that have nodes on the element boundary

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators

**Returns**

* `sbp`: an SBP-Gamma operator of the appropriate degree

"""->
function getTetSBPGamma(;degree::Int=1, Tsbp::Type=Float64)
  
  cub, vtx = getTetCubatureGamma(2*degree-1, Tsbp)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 3))
  w, Q = SummationByParts.buildoperators(cub, vtx, degree)
  return TetSBP{Tsbp}(degree, cub, vtx, w, Q)
end

@doc """
### SBP.getTetSBPOmega

Returns SBP-Omega type elements, that have no nodes on the element boundary

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators

**Returns**

* `sbp`: an SBP-Omega operator of the appropriate degree

"""->
function getTetSBPOmega(;degree::Int=1, Tsbp::Type=Float64)
  cub, vtx = getTetCubatureOmega(2*degree-1, Tsbp)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 3))
  w, Q = SummationByParts.buildoperators(cub, vtx, degree)
  return TetSBP{Tsbp}(degree, cub, vtx, w, Q)  
end

@doc """
### SBP.getTetSBPDiagE

Returns SBP-DiagE type elements, whose boundary nodes are positioned at cubature
points

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators

**Returns**

* `sbp`: an SBP-DiagE operator of the appropriate degree

"""->
function getTetSBPDiagE(;degree::Int=1, Tsbp::Type=Float64,
                        edges::Bool=false, vertices::Bool=false)
  @assert( degree >= 1 && degree <= 4 )
  cub, vtx = getTetCubatureDiagE(2*degree, Tsbp, vertices=vertices)
  w = zeros(Tsbp, (cub.numnodes))
  w = SymCubatures.calcweights(cub)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 3))
  if degree == 1
    Q = reshape(readdlm(dirname(@__FILE__)"/tet_diage_p1.dat", '\t', Tsbp),
                size(Q))
  elseif degree == 2
    Q = reshape(readdlm(dirname(@__FILE__)"/tet_diage_p2.dat", '\t', Tsbp),
                size(Q))
  elseif degree == 3
    Q = reshape(readdlm(dirname(@__FILE__)"/tet_diage_p3.dat", '\t', Tsbp),
                size(Q))
  elseif degree == 4
    Q = reshape(readdlm(dirname(@__FILE__)"/tet_diage_p4.dat", '\t', Tsbp),
                size(Q))
  end
  return TetSBP{Tsbp}(degree, cub, vtx, w, Q)
end

@doc """
### SBP.TriFace

Outer constructor for backward compatibility

**Inputs**

* `degree`: face integration is exact for polys of degree 2*`degree`
* `volcub`: cubature associated with the volume domain this face is part of
* `vtx`: vertices of the volume domain
* `vertices`: if true, element vertices are included in the nodes

**Returns**

* `sbpface`: an SBP face type for triangle elements

"""->
function call{T}(::Type{TriFace{T}}, degree::Int, volcub::TriSymCub{T},
                 vtx::Array{T,2}; vertices::Bool=false)
  @assert( degree >= 1 && degree <= 5 )
  R::Array{T,2}
  perm::Array{Int,2}
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

@doc """
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

"""->
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

@doc """
### SBP.TetFace

Outer constructor for backward compatibility

**Inputs**

* `degree`: face integration is exact for polys of degree 2*`degree`
* `volcub`: cubature associated with the volume domain this face is part of
* `vtx`: vertices of the volume domain

**Returns**

* `sbpface`: an SBP face type for tetrahedral elements

"""->
function call{T}(::Type{TetFace{T}}, degree::Int, volcub::TetSymCub{T},
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

@doc """
### SBP.getTetFaceForDiagE

Returns a quadrature that can be used to construct SBP operators with diagonal
boundary operators, E.

**Inputs**

* `degree`: face integration is exact for polys of degree 2*`degree`+1
* `volcub`: cubature rule for the associated volume
* `vtx`: vertices of the triangle

**Returns**

* `sbpface`: an SBP face type for triangle elements

"""->
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
