# This file gathers together outer constructors for the SBP operators

@doc """
### SBP.TriSBP

Outer constructor for backward compatibility

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `internal`: if true, all element nodes are strictly internal
* `vertices`: if true, element vertices are included in the nodes

**Returns**

* `sbp`: an SBP type for triangular elements

"""->
function call{T}(::Type{TriSBP{T}}; degree::Int=1, internal::Bool=false,
                 vertices::Bool=true)
  cub, vtx = tricubature(2*degree-1, T, internal=internal,
                         vertices=vertices)
  Q = zeros(T, (cub.numnodes, cub.numnodes, 2))
  w, Q = SummationByParts.buildoperators(cub, vtx, degree)
  TriSBP{T}(degree, cub, vtx, w, Q)
end

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
  return TriSBP{Tsbp}(degree=degree, internal=false, vertices=true)
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
  return TriSBP{Tsbp}(degree=degree, internal=true, vertices=false)
end

@doc """
### SBP.TetSBP

Outer constructor for backward compatibility

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `internal`: if true, all element nodes are strictly internal
* `vertices`: if true, element vertices are included in the nodes

**Returns**

* `sbp`: an SBP type for tetrahedral elements

"""->
function call{T}(::Type{TetSBP{T}}; degree::Int=1, internal::Bool=false)
  cub, vtx = tetcubature(2*degree-1, T, internal=internal)
  Q = zeros(T, (cub.numnodes, cub.numnodes, 3))
  w, Q = SummationByParts.buildoperators(cub, vtx, degree)
  TetSBP{T}(degree, cub, vtx, w, Q)
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
  return TetSBP{Tsbp}(degree=degree, internal=false, vertices=true)
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
  return TetSBP{Tsbp}(degree=degree, internal=true, vertices=false)
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
  if vertices
    facecub, facevtx = quadrature(2*degree, T, internal=false)
  else
    facecub, facevtx = quadrature(2*degree, T, internal=true)
  end
  R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx,
                                                     degree)
  #R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx,
  #                                                   degree+1)
  D, Dperm = SummationByParts.buildfacederivatives(facecub, volcub, vtx,
                                                   degree)
  nbrperm = SymCubatures.getneighbourpermutation(facecub)
  wface = SymCubatures.calcweights(facecub)
  stencilsize = size(R,2)
  dstencilsize = size(D,1)
  TriFace{T}(degree, facecub, facevtx, R.', perm, D, Dperm)
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
  facecub, facevtx = tricubature(2*degree, T, internal=true)
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
