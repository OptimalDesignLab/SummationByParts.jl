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
  cub, vtx = SummationByParts.Cubature.quadrature(2*degree-1, Tsbp, internal=false)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 2))
  w, Q, E = SummationByParts.buildoperators(cub, vtx, degree)
  return LineSegSBP{Tsbp}(degree, cub, vtx, w, Q, E)
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
  cub, vtx = SummationByParts.Cubature.quadrature(2*degree, Tsbp, internal=true)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 2))
  w, Q, E = SummationByParts.buildoperators(cub, vtx, degree)
  return LineSegSBP{Tsbp}(degree, cub, vtx, w, Q, E)
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
  cub, vtx = SummationByParts.Cubature.getTriCubatureGamma(2*degree-1, Tsbp)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 2))
  w, Q, E = SummationByParts.buildoperators(cub, vtx, degree)
  return TriSBP{Tsbp}(degree, cub, vtx, w, Q, E)
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
  cub, vtx = SummationByParts.Cubature.getTriCubatureOmega(2*degree, Tsbp)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 2))
  w, Q, E = SummationByParts.buildoperators(cub, vtx, degree)
  return TriSBP{Tsbp}(degree, cub, vtx, w, Q, E)
end

"""
### SBP.getTriSBPDiagE

Returns SBP-DiagE type elements, whose boundary nodes are positioned at cubature
points

**Inputs**

* `degree`: maximum polynomial degree for which the derivatives are exact
* `Tsbp`: floating point type used for the operators

**Returns**

* `sbp`: an SBP-DiagE operator of the appropriate degree

"""
function getTriSBPDiagE(;degree::Int=1, Tsbp::Type=Float64,
                        vertices::Bool=true, quad_degree::Int=-1)
  if quad_degree==-1
    quad_degree = 2*degree-1
  end
  cub, vtx = SummationByParts.Cubature.getTriCubatureDiagE(quad_degree, Tsbp, vertices=vertices)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 2))
  w, Q, E = SummationByParts.buildoperators(cub, vtx, degree,vertices=vertices)
  #-------
  # D,H,E,Q,S = SummationByParts.buildoperators_pocs(cub, vtx, degree,vertices=vertices)
  # w = diag(H)
  #-------
  # w, Q = SummationByParts.buildMinConditionOperators(cub, vtx, degree, vertices=vertices)
  return TriSBP{Tsbp}(degree, cub, vtx, w, Q, E)
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
  
  cub, vtx = SummationByParts.Cubature.getTetCubatureGamma(2*degree-1, Tsbp)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 3))
  w, Q, E = SummationByParts.buildoperators(cub, vtx, degree)
  return TetSBP{Tsbp}(degree, cub, vtx, w, Q, E)
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
  cub, vtx = SummationByParts.Cubature.getTetCubatureOmega(2*degree-1, Tsbp)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 3))
  w, Q, E = SummationByParts.buildoperators(cub, vtx, degree)
  return TetSBP{Tsbp}(degree, cub, vtx, w, Q, E)  
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
                        edges::Bool=false, faceopertype::Symbol = :DiagE, cubdegree=nothing)
  @assert( degree >= 1 && degree <= 5 )
  @assert isa(cubdegree, Union{Integer, Nothing})
  if isnothing(cubdegree)
    cubdegree = 2*degree-1
  end
  cub, vtx = SummationByParts.Cubature.getTetCubatureDiagE(cubdegree, Tsbp, faceopertype=faceopertype)
  w = zeros(Tsbp, (cub.numnodes))
  w = SymCubatures.calcweights(cub)
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 3))
  # if degree == 1
  #   Q = reshape(readdlm(joinpath(dirname(@__FILE__), "tet_diage_p1.dat"), '\t', Tsbp),
  #               size(Q))
  # elseif degree == 2
  #   Q = reshape(readdlm(joinpath(dirname(@__FILE__), "tet_diage_p2.dat"), '\t', Tsbp),
  #               size(Q))
  # elseif degree == 3
  #   Q = reshape(readdlm(joinpath(dirname(@__FILE__), "tet_diage_p3.dat"), '\t', Tsbp),
  #               size(Q))
  # elseif degree == 4
  #   Q = reshape(readdlm(joinpath(dirname(@__FILE__), "tet_diage_p4.dat"), '\t', Tsbp),
  #               size(Q))
  # end
  w, Q, E = SummationByParts.buildoperators(cub, vtx, degree, faceopertype=faceopertype)
  return TetSBP{Tsbp}(degree, cub, vtx, w, Q, E)
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
function getLineSegFace(degree::Int, volcub::LineSymCub{T}, vtx::Array{T,2}) where {T}
  facecub, facevtx = SummationByParts.Cubature.pointCubature()
  R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx,
                                                     degree)
  D, Dperm = SummationByParts.buildfacederivatives(facecub, volcub, vtx, degree)
  wface = SymCubatures.calcweights(facecub)
  stencilsize = size(R,2)
  dstencilsize = size(D,1)
  return LineSegFace{T}(degree, facecub, facevtx, Matrix(R'), perm, D, Dperm)
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
function (::Type{TriFace{T}})(degree::Int, volcub::TriSymCub{T},
                                 vtx::Array{T,2}; vertices::Bool=false) where {T}
  @assert( degree >= 1 && degree <= 10 )
  local R::Array{T,2}
  local perm::Array{Int,2}
  if vertices
    facecub, facevtx = SummationByParts.Cubature.quadrature(2*degree, T, internal=false)
    R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx, degree+1)
  else
    facecub, facevtx = SummationByParts.Cubature.quadrature(2*degree, T, internal=true)
    R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx, degree)    
  end
  D, Dperm = SummationByParts.buildfacederivatives(facecub, volcub, vtx,
                                                   degree)
  nbrperm = SymCubatures.getneighbourpermutation(facecub)
  wface = SymCubatures.calcweights(facecub)
  stencilsize = size(R,2)
  dstencilsize = size(D,1)
  TriFace{T}(degree, facecub, facevtx, Matrix(R'), perm, D, Dperm)
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
function getTriFaceForDiagE(degree::Int, volcub::TriSymCub{T},
                               vtx::Array{T,2}; vertices::Bool=true) where {T}
  #@assert( degree >= 1 && degree <= 4 )
  if vertices
    facecub, facevtx = SummationByParts.Cubature.quadrature(2*degree, T, internal=false)
    R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx, degree+1)
  else
    facecub, facevtx = SummationByParts.Cubature.quadrature(2*degree, T, internal=true)
    R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx, degree)
  end
  perm_red = zeros(Int, (facecub.numnodes,3))
  for i = 1:facecub.numnodes
    node = argmax(R[i,:])
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
function (::Type{TetFace{T}})(degree::Int, volcub::TetSymCub{T},
                                 vtx::Array{T,2}; faceopertype::Symbol = :Omega,
                                 facecub::Union{TriSymCub{T}, Nothing}=nothing, 
                                 facevtx::Union{Array{T,2}, Nothing}=nothing) where {T}
  @assert( degree >= 1 && degree <= 5 )
  # facecub, facevtx = SummationByParts.Cubature.getTriCubatureOmega(2*degree, T)
  if (facecub === nothing || facevtx === nothing)
    if faceopertype == :DiagE
      facecub, facevtx = SummationByParts.Cubature.getTriCubatureForTetFaceDiagE(2*degree, T)
      # facecub, facevtx = SummationByParts.Cubature.getTriCubatureOmega(2*degree, T)
    else
      facecub, facevtx = SummationByParts.Cubature.getTriCubatureOmega(2*degree, T)
    end
  end
  
  normal = T[0 0 -1; 0 -1 0; 1 1 1; -1 0 0]'
  R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx, degree, faceopertype=faceopertype)
  #D, Dperm = SummationByParts.buildfacederivatives(facecub, volcub, vtx,
  #                                                 degree)
  nbrperm = SymCubatures.getneighbourpermutation(facecub)
  wface = SymCubatures.calcweights(facecub)
  stencilsize = size(R,2)
  #dstencilsize = size(D,1)
  #new(degree, facecub.numnodes, stencilsize, dstencilsize, facecub, wface,
  #    normal, R', perm, D, Dperm, nbrperm)
  TetFace{T}(degree, facecub, facevtx, Matrix(R'), perm)
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
function getTetFaceForDiagE(degree::Int, volcub::TetSymCub{T},
                               vtx::Array{T,2}; faceopertype::Symbol = :DiagE) where {T}
  @assert( degree >= 1 && degree <= 4 )
  
  # facecub, facevtx = SummationByParts.Cubature.getTriCubatureOmega(2*degree, T)
  if faceopertype == :DiagE
    facecub, facevtx = SummationByParts.Cubature.getTriCubatureDiagE(2*degree, T)
  else
    facecub, facevtx = SummationByParts.Cubature.getTriCubatureOmega(2*degree, T)
  end
  
  #facecub = SymCubatures.TriSymCub{T}(vertices=true, numS21=1)
  #SymCubatures.setparams!(facecub, T[0.9055050463303657])
  #SymCubatures.setweights!(facecub, T[(1+sqrt(21))/60;(39+sqrt(21))/60])
  #facevtx = T[-1 -1; 1 -1; -1 1]
  R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx, degree, faceopertype=faceopertype)
  #println("TEMP change in getTriFaceForDiagE!!!!")
  #facecub, facevtx = quadrature(2*degree, T, internal=true)
  #R, perm = SummationByParts.buildfacereconstruction(facecub, volcub, vtx,
  #                                                   degree)

  perm_red = zeros(Int, (facecub.numnodes,4))
  for i = 1:facecub.numnodes
    node = argmax(R[i,:])
    perm_red[i,:] = perm[node,:]
  end    
  nbrperm = SymCubatures.getneighbourpermutation(facecub)
  wface = SymCubatures.calcweights(facecub)
  return TetSparseFace{T}(degree, facecub, facevtx, perm_red)
end
