# SBP abstract and concrete type definitions

"""
### SBP.AbstractSBP

`AbstractSBP` is a parametric abstract type that defines summation-by-parts
finite-difference operators.

"""
abstract AbstractSBP{T<:Number}
#abstract AbstractSBP{T<:AbstractFloat}

"""
### SBP.TriSBP

Defines diagonal-norm SBP first-derivative operators on a line segment.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes on line segment required for these operators
* `cub` : a symmetric cubature type for line segments (usually LG or LGL)
* `vtx` : vertices of the reference element in computational space
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,1]` : discrete stiffness matrix operator
  """
immutable LineSegSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  cub::LineSymCub{T}
  vtx::Array{T,2}
  w::Array{T,1}
  Q::Array{T,3}

  # inner constructor
  function LineSegSBP(degree::Int, cub::LineSymCub{T}, vtx::Array{T,2},
                      w::Array{T,1}, Q::Array{T,3})
    numnodes = cub.numnodes
    @assert( size(Q,1) == size(Q,2) == size(w,1) == numnodes )
    @assert( size(Q,3) == 1 )
    new(degree, numnodes, cub, vtx, w, Q)
  end
end


"""
### SBP.TriSBP

Defines diagonal-norm SBP first-derivative operators on a right-triangle.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the triangle required for these operators
* `cub` : a symmetric cubature type for triangles
* `vtx` : vertices of the reference element in computational space
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""
immutable TriSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  cub::TriSymCub{T}
  vtx::Array{T,2}
  w::Array{T,1}
  Q::Array{T,3}

  # inner constructor
  function TriSBP(degree::Int, cub::TriSymCub{T}, vtx::Array{T,2},
                  w::Array{T,1}, Q::Array{T,3})
    @assert( degree >= 1 && degree <= 4)
    numnodes = cub.numnodes
    @assert( size(Q,1) == size(Q,2) == size(w,1) == numnodes )
    @assert( size(Q,3) == 2 )
    new(degree, numnodes, cub, vtx, w, Q)
  end
end

"""
### SBP.SparseTriSBP

Defines diagonal-norm SBP first-derivative operators on a right-triangle using a
cubature rule that is greater than 2*p-1.  This provides additional flexiblity
in the SBP operator that is used to make a sparse S.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the triangle required for these operators
* `cub` : a symmetric cubature type for triangles
* `vtx` : vertices of the reference element in computational space
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""
immutable SparseTriSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  cub::TriSymCub{T}
  vtx::Array{T,2}
  w::Array{T,1}
  Q::Array{T,3}

  function SparseTriSBP(;degree::Int=1, faceorder::Array{Int,1}=[1;2;3], 
                        internal=false, cubdegree::Int=2*degree+1)
    @assert( degree >= 1 && degree <= 4 )
    if internal
      cub, vtx = getTriCubatureOmega(cubdegree, T)
    else
      cub, vtx = getTriCubatureGamma(cubdegree, T)
    end
    numnodes = cub.numnodes
    Q = zeros(T, (numnodes, numnodes, 2))
    w, Q = SummationByParts.buildsparseoperators(cub, vtx, degree)
    new(degree, numnodes, cub, vtx, w, Q)
  end
end

"""
### SBP.TetSBP

Defines diagonal-norm SBP first-derivative operators on a right-tetrahedron.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the tetrahedron required for these operators
* `cub` : a symmetric cubature type for tetrahedra
* `vtx` : vertices of the reference element in computational space
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""
immutable TetSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  cub::TetSymCub{T}
  vtx::Array{T,2}
  w::Array{T,1}
  Q::Array{T,3}
  
  # inner constructor
  function TetSBP(degree::Int, cub::TetSymCub{T}, vtx::Array{T,2},
                  w::Array{T,1}, Q::Array{T,3})
    @assert( degree >= 1 && degree <= 4)
    numnodes = cub.numnodes
    @assert( size(Q,1) == size(Q,2) == size(w,1) == numnodes )
    @assert( size(Q,3) == 3 )
    new(degree, numnodes, cub, vtx, w, Q)
  end
end

"""
### SBP.SparseTetSBP

Defines diagonal-norm SBP first-derivative operators on a right-tetrahedron
using a cubature rule that is greater than 2*p-1.  This provides additional
flexiblity in the SBP operator that is used to make a sparse S.

**Fields**

* `degree` : maximum polynomial degree for which the derivatives are exact
* `numnodes` : number of nodes in the tetrahedron required for these operators
* `cub` : a symmetric cubature type for tetrahedra
* `vtx` : vertices of the reference element in computational space
* `w` : cubature weights, i.e. the diagonal SBP norm, stored as an array
* `Q[:,:,i]` : discrete stiffness matrix operator in ith coordinate direction

"""
immutable SparseTetSBP{T} <: AbstractSBP{T}
  degree::Int
  numnodes::Int
  cub::TetSymCub{T}
  vtx::Array{T,2}
  w::Array{T,1}
  Q::Array{T,3}

  function SparseTetSBP(;degree::Int=1, faceorder::Array{Int,1}=[1;2;3;4],
                        internal=false, cubdegree::Int=2*degree-1)
    @assert( degree >= 1 && degree <= 3 )
    if internal
      cub, vtx = getTetCubatureOmega(cubdegree, T)
    else
      cub, vtx = getTetCubatureGamma(cubdegree, T)
    end
    numnodes = cub.numnodes
    Q = zeros(T, (numnodes, numnodes, 3))
    w, Q = SummationByParts.buildsparseoperators(cub, vtx, degree)
    new(degree, numnodes, cub, vtx, w, Q)
  end
end

"""
### SBP.AbstractFace

`AbstractFace` is a parametric abstract type that defines face-based data and
operations (e.g. volume-to-face reconstruction, face integration, etc) for
summation-by-parts finite-difference operators.

"""
abstract AbstractFace{T<:Number}

"""
### SBP.DenseFace

`DenseFace` is a parametric abstract type that defines face-based data and
operations (e.g. volume-to-face reconstruction, face integration, etc) for
summation-by-parts finite-difference operators.  This is a subtype for which
interpolation is a dense matrix.

"""
abstract DenseFace{T} <: AbstractFace{T} 

"""
### SBP.LineSegFace

Defines a "face" between two LineSegSBP operators with the same cubature nodes.

**Fields**

* `degree` : face integration is exact for polys of degree 2*`degree`
* `numnodes` : number of cubature nodes (always 1)
* `stencilsize` : number of nodes in the reconstruction stencil
* `dstencilsize` : number of nodes in the derivative operator stencils
* `cub` : a symmetric cubature type for line-segment faces (i.e. points)
* `vtx` : the vertices of the face in reference space, [-1]
* `wface` : mass matrix (quadrature) for the face (always 1.0)
* `interp[:,:]` : volume-to-face-nodes reconstruction operator
* `perm[:,:]` : permutation for volume nodes so `interp` can be used on both sides
* `deriv[:,:]` : derivative operators for face-based coordinate system
* `dperm[:,:]` : permutation for volume nodes so `deriv` can be used on both sides
* `nbrperm[:,:]` : permutation for face nodes on neighbour element

"""
immutable LineSegFace{T} <: DenseFace{T}
  degree::Int
  numnodes::Int
  stencilsize::Int
  dstencilsize::Int
  cub::PointSymCub{T}
  vtx::Array{T,2}
  wface::Array{T,1}
  normal::Array{T,2}
  interp::Array{T,2}
  perm::Array{Int,2}
  deriv::Array{T,3}
  dperm::Array{Int,2}
  nbrperm::Array{Int,2}

  # inner constructor
  function LineSegFace(degree::Int, facecub::PointSymCub{T}, facevtx::Array{T,2},
                       interp::Array{T,2}, perm::Array{Int,2},
                       deriv::Array{T,3}, dperm::Array{Int,2})
    @assert( degree >= 1 )
    numnodes = facecub.numnodes
    @assert( size(interp,2) == size(deriv,2) == numnodes )
    normal = T[-1; 1].'
    nbrperm = SymCubatures.getneighbourpermutation(facecub)
    wface = SymCubatures.calcweights(facecub)
    stencilsize = size(interp,1)
    dstencilsize = size(deriv,1)
    new(degree, facecub.numnodes, stencilsize, dstencilsize, facecub, facevtx, 
        wface, normal, interp, perm, deriv, dperm, nbrperm)
  end
end

"""
### SBP.TriFace

Defines a face between two TriSBP operators with the same cubature nodes

**Fields**

* `degree` : face integration is exact for polys of degree 2*`degree`
* `numnodes` : number of cubature nodes
* `stencilsize` : number of nodes in the reconstruction stencil
* `dstencilsize` : number of nodes in the derivative operator stencils
* `cub` : a symmetric cubature type for triangle faces (i.e. edges)
* `vtx` : the vertices of the face in reference space, [-1,1]
* `wface` : mass matrix (quadrature) for the face
* `interp[:,:]` : volume-to-face-nodes reconstruction operator
* `perm[:,:]` : permutation for volume nodes so `interp` can be used on all sides
* `deriv[:,:]` : derivative operators for face-based coordinate system
* `dperm[:,:]` : permutation for volume nodes so `deriv` can be used on all sides
* `nbrperm[:,:]` : permutation for face nodes on neighbour element

"""
immutable TriFace{T} <: DenseFace{T}
  degree::Int
  numnodes::Int
  stencilsize::Int
  dstencilsize::Int
  cub::LineSymCub{T}
  vtx::Array{T,2}
  wface::Array{T,1}
  normal::Array{T,2}
  interp::Array{T,2}
  perm::Array{Int,2}
  deriv::Array{T,3}
  dperm::Array{Int,2}
  nbrperm::Array{Int,2}

  # inner constructor
  function TriFace(degree::Int, facecub::LineSymCub{T}, facevtx::Array{T,2},
                   interp::Array{T,2}, perm::Array{Int,2},
                   deriv::Array{T,3}, dperm::Array{Int,2})
    @assert( degree >= 1 && degree <= 5 )
    numnodes = facecub.numnodes
    @assert( size(interp,2) == size(deriv,2) == numnodes )
    normal = T[0 -1; 1 1; -1 0].'
    nbrperm = SymCubatures.getneighbourpermutation(facecub)
    wface = SymCubatures.calcweights(facecub)
    stencilsize = size(interp,1)
    dstencilsize = size(deriv,1)
    new(degree, facecub.numnodes, stencilsize, dstencilsize, facecub, facevtx, 
        wface, normal, interp, perm, deriv, dperm, nbrperm)
  end
end

"""
### SBP.TetFace

Defines a face between two TetSBP operators with the same cubature nodes

**Fields**

* `degree` : face integration is exact for polys of degree 2*`degree`
* `numnodes` : number of cubature nodes
* `stencilsize` : number of nodes in the reconstruction stencil
* `dstencilsize` : number of nodes in the derivative operator stencils
* `cub` : a symmetric cubature type for tetrahedral faces (i.e. triangles)
* `vtx` : the vertices of the face in the reference space of the face
* `wface` : mass matrix (quadrature) for the face
* `interp[:,:]` : volume-to-face-nodes reconstruction operator
* `perm[:,:]` : permutation for volume nodes so `interp` can be used on all sides
* `deriv[:,:]` : derivative operators for face-based coordinate system
* `dperm[:,:]` : permutation for volume nodes so `deriv` can be used on all sides
* `nbrperm[:,:]` : permutation for face nodes on neighbour element

"""
immutable TetFace{T} <: DenseFace{T}
  degree::Int
  numnodes::Int
  stencilsize::Int
  #dstencilsize::Int
  cub::TriSymCub{T}
  vtx::Array{T,2}
  wface::Array{T,1}
  normal::Array{T,2}
  interp::Array{T,2}
  perm::Array{Int,2}
  #deriv::Array{T,3}
  #dperm::Array{Int,2}
  nbrperm::Array{Int,2}

  # inner constructor
  function TetFace(degree::Int, facecub::TriSymCub{T}, facevtx::Array{T,2},
                   interp::Array{T,2}, perm::Array{Int,2})
    @assert( degree >= 1 && degree <= 4 )
    numnodes = facecub.numnodes
    @assert( size(interp,2) == numnodes )
    normal = T[0 0 -1; 0 -1 0; 1 1 1; -1 0 0].'
    nbrperm = SymCubatures.getneighbourpermutation(facecub)
    wface = SymCubatures.calcweights(facecub)
    stencilsize = size(interp,1)
    new(degree, numnodes, stencilsize, facecub, facevtx, wface, normal, interp,
        perm, nbrperm)
  end
end

"""
### SBP.SparseFace

`SparseFace` is a parametric abstract type that defines face-based data and
operations (e.g. volume-to-face reconstruction, face integration, etc) for
summation-by-parts finite-difference operators in the case where the
face-cubature nodes and volume nodes coincide (i.e. diagonal E operators).

"""
abstract SparseFace{T} <: AbstractFace{T}

"""
### SBP.TriSparseFace

Defines a face between two TriSBP operators with the same cubature nodes, in
which the face-cubature nodes and volume nodes coincide (i.e. diagonal E
operators).

**Fields**

* `degree` : face integration is exact for polys of degree 2*`degree`
* `numnodes` : number of cubature nodes
* `dstencilsize` : number of nodes in the derivative operator stencils
* `cub` : a symmetric cubature type for triangle faces (i.e. edges)
* `vtx` : the vertices of the face in reference space, [-1,1]
* `wface` : mass matrix (quadrature) for the face
* `perm[:,:]` : maps volume nodes to face nodes on each side
* `deriv[:,:]` : derivative operators for face-based coordinate system
* `dperm[:,:]` : permutation for volume nodes so `deriv` can be used on all sides
* `nbrperm[:,:]` : permutation for face nodes on neighbour element

"""
immutable TriSparseFace{T} <: SparseFace{T}
  degree::Int
  numnodes::Int
  dstencilsize::Int
  cub::LineSymCub{T}
  vtx::Array{T,2}
  wface::Array{T,1}
  normal::Array{T,2}
  perm::Array{Int,2}
  deriv::Array{T,3}
  dperm::Array{Int,2}
  nbrperm::Array{Int,2}

  # inner constructor
  function TriSparseFace(degree::Int, facecub::LineSymCub{T}, facevtx::Array{T,2},
                         perm::Array{Int,2}, deriv::Array{T,3}, dperm::Array{Int,2})
    # @assert( degree >= 1 && degree <= 5 )
    numnodes = facecub.numnodes
    @assert( size(deriv,2) == numnodes )
    normal = T[0 -1; 1 1; -1 0].'
    nbrperm = SymCubatures.getneighbourpermutation(facecub)
    wface = SymCubatures.calcweights(facecub)
    dstencilsize = size(deriv,1)
    new(degree, facecub.numnodes, dstencilsize, facecub, facevtx, wface, normal,
        perm, deriv, dperm, nbrperm)
  end
end

"""
### SBP.TetSparseFace

Defines a face between two TetSBP operators with the same cubature nodes, in
which the face-cubature nodes and volume nodes coincide (i.e. diagonal E
operators).

**Fields**

* `degree` : face integration is exact for polys of degree 2*`degree`
* `numnodes` : number of cubature nodes
* `dstencilsize` : number of nodes in the derivative operator stencils
* `cub` : a symmetric cubature type for tetrahedral faces (i.e. triangles)
* `vtx` : the vertices of the face in the reference space of the face
* `wface` : mass matrix (quadrature) for the face
* `perm[:,:]` : permutation for volume nodes to face nodes on each side
* `deriv[:,:]` : derivative operators for face-based coordinate system
* `dperm[:,:]` : permutation for volume nodes so `deriv` can be used on all sides
* `nbrperm[:,:]` : permutation for face nodes on neighbour element

"""
immutable TetSparseFace{T} <: SparseFace{T}
  degree::Int
  numnodes::Int
  #dstencilsize::Int
  cub::TriSymCub{T}
  vtx::Array{T,2}
  wface::Array{T,1}
  normal::Array{T,2}
  perm::Array{Int,2}
  #deriv::Array{T,3}
  #dperm::Array{Int,2}
  nbrperm::Array{Int,2}

  # inner constructor
  function TetSparseFace(degree::Int, facecub::TriSymCub{T}, facevtx::Array{T,2},
                         perm::Array{Int,2})
    @assert( degree >= 1 && degree <= 4 )
    numnodes = facecub.numnodes
    normal = T[0 0 -1; 0 -1 0; 1 1 1; -1 0 0].'
    nbrperm = SymCubatures.getneighbourpermutation(facecub)
    wface = SymCubatures.calcweights(facecub)
    new(degree, numnodes, facecub, facevtx, wface, normal, perm, nbrperm)
  end
end
