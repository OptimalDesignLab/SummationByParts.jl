module SymCubatures
# types and methods for mapping between symmetry groups and nodes for cubatures
# on various domains

export SymCub, PointSymCub, LineSymCub, TriSymCub, TetSymCub

"""
### SymCubatures.SymCub

`SymCub` is an parametric abstract type that defines cubatures for symmetric
nodal distributions.  It is parameterized on `T` in order to allow for future
implementations of arbitrary precision types.  The parameterization also permits
the use of the complex-step method for verification.

""" abstract type SymCub{T<:Number} end

"""
### SymCubatures.PointSymCub

Defines the trivial point cubature for uniformity across methods.

**Fields**

* `numparams` : total number of nodal degrees of freedom
* `numweights` : total number of unique weights
* `numnodes` : total number of nodes
* `vertices` : if true, vertices (ends of interval) are in the set of nodes
* `centroid` : if true, centroid is present in set of nodes
* `numedge` : number of unique edge parameters
* `numsym` : number of node sets in each symmetry group (0 for [1])
* `params` : the actual values of the orbit nodal parameters
* `weights` : values of the unique weights

"""
type PointSymCub{T} <: SymCub{T}
  numparams::Int
  numweights::Int
  numnodes::Int
  vertices::Bool
  centroid::Bool
  numedge::Int
  numsym::Array{Int,1}
  params::Array{T,1}
  weights::Array{T,1}

  function PointSymCub{T}() where T
    numparams = 0
    numweights = 1
    numnodes = 1
    vertices = true
    centroid = false
    numedge = 0
    numsym = [1;]
    # initialize parameter arrays
    @assert(numweights == sum(numsym))
    params = zeros(T, numparams)
    weights = zeros(T, numweights)
    new(numparams, numweights, numnodes, vertices, centroid, numedge, numsym,
        params, weights)
  end
end

"""
### SymCubatures.LineSymCub
  
Defines a symmetric quadrature rule on the interval [-1,1].  Current choices are
Legendre-Gauss-Lobatto (LGL) or Legendre-Gauss (LG) rules.

**Fields**

* `numparams` : total number of nodal degrees of freedom
* `numweights` : total number of unique weights
* `numnodes` : total number of nodes
* `vertices` : if true, vertices (ends of interval) are in the set of nodes
* `centroid` : if true, centroid is present in set of nodes
* `numedge` : number of unique edge parameters
* `numsym` : number of node sets in each symmetry group (2 for [-1,1])
* `params` : the actual values of the orbit nodal parameters
* `weights` : values of the unique weights

"""
type LineSymCub{T} <: SymCub{T}
  numparams::Int
  numweights::Int
  numnodes::Int
  vertices::Bool
  centroid::Bool
  numedge::Int
  numsym::Array{Int,1}
  params::Array{T,1}
  weights::Array{T,1}
    
  function LineSymCub{T}(;numedge::Int=0, vertices::Bool=true,
                         centroid::Bool=false) where T
    @assert(numedge >= 0)
    # compute the number of degrees of freedom and unique weights
    numparams = 0
    numweights = 0
    numparams += numedge
    # compute the number of nodes and number of symmetries
    numsym = zeros(Int, 2)
    numnodes = 0
    if vertices
      numnodes += 2
      numweights += 1
      numsym[2] += 1
    end
    if centroid
      numnodes += 1
      numweights += 1
      numsym[1] += 1
    end
    numnodes += 2*numedge # 2 * (3 edges) X (number of edge parameters)
    numsym[2] += numedge
    numweights += numedge
    # initialize parameter arrays
    @assert(numweights == sum(numsym))
    params = zeros(T, numparams)
    weights = zeros(T, numweights)
    new(numparams, numweights, numnodes, vertices, centroid, numedge, numsym,
        params, weights)
  end
end

"""
### SymCubatures.TriSymCub

Used to define symmetric cubature rules on the triangle.  The `params` array
determines the position of the parameterized nodes, and the `weights` array
determines the value of the weight for each symmetric orbit.  Note that boolean
fields are used to activate some degenerate orbits.  For example, vertices are a
special case of several of the orbits, and should be activated by setting
vertices=true rather than relying on a specific value of a parameter.

**Fields**

* `numparams` : total number of nodal degrees of freedom
* `numweights` : total number of unique weights
* `numnodes` : total number of nodes
* `vertices` : if true, vertices are present in the set of nodes
* `midedges` : if true, edge midpoints are present in set of nodes
* `centroid` : if true, centroid is present in set of nodes
* `numedge` : number of unique edge parameters
* `numS21` : number of S21 orbits (vertex to opposite face)
* `numS111` : number of S111 orbits
* `numsym` : number of node sets in each symmetry group (3 groups for Tri)
* `params` : the actual values of the orbit nodal parameters
* `weights` : values of the unique weights

"""
type TriSymCub{T} <: SymCub{T}
  numparams::Int
  numweights::Int
  numnodes::Int
  vertices::Bool
  midedges::Bool
  centroid::Bool
  numedge::Int
  numS21::Int
  numS111::Int
  numsym::Array{Int,1}
  params::Array{T,1}
  weights::Array{T,1}

  function TriSymCub{T}(;numedge::Int=0, numS21::Int=0, numS111::Int=0,
                        vertices::Bool=true, midedges::Bool=false,
                        centroid::Bool=false) where T
    @assert(numedge >= 0)
    @assert(numS21 >= 0)
    @assert(numS111 >= 0)
    # compute the number of degrees of freedom and unique weights
    numparams = 0
    numweights = 0
    numparams += numedge + numS21 + 2*numS111
    numsym = zeros(Int, 3)
    # compute the number of nodes and symmetries
    numnodes = 0
    if vertices
      numnodes += 3
      numweights += 1
      numsym[2] += 1
    end
    if midedges
      numnodes += 3
      numweights += 1
      numsym[2] += 1
    end
    if centroid
      numnodes += 1
      numweights += 1
      numsym[1] += 1
    end
    numnodes += 6*numedge # 2 * (3 edges) X (number of edge parameters)
    numsym[3] += numedge
    numweights += numedge
    numnodes += 3*numS21 # (3 permutations) X (number of S21 orbits)
    numsym[2] += numS21
    numweights += numS21
    numnodes += 6*numS111 # (6 permutations) X (number of S111 orbits)
    numsym[3] += numS111
    numweights += numS111
    # initialize parameter arrays
    @assert(numweights == sum(numsym))
    params = zeros(T, numparams)
    weights = zeros(T, numweights)
    new(numparams, numweights, numnodes, vertices, midedges, centroid,
        numedge, numS21, numS111, numsym, params, weights)
  end
end

"""
### SymCubatures.TetSymCub

Used to define symmetric cubature rules on the tetrahedron.  The `params` array
determines the position of the parameterized nodes, and the `weights` array
determines the value of the weight for each symmetric orbit.  Note that boolean
fields are used to activate some degenerate symmetry orbits.  For example,
vertices are a special case of several of the orbits, and should be activated by
setting vertices=true rather than relying on a specific value of a parameter.

**Fields**

* `numparams` : total number of nodal degrees of freedom
* `numweights` : total number of unique weights
* `numnodes` : total number of nodes
* `vertices` : if true, vertices are present in the set of nodes
* `midedges` : if true, edge midpoints are present in set of nodes
* `centroid` : if true, centroid is present in set of nodes
* `facecentroid` : if true, face centroids are present in the set of nodes
* `numedge` : number of unique edge parameters
* `numfaceS21` : number of S21 face orbits (same tri orbit on face)
* `numfaceS111` : number of S111 face orbits (same tri orbit on face)
* `numS31` : number of S31 orbits (vertex to opposite face)
* `numS22` : number of S22 orbits
* `numS211`: number of S211 orbits
* `numS1111`: number of S1111 orbits
* `numsym` : number of node sets in each symmetry group (5 groups for Tet)
* `params` : the actual values of the orbit nodal parameters
* `weights` : values of the unique weights

"""
type TetSymCub{T} <: SymCub{T}
  numparams::Int
  numweights::Int
  numnodes::Int
  vertices::Bool
  midedges::Bool
  centroid::Bool
  facecentroid::Bool
  numedge::Int
  numfaceS21::Int
  numfaceS111::Int
  numS31::Int
  numS22::Int
  numS211::Int
  numS1111::Int
  numsym::Array{Int,1}
  params::Array{T,1}
  weights::Array{T,1}

  function TetSymCub{T}(;numedge::Int=0, numfaceS21::Int=0, numfaceS111::Int=0,
                        numS31::Int=0, numS22::Int=0, numS211::Int=0,
                        numS1111::Int=0, vertices::Bool=true,
                        midedges::Bool=false, centroid::Bool=false,
                        facecentroid::Bool=false) where T
    @assert(numedge >= 0)
    @assert(numfaceS21 >= 0)
    @assert(numfaceS111 >= 0)
    @assert(numS31 >= 0)
    @assert(numS22 >= 0)
    @assert(numS211 >= 0)
    @assert(numS1111 >= 0)
    # compute the number of degrees of freedom and unique weights
    numparams = 0
    numweights = 0
    numparams += (numedge + numfaceS21 + 2*numfaceS111 + numS31 + numS22
                  + 2*numS211 + 3*numS1111)
    numsym = zeros(Int, 5)
    # compute the number of nodes and symmetries
    numnodes = 0
    if vertices
      numnodes += 4
      numweights += 1
      numsym[2] += 1
    end
    if midedges
      numnodes += 6
      numweights += 1
      numsym[3] += 1
    end
    if centroid
      numnodes += 1
      numweights += 1
      numsym[1] += 1
    end
    if facecentroid
      numnodes += 4
      numweights += 1
      numsym[2] += 1
    end
    numnodes += 12*numedge # 2 X (6 edges) X (number of edge parameters)
    numsym[4] += numedge
    numweights += numedge
    numnodes += 12*numfaceS21 # 3 X (4 faces) X (number of face S21 orbits)
    numsym[4] += numfaceS21
    numweights += numfaceS21
    numnodes += 24*numfaceS111 # 6 X (4 faces) X (number of face S111 orbits)
    numsym[5] += numfaceS111
    numweights += numfaceS111
    numnodes += 4*numS31 # (4 permutations) X (number of S31 orbits)
    numsym[2] += numS31
    numweights += numS31
    numnodes += 6*numS22 # (6 permutations) X (number of S22 orbits)
    numsym[3] += numS22
    numweights += numS22
    numnodes += 12*numS211 # (12 permutations) X (number of S211 orbits)
    numsym[4] += numS211
    numweights += numS211
    numnodes += 24*numS1111 # (24 permutations) X (number of S1111 orbits)
    numsym[5] += numS1111
    numweights += numS1111
    # initialize parameter arrays
    @assert(numweights == sum(numsym))
    params = zeros(T, numparams)
    weights = zeros(T, numweights)
    new(numparams, numweights, numnodes, vertices, midedges, centroid,
        facecentroid, numedge, numfaceS21, numfaceS111, numS31, numS22, numS211,
        numS1111, numsym, params, weights)
  end
end

"""
### SymCubatures.getnumboundarynodes

Returns the number of (explicit) boundary nodes

*Notes*: if the parameter value for an internal orbit is such that the
corresponding node lies on the boundary, this node is **NOT** included in the
boundary-node count returned.

**Inputs**

* `cub`: symmetric cubature rule

**Outputs**

* `numboundary`: number of boundary nodes

"""
function getnumboundarynodes{T}(cub::PointSymCub{T})
  return 1
end

function getnumboundarynodes{T}(cub::LineSymCub{T})
  cub.vertices ? (return 2) : (return 0)
end

function getnumboundarynodes{T}(cub::TriSymCub{T})
  numboundary = 0
  cub.vertices ? numboundary += 3 : nothing
  cub.midedges ? numboundary += 3 : nothing
  numboundary += 6*cub.numedge
  return numboundary
end

function getnumboundarynodes{T}(cub::TetSymCub{T})
  numboundary = 0
  cub.vertices ? numboundary += 4 : nothing
  cub.midedges ? numboundary += 6 : nothing
  cub.facecentroid ? numboundary += 4 : nothing
  numboundary += 12*cub.numedge
  numboundary += 12*cub.numfaceS21
  numboundary += 24*cub.numfaceS111
  return numboundary
end

"""
### SymCubatures.getnumfacenodes

Returns the number of nodes on an individual face of the element.

**Inputs**

* `cub`: symmetric cubature rule

**Outputs**

* `numfacenodes`: number of nodes on a face

"""
function getnumfacenodes{T}(cub::PointSymCub{T})
  return 1
end

function getnumfacenodes{T}(cub::LineSymCub{T})
  cub.vertices ? (return 1) : (return 0)
end

function getnumfacenodes{T}(cub::TriSymCub{T})
  numfacenodes = 0
  cub.vertices ? numfacenodes += 2 : nothing
  cub.midedges ? numfacenodes += 1 : nothing
  numfacenodes += 2*cub.numedge
  return numfacenodes
end

function getnumfacenodes{T}(cub::TetSymCub{T})
  numfacenodes = 0
  cub.vertices ? numfacenodes += 3 : nothing
  cub.midedges ? numfacenodes += 3 : nothing
  cub.facecentroid ? numfacenodes += 1 : nothing
  numfacenodes += 6*cub.numedge
  numfacenodes += 3*cub.numfaceS21
  numfacenodes += 6*cub.numfaceS111
  return numfacenodes
end

"""
### SymCubatures.getbndrynodeindices

Returns the indices of the nodes that lie on the boundary, in their natural
order.  See getfacenodeindices for a method returns node indices for each face.

**Inputs**

* `cub`: a symmetric cubature rule whose boundary-node indices are sought

**Outputs**

* `bndryindices`: indicies of nodes that lie on boundary

"""
function getbndrynodeindices{T}(cub::PointSymCub{T})
  bndryindices = zeros(Int, getnumboundarynodes(cub))
  bndryindices[:] = [1;]
  return bndryindices
end

function getbndrynodeindices{T}(cub::LineSymCub{T})
  bndryindices = zeros(Int, getnumboundarynodes(cub))
  # add vertices to the indices
  if cub.vertices
    bndryindices[1:2] = [1;2]
  end
  return bndryindices
end

function getbndrynodeindices{T}(cub::TriSymCub{T})
  bndryindices = zeros(Int, getnumboundarynodes(cub))
  ptr = 0
  idxptr = 0
  # add vertices to indices
  if cub.vertices
    bndryindices[idxptr+1:idxptr+3] = [ptr+1; ptr+2; ptr+3]
    ptr += 3
    idxptr += 3
  end
  # add midedge nodes to indices
  if cub.midedges
    bndryindices[idxptr+1:idxptr+3] = [ptr+1; ptr+2; ptr+3]
    ptr += 3
    idxptr += 3
  end
  # account for S21 orbits
  ptr += 3*cub.numS21
  # add remaining edge nodes to indices
  for i = 1:cub.numedge
    bndryindices[idxptr+1:idxptr+6] = [ptr+1; ptr+2; ptr+3;
                                       ptr+4; ptr+5; ptr+6]
    ptr += 6
    idxptr += 6
  end
  return bndryindices
end

function getbndrynodeindices{T}(cub::TetSymCub{T})  
  bndryindices = zeros(Int, getnumboundarynodes(cub))
  ptr = 0
  idxptr = 0
  # add vertices to indices
  if cub.vertices
    bndryindices[idxptr+1:idxptr+4] = [ptr+1; ptr+2; ptr+3; ptr+4]
    ptr += 4
    idxptr += 4
  end
  # add face centroids
  if cub.facecentroid
    bndryindices[idxptr+1:idxptr+4] = [ptr+1; ptr+2; ptr+3; ptr+4]
    ptr += 4
    idxptr += 4
  end
  # account for S31 orbits
  ptr += 4*cub.numS31
  # add mid-edge nodes
  if cub.midedges
    bndryindices[idxptr+1:idxptr+6] = [ptr+1; ptr+2; ptr+3;
                                       ptr+4; ptr+5; ptr+6]
    ptr += 6
    idxptr += 6
  end
  # account for S22 nodes
  ptr += 6*cub.numS22
  # add remaining edge nodes to indices
  for i = 1:cub.numedge
    bndryindices[idxptr+1:idxptr+12] = [ptr+1; ptr+2; ptr+3;
                                        ptr+4; ptr+5; ptr+6;
                                        ptr+7; ptr+8; ptr+9;
                                        ptr+10; ptr+11; ptr+12]
    ptr += 12
    idxptr += 12
  end
  # add face nodes corresponding to S21 orbit
  for i = 1:cub.numfaceS21
    bndryindices[idxptr+1:idxptr+12] = [ptr+1; ptr+2; ptr+3;
                                        ptr+4; ptr+5; ptr+6;
                                        ptr+7; ptr+8; ptr+9;
                                        ptr+10; ptr+11; ptr+12]
    ptr += 12
    idxptr += 12
  end
  # account for S211 nodes
  ptr += 12*cub.numS211
  # add face nodes corresponding to S111 orbit
  for i = 1:cub.numfaceS111
    bndryindices[idxptr+1:idxptr+24] = (ptr+1):(ptr+24)
    ptr += 24
    idxptr += 24
  end  
  return bndryindices
end

"""
### SymCubatures.getinteriornodeindices

Returns the indices of the nodes that are strictly interior.

**Inputs**

* `cub`: a symmetric cubature rule whose interior-node indices are sought

**Outputs**

* `indices`: indicies of nodes that are strictly interior.

"""
function getinteriornodeindices{T}(cub::PointSymCub{T})
  numinterior = cub.numnodes - getnumboundarynodes(cub)
  @assert( numinterior == 0)
  indices = zeros(Int, numinterior)
  return indices
end

function getinteriornodeindices{T}(cub::LineSymCub{T})
  numinterior = cub.numnodes - getnumboundarynodes(cub)
  indices = zeros(Int, numinterior)
  ptr = 0
  if cub.vertices
    ptr += 2
  end
  indices = [ptr+1:ptr+numinterior;]
  return indices
end

function getinteriornodeindices{T}(cub::TriSymCub{T})
  numinterior = cub.numnodes - getnumboundarynodes(cub)
  indices = zeros(Int, numinterior)
  ptr = 0
  idxptr = 0
  # account for vertices
  if cub.vertices
    ptr += 3
  end
  # account for mid-edge nodes
  if cub.midedges
    ptr += 3
  end
  # set S21 orbit node indices
  for i = 1:cub.numS21
    indices[idxptr+1:idxptr+3] = [ptr+1:ptr+3;]
    idxptr += 3
    ptr += 3
  end
  # account for edge nodes
  ptr += 6*cub.numedge
  # set S111 orbit node indices
  for i = 1:cub.numS111
    indices[idxptr+1:idxptr+6] = [ptr+1:ptr+6;]
    ptr += 6
    idxptr += 6
  end
  # set centroid index
  if cub.centroid
    indices[idxptr+1] = ptr+1
    ptr += 1
    idxptr += 1
  end
  return indices
end

function getinteriornodeindices{T}(cub::TetSymCub{T})
  numinterior = cub.numnodes - getnumboundarynodes(cub)
  indices = zeros(Int, numinterior)
  ptr = 0
  idxptr = 0
  # account for vertices
  if cub.vertices
    ptr += 4
  end
  # account for face centroids
  if cub.facecentroid
    ptr += 4
  end
  # set S31 orbit node indices
  for i = 1:cub.numS31
    indices[idxptr+1:idxptr+4] = [ptr+1:ptr+4;]
    ptr += 4
    idxptr += 4
  end
  # account for mid-edge nodes
  if cub.midedges
    ptr += 6
  end
  # set S22 orbit node indices
  for i = 1:cub.numS22
    indices[idxptr+1:idxptr+6] = [ptr+1:ptr+6;]
    ptr += 6
    idxptr += 6
  end
  # account for edge nodes
  ptr += 12*cub.numedge
  # account for face nodes corresponding to S21 orbit
  ptr += 12*cub.numfaceS21
  # set S211 orbit node indices
  for i = 1:cub.numS211
    indices[idxptr+1:idxptr+12] = [ptr+1:ptr+12;]
    ptr += 12
    idxptr += 12
  end
  # account for face nodes corresponding to S111 orbit
  ptr += 24*cub.numfaceS111
  # set S1111 oribt node indices
  for i = 1:cub.numS1111
    indices[idxptr+1:idxptr+24] = [ptr+1:ptr+24;]
    ptr += 24
    idxptr += 24
  end
  # set centroid node
  if cub.centroid
    indices[idxptr+1] = ptr+1
    ptr += 1
    idxptr += 1
  end
  return indices
end

"""
### SymCubatures.getfacevertexindices

Returns the indices of the vertices that make up each face.  This is useful when
building nodes on a given face using Barycentric coordinates

**Inputs**

* `cub`: a symmetric cubature rule whose face vertices are sought

**Outputs**

* `facevtx`: subarray `facevtx[:,f]` lists the vertices of face `f` 

"""
function getfacevertexindices{T}(cub::PointSymCub{T})
  return [1]
end

function getfacevertexindices{T}(cub::LineSymCub{T})
  return [1 2]
end

function getfacevertexindices{T}(cub::TriSymCub{T})
  return [1 2 3; 2 3 1]
end

function getfacevertexindices{T}(cub::TetSymCub{T})
  return [1 1 2 1; 2 4 4 3; 3 2 3 4]
end

"""
### SymCubatures.getfacenodeindices

Returns the indices of the nodes that lie on each face.  See getbndrynodeindices
for a method that returns a single array of boundary nodes.

**Inputs**

* `cub`: a symmetric cubature rule whose boundary-node indices are sought

**Outputs**

* `bndryindices`: indicies of nodes that lie on boundary; there is a separate
  column of indices for each edge/face.

"""
function getfacenodeindices{T}(cub::PointSymCub{T})
  numedge = getnumfacenodes(cub)
  bndryindices = zeros(Int, (numedge,1) )
  bndryindices[:,:] = 1
  return bndryindices
end

function getfacenodeindices{T}(cub::LineSymCub{T})
  numedge = getnumfacenodes(cub)
  bndryindices = zeros(Int, (numedge,2) )
  # add vertices to indices
  if cub.vertices
    bndryindices[:,:] = [1 2]
  end
  return bndryindices
end

function getfacenodeindices{T}(cub::TriSymCub{T})
  # get the number of nodes on one edge
  numedge = getnumfacenodes(cub)
  bndryindices = zeros(Int, (numedge,3) )
  ptr = 0
  idxptr = 0
  # add vertices to indices
  if cub.vertices
    bndryindices[idxptr+1:idxptr+2,:] = [ptr+1 ptr+2 ptr+3;
                                         ptr+2 ptr+3 ptr+1] 
    ptr += 3
    idxptr += 2
  end
  # add midedge nodes to indices
  if cub.midedges
    bndryindices[idxptr+1,:] = [ptr+1 ptr+2 ptr+3]
    ptr += 3
    idxptr += 1
  end
  # account for S21 orbits
  ptr += 3*cub.numS21
  # add remaining edge nodes to indices
  for i = 1:cub.numedge
    bndryindices[idxptr+1:idxptr+2,:] = [ptr+1 ptr+3 ptr+5;
                                         ptr+2 ptr+4 ptr+6]
    ptr += 6
    idxptr += 2
  end
  return bndryindices
end

function getfacenodeindices{T}(cub::TetSymCub{T})
  # get the number of nodes on one face
  numface = getnumfacenodes(cub)
  bndryindices = zeros(Int, (numface,4) )
  ptr = 0
  idxptr = 0
  # add vertices to indices
  if cub.vertices
    bndryindices[idxptr+1:idxptr+3,:] = [ptr+1 ptr+1 ptr+2 ptr+1;
                                         ptr+2 ptr+4 ptr+4 ptr+3;
                                         ptr+3 ptr+2 ptr+3 ptr+4] 
    ptr += 4
    idxptr += 3
  end
  # add face centroids to indices
  if cub.facecentroid
    bndryindices[idxptr+1,:] = [ptr+1 ptr+2 ptr+3 ptr+4]
    ptr += 4
    idxptr += 1
  end
  # account for S31 orbits
  ptr += 4*cub.numS31
  # add mid-edge to indices
  if cub.midedges
    bndryindices[idxptr+1:idxptr+3,:] = [ptr+1 ptr+4 ptr+5 ptr+3;
                                         ptr+2 ptr+5 ptr+6 ptr+6;
                                         ptr+3 ptr+1 ptr+2 ptr+4]
    ptr += 6
    idxptr += 3
  end
  # account for S22 orbits
  ptr += 6*cub.numS22
  # add edge nodes to indices
  for i = 1:cub.numedge
    bndryindices[idxptr+1:idxptr+6,:] = [ptr+1 ptr+7  ptr+9  ptr+6;
                                         ptr+2 ptr+8  ptr+10 ptr+5;
                                         ptr+3 ptr+10 ptr+12 ptr+11;
                                         ptr+4 ptr+9  ptr+11 ptr+12;
                                         ptr+5 ptr+2  ptr+4  ptr+8;
                                         ptr+6 ptr+1  ptr+3  ptr+7]
    ptr += 12
    idxptr += 6
  end
  # add face S21 orbits to indices
  for i = 1:cub.numfaceS21
    bndryindices[idxptr+1:idxptr+3,:] = [ptr+1 ptr+4 ptr+7 ptr+10;
                                         ptr+2 ptr+5 ptr+8 ptr+11;
                                         ptr+3 ptr+6 ptr+9 ptr+12]
    ptr += 12
    idxptr += 3
  end
  # account for S211 oribts
  ptr += 12*cub.numS211
  # add face S111 orbits to indices
  for i = 1:cub.numfaceS111
    bndryindices[idxptr+1:idxptr+6,:] = [ptr+1 ptr+7  ptr+13 ptr+19;
                                         ptr+2 ptr+8  ptr+14 ptr+20;
                                         ptr+3 ptr+9  ptr+15 ptr+21;
                                         ptr+4 ptr+10 ptr+16 ptr+22;
                                         ptr+5 ptr+11 ptr+17 ptr+23;
                                         ptr+6 ptr+12 ptr+18 ptr+24]
    ptr += 24
    idxptr += 6
  end    
  return bndryindices
end

"""
### SymCubatures.findleftperm!

For a matrix `A`, we are given a right permutation of the columns, `A[:,permR]`.
This function attempts to find the left permultation of rows such that
                      `A[permL,:] = A[:,permR]`

**Inputs**

* `A`: a rectangular matrix for which the left permutation is sought
* `permR`: the given right permutation of the columns

**Outputs**

* `permL`: the left permutation of the rows, if it exists

**Returns**

* `true` if the permutation exists, `false` otherwise

"""
function findleftperm!{T}(A::AbstractArray{T,2}, permR::AbstractVector{Int},
                          permL::AbstractVector{Int})
  @assert( size(A,1) == length(permL) )
  @assert( size(A,2) == length(permR) )
  rows = [ view(A, i, 1:size(A,2)) for i=1:size(A,1) ]
  permA = sortperm(rows; order=Base.Lexicographic)
  AR = A[:,permR]
  rows = [ view(AR, i, 1:size(AR,2)) for i=1:size(AR,1) ]
  permAR = sortperm(rows, order=Base.Lexicographic)
  invpermAR = invperm(permAR)
  permL[:] = permA[invpermAR]
  if A[permL,:] == AR
    return true
  else
    return false
  end
end

"""
### SymCubatures.getpermutation

Returns a permutation of the cubature nodes based on a reordering of the
vertices from their canonical orientation.  That is, finds the node ordering
such that the new nodes have the same order as the old nodes, but relative to
`vtxperm`.  This is useful for face-based operations.

**Inputs**

* `cub`: a symmetric cubature rule for which the permutation is sought
* `vtxperm`: the permutation that is applied to the vertices

**Outputs**

* `perm`: permutation of the cubature nodes

"""
function getpermutation{T}(cub::PointSymCub{T}, vtxperm::Array{Int,1})
  @assert( length(vtxperm) == 1 )
  perm = zeros(Int, (cub.numnodes))
  perm[1] = 1
  return perm
end

function getpermutation{T}(cub::LineSymCub{T}, vtxperm::Array{Int,1})
  @assert( length(vtxperm) == 2 )
  perm = zeros(Int, (cub.numnodes))
  ptr = 0
  # set permutation for nodes with 2-symmetries
  # set vertices
  if cub.vertices
    perm[ptr+1:ptr+2] = invperm(vtxperm) + ptr
    ptr += 2
  end
  # set edge nodes
  A = T[1 -1; -1 1]
  permL = zeros(Int, 2)
  @assert( findleftperm!(A, vtxperm, permL) )
  for i = 1:cub.numedge
    perm[ptr+1:ptr+2] = permL + ptr
    ptr += 2
  end
  # set permutation for node with 1-symmetry
  if cub.centroid
    perm[ptr+1] = ptr+1
    ptr += 1
  end
  @assert( ptr == cub.numnodes )
  return perm  
end

function getpermutation{T}(cub::TriSymCub{T}, vtxperm::AbstractArray{Int,1})
  @assert( length(vtxperm) == 3 )
  perm = zeros(Int, (cub.numnodes))
  ptr = 0
  # set permutation for nodes with 3-symmetries
  # set vertices
  if cub.vertices
    perm[ptr+1:ptr+3] = invperm(vtxperm) + ptr
    ptr += 3
  end
  # mid-edge nodes
  A = T[0.5 0.5 0; 0 0.5 0.5; 0.5 0 0.5]
  permL = zeros(Int, 3)
  @assert( findleftperm!(A, vtxperm, permL) )
  if cub.midedges
    perm[ptr+1:ptr+3] = permL + ptr
    ptr += 3
  end
  # set S21 orbit nodes
  A = T[0.5 0.5 -1; -1 0.5 0.5; 0.5 -1 0.5]
  permL = zeros(Int, 3)
  @assert( findleftperm!(A, vtxperm, permL) )
  for i = 1:cub.numS21
    perm[ptr+1:ptr+3] = permL + ptr
    ptr += 3
  end
  # set permutation for nodes with 6-symmetries
  # set edge nodes
  A = T[1 -1 0; -1 1 0; 0 1 -1; 0 -1 1; -1 0 1; 1 0 -1]
  permL = zeros(Int, 6)
  @assert( findleftperm!(A, vtxperm, permL) )
  for i = 1:cub.numedge
    perm[ptr+1:ptr+6] = permL + ptr
    ptr += 6
  end
  # set S111 orbit nodes
  alpha = 0.1
  beta = 0.4
  A = T[alpha beta (1-alpha-beta);
        beta alpha (1-alpha-beta);
        (1-alpha-beta) alpha beta;
        (1-alpha-beta) beta alpha;
        beta (1-alpha-beta) alpha;
        alpha (1-alpha-beta) beta]
  permL = zeros(Int, 6)
  @assert( findleftperm!(A, vtxperm, permL) )
  for i = 1:cub.numS111
    perm[ptr+1:ptr+6] = permL + ptr
    ptr += 6
  end
  # set permutation for node with 1-symmetry
  if cub.centroid
    perm[ptr+1:ptr+1] = ptr+1
    ptr += 1
  end
  @assert( ptr == cub.numnodes )  
  return perm
end

function getpermutation{T}(cub::TetSymCub{T}, vtxperm::AbstractArray{Int,1})
  @assert( length(vtxperm) == 4 )
  perm = zeros(Int, (cub.numnodes))
  ptr = 0
  # set permutation for nodes with 4-symmetries
  # set vertices
  if cub.vertices
    perm[ptr+1:ptr+4] = invperm(vtxperm) + ptr
    ptr += 4
  end
  # set face centroid
  A = T[1 1 1 0; 1 1 0 1; 0 1 1 1; 1 0 1 1]
  permL = zeros(Int, 4)
  @assert( findleftperm!(A, vtxperm, permL) )
  if cub.facecentroid
    perm[ptr+1:ptr+4] = permL + ptr
    ptr += 4
  end
  # set S31 orbits
  A = T[1 1 1 -3; 1 1 -3 1; -3 1 1 1; 1 -3 1 1]
  permL = zeros(Int, 4)
  @assert( findleftperm!(A, vtxperm, permL) )
  for i = 1:cub.numS31
    perm[ptr+1:ptr+4] = permL + ptr
    ptr += 4
  end
  # set permutation for nodes with 6-symmetries
  # set mid-edge nodes
  A = T[0.5 0.5 0 0;
        0 0.5 0.5 0;
        0.5 0 0.5 0;
        0.5 0 0 0.5;
        0 0.5 0 0.5;
        0 0 0.5 0.5]
  permL = zeros(Int, 6)
  @assert( findleftperm!(A, vtxperm, permL) )
  if cub.midedges
    perm[ptr+1:ptr+6] = permL + ptr
    ptr += 6
  end
  # set S22 oribt nodes
  A = T[1 1 -1 -1; 1 -1 1 -1; 1 -1 -1 1; -1 1 1 -1; -1 1 -1 1; -1 -1 1 1]
  permL = zeros(Int, 6)
  findleftperm!(A, vtxperm, permL)
  for i = 1:cub.numS22
    perm[ptr+1:ptr+6] = permL + ptr
    ptr += 6
  end
  # set all nodes with 12-symmetries
  # set edge nodes
  A = T[1 -1 0 0; # edge 1
        -1 1 0 0; # edge 1
        0 1 -1 0; # edge 2
        0 -1 1 0; # edge 2
        -1 0 1 0; # edge 3
        1 0 -1 0; # edge 3
        1 0 0 -1; # edge 4
        -1 0 0 1; # edge 4
        0 1 0 -1; # edge 5
        0 -1 0 1; # edge 5
        0 0 1 -1; # edge 6
        0 0 -1 1] # edge 6
  permL = zeros(Int, 12)
  findleftperm!(A, vtxperm, permL)
  for i = 1:cub.numedge
    perm[ptr+1:ptr+12] = permL + ptr
    ptr += 12
  end
  # set face nodes corresponding to S21 orbit
  A = T[0.5 0.5 -1 0; # face 1
        -1 0.5 0.5 0;
        0.5 -1 0.5 0;
        0.5 -1 0 0.5; # face 2
        -1 0.5 0 0.5;
        0.5 0.5 0 -1;
        0 0.5 -1 0.5; # face 3
        0 -1 0.5 0.5;
        0 0.5 0.5 -1;
        0.5 0 0.5 -1; # face 4
        -1 0 0.5 0.5;
        0.5 0 -1 0.5]
  permL = zeros(Int, 12)
  findleftperm!(A, vtxperm, permL)
  for i = 1:cub.numfaceS21
    perm[ptr+1:ptr+12] = permL + ptr
    ptr += 12
  end
  # set S211 oribt nodes
  alpha = 0.1
  beta = 0.3
  Aface = [alpha alpha (1-2*alpha-beta) beta;
           (1-2*alpha-beta) alpha alpha beta;
           alpha (1-2*alpha-beta) alpha beta]
  facevtx = [1 1 2 1;
             2 4 4 3;
             3 2 3 4;
             4 3 1 2]
  A = zeros(12, 4)
  for face = 1:4
    A[3*(face-1)+1:3*face,:] = Aface[:,invperm(facevtx[:,face])]
  end
  permL = zeros(Int, 12)
  findleftperm!(A, vtxperm, permL)
  for i = 1:cub.numS211
    perm[ptr+1:ptr+12] = permL + ptr
    ptr += 12
  end
  # set all nodes with 24-symmetries
  # set face nodes corresponding to S111 orbit
  alpha = 0.1
  beta = 0.3
  gamma = 0.25 # these have gamma = 0, but for permL this is not important
  Aface = [alpha beta (1-alpha-beta-gamma) gamma;
           beta alpha (1-alpha-beta-gamma) gamma;
           (1-alpha-beta-gamma) alpha beta gamma;
           (1-alpha-beta-gamma) beta alpha gamma;
           beta (1-alpha-beta-gamma) alpha gamma;
           alpha (1-alpha-beta-gamma) beta gamma]  
  A = zeros(24, 4)
  for face = 1:4
    A[6*(face-1)+1:6*face,:] = Aface[:,invperm(facevtx[:,face])]
  end
  permL = zeros(Int, 24)
  findleftperm!(A, vtxperm, permL)
  for i = 1:cub.numfaceS111
    perm[ptr+1:ptr+24] = permL + ptr
    ptr += 24
  end
  # set S1111 orbit nodes
  # ...these use the same permL as the face S111 orbits above
  for i = 1:cub.numS1111
    perm[ptr+1:ptr+24] = permL + ptr
    ptr += 24
  end
  # set node with 1 symmetry (i.e. the centroid)
  if cub.centroid
    perm[ptr+1:ptr+1] = ptr+1
    ptr += 1
  end
  @assert( ptr == cub.numnodes )
  return perm
end

"""
### SymCubatures.getfacebasedpermutation

Returns a permutation of the volume nodes (or a subset of them) for each face,
such that the same face operator can be applied to all faces.  This is useful
for volume-to-face interpolation or differentiation.

**Inputs**

* `cub`: a symmetric cubature rule for which a face-based permutation is sought
* `faceonly`: if true, only face nodes are used in the permutation.

**Outputs**

* `perm`: permutation of the volume nodes for each face

"""
function getfacebasedpermutation{T}(cub::PointSymCub{T}; faceonly::Bool=false)
  perm = zeros(Int, (cub.numnodes, 1))
  perm[1,1] = 1
  return perm
end

function getfacebasedpermutation{T}(cub::LineSymCub{T}; faceonly::Bool=false)
  if faceonly
    @assert(cub.vertices) # vertices must be active
    perm = zeros(Int, (getnumfacenodes(cub), 2))
    perm = getfacenodeindices(cub)
  else
    perm = zeros(Int, (cub.numnodes, 2))
    perm[:,1] = [1:cub.numnodes;]
    perm[:,2] = getpermutation(cub, [2;1])
  end
  return perm
end

function getfacebasedpermutation{T}(cub::TriSymCub{T}; faceonly::Bool=false)
  if faceonly
    perm = zeros(Int, (getnumfacenodes(cub), 3))
    perm = getfacenodeindices(cub)
  else
    perm = zeros(Int, (cub.numnodes, 3))
    perm[:,1] = [1:cub.numnodes;] # no permutation on face 1
    perm[:,2] = getpermutation(cub, invperm([2;3;1]))
    perm[:,3] = getpermutation(cub, invperm([3;1;2]))
  end
  return perm
end

function getfacebasedpermutation{T}(cub::TetSymCub{T}; faceonly::Bool=false)
  if faceonly
    perm = zeros(Int, (getnumfacenodes(cub), 4))
    perm = getfacenodeindices(cub)
  else
    perm = zeros(Int, (cub.numnodes, 4))
    perm[:,1] = [1:cub.numnodes;] # no permutation on face 1
    perm[:,2] = getpermutation(cub, invperm([1;4;2;3]))
    perm[:,3] = getpermutation(cub, invperm([2;4;3;1]))
    perm[:,4] = getpermutation(cub, invperm([1;3;4;2]))
  end
  return perm
end

"""
### SymCubatures.getneighbourpermutation

At element interfaces, the cubature nodes of the common face will not match when
natural ordering is provided.  This routine produces the permutation that makes
the 'right' element's nodes match the 'left' element's nodes.  The permutation
depends on the face dimension:

* For line segment faces, i.e. points, a trival permutation is returned
* For triangle faces, i.e. line segments, there is only one possible orientation.
* For tetrahedral faces, i.e. triangles, there are three orientations:

1. The \"1\" vertex from each element's face is coincident;
2. The \"1\" vertex from face 1 coincides with \"2\" vertex from face 2;
3. The \"1\" vertex from face 1 coincides with \"3\" vertex from face 2.

**Inputs**

* `cub`: a symmetric cubature rule for which the permutation is sought

**Outputs**

* `perm`: permutation of the interface nodes for each possible orientation

"""
function getneighbourpermutation{T}(cub::PointSymCub{T})
  perm = zeros(Int, (1,1))
  perm[1,1] = 1
  return perm
end

function getneighbourpermutation{T}(cub::LineSymCub{T})
  perm = zeros(Int, (cub.numnodes, 1))
  perm[:,1] = getpermutation(cub, [2;1])
  return perm
end

function getneighbourpermutation{T}(cub::TriSymCub{T})
  perm = zeros(Int, (cub.numnodes, 3))
  perm[:,1] = getpermutation(cub, invperm([1;3;2]))
  perm[:,2] = getpermutation(cub, invperm([2;1;3]))
  perm[:,3] = getpermutation(cub, invperm([3;2;1]))
  return perm
end

""" 
### SymCubatures.getnodevalences

Returns 1 for interior, 2 for faces/edges, 6 for vertices/edges

**Inputs**: 

* `cub`: symmetric cubature rule

**Returns**:
 
* `val`: appropriate valences

This function was designed for the buildminfrobeniusoperator().

"""
function getnodevalences(cub::TriSymCub{T}) where {T}
  val = zeros(T, (cub.numnodes))
  ptr = 0
  if cub.vertices
    val[1:3] = 6.0
    ptr = 3
  end
  # set mid-edge nodes
  if cub.midedges
    val[ptr+1:ptr+3] = 2.0
    ptr += 3
  end
  # set S21 orbit nodes
  for i = 1:cub.numS21
    val[ptr+1:ptr+3] = 1.0
    ptr += 3
  end
  # set all nodes with 6-symmetries
  # set edge nodes
  for i = 1:cub.numedge
    val[ptr+1:ptr+6] = 2.0
    ptr += 6
  end
  # set S111 orbit nodes
  for i = 1:cub.numS111
    val[ptr+1:ptr+6] = 1.0
    ptr += 6
  end
  # set node with 1-symmetry (i.e. centroid)
  if cub.centroid
    val[ptr+1] = 1.0
    ptr += 1
  end
  @assert( ptr == cub.numnodes )
  return val
end


"""
### SymCubatures.setparams!

Sets the nodal parameters for any parameterized symmetry orbits in the cubature.

**Inputs**

* `params`: parameter values

**In/Outs**

* `cub`: symmetric cubature rule whose nodal parameters are being updated

"""
function setparams!{T}(cub::SymCub{T}, params::Array{T})
  @assert( length(params) == cub.numparams )
  cub.params = vec(params)
end

"""
### SymCubatures.setweights!

Sets a cubature's (unique) weights.

**Inputs**

* `weights`: cubature weights grouped by orbit

**In/Outs**

* `cub`: symmetric cubature rule whose weights are being updated

"""
function setweights!{T}(cub::SymCub{T}, weights::Array{T})
  @assert( length(weights) == cub.numweights )
  cub.weights = vec(weights)
end

"""
### SymCubatures.calcnodes

Use the orbital parameter values to compute a cubature's nodal coordinates.  The
second dimension of the `vtx` array of vertices does not need to match the
dimension as the cubature; for example, a line quadrature can be over a line in
1D, 2D, 3D or ND, and a triangle cubature can be in 2D, 3D, or ND, etc.

**Inputs**

* `cub`: symmetric cubature rule
* `vtx`: vertices that define the domain

**Outputs**

* `x`: cubature's nodal coordinates (potentially in a subspace)

"""
function calcnodes{T}(cub::PointSymCub, vtx::Array{T,2})
  @assert(cub.numparams == 0)
  @assert(cub.numnodes == 1)
  @assert(size(vtx,1) == 1)
  x = zeros(T, (size(vtx,2),cub.numnodes))
  x[:,1] = vtx[:,:].'
  return x
end

function calcnodes{T}(cub::LineSymCub, vtx::Array{T,2})
  @assert(cub.numparams >= 0)
  @assert(cub.numnodes >= 1)
  @assert(size(vtx,1) == 2)
  x = zeros(T, (size(vtx,2),cub.numnodes))
  ptr = 0
  paramptr = 0
  # set all nodes with 2-symmetries
  # set vertices
  if cub.vertices
    x[:,1:2] = vtx[:,:].'
    ptr = 2
  end
  # set edge nodes
  for i = 1:cub.numedge
    alpha = cub.params[paramptr+1]
    A = T[alpha (1-alpha);
          (1-alpha) alpha]
    x[:,ptr+1:ptr+2] = (A*vtx[:,:]).'
    ptr += 2
    paramptr += 1
  end
  # set node with 1-symmetry (i.e. centroid)
  if cub.centroid
    A = T[1/2 1/2]
    x[:,ptr+1] = A*vtx[:,:]
    ptr += 1
  end
  return x
end

function calcnodes{T}(cub::TriSymCub, vtx::Array{T,2})
  @assert(cub.numparams >= 0)
  @assert(cub.numnodes >= 1)
  @assert(size(vtx,1) == 3 && size(vtx,2) >= 2)
  x = zeros(T, (size(vtx,2),cub.numnodes))
  ptr = 0
  paramptr = 0
  # set all nodes with 3-symmetries
  # set vertices
  if cub.vertices
    x[:,1:3] = vtx.'
    ptr = 3
  end
  # set mid-edge nodes
  if cub.midedges
    A = T[0.5 0.5 0;
          0 0.5 0.5;
          0.5 0 0.5]
    x[:,ptr+1:ptr+3] = (A*vtx).'
    ptr += 3
  end
  # set S21 orbit nodes
  for i = 1:cub.numS21
    alpha = 0.5*cub.params[paramptr+1]
    A = T[alpha alpha (1-2*alpha);
          (1-2*alpha) alpha alpha;
          alpha (1-2*alpha) alpha]
    x[:,ptr+1:ptr+3] = (A*vtx).'
    ptr += 3
    paramptr += 1
  end
  # set all nodes with 6-symmetries
  # set edge nodes
  for i = 1:cub.numedge
    alpha = cub.params[paramptr+1]
    A = T[alpha (1-alpha) 0;
          (1-alpha) alpha 0;
          0 alpha (1-alpha);
          0 (1-alpha) alpha;
          (1-alpha) 0 alpha;
          alpha 0 (1-alpha)]
    x[:,ptr+1:ptr+6] = (A*vtx).'
    ptr += 6
    paramptr += 1
  end
  # set S111 orbit nodes
  for i = 1:cub.numS111
    alpha = 0.5*cub.params[paramptr+1]
    beta = 0.5*cub.params[paramptr+2]
    A = T[alpha beta (1-alpha-beta);
          beta alpha (1-alpha-beta);
          (1-alpha-beta) alpha beta;
          (1-alpha-beta) beta alpha;
          beta (1-alpha-beta) alpha;
          alpha (1-alpha-beta) beta]
    x[:,ptr+1:ptr+6] = (A*vtx).'
    ptr += 6
    paramptr += 2
  end
  # set node with 1-symmetry (i.e. centroid)
  if cub.centroid
    A = T[1/3 1/3 1/3]
    x[:,ptr+1] = (A*vtx).'
    ptr += 1
  end
  return x
end

function calcnodes{T}(cub::TetSymCub, vtx::Array{T,2})
  @assert(cub.numparams >= 0)
  @assert(cub.numnodes >= 1)
  @assert(size(vtx,1) == 4 && size(vtx,2) >= 3)
  x = zeros(T, (size(vtx,2),cub.numnodes))
  ptr = 0
  paramptr = 0
  # set all nodes with 4-symmetries
  # set vertices
  if cub.vertices
    x[:,ptr+1:ptr+4] = vtx.'
    ptr = 4
  end
  # set face centroids
  if cub.facecentroid
    alpha = (T)(1/3)
    A = T[alpha alpha alpha 0;
          alpha alpha 0 alpha;
          0 alpha alpha alpha;
          alpha 0 alpha alpha]
    x[:,ptr+1:ptr+4] = (A*vtx).'
    ptr += 4
  end
  # set S31 orbit nodes
  for i = 1:cub.numS31
    alpha = cub.params[paramptr+1]/3
    A = T[alpha alpha alpha (1-3*alpha);
          alpha alpha (1-3*alpha) alpha;
          (1-3*alpha) alpha alpha alpha;
          alpha (1-3*alpha) alpha alpha]
    x[:,ptr+1:ptr+4] = (A*vtx).'
    ptr += 4
    paramptr += 1
  end
  # set all nodes with 6-symmetries
  # set mid-edge nodes
  if cub.midedges
    A = T[0.5 0.5 0 0;
          0 0.5 0.5 0;
          0.5 0 0.5 0;
          0.5 0 0 0.5;
          0 0.5 0 0.5;
          0 0 0.5 0.5]
    x[:,ptr+1:ptr+6] = (A*vtx).'
    ptr += 6
  end
  # set S22 oribt nodes
  for i = 1:cub.numS22
    alpha = 0.5*cub.params[paramptr+1]
    A = T[alpha alpha (0.5-alpha) (0.5-alpha);
          alpha (0.5-alpha) alpha (0.5-alpha);
          alpha (0.5-alpha) (0.5-alpha) alpha;
          (0.5-alpha) alpha alpha (0.5-alpha);
          (0.5-alpha) alpha (0.5-alpha) alpha;
          (0.5-alpha) (0.5-alpha) alpha alpha]
    x[:,ptr+1:ptr+6] = (A*vtx).'
    ptr += 6
    paramptr += 1
  end
  # set all nodes with 12-symmetries
  # set edge nodes
  for i = 1:cub.numedge
    alpha = cub.params[paramptr+1]
    A = T[alpha (1-alpha) 0 0; # edge 1
          (1-alpha) alpha 0 0; # edge 1
          0 alpha (1-alpha) 0; # edge 2
          0 (1-alpha) alpha 0; # edge 2
          (1-alpha) 0 alpha 0; # edge 3
          alpha 0 (1-alpha) 0; # edge 3
          alpha 0 0 (1-alpha); # edge 4
          (1-alpha) 0 0 alpha; # edge 4
          0 alpha 0 (1-alpha); # edge 5
          0 (1-alpha) 0 alpha; # edge 5
          0 0 alpha (1-alpha); # edge 6
          0 0 (1-alpha) alpha] # edge 6
    x[:,ptr+1:ptr+12] = (A*vtx).'
    ptr += 12
    paramptr += 1
  end
  # set face nodes corresponding to S21 orbit
  for i = 1:cub.numfaceS21
    alpha = 0.5*cub.params[paramptr+1]
    A = T[alpha alpha (1-2*alpha);
          (1-2*alpha) alpha alpha;
          alpha (1-2*alpha) alpha]
    facevtx = [1 1 2 1;
               2 4 4 3;
               3 2 3 4]
    for face = 1:4
      x[:,ptr+1:ptr+3] = (A*vtx[facevtx[:,face],:]).'
      ptr += 3
    end
    paramptr += 1
  end
  # set S211 orbit nodes
  for i = 1:cub.numS211
    alpha = 0.5*cub.params[paramptr+1]
    beta = 0.5*cub.params[paramptr+2]
    A = T[alpha alpha (1-2*alpha-beta) beta;
          (1-2*alpha-beta) alpha alpha beta;
          alpha (1-2*alpha-beta) alpha beta]
    facevtx = [1 1 2 1;
               2 4 4 3;
               3 2 3 4;
               4 3 1 2]
    for face = 1:4
      x[:,ptr+1:ptr+3] = (A*vtx[facevtx[:,face],:]).'
      ptr += 3
    end
    paramptr += 2
  end
  # set all nodes with 24-symmetries
  # set face nodes corresponding to S111 orbit
  for i = 1:cub.numfaceS111
    alpha = 0.5*cub.params[paramptr+1]
    beta = 0.5*cub.params[paramptr+2]
    A = T[alpha beta (1-alpha-beta);
          beta alpha (1-alpha-beta);
          (1-alpha-beta) alpha beta;
          (1-alpha-beta) beta alpha;
          beta (1-alpha-beta) alpha;
          alpha (1-alpha-beta) beta]
    facevtx = [1 1 2 1;
               2 4 4 3;
               3 2 3 4]
    for face = 1:4
      x[:,ptr+1:ptr+6] = (A*vtx[facevtx[:,face],:]).'
      ptr += 6
    end
    paramptr += 2
  end
  # set S1111 orbit nodes
  for i = 1:cub.numS1111
    alpha = 0.5*cub.params[paramptr+1]
    beta = 0.5*cub.params[paramptr+2]
    gamma = 0.5*cub.params[paramptr+3]
    A = T[alpha beta (1-alpha-beta-gamma) gamma;
          beta alpha (1-alpha-beta-gamma) gamma;
          (1-alpha-beta-gamma) alpha beta gamma;
          (1-alpha-beta-gamma) beta alpha gamma;
          beta (1-alpha-beta-gamma) alpha gamma;
          alpha (1-alpha-beta-gamma) beta gamma]
    facevtx = [1 1 2 1;
               2 4 4 3;
               3 2 3 4;
               4 3 1 2]
    for face = 1:4
      x[:,ptr+1:ptr+6] = (A*vtx[facevtx[:,face],:]).'
      ptr += 6
    end
    paramptr += 3
  end
  # set node with 1 symmetry (i.e. the centroid)
  if cub.centroid
    A = T[0.25 0.25 0.25 0.25]
    x[:,ptr+1] = (A*vtx).'
    ptr += 1
  end
  return x
end

"""
### SymCubatures.calcjacobianofnodes

Returns the Jacobian of the nodes with respect to the orbit parameters.

*Notes*: Jac stores all the x-coordinate Jacobians first, then y (then z)

**Inputs**

* `cub`: symmetric cubature rule
* `vtx`: vertices of the cubature domain

**Outputs**

* `Jac`: Jacobian of the mapping from node parameters to nodes

"""
function calcjacobianofnodes{T}(cub::PointSymCub{T}, vtx::Array{T,2})
  error("SymCubatures.calcjacobianofnodes called with PointSymCub")
end

function calcjacobianofnodes{T}(cub::TriSymCub{T}, vtx::Array{T,2})
  @assert(cub.numparams >= 0)
  @assert(cub.numnodes >= 1)
  @assert( size(vtx,1) == 3 && size(vtx,2) == 2)

  Jac = zeros(T, (2*cub.numnodes, cub.numparams) )
  ptr = 0
  paramptr = 0
  # set Jacobian for all nodes with 3-symmetries
  if cub.vertices
    ptr += 3 # block of zeros, because the vertices are not parameterized
  end
  if cub.midedges
    ptr += 3 # block of zeros, because the midedges are not parameterized
  end
  # set Jacobian for S21 nodes
  A = T[0.5 0.5 -1;
        -1 0.5 0.5;
        0.5 -1 0.5]
  for i = 1:cub.numS21
    for j = 1:2
      Jac[(j-1)*cub.numnodes+ptr+1:(j-1)*cub.numnodes+ptr+3,
          paramptr+1] = A*vtx[:,j]
    end
    ptr += 3
    paramptr += 1
  end  
  # set Jacobian for all nodes with 6-symmetries
  # set Jacobian for edge nodes
  A = T[1 -1 0;
        -1 1 0;
        0 1 -1;
        0 -1 1;
        -1 0 1;
        1 0 -1]
  for i = 1:cub.numedge
    for j = 1:2
      Jac[(j-1)*cub.numnodes+ptr+1:(j-1)*cub.numnodes+ptr+6,
          paramptr+1] = A*vtx[:,j]
    end
    ptr += 6
    paramptr += 1
  end
  # set S111 orbit nodes
  Aalpha = T[0.5 0 -0.5; 0 0.5 -0.5; -0.5 0.5 0;
             -0.5 0 0.5; 0 -0.5 0.5; 0.5 -0.5 0]
  Abeta =  T[0 0.5 -0.5; 0.5 0 -0.5; -0.5 0 0.5;
             -0.5 0.5 0; 0.5 -0.5 0; 0 -0.5 0.5]
  for i = 1:cub.numS111
    for j = 1:2
      Jac[(j-1)*cub.numnodes+ptr+1:(j-1)*cub.numnodes+ptr+6,
          paramptr+1] = Aalpha*vtx[:,j]
      Jac[(j-1)*cub.numnodes+ptr+1:(j-1)*cub.numnodes+ptr+6,
          paramptr+2] = Abeta*vtx[:,j]
    end
    ptr += 6
    paramptr += 2
  end
  # set Jacobian for node with 1-symmetry
  if cub.centroid
    ptr += 1 # block of zeros, because the centroid is not parameterized
  end
  return Jac
end

function calcjacobianofnodes{T}(cub::TetSymCub{T}, vtx::Array{T,2})
  @assert(cub.numparams >= 0)
  @assert(cub.numnodes >= 1)
  @assert( size(vtx,1) == 4 && size(vtx,2) == 3)

  Jac = zeros(T, (3*cub.numnodes, cub.numparams) )
  ptr = 0
  paramptr = 0
  # set Jacobian for all nodes with 4-symmetries
  if cub.vertices
    ptr += 4 # block of zeros, because the vertices are not parameterized
  end
  if cub.facecentroid
    ptr += 4 # block of zeros, because the face centroids are not parameterized
  end
  # set Jacobian for S31 nodes
  A = T[1 1 1 -3;
        1 1 -3 1;
        -3 1 1 1;
        1 -3 1 1]./3
  for i = 1:cub.numS31
    for j = 1:3
      Jac[(j-1)*cub.numnodes+ptr+1:(j-1)*cub.numnodes+ptr+4,
          paramptr+1] = A*vtx[:,j]
    end
    ptr += 4
    paramptr += 1
  end
  # set Jacobian for all nodes with 6-symmetries  
  if cub.midedges
    ptr += 6 # block of zeros, because the midedges are not parameterized
  end
  # set Jacobian for S22 oribt nodes
  A = T[1 1 -1 -1; 1 -1 1 -1; 1 -1 -1 1; -1 1 1 -1; -1 1 -1 1; -1 -1 1 1]./2
  for i = 1:cub.numS22
    for j = 1:3
      Jac[(j-1)*cub.numnodes+ptr+1:(j-1)*cub.numnodes+ptr+6,
          paramptr+1] = A*vtx[:,j]
    end
    ptr += 6
    paramptr += 1
  end
  # set Jacobian for all nodes with 12-symmetries  
  # set Jacobian for edge nodes
  A = T[1 -1 0 0; # edge 1
        -1 1 0 0; # edge 1
        0 1 -1 0; # edge 2
        0 -1 1 0; # edge 2
        -1 0 1 0; # edge 3
        1 0 -1 0; # edge 3
        1 0 0 -1; # edge 4
        -1 0 0 1; # edge 4
        0 1 0 -1; # edge 5
        0 -1 0 1; # edge 5
        0 0 1 -1; # edge 6
        0 0 -1 1] # edge 6
  for i = 1:cub.numedge
    for j = 1:3
      Jac[(j-1)*cub.numnodes+ptr+1:(j-1)*cub.numnodes+ptr+12,
          paramptr+1] = A*vtx[:,j]
    end
    ptr += 12
    paramptr += 1
  end
  # set Jacobian for face nodes corresponding to S21 orbit
  A = T[0.5 0.5 -1;
        -1 0.5 0.5;
        0.5 -1 0.5]
  facevtx = [1 1 2 1;
             2 4 4 3;
             3 2 3 4]
  for i = 1:cub.numfaceS21
    for face = 1:4
      for j = 1:3
        Jac[(j-1)*cub.numnodes+ptr+1:(j-1)*cub.numnodes+ptr+3,
            paramptr+1] = A*vtx[facevtx[:,face],j]
      end
      ptr += 3
    end
    paramptr += 1
  end
  # set Jacobian for S211 oribt nodes
  Aalpha = T[ 0.5  0.5 -1.0 0.0;
             -1.0  0.5  0.5 0.0;
              0.5 -1.0  0.5 0.0]
  Abeta = T[ 0.0  0.0 -0.5 0.5;
            -0.5  0.0  0.0 0.5;
             0.0 -0.5  0.0 0.5]
  facevtx = [1 1 2 1;
             2 4 4 3;
             3 2 3 4;
             4 3 1 2]
  for i = 1:cub.numS211
    for face = 1:4
      for j = 1:3
        Jac[(j-1)*cub.numnodes+ptr+1:(j-1)*cub.numnodes+ptr+3,
            paramptr+1] = Aalpha*vtx[facevtx[:,face],j]
        Jac[(j-1)*cub.numnodes+ptr+1:(j-1)*cub.numnodes+ptr+3,
            paramptr+2] = Abeta*vtx[facevtx[:,face],j]
      end
      ptr += 3
    end
    paramptr += 2
  end
  # set Jacobian for all nodes with 24-symmetries
  # set Jacobian for face nodes corresponding to S111 orbit
  Aalpha = T[0.5 0 -0.5; 0 0.5 -0.5; -0.5 0.5 0;
             -0.5 0 0.5; 0 -0.5 0.5; 0.5 -0.5 0]
  Abeta =  T[0 0.5 -0.5; 0.5 0 -0.5; -0.5 0 0.5;
             -0.5 0.5 0; 0.5 -0.5 0; 0 -0.5 0.5]
  facevtx = [1 1 2 1;
             2 4 4 3;
             3 2 3 4]   
  for i = 1:cub.numfaceS111
    for face = 1:4
      for j = 1:3
        Jac[(j-1)*cub.numnodes+ptr+1:(j-1)*cub.numnodes+ptr+6,
            paramptr+1] = Aalpha*vtx[facevtx[:,face],j]
        Jac[(j-1)*cub.numnodes+ptr+1:(j-1)*cub.numnodes+ptr+6,
            paramptr+2] = Abeta*vtx[facevtx[:,face],j]
      end
      ptr += 6
    end
    paramptr += 2
  end
  # set Jacobian for S1111 node orbits 
  Aalpha = T[0.5 0 -0.5 0; 0 0.5 -0.5 0; -0.5 0.5 0 0;
             -0.5 0 0.5 0; 0 -0.5 0.5 0; 0.5 -0.5 0 0]
  Abeta = T[0 0.5 -0.5 0; 0.5 0 -0.5 0; -0.5 0 0.5 0;
            -0.5 0.5 0 0; 0.5 -0.5 0 0; 0 -0.5 0.5 0]
  Agamma = T[0 0 -0.5 0.5; 0 0 -0.5 0.5; -0.5 0 0 0.5;
             -0.5 0 0 0.5; 0 -0.5 0 0.5; 0 -0.5 0 0.5]
  facevtx = [1 1 2 1;
             2 4 4 3;
             3 2 3 4;
             4 3 1 2]
  for i = 1:cub.numS1111
    for face = 1:4
      for j = 1:3
        Jac[(j-1)*cub.numnodes+ptr+1:(j-1)*cub.numnodes+ptr+6,
            paramptr+1] = Aalpha*vtx[facevtx[:,face],j]
        Jac[(j-1)*cub.numnodes+ptr+1:(j-1)*cub.numnodes+ptr+6,
            paramptr+2] = Abeta*vtx[facevtx[:,face],j]
        Jac[(j-1)*cub.numnodes+ptr+1:(j-1)*cub.numnodes+ptr+6,
            paramptr+3] = Agamma*vtx[facevtx[:,face],j]
      end
      ptr += 6
    end
    paramptr += 3
  end
  # set Jacobian for node with 1-symmetry
  if cub.centroid
    ptr += 1 # block of zeros, because the centroid is not parameterized
  end
  return Jac
end

"""
### SymCubatures.calcweights

Map the unique cubature weights to the weights of all nodes.

**Inputs**

* `cub`: symmetric cubature rule

**Outputs**

* `w`: cubature's weights at all nodes

"""
function calcweights{T}(cub::PointSymCub{T})
  @assert(cub.numweights == 1)
  @assert(cub.numnodes == 1)
  w = zeros(T, (cub.numnodes))
  w[1] = cub.weights[1]
  return w
end

function calcweights{T}(cub::LineSymCub{T})
  @assert(cub.numweights >= 0)
  @assert(cub.numnodes >= 1)
  w = zeros(T, (cub.numnodes))
  ptr = 0
  wptr = 0
  # set weights for all nodes with 2-symmetries
  # set vertices weights
  if cub.vertices
    w[1:2] = cub.weights[wptr+1]
    ptr += 2
    wptr += 1
  end
  # set edge weights
  for i = 1:cub.numedge
    w[ptr+1:ptr+2] = cub.weights[wptr+1]
    ptr += 2
    wptr += 1
  end
  # set weight for node with 1-symmetry (i.e. centroid)
  if cub.centroid
    w[ptr+1:ptr+1] = cub.weights[wptr+1]
    ptr += 1
    wptr += 1
  end
  return w
end

function calcweights{T}(cub::TriSymCub{T})
  @assert(cub.numweights >= 0)
  @assert(cub.numnodes >= 1)

  w = zeros(T, (cub.numnodes))
  ptr = 0
  wptr = 0
  # set weights for all nodes with 3-symmetries
  # set vertices weights
  if cub.vertices
    w[1:3] = cub.weights[wptr+1]
    ptr = 3
    wptr += 1
  end
  # set mid-edge weights
  if cub.midedges
    w[ptr+1:ptr+3] = cub.weights[wptr+1]
    ptr += 3
    wptr += 1
  end
  # set S21 orbit weights
  for i = 1:cub.numS21
    w[ptr+1:ptr+3] = cub.weights[wptr+1]
    ptr += 3
    wptr += 1
  end
  # set weights for all nodes with 6-symmetries
  # set edge weights
  for i = 1:cub.numedge
    w[ptr+1:ptr+6] = cub.weights[wptr+1]
    ptr += 6
    wptr += 1
  end
  # set S111 orbit weights
  for i = 1:cub.numS111
    w[ptr+1:ptr+6] = cub.weights[wptr+1]
    ptr += 6
    wptr += 1
  end
  # set weight for node with 1-symmetry (i.e. centroid)
  if cub.centroid
    w[ptr+1:ptr+1] = cub.weights[wptr+1]
    ptr += 1
    wptr += 1
  end
  return w
end

function calcweights{T}(cub::TetSymCub{T})
  @assert(cub.numweights >= 0)
  @assert(cub.numnodes >= 1)

  w = zeros(T, (cub.numnodes))
  ptr = 0
  wptr = 0
  # set weights for all nodes with 4-symmetries
  # set vertices weights
  if cub.vertices
    w[1:4] = cub.weights[wptr+1]
    ptr = 4
    wptr += 1
  end
  # set face centroid weights
  if cub.facecentroid
    w[ptr+1:ptr+4] = cub.weights[wptr+1]
    ptr += 4
    wptr += 1
  end
  # set S31 orbit weights
  for i = 1:cub.numS31
    w[ptr+1:ptr+4] = cub.weights[wptr+1]
    ptr += 4
    wptr += 1
  end
  # set weights for all nodes with 6-symmetries
  # set mid-edge weights
  if cub.midedges
    w[ptr+1:ptr+6] = cub.weights[wptr+1]
    ptr += 6
    wptr += 1
  end
  # set S22 orbit weights
  for i = 1:cub.numS22
    w[ptr+1:ptr+6] = cub.weights[wptr+1]
    ptr += 6
    wptr += 1
  end
  # set weights for all nodes with 12-symmetries
  # set edge weights
  for i = 1:cub.numedge
    w[ptr+1:ptr+12] = cub.weights[wptr+1]
    ptr += 12
    wptr += 1
  end
  # set face S21 weights
  for i = 1:cub.numfaceS21
    w[ptr+1:ptr+12] = cub.weights[wptr+1]
    ptr += 12
    wptr += 1
  end
  # set S211 orbit weights
  for i = 1:cub.numS211
    w[ptr+1:ptr+12] = cub.weights[wptr+1]
    ptr += 12
    wptr += 1
  end
  # set weights for all nodes with 24-symmetries
  # set face S111 weights
  for i = 1:cub.numfaceS111
    w[ptr+1:ptr+24] = cub.weights[wptr+1]
    ptr += 24
    wptr += 1
  end
  # set S1111 orbit weights
  for i = 1:cub.numS1111
    w[ptr+1:ptr+24] = cub.weights[wptr+1]
    ptr += 24
    wptr += 1
  end
  # set weight for node with 1-symmetry
  # set centroid weight
  if cub.centroid
    w[ptr+1:ptr+1] = cub.weights[wptr+1]
    ptr += 1
    wptr += 1
  end
  return w
end

"""
### SymCubatures.calcjacobianofweights

Returns the Jacobian of the nodal weights with respect to the unique weights.
The resulting Jacobian is a rectangular matrix of ones and zeros that indicates
the mapping from the unique weights to the nodal weights.

**Inputs**

* `cub`: symmetric cubature rule

**Outputs**

* `Jac`: Jacobian of the mapping from (unique) weights to nodal weights

"""
function calcjacobianofweights{T}(cub::PointSymCub{T})
  error("SymCubatures.calcjacobianofweights called for PointSymCub")
end

function calcjacobianofweights{T}(quad::LineSymCub{T})
  @assert(quad.numweights >= 0)
  @assert(quad.numnodes >= 1)

  Jac = zeros(T, (quad.numnodes, quad.numweights) )
  ptr = 0
  wptr = 0
  # set Jacobian for all nodes with 2-symmetries
  # set Jacobian of vertex weights
  if quad.vertices
    Jac[1:2,wptr+1] = ones(T, (2,1))
    ptr = 2
    wptr += 1
  end
  # set Jacobian of edge weights
  for i = 1:quad.numedge
    Jac[ptr+1:ptr+2,wptr+1] = ones(T, (2,1))
    ptr += 2
    wptr += 1
  end
  # set Jacobian for all nodes with 1-symmetries
  if quad.centroid
    Jac[ptr+1,wptr+1] = one(T)
    ptr += 1
    wptr += 1
  end
  return Jac
end

function calcjacobianofweights{T}(cub::TriSymCub{T})
  @assert(cub.numweights >= 0)
  @assert(cub.numnodes >= 1)

  Jac = zeros(T, (cub.numnodes, cub.numweights) )
  ptr = 0
  wptr = 0
  # set Jacobian for all nodes with 3-symmetries
  # set Jacobian of vertex weights
  if cub.vertices
    Jac[1:3,wptr+1] = ones(T, (3,1))
    ptr = 3
    wptr += 1
  end
  # set Jacobian of mid-edge weights
  if cub.midedges
    Jac[ptr+1:ptr+3,wptr+1] = ones(T, (3,1))
    ptr += 3
    wptr += 1
  end
  # set S21 orbit Jacobian weights
  for i = 1:cub.numS21
    Jac[ptr+1:ptr+3,wptr+1] = ones(T, (3,1))
    ptr += 3
    wptr += 1
  end
  # set Jacobian of all nodes with 6-symmetries
  # set Jacobian of edge nodes
  for i = 1:cub.numedge
    Jac[ptr+1:ptr+6,wptr+1] = ones(T, (6,1))
    ptr += 6
    wptr += 1
  end
  # set S111 orbit weights
  for i = 1:cub.numS111
    Jac[ptr+1:ptr+6,wptr+1] = ones(T, (6,1))
    ptr += 6
    wptr += 1
  end
  # set Jacobian of node with 1-symmetry (i.e. centroid)
  if cub.centroid
    Jac[ptr+1:ptr+1,wptr+1] = ones(T, (1,1))
    ptr += 1
    wptr += 1
  end
  return Jac
end

function calcjacobianofweights{T}(cub::TetSymCub{T})
  @assert(cub.numweights >= 0)
  @assert(cub.numnodes >= 1)

  Jac = zeros(T, (cub.numnodes, cub.numweights) )
  ptr = 0
  wptr = 0
  # set Jacobian for all nodes with 4-symmetries
  # set Jacobian of vertex weights
  if cub.vertices
    Jac[1:4,wptr+1] = ones(T, (4,1))
    ptr = 4
    wptr += 1
  end
  # set Jacobian of face centroid weights
  if cub.facecentroid
    Jac[ptr+1:ptr+4,wptr+1] = ones(T, (4,1))
    ptr += 4
    wptr += 1
  end
  # set Jacobian of S31 orbit weights
  for i = 1:cub.numS31
    Jac[ptr+1:ptr+4,wptr+1] = ones(T, (4,1))
    ptr += 4
    wptr += 1
  end
  # set Jacobian for all nodes with 6-symmetries
  # set Jacobian of mid-edge weights
  if cub.midedges
    Jac[ptr+1:ptr+6,wptr+1] = ones(T, (6,1))
    ptr += 6
    wptr += 1
  end
  # set Jacobian of S22 orbit weights
  for i = 1:cub.numS22
    Jac[ptr+1:ptr+6,wptr+1] = ones(T, (6,1))
    ptr += 6
    wptr += 1
  end
  # set Jacobian for all nodes with 12-symmetries
  # set Jacobian of edge nodes
  for i = 1:cub.numedge
    Jac[ptr+1:ptr+12,wptr+1] = ones(T, (12,1))
    ptr += 12
    wptr += 1
  end
  # set Jacobian of face S21 nodes
  for i = 1:cub.numfaceS21
    Jac[ptr+1:ptr+12,wptr+1] = ones(T, (12,1))
    ptr += 12
    wptr += 1
  end
  # set Jacobian of S211 orbit weights
  for i = 1:cub.numS211
    Jac[ptr+1:ptr+12,wptr+1] = ones(T, (12,1))
    ptr += 12
    wptr += 1
  end
  # set Jacobian for all nodes with 24-symmetries
  # set Jacobian of face S111 nodes
  for i = 1:cub.numfaceS111
    Jac[ptr+1:ptr+24,wptr+1] = ones(T, (24,1))
    ptr += 24
    wptr += 1
  end
  # set Jacobian of S1111 nodes
  for i = 1:cub.numS1111
    Jac[ptr+1:ptr+24,wptr+1] = ones(T, (24,1))
    ptr += 24
    wptr += 1
  end
  # set Jacobian for node with 1-symmetry
  # set Jacobian of centroid weight
  if cub.centroid
    Jac[ptr+1:ptr+1,wptr+1] = ones(T, (1,1))
    ptr += 1
    wptr += 1
  end
  return Jac
end

"""
### SymCubatures.calcjacobian

Returns the Jacobian of the nodal coordinates and weights with respect to their
parameters.  In other words, returns the block-rectangular matrix [Jcoords, 0;
0, Jweights], where Jcoords is the Jacobian of the coordinates, Jweights is the
Jacobian of the weights, and the 0s indicate zero blocks of the appropriate
size.

**Inputs**

* `cub`: symmetric cubature rule
* `vtx`: vertices of the cubature domain

**Outputs**

* `Jac`: Jacobian of the nodal coordinates and weights

"""
function calcjacobian{T}(cub::TriSymCub{T},
                         vtx::Array{T,2}=T[-1 -1; 1 -1; -1 1])
  Jac = zeros(T, (3*cub.numnodes, cub.numparams + cub.numweights) )
  Jac[1:2*cub.numnodes, 1:cub.numparams] =
  SymCubatures.calcjacobianofnodes(cub, vtx)
  Jac[2*cub.numnodes+1:end, cub.numparams+1:end] =
  SymCubatures.calcjacobianofweights(cub)
  return Jac
end

function calcjacobian{T}(cub::TetSymCub{T},
                         vtx::Array{T,2}=T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1])
  Jac = zeros(T, (4*cub.numnodes, cub.numparams + cub.numweights) )
  Jac[1:3*cub.numnodes, 1:cub.numparams] =
  SymCubatures.calcjacobianofnodes(cub, vtx)
  Jac[3*cub.numnodes+1:end, cub.numparams+1:end] =
  SymCubatures.calcjacobianofweights(cub)
  return Jac
end

"""
### SymCubatures.getInternalParamMask

Returns the set of parameter indices corresponding to internal nodes; this is
useful when finding cubature rules for which we wish to fix the boundary nodes
and only allow the internal nodes to move.

**Inputs**

* `cub`: symmetric cubature rule

**Returns**

* `mask`: integer array of parameter indices associated with internal nodes

"""
function getInternalParamMask{T}(cub::TriSymCub{T})
  mask = Array{Int64}(0)
  paramptr = 1
  # include S21 orbit parameters
  for i = 1:cub.numS21
    push!(mask, paramptr)
    paramptr += 1
  end
  # account for edge node parameters (not included in mask)
  paramptr += cub.numedge
  # include S111 orbit nodes
  for i = 1:cub.numS111
    push!(mask, paramptr)
    paramptr += 1
    push!(mask, paramptr)
    paramptr += 1
  end
  return mask
end

function getInternalParamMask{T}(cub::TetSymCub{T})
  mask = Array{Int64}(0)
  paramptr = 1
  # include S31 orbit parameters
  for i = 1:cub.numS31
    push!(mask, paramptr)
    paramptr += 1
  end
  # include S22 orbit parameters
  for i = 1:cub.numS22
    push!(mask, paramptr)
    paramptr += 1
  end
  # account for edge node parameters (not included in mask)
  paramptr += cub.numedge
  # account for face S21 orbits (not included in mask)
  paramptr += cub.numfaceS21
  # include S211 orbit parameters
  for i = 1:cub.numS211
    push!(mask, paramptr)
    paramptr += 1
    push!(mask, paramptr)
    paramptr += 1
  end
  # account for face S111 orbits (not included in mask)
  paramptr += 2*cub.numfaceS111
  # include S1111 orbit parameters
  for i = 1:cub.numS1111
    push!(mask, paramptr)
    paramptr += 1
    push!(mask, paramptr)
    paramptr += 1
    push!(mask, paramptr)
    paramptr += 1
  end
  return mask
end
                                 
end
