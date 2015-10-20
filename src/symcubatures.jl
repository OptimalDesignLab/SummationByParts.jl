module SymCubatures
# types and methods for mapping between symmetry groups and nodes for cubatures
# on various domains

export SymCub, LineSymCub, TriSymCub, TetSymCub

@doc """
### SymCubatures.SymCub

`SymCub` is an parametric abstract type that defines cubatures for symmetric
nodal distributions.  It is parameterized on `T` in order to allow for future
implementations of arbitrary precision types.  The parameterization also permits
the use of the complex-step method for verification.

"""-> abstract SymCub{T<:Number}

@doc """
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

"""->
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

  function LineSymCub(;numedge::Int=0, vertices::Bool=true, centroid::Bool=false)
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

@doc """
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

"""->
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

  function TriSymCub(;numedge::Int=0, numS21::Int=0, numS111::Int=0,
                     vertices::Bool=true, midedges::Bool=false,
                     centroid::Bool=false)
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

@doc """
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
* `numS31` : number of S31 orbits (vertex to opposite face)
* `numS22` : number of S22 orbits
* `numsym` : number of node sets in each symmetry group (5 groups for Tet)
* `params` : the actual values of the orbit nodal parameters
* `weights` : values of the unique weights

"""->
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
  numS31::Int
  numS22::Int
  numsym::Array{Int,1}
  params::Array{T,1}
  weights::Array{T,1}

  function TetSymCub(;numedge::Int=0, numfaceS21::Int=0, numS31::Int=0,
                     numS22::Int=0, vertices::Bool=true, midedges::Bool=false,
                     centroid::Bool=false, facecentroid::Bool=false)
    @assert(numedge >= 0)
    @assert(numfaceS21 >= 0)
    @assert(numS31 >= 0)
    @assert(numS22 >= 0)
    # compute the number of degrees of freedom and unique weights
    numparams = 0
    numweights = 0
    numparams += numedge + numfaceS21 + numS31 + numS22
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
    numnodes += 12*numedge # 2 * (6 edges) X (number of edge parameters)
    numsym[4] += numedge
    numweights += numedge
    numnodes += 12*numfaceS21 # 3 * (4 faces) X (number of S21 orbits)
    numsym[4] += numfaceS21
    numweights += numfaceS21
    numnodes += 4*numS31 # (4 permutations) X (number of S31 orbits)
    numsym[2] += numS31
    numweights += numS31
    numnodes += 6*numS22 # (6 permutations) X (number of S22 orbits)
    numsym[3] += numS22
    numweights += numS22
    # initialize parameter arrays
    @assert(numweights == sum(numsym))
    params = zeros(T, numparams)
    weights = zeros(T, numweights)
    new(numparams, numweights, numnodes, vertices, midedges, centroid,
        facecentroid, numedge, numfaceS21, numS31, numS22, numsym, params,
        weights)
  end
end

@doc """
### SymCubatures.getnumboundarynodes

Returns the number of (explicit) boundary nodes

*Notes*: if the parameter value for an internal orbit is such that the
corresponding node lies on the boundary, this node is **NOT** included in the
boundary-node count returned.

**Inputs**

* `cub`: symmetric cubature rule

**Outputs**

* `numboundary`: number of boundary nodes

"""->
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
  return numboundary
end

@doc """
### SymCubatures.getnumfacenodes

Returns the number of nodes on an individual face of the element.

**Inputs**

* `cub`: symmetric cubature rule

**Outputs**

* `numfacenodes`: number of nodes on a face

"""->
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
  return numfacenodes
end

@doc """
### SymCubatures.getbndrynodeindices

Returns the indices of the nodes that on the boundary, in their natural order.
See getfacenodeindices for a method returns node indices for each face.

**Inputs**

* `cub`: a symmetric cubature rule whose boundary-node indices are sought

**Outputs**

* `bndryindices`: indicies of nodes that lie on boundary

"""->
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
  # more to come...

  return bndryindices
end

@doc """
### SymCubatures.getinteriornodeindices

Returns the indices of the nodes that are strictly interior.

**Inputs**

* `cub`: a symmetric cubature rule whose interior-node indices are sought

**Outputs**

* `indices`: indicies of nodes that are strictly interior.

"""->
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
  # ...
  if cub.centroid
    indices[idxptr+1] = ptr+1
    ptr += 1
    idxptr += 1
  end
  return indices
end

@doc """
### SymCubatures.getfacenodeindices

Returns the indices of the nodes that lie on each face.  See getbndrynodeindices
for a method that returns a single array of boundary nodes.

**Inputs**

* `cub`: a symmetric cubature rule whose boundary-node indices are sought

**Outputs**

* `bndryindices`: indicies of nodes that lie on boundary; there is a separate
  column of indices for each edge/face.

"""->
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
    idxptr += 1
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
    bndryindices[idxptr+1:idxptr+3,:] = [ptr+1 ptr+2 ptr+3 ptr+4;
                                         ptr+3 ptr+3 ptr+1 ptr+1;
                                         ptr+2 ptr+4 ptr+4 ptr+2] 
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
    bndryindices[idxptr+1:idxptr+3,:] = [ptr+5 ptr+2 ptr+5 ptr+4;
                                         ptr+2 ptr+3 ptr+4 ptr+1;
                                         ptr+1 ptr+6 ptr+3 ptr+6]
    ptr += 6
    idxptr += 3
  end
  # account for S22 orbits
  ptr += 6*cub.numS22
  # add edge nodes to indices
  for i = 1:cub.numedge
    bndryindices[idxptr+1:idxptr+6,:] = [ptr+9  ptr+3  ptr+10 ptr+8;
                                         ptr+10 ptr+4  ptr+9  ptr+7;
                                         ptr+4  ptr+5  ptr+7  ptr+1;
                                         ptr+3  ptr+6  ptr+8  ptr+2;
                                         ptr+2  ptr+12 ptr+6  ptr+11;
                                         ptr+1  ptr+11 ptr+5  ptr+12]
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
  return bndryindices
end

@doc """
### SymCubatures.getfacebasedpermutation

Returns a permutation of the volume nodes (or a subset of them) for each face,
such that the same face operator can be applied to all faces.  This is useful
for volume-to-face interpolation or differentiation.

**Inputs**

* `cub`: a symmetric cubature rule for which a face-based permutation is sought
* `faceonly`: if true, only face nodes are used in the permutation.

**Outputs**

* `perm`: permutation of the volume nodes for each face

"""->
function getfacebasedpermutation{T}(cub::LineSymCub{T}; faceonly::Bool=false)
  if faceonly
    @assert(cub.vertices) # vertices must be active
    perm = zeros(Int, (getnumfacenodes(cub), 2))
    perm = getfacenodeindices(cub)
  else
    perm = zeros(Int, (cub.numnodes, 2))
    perm[:,1] = [1:cub.numnodes;]
    ptr = 0
    # set permutation for nodes with 2-symmetries
    # set vertices
    if cub.vertices
      perm[ptr+1:ptr+2,2] = [2; 1]
      ptr += 2
    end
    # set edge nodes
    for i = 1:cub.numedge
      perm[ptr+1:ptr+2,2] = [ptr+2; ptr+1]
      ptr += 2
    end
    # set permutation for node with 1-symmetry
    if cub.centroid
      perm[ptr+1,2] = ptr+1
      ptr =+ 1
    end
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
    ptr = 0
    # set permutation for nodes with 3-symmetries
    # set vertices
    if cub.vertices
      perm[ptr+1:ptr+3,2] = [ptr+2; ptr+3; ptr+1]
      perm[ptr+1:ptr+3,3] = [ptr+3; ptr+1; ptr+2]
      ptr += 3
    end
    # mid-edge nodes
    if cub.midedges
      perm[ptr+1:ptr+3,2] = [ptr+2; ptr+3; ptr+1]
      perm[ptr+1:ptr+3,3] = [ptr+3; ptr+1; ptr+2]
      ptr += 3
    end
    # set S21 orbit nodes
    for i = 1:cub.numS21
      perm[ptr+1:ptr+3,2] = [ptr+2; ptr+3; ptr+1]
      perm[ptr+1:ptr+3,3] = [ptr+3; ptr+1; ptr+2]
      ptr += 3
    end
    # set permutation for nodes with 6-symmetries
    # set edge nodes
    for i = 1:cub.numedge
      perm[ptr+1:ptr+6,2] = [ptr+3; ptr+4; ptr+5; ptr+6; ptr+1; ptr+2]
      perm[ptr+1:ptr+6,3] = [ptr+5; ptr+6; ptr+1; ptr+2; ptr+3; ptr+4]
      ptr += 6
    end
    # set S111 orbit nodes
    for i = 1:cub.numS111
      perm[ptr+1:ptr+6,2] = [ptr+3; ptr+4; ptr+5; ptr+6; ptr+1; ptr+2]
      perm[ptr+1:ptr+6,3] = [ptr+5; ptr+6; ptr+1; ptr+2; ptr+3; ptr+4]
      ptr += 6
    end
    # set permutation for node with 1-symmetry
    if cub.centroid
      perm[ptr+1,2] = ptr+1
      perm[ptr+1,3] = ptr+1
      ptr += 1
    end
  end
  return perm
end

function getfacebasedpermutation{T}(cub::TetSymCub{T}; faceonly::Bool=false)
  if faceonly
    perm = zeros(Int, (getnumfacenodes(cub), 3))
    perm = getfacenodeindices(cub)
  else
    perm = zeros(Int, (cub.numnodes, 4))
    perm[:,1] = [1:cub.numnodes;] # no permutation on face 1
    ptr = 0
    # set permutation for nodes with 4-symmetries
    # set vertices
    if cub.vertices
      perm[ptr+1:ptr+4,2] = [ptr+2; ptr+3; ptr+1; ptr+4]
      perm[ptr+1:ptr+4,3] = [ptr+3; ptr+1; ptr+2; ptr+4]
      perm[ptr+1:ptr+4,4] = [ptr+4; ptr+2; ptr+1; ptr+3]
      ptr += 4
    end
    if cub.facecentroid
      perm[ptr+1:ptr+4,2] = [ptr+2; ptr+3; ptr+1; ptr+4]
      perm[ptr+1:ptr+4,3] = [ptr+3; ptr+1; ptr+2; ptr+4]
      perm[ptr+1:ptr+4,4] = [ptr+4; ptr+2; ptr+1; ptr+3]
      ptr += 4
    end
    for i = 1:cub.numS31
      perm[ptr+1:ptr+4,2] = [ptr+2; ptr+3; ptr+1; ptr+4]
      perm[ptr+1:ptr+4,3] = [ptr+3; ptr+1; ptr+2; ptr+4]
      perm[ptr+1:ptr+4,4] = [ptr+4; ptr+2; ptr+1; ptr+3]
      ptr += 4
    end
    # set permutation for nodes with 6-symmetries
    # set mid-edge nodes
    if cub.midedges
      # ...
    end
    # To be continued...
  end
  return perm
end

@doc """
### SymCubatures.setparams!

Sets the nodal parameters for any parameterized symmetry orbits in the cubature.

**Inputs**

* `params`: parameter values

**In/Outs**

* `cub`: symmetric cubature rule whose nodal parameters are being updated

"""->
function setparams!{T}(cub::SymCub{T}, params::Array{T})
  @assert( length(params) == cub.numparams )
  cub.params = vec(params)
end

@doc """
### SymCubatures.setweights!

Sets a cubature's (unique) weights.

**Inputs**

* `weights`: cubature weights grouped by orbit

**In/Outs**

* `cub`: symmetric cubature rule whose weights are being updated

"""->
function setweights!{T}(cub::SymCub{T}, weights::Array{T})
  @assert( length(weights) == cub.numweights )
  cub.weights = vec(weights)
end

@doc """
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

"""->
function calcnodes{T}(cub::LineSymCub{T}, vtx::Array{T,2})
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

function calcnodes{T}(cub::TriSymCub{T}, vtx::Array{T,2})
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

function calcnodes{T}(cub::TetSymCub{T}, vtx::Array{T,2})
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
          0 alpha alpha alpha;
          alpha 0 alpha alpha;
          alpha alpha 0 alpha]
    x[:,ptr+1:ptr+4] = (A*vtx).'
    ptr += 4
  end
  # set S31 orbit nodes
  for i = 1:cub.numS31
    alpha = cub.params[paramptr+1]/3
    A = T[alpha alpha alpha (1-3*alpha);
          (1-3*alpha) alpha alpha alpha;
          alpha (1-3*alpha) alpha alpha;
          alpha alpha (1-3*alpha) alpha]
    x[:,ptr+1:ptr+4] = (A*vtx).'
    ptr += 4
    paramptr += 1
  end
  # set all nodes with 6-symmetries
  # set mid-edge nodes
  if cub.midedges
    A = T[0.5 0.5 0 0;
          0 0.5 0.5 0;
          0 0 0.5 0.5;
          0.5 0 0 0.5;
          0.5 0 0.5 0;
          0 0.5 0 0.5]
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
    A = T[alpha (1-alpha) 0 0;
          (1-alpha) alpha 0 0;
          0 alpha (1-alpha) 0;
          0 (1-alpha) alpha 0;
          0 0 alpha (1-alpha);
          0 0 (1-alpha) alpha;
          alpha 0 0 (1-alpha);
          (1-alpha) 0 0 alpha;
          alpha 0 (1-alpha) 0;
          (1-alpha) 0 alpha 0;
          0 alpha 0 (1-alpha);
          0 (1-alpha) 0 alpha]
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
    facevtx = [1 2 3 4;
               3 3 1 1;
               2 4 4 2]
    for face = 1:4
      x[:,ptr+1:ptr+3] = (A*vtx[facevtx[:,face],:]).'
      ptr += 3
    end
    paramptr += 1
  end
  # set all nodes with 24-symmetries
  # ... face nodes with 6 nodes and volume nodes

  # set node with 1 symmetry (i.e. the centroid)
  if cub.centroid
    A = T[0.25 0.25 0.25 0.25]
    x[:,ptr+1] = (A*vtx).'
    ptr += 1
  end

  return x
end

@doc """
### SymCubatures.calcjacobianofnodes

Returns the Jacobian of the nodes with respect to the orbit parameters.

*Notes*: Jac stores all the x-coordinate Jacobians first, then y (then z)

**Inputs**

* `cub`: symmetric cubature rule
* `vtx`: vertices of the cubature domain

**Outputs**

* `Jac`: Jacobian of the mapping from node parameters to nodes

"""->
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
        -3 1 1 1;
        1 -3 1 1;
        1 1 -3 1]./3
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
  A = T[1 -1 0 0;
        -1 1 0 0;
        0 1 -1 0;
        0 -1 1 0;
        0 0 1 -1;
        0 0 -1 1;
        1 0 0 -1;
        -1 0 0 1;
        1 0 -1 0;
        -1 0 1 0;
        0 1 0 -1;
        0 -1 0 1]
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
  facevtx = [1 2 3 4;
             3 3 1 1;
             2 4 4 2]
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
  # set Jacobian for all nodes with 24-symmetries
  # ...

  # set Jacobian for node with 1-symmetry
  if cub.centroid
    ptr += 1 # block of zeros, because the centroid is not parameterized
  end
  return Jac
end

@doc """
### SymCubatures.calcweights

Map the unique cubature weights to the weights of all nodes.

**Inputs**

* `cub`: symmetric cubature rule

**Outputs**

* `w`: cubature's weights at all nodes

"""->
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
  # set weights for all nodes with 24-symmetries
  # ...

  # set weight for node with 1-symmetry
  # set centroid weight
  if cub.centroid
    w[ptr+1:ptr+1] = cub.weights[wptr+1]
    ptr += 1
    wptr += 1
  end
  return w
end

@doc """
### SymCubatures.calcjacobianofweights

Returns the Jacobian of the nodal weights with respect to the unique weights.
The resulting Jacobian is a rectangular matrix of ones and zeros that indicates
the mapping from the unique weights to the nodal weights.

**Inputs**

* `cub`: symmetric cubature rule

**Outputs**

* `Jac`: Jacobian of the mapping from (unique) weights to nodal weights

"""->
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
  # set Jacobian for all nodes with 24-symmetries
  # ...
  # set Jacobian for node with 1-symmetry
  # set Jacobian of centroid weight
  if cub.centroid
    Jac[ptr+1:ptr+1,wptr+1] = ones(T, (1,1))
    ptr += 1
    wptr += 1
  end
  return Jac
end

@doc """
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

"""->
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

end
