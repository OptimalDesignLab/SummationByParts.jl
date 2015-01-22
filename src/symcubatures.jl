module SymCubatures
# types and methods for mapping between symmetry groups and nodes for cubatures
# on various domains

export SymCub, TriSymCub, TetSymCub

@doc """
### SymCubatures.SymCub

`SymCub` is an parametric abstract type that defines cubatures for symmetric
nodal distributions.  It is parameterized on `T` in order to allow for future
implementations of arbitrary precision types.  The parameterization also permits
the use of the complex-step method for verification.

"""-> abstract SymCub{T<:Number}

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
  params::Array{T}
  weights::Array{T}

  function TriSymCub(;numedge::Int=0, numS21::Int=0, vertices::Bool=true,
                     midedges::Bool=false, centroid::Bool=false)
    @assert(numedge >= 0)
    @assert(numS21 >= 0)
    # compute the number of degrees of freedom and unique weights
    numparams = 0
    numweights = 0
    numparams += numedge + numS21
    # compute the number of nodes
    numnodes = 0
    if vertices
      numnodes += 3
      numweights += 1
    end
    if midedges
      numnodes += 3
      numweights += 1
    end
    if centroid
      numnodes += 1
      numweights += 1
    end
    numnodes += 6*numedge # 2 * (3 edges) X (number of edge parameters)
    numweights += numedge
    numnodes += 3*numS21 # (3 permutations) X (number of S21 orbits)
    numweights += numS21
    # initialize parameter arrays
    params = zeros(T, (numparams))
    weights = zeros(T, (numweights))
    new(numparams, numweights, numnodes, vertices, midedges, centroid,
        numedge, numS21, params, weights)
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
* `numS31` : number of S31 orbits (vertex to opposite face)
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
  numS31::Int
  params::Array{T}
  weights::Array{T}

  function TetSymCub(;numedge::Int=0, numS31::Int=0, vertices::Bool=true,
                        midedges::Bool=false, centroid::Bool=false,
                        facecentroid::Bool=false)
    @assert(numedge >= 0)
    @assert(numS31 >= 0)
    # compute the number of degrees of freedom and unique weights
    numparams = 0
    numweights = 0
    numparams += numedge + numS31
    # compute the number of nodes
    numnodes = 0
    if vertices
      numnodes += 4
      numweights += 1
    end
    if midedges
      numnodes += 6
      numweights += 1
    end
    if centroid
      numnodes += 1
      numweights += 1
    end
    if facecentroid
      numnodes += 4
      numweights += 1
    end
    numnodes += 12*numedge # 2 * (6 edges) X (number of edge parameters)
    numweights += numedge
    numnodes += 4*numS31 # (4 permutations) X (number of S31 orbits)
    numweights += numS31
    # initialize parameter arrays
    params = zeros(T, (numparams))
    weights = zeros(T, (numweights))
    new(numparams, numweights, numnodes, vertices, midedges, centroid,
        facecentroid, numedge, numS31, params, weights)
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
  return numboundary
end

@doc """
### SymCubatures.setparams!

Sets the nodal parameters for any parameterized symmetry orbits in the cubature.

**Inputs**

* `params`: parameter values

**In/Outs**

* `cub`: symmetric cubature rule whose nodal parameters are being updated

"""->
function setparams!{T}(cub::TriSymCub{T}, params::Array{T})
  @assert( length(params) == cub.numparams )
  size(params) == size(cub.params) ? cub.params = params :
  cub.params = params.'
end

function setparams!{T}(cub::TetSymCub{T}, params::Array{T})
  @assert( length(params) == cub.numparams )
  size(params) == size(cub.params) ? cub.params = params :
  cub.params = params.'
end

@doc """
### SymCubatures.setweights!

Sets a cubature's (unique) weights.

**Inputs**

* `weights`: cubature weights grouped by orbit

**In/Outs**

* `cub`: symmetric cubature rule whose weights are being updated

"""->
function setweights!{T}(cub::TriSymCub{T}, weights::Array{T})
  @assert( length(weights) == cub.numweights )
  size(weights) == size(cub.weights) ? cub.weights = weights : 
  cub.weights = weights.'
end

function setweights!{T}(cub::TetSymCub{T}, weights::Array{T})
  @assert( length(weights) == cub.numweights )
  size(weights) == size(cub.weights) ? cub.weights = weights :
  cub.weights = weights.'
end

@doc """
### SymCubatures.calcnodes

Use the orbital parameter values to compute a cubature's nodal coordinates.

**Inputs**

* `cub`: symmetric cubature rule
* `vtx`: vertices that define the domain

**Outputs**

* `x`,`y` (`z`): cubature's nodal coordinates

"""->
function calcnodes{T}(cub::TriSymCub{T}, vtx::Array{T,2})
  @assert(cub.numparams >= 0)
  @assert(cub.numnodes >= 1)
  @assert( size(vtx,1) == 3 && size(vtx,2) == 2)
  x = zeros(T, (cub.numnodes))
  y = zeros(T, (cub.numnodes))
  ptr = 0
  paramptr = 0
  # set vertices
  if cub.vertices
    x[1:3] = vtx[:,1]
    y[1:3] = vtx[:,2]
    ptr = 3
  end
  # set mid-edge nodes
  if cub.midedges
    A = T[0.5 0.5 0;
          0 0.5 0.5;
          0.5 0 0.5]
    x[ptr+1:ptr+3] = A*vtx[:,1]
    y[ptr+1:ptr+3] = A*vtx[:,2]
    ptr += 3
  end
  # set centroid !!!! This is a volume node!! move down ?!
  if cub.centroid
    A = T[1/3 1/3 1/3]
    x[ptr+1:ptr+1] = A*vtx[:,1]
    y[ptr+1:ptr+1] = A*vtx[:,2]
    ptr += 1
  end
  # set edge nodes
  for i = 1:cub.numedge
    alpha = cub.params[paramptr+1]
    A = T[alpha (1-alpha) 0;
          (1-alpha) alpha 0;
          0 alpha (1-alpha);
          0 (1-alpha) alpha;
          (1-alpha) 0 alpha;
          alpha 0 (1-alpha)]
    x[ptr+1:ptr+6] = A*vtx[:,1]
    y[ptr+1:ptr+6] = A*vtx[:,2]
    ptr += 6
    paramptr += 1
  end
  # Now set internal nodes, if any,...
  # set S21 orbit nodes
  for i = 1:cub.numS21
    alpha = 0.5*cub.params[paramptr+1]
    A = T[alpha alpha (1-2*alpha);
          (1-2*alpha) alpha alpha;
          alpha (1-2*alpha) alpha]
    x[ptr+1:ptr+3] = A*vtx[:,1]
    y[ptr+1:ptr+3] = A*vtx[:,2]
    ptr += 3
    paramptr += 1
  end
  # etc.
  return x, y
end

function calcnodes{T}(cub::TetSymCub{T}, vtx::Array{T,2})
  @assert(cub.numparams >= 0)
  @assert(cub.numnodes >= 1)
  @assert( size(vtx,1) == 4 && size(vtx,2) == 3)
  x = zeros(T, (cub.numnodes))
  y = zeros(T, (cub.numnodes))
  z = zeros(T, (cub.numnodes))
  ptr = 0
  paramptr = 0
  # set vertices
  if cub.vertices
    x[1:4] = vtx[:,1]
    y[1:4] = vtx[:,2]
    z[1:4] = vtx[:,3]
    ptr = 4
  end
  # set mid-edge nodes
  if cub.midedges
    A = T[0.5 0.5 0 0;
          0 0.5 0.5 0;
          0 0 0.5 0.5;
          0.5 0 0 0.5;
          0.5 0 0.5 0;
          0 0.5 0 0.5]
    x[ptr+1:ptr+6] = A*vtx[:,1]
    y[ptr+1:ptr+6] = A*vtx[:,2]
    z[ptr+1:ptr+6] = A*vtx[:,3]
    ptr += 6
  end
  # set centroid !!!! This is a volume node!! move down ?!
  if cub.centroid
    A = T[0.25 0.25 0.25 0.25]
    x[ptr+1:ptr+1] = A*vtx[:,1]
    y[ptr+1:ptr+1] = A*vtx[:,2]
    z[ptr+1:ptr+1] = A*vtx[:,3]
    ptr += 1
  end
  # set face centroids
  if cub.facecentroid
    alpha = (T)(1/3)
    A = T[alpha alpha alpha 0;
          0 alpha alpha alpha;
          alpha 0 alpha alpha;
          alpha alpha 0 alpha]
    x[ptr+1:ptr+4] = A*vtx[:,1]
    y[ptr+1:ptr+4] = A*vtx[:,2]
    z[ptr+1:ptr+4] = A*vtx[:,3]
    ptr += 4
  end
  # set edge nodes
  for i = 1:cub.numedge
    alpha = cub.params[paramptr+1]
    A = T[alpha (1-alpha) 0 0;
          (1-alpha) alpha 0 0;
          0 alpha (1-alpha) 0;
          0 (1-alpha) alpha 0;
          0 0 alpha (1-alpha);
          0 0 (1-alpha) alpha;
          (1-alpha) 0 0 alpha;
          alpha 0 0 (1-alpha);
          alpha 0 (1-alpha) 0;
          (1-alpha) 0 alpha 0;
          0 alpha 0 (1-alpha)
          0 (1-alpha) 0 alpha]
    x[ptr+1:ptr+12] = A*vtx[:,1]
    y[ptr+1:ptr+12] = A*vtx[:,2]
    z[ptr+1:ptr+12] = A*vtx[:,3]
    ptr += 12
    paramptr += 1
  end
  # Now set volume nodes, if any,...
  # set S31 orbit nodes
  for i = 1:cub.numS31
    alpha = cub.params[paramptr+1]/3
    A = T[alpha alpha alpha (1-3*alpha);
          (1-3*alpha) alpha alpha alpha;
          alpha (1-3*alpha) alpha alpha;
          alpha alpha (1-3*alpha) alpha]
    x[ptr+1:ptr+4] = A*vtx[:,1]
    y[ptr+1:ptr+4] = A*vtx[:,2]
    z[ptr+1:ptr+4] = A*vtx[:,3]
    ptr += 4
    paramptr += 1
  end
  # etc.
  return x, y, z
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
  if cub.vertices
    ptr += 3 # block of zeros, because the vertices are not parameterized
  end
  if cub.midedges
    ptr += 3 # block of zeros, because the midedges are not parameterized
  end
  if cub.centroid
    ptr += 1 # block of zeros, because the midedges are not parameterized
  end
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
  # etc.
  return Jac
end

function calcjacobianofnodes{T}(cub::TetSymCub{T}, vtx::Array{T,2})
  @assert(cub.numparams >= 0)
  @assert(cub.numnodes >= 1)
  @assert( size(vtx,1) == 4 && size(vtx,2) == 3)

  Jac = zeros(T, (3*cub.numnodes, cub.numparams) )
  ptr = 0
  paramptr = 0
  if cub.vertices
    ptr += 4 # block of zeros, because the vertices are not parameterized
  end
  if cub.midedges
    ptr += 6 # block of zeros, because the midedges are not parameterized
  end
  if cub.centroid
    ptr += 1 # block of zeros, because the midedges are not parameterized
  end
  if cub.facecentroid
    ptr += 4 # block of zeros, because the midedges are not parameterized
  end
  # set Jacobian for edge nodes
  A = T[1 -1 0 0;
        -1 1 0 0;
        0 1 -1 0;
        0 -1 1 0;
        0 0 1 -1;
        0 0 -1 1;
        -1 0 0 1;
        1 0 0 -1;
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
  # etc.
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
function calcweights{T}(cub::TriSymCub{T})
  @assert(cub.numweights >= 0)
  @assert(cub.numnodes >= 1)

  w = zeros(T, (cub.numnodes))
  ptr = 0
  wptr = 0
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
  # set centroid weight
  if cub.centroid
    w[ptr+1:ptr+1] = cub.weights[wptr+1]
    ptr += 1
    wptr += 1
  end
  # set edge weights
  for i = 1:cub.numedge
    w[ptr+1:ptr+6] = cub.weights[wptr+1]
    ptr += 6
    wptr += 1
  end

  # Now set volume weights, if any,...
  # set S21 orbit nodes
  for i = 1:cub.numS21
    w[ptr+1:ptr+3] = cub.weights[wptr+1]
    ptr += 3
    wptr += 1
  end
  # etc.
  return w
end

function calcweights{T}(cub::TetSymCub{T})
  @assert(cub.numweights >= 0)
  @assert(cub.numnodes >= 1)

  w = zeros(T, (cub.numnodes))
  ptr = 0
  wptr = 0
  # set vertices weights
  if cub.vertices
    w[1:4] = cub.weights[wptr+1]
    ptr = 4
    wptr += 1
  end
  # set mid-edge weights
  if cub.midedges
    w[ptr+1:ptr+6] = cub.weights[wptr+1]
    ptr += 6
    wptr += 1
  end
  # set centroid weight
  if cub.centroid
    w[ptr+1:ptr+1] = cub.weights[wptr+1]
    ptr += 1
    wptr += 1
  end
  # set face centroid weights
  if cub.facecentroid
    w[ptr+1:ptr+4] = cub.weights[wptr+1]
    ptr += 4
    wptr += 1
  end
  # set edge weights
  for i = 1:cub.numedge
    w[ptr+1:ptr+12] = cub.weights[wptr+1]
    ptr += 12
    wptr += 1
  end

  # Now set volume weights, if any,...
  # set S31 orbit nodes
  for i = 1:cub.numS31
    w[ptr+1:ptr+4] = cub.weights[wptr+1]
    ptr += 4
    wptr += 1
  end
  # etc.
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
  # set Jacobian of centroid weight
  if cub.centroid
    Jac[ptr+1:ptr+1,wptr+1] = ones(T, (1,1))
    ptr += 1
    wptr += 1
  end
  # set Jacobian of edge nodes
  for i = 1:cub.numedge
    Jac[ptr+1:ptr+6,wptr+1] = ones(T, (6,1))
    ptr += 6
    wptr += 1
  end

  # Now set Jacobian of volume weights, if any,...
  # set S21 orbit Jacobian weights
  for i = 1:cub.numS21
    Jac[ptr+1:ptr+3,wptr+1] = ones(T, (3,1))
    ptr += 3
    wptr += 1
  end
  # etc.
  return Jac
end

function calcjacobianofweights{T}(cub::TetSymCub{T})
  @assert(cub.numweights >= 0)
  @assert(cub.numnodes >= 1)

  Jac = zeros(T, (cub.numnodes, cub.numweights) )
  ptr = 0
  wptr = 0
  # set Jacobian of vertex weights
  if cub.vertices
    Jac[1:4,wptr+1] = ones(T, (4,1))
    ptr = 4
    wptr += 1
  end
  # set Jacobian of mid-edge weights
  if cub.midedges
    Jac[ptr+1:ptr+6,wptr+1] = ones(T, (6,1))
    ptr += 6
    wptr += 1
  end
  # set Jacobian of centroid weight
  if cub.centroid
    Jac[ptr+1:ptr+1,wptr+1] = ones(T, (1,1))
    ptr += 1
    wptr += 1
  end
  # set Jacobian of face centroid weights
  if cub.facecentroid
    Jac[ptr+1:ptr+4,wptr+1] = ones(T, (4,1))
    ptr += 4
    wptr += 1
  end
  # set Jacobian of edge nodes
  for i = 1:cub.numedge
    Jac[ptr+1:ptr+12,wptr+1] = ones(T, (12,1))
    ptr += 12
    wptr += 1
  end

  # Now set Jacobian of volume weights, if any,...
  # set S31 orbit Jacobian weights
  for i = 1:cub.numS31
    Jac[ptr+1:ptr+4,wptr+1] = ones(T, (4,1))
    ptr += 4
    wptr += 1
  end
  # etc.
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
