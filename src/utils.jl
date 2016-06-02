# This file gathers together a hodge-podge of functions that are not easily
# categorized

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
### SummationByParts.findleftperm

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

"""->
function findleftperm!{T}(A::AbstractArray{T,2}, permR::AbstractVector{Int},
                          permL::AbstractVector{Int})
  @assert( size(A,1) == length(permL) )
  @assert( size(A,2) == length(permR) )
  rows = [ sub(A,i,1:size(A,2)) for i=1:size(A,1) ]
  permA = sortperm(rows; order=Base.Lexicographic)
  AR = A[:,permR]
  rows = [ sub(AR,i,1:size(AR,2)) for i=1:size(AR,1) ]
  permAR = sortperm(rows, order=Base.Lexicographic)
  invpermAR = invperm(permAR)
  permL = permA[invpermAR]
  if A[permL,:] == AR
    return true
  else
    return false
  end
end
