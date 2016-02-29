# This file gathers together a hodge-podge of functions that are not easily
# categorized

@doc """
### SummationByParts.buildinterpolation

Builds a matrix operator that can reconstruct a field located at the sbp nodes
to an auxlliary set of nodes.

**Inputs**

* `sbp`: an SBP operator
* `xinterp`: location of points to interpolate to, size = [ndim,numpoints]

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