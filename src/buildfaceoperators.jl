# This file gathers together functions used to build the SBP face operators

function reconstructionoperators{T}(facecub::LineSymCub{T}, cub::TriSymCub{T},
                                    vtx::Array{T,2}, d::Int; faceonly::Bool=false)
  perm = SymCubatures.getfacebasedpermutation(cub)
  # evaluate the basis at the volume and face cubature points
  N = convert(Int, (d+1)*(d+2)/2 )
  Pv = zeros(T, (size(perm,1),N) )  
  Pf = zeros(T, (facecub.numnodes,N) ) 
  xv, yv = SymCubatures.calcnodes(cub, vtx)
  xf, yf = SymCubatures.calcnodes(facecub, vtx[[3;1],:])
  ptr = 1
  for r = 0:d
    for j = 0:r
      i = r-j
      Pv[:,ptr] = OrthoPoly.proriolpoly(xv[perm[:,1]], yv[perm[:,1]], i, j)
      Pf[:,ptr] = OrthoPoly.proriolpoly(xf, yf, i, j)
      ptr += 1
    end
  end
  R = Pb/P
  return R, perm
end