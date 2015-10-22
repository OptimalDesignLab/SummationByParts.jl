# This file gathers together functions used to build the SBP face operators

function buildfacereconstruction{T}(facecub::LineSymCub{T}, cub::TriSymCub{T},
                                    vtx::Array{T,2}, d::Int; faceonly::Bool=false)
  perm = SymCubatures.getfacebasedpermutation(cub, faceonly=faceonly)
  # evaluate the basis at the volume and face cubature points
  N = convert(Int, (d+1)*(d+2)/2 )
  Pv = zeros(T, (size(perm,1),N) )  
  Pf = zeros(T, (facecub.numnodes,N) ) 
  xv = SymCubatures.calcnodes(cub, vtx)
  xf = SymCubatures.calcnodes(facecub, vtx[[1;2],:])
  ptr = 1
  for r = 0:d
    for j = 0:r
      i = r-j
      Pv[:,ptr] = OrthoPoly.proriolpoly(vec(xv[1,perm[:,1]]),
                                        vec(xv[2,perm[:,1]]), i, j)
      Pf[:,ptr] = OrthoPoly.proriolpoly(vec(xf[1,:]), vec(xf[2,:]), i, j)
      ptr += 1
    end
  end
  R = Pf/Pv
  return R, perm
end