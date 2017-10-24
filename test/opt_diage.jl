using SummationByParts
using SummationByParts.Cubature
using SummationByParts.SymCubatures

function getOptTetDiagE(degree::Int=1, Tsbp::Type=Float64;
                        vertices::Bool=false, opt_tol::Float64=1e-10)
  cub, vtx = getTetCubatureDiagE(2*degree, Tsbp, vertices=vertices)
  w = zeros(Tsbp, (cub.numnodes))
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 3))
  w, Q = SummationByParts.buildMinConditionOperators(cub, vtx, degree,
                                                     tol=opt_tol,
                                                     vertices=vertices,
                                                     opthist=true)
  f = open("tet_diage_p$degree.jl", "w")
  print(f, "w[:] = ")
  println(f, w)
  println(f)
  print(f, "Q[:,:,1] = ")
  println(f, Q[:,:,1])
  println(f)
  print(f, "Q[:,:,2] = ")
  println(f, Q[:,:,2])
  println(f)
  print(f, "Q[:,:,3] = ")
  println(f, Q[:,:,3])
  close(f)
end

