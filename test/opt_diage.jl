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
  writedlm("tet_diage_p$degree.dat", Q)
end

