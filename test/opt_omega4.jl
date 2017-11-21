using SummationByParts
using SummationByParts.Cubature
using SummationByParts.SymCubatures

function getOptOmega4(degree::Int=1, Tsbp::Type=Float64;
                        vertices::Bool=false, opt_tol::Float64=1e-10)
  cub, vtx = getTriCubatureOmega(2*degree, Tsbp)
  w = zeros(Tsbp, (cub.numnodes))
  Q = zeros(Tsbp, (cub.numnodes, cub.numnodes, 3))
  w, Q = SummationByParts.buildMinConditionOperators(cub, vtx, degree,
                                                     tol=opt_tol,
                                                     vertices=vertices,
                                                     opthist=true)
  writedlm("tri_omega4_p$degree.dat", Q)
end

getOptOmega4(3, opt_tol=1e-60)
