# SummationByParts.jl uses the following (very small) Package, which is not
# listed in the METADATA; therefore, it must be added explicitly
pdesolvercommon_git = "https://github.com/OptimalDesignLab/PDESolverCommon.jl.git"
pkg_dict = Pkg.installed()
if !haskey(pkg_dict, "PDESolverCommon")
  Pkg.clone(pdesolvercommon_git)
  Pkg.build("PDESolverCommon")
end