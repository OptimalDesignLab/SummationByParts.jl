# SummationByParts.jl uses the following (very small) Package, which is not
# listed in the METADATA; therefore, it must be added explicitly
odlcommontools_git = "https://github.com/OptimalDesignLab/ODLCommonTools.jl.git"
pkg_dict = Pkg.installed()
if !haskey(pkg_dict, "ODLCommonTools")
  Pkg.clone(odlcommontools_git)
  Pkg.build("ODLCommonTools")
end