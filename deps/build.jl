# SummationByParts.jl uses the following (very small) Package, which is not
# listed in the METADATA; therefore, it must be added explicitly

# install PkgFix if not present
if !isdir(joinpath(Pkg.dir(), "PkgFix"))
  Pkg.clone("https://github.com/OptimalDesignLab/PkgFix.jl.git")
end
Pkg.checkout("PkgFix", "upgrade_0.6")

using PkgFix  # from now on, use PkgFix instead of Pkg for everything

ODL_URL = "https://github.com/OptimalDesignLab/ODLCommonTools.jl.git"
ODL_VER = "upgrade_0.6"
ARRAYVIEWS_VER = "master"


pkg_dict = PkgFix.installed()
if !haskey(pkg_dict, "ODLCommonTools")
  PkgFix.add(ODL_URL, branch_ish=ODL_VER)
end
# get the right commit of ArrayViews (after the aview change)
PkgFix.checkout("ArrayViews", ARRAYVIEWS_VER)
