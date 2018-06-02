# SummationByParts.jl uses the following (very small) Package, which is not
# listed in the METADATA; therefore, it must be added explicitly

# install PkgFix if not present
if !isdir(joinpath(Pkg.dir(), "PkgFix"))
  Pkg.clone("https://github.com/OptimalDesignLab/PkgFix.jl.git")
end
Pkg.checkout("PkgFix", "upgrade_0.5")

using PkgFix  # from now on, use PkgFix instead of Pkg for everything

ODL_URL = "https://github.com/OptimalDesignLab/ODLCommonTools.jl.git"
ODL_VER = "update_0.5"
ARRAYVIEWS_VER = "93e80390aeedb1dbcd90281b6dff7f760f430bc8"


pkg_dict = PkgFix.installed()
if !haskey(pkg_dict, "ODLCommonTools")
  PkgFix.add(ODL_URL, branch_ish=ODL_VER)
end
# get the right commit of ArrayViews (before the aview change)
PkgFix.checkout("ArrayViews", ARRAYVIEWS_VER)
