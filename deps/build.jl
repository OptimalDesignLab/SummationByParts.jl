# SummationByParts.jl uses the following (very small) Package, which is not
# listed in the METADATA; therefore, it must be added explicitly
pkg_dict = Pkg.installed()
start_dir = pwd()
deps_path = joinpath(Pkg.dir("SummationByParts"), "deps")
cd(deps_path)
if !haskey(pkg_dict, "ODLCommonTools")
  run(`./download.sh`)
  start_dir2 = pwd()
  cd(Pkg.dir("ODLCommonTools"))
  run(`git checkout new_parallel`)
  cd(start_dir2)
  Pkg.build("ODLCommonTools")
end
cd(start_dir)

# get the right commit of ArrayViews (before the aview change)
dir = Pkg.dir("ArrayViews")
cd(dir)
run(`git checkout 93e80390aeedb1dbcd90281b6dff7f760f430bc8`)
cd(start_dir)
