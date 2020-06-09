#include("/home/jehicken/.julia/v0.4/SummationByParts/src/SummationByParts.jl")
using SummationByParts

# function to print triangle SBP nodes to file
function printTriNodes(degree::Int, internal::Bool=false,
                       filename::AbstractString="nodes.dat")
  @assert( degree >= 1 && degree <= 4)
  sbp = TriSBP{Float64}(degree=degree, reorder=false, internal=internal)
  sbpface = TriFace{Float64}(degree, sbp.cub, [-1. -1; 1 -1; -1 1])
  # compute the volume nodes
  vtx = [0.0 0.0; 1 0.0; 0.5 sqrt(3)*0.5]
  x = SummationByParts.SymCubatures.calcnodes(sbp.cub, vtx)
  xf = zeros(2, sbpface.numnodes, 3)
  for f = 1:3
    xf[:,:,f] = SummationByParts.SymCubatures.calcnodes(sbpface.cub,
                                                        vtx[[f;mod(f,3)+1],:])
  end
  f = open(filename, "w")
  for xelem in x[1,:]
    print(f, xelem, " ")
  end
  println(f)
  for yelem in x[2,:]
    print(f, yelem, " ")
  end
  println(f)
  for xface in xf[1,:,:]
    print(f, xface, " ")
  end
  println(f)
  for yface in xf[2,:,:]
    print(f, yface, " ")
  end
  println(f)
  close(f)
end

# function to print triangle SBP nodes to file
function printTriOmegaNodes(degree::Int,
                            filename::AbstractString="nodes.dat")
  @assert( degree >= 1 && degree <= 4)
  sbp = getTriSBPOmega(degree=degree)
  sbpface = TriFace{Float64}(degree, sbp.cub, [-1. -1; 1 -1; -1 1])
  # compute the volume nodes
  vtx = [0.0 0.0; 1 0.0; 0.5 sqrt(3)*0.5]
  x = SummationByParts.SymCubatures.calcnodes(sbp.cub, vtx)
  xf = zeros(2, sbpface.numnodes, 3)
  for f = 1:3
    xf[:,:,f] = SummationByParts.SymCubatures.calcnodes(sbpface.cub,
                                                        vtx[[f;mod(f,3)+1],:])
  end
  f = open(filename, "w")
  for xelem in x[1,:]
    print(f, xelem, " ")
  end
  println(f)
  for yelem in x[2,:]
    print(f, yelem, " ")
  end
  println(f)
  for xface in xf[1,:,:]
    print(f, xface, " ")
  end
  println(f)
  for yface in xf[2,:,:]
    print(f, yface, " ")
  end
  println(f)
  close(f)
end

# function to print triangle SBP nodes to file
function printTriDiagENodes(degree::Int,
                            filename::AbstractString="nodes.dat")
  @assert( degree >= 1 && degree <= 4)
  sbp = getTriSBPDiagE(degree=degree, mincond=true)
  sbpface = getTriFaceForDiagE(degree, sbp.cub, sbp.vtx)
  # compute the volume nodes
  vtx = [0.0 0.0; 1 0.0; 0.5 sqrt(3)*0.5]
  #vtx = sbp.vtx
  x = SummationByParts.SymCubatures.calcnodes(sbp.cub, vtx)
  xf = zeros(2, sbpface.numnodes, 3)
  for f = 1:3
    xf[:,:,f] = SummationByParts.SymCubatures.calcnodes(sbpface.cub,
                                                        vtx[[f;mod(f,3)+1],:])
  end
  f = open(filename, "w")
  for xelem in x[1,:]
    print(f, xelem, " ")
  end
  println(f)
  for yelem in x[2,:]
    print(f, yelem, " ")
  end
  println(f)
  for xface in xf[1,:,:]
    print(f, xface, " ")
  end
  println(f)
  for yface in xf[2,:,:]
    print(f, yface, " ")
  end
  println(f)
  close(f)
end

# function to print tetrahedral SBP nodes to file
function printTetNodes(degree::Int, internal::Bool=false,
                       filename::AbstractString="nodes.dat")
  @assert( degree >= 1 && degree <= 4)
  sbp = getTetSBPGamma(degree=degree) #, reorder=false, internal=internal)
  sbpface = TetFace{Float64}(degree, sbp.cub, sbp.vtx)
  # compute the volume nodes
  vtx = [0.0 0.0 0.0;
         1.0 0.0 0.0;
         0.5 sqrt(3)/2 0.0;
         0.5 sqrt(3)/6 sqrt(2/3)]  
  x = SummationByParts.SymCubatures.calcnodes(sbp.cub, vtx)
  xf = zeros(3, sbpface.numnodes, 4)
  facevtx = [1 1 2 1;
             2 4 4 3;
             3 2 3 4]
  for f = 1:4
    xf[:,:,f] = SummationByParts.SymCubatures.calcnodes(sbpface.cub,
                                                        vtx[facevtx[:,f],:])
  end
  f = open(filename, "w")
  for xelem in x[1,:]
    print(f, xelem, " ")
  end
  println(f)
  for yelem in x[2,:]
    print(f, yelem, " ")
  end
  println(f)
  for zelem in x[3,:]
    print(f, zelem, " ")
  end
  println(f)
  for xface in xf[1,:,:]
    print(f, xface, " ")
  end
  println(f)
  for yface in xf[2,:,:]
    print(f, yface, " ")
  end
  println(f)
  for zface in xf[3,:,:]
    print(f, zface, " ")
  end
  println(f)
  close(f)
end
