using SummationByParts
using SummationByParts.OrthoPoly
using SummationByParts.Cubature
using SummationByParts.SymCubatures
using FactCheck
using ArrayViews

facts("Checking run time of various operations...") do

if false
  context("Timing SummationByParts.weakdifferentiate! (TriSBP, scalar field method)") do
    # warm-up 
    sbp = TriSBP{Float64}(degree=1)
    u = ones(Float64, (sbp.numnodes,2))
    di = 1
    res = zeros(u)
    print("\tIgnore this --->")
    @time weakdifferentiate!(sbp, di, u, res, trans=true)
    numelem = 10000
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,numelem))
      di = 1
      res = zeros(u)
      print("p = ",p,": ")
      @time weakdifferentiate!(sbp, di, u, res, trans=true)
    end
    println()
  end

  context("Timing SummationByParts.weakdifferentiate! (TetSBP, scalar field method)") do
    # warm-up 
    sbp = TetSBP{Float64}(degree=1)
    u = ones(Float64, (sbp.numnodes,2))
    di = 1
    res = zeros(u)
    print("\tIgnore this --->")
    @time weakdifferentiate!(sbp, di, u, res, trans=true)
    numelem = 10000
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,numelem))
      di = 1
      res = zeros(u)
      print("p = ",p,": ")
      @time weakdifferentiate!(sbp, di, u, res, trans=true)
    end
    println()
  end

  context("Timing SummationByParts.weakdifferentiate! (TriSBP, vector field method)") do
    # warm-up 
    sbp = TriSBP{Float64}(degree=1)
    u = ones(Float64, (5,sbp.numnodes,2))
    di = 1
    res = zeros(u)
    print("\tIgnore this --->")
    @time weakdifferentiate!(sbp, di, u, res, trans=true)
    numelem = 10000
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      u = ones(Float64, (5,sbp.numnodes,numelem))
      di = 1
      res = zeros(u)
      print("p = ",p,": ")
      @time weakdifferentiate!(sbp, di, u, res, trans=true)
    end
    println()
  end
  
  context("Timing SummationByParts.weakdifferentiate! (TetSBP, vector field method)") do
    # warm-up 
    sbp = TetSBP{Float64}(degree=1)
    u = ones(Float64, (5,sbp.numnodes,2))
    di = 1
    res = zeros(u)
    print("\tIgnore this --->")
    @time weakdifferentiate!(sbp, di, u, res, trans=true)
    numelem = 10000
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (5,sbp.numnodes,numelem))
      di = 1
      res = zeros(u)
      print("p = ",p,": ")
      @time weakdifferentiate!(sbp, di, u, res, trans=true)
    end
    println()
  end

  context("Timing SummationByParts.differentiate! (TriSBP, scalar field method)") do
    # warm-up 
    sbp = TriSBP{Float64}(degree=1)
    u = ones(Float64, (sbp.numnodes,2))
    di = 1
    res = zeros(u)
    print("\tIgnore this --->")
    @time differentiate!(sbp, di, u, res)
    numelem = 10000
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,numelem))
      di = 1
      res = zeros(u)
      print("p = ",p,": ")
      @time differentiate!(sbp, di, u, res)
    end
    println()
  end

  context("Timing SummationByParts.differentiate! (TetSBP, scalar field method)") do
    # warm-up 
    sbp = TetSBP{Float64}(degree=1)
    u = ones(Float64, (sbp.numnodes,2))
    di = 1
    res = zeros(u)
    print("\tIgnore this --->")
    @time differentiate!(sbp, di, u, res)
    numelem = 10000
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,numelem))
      di = 1
      res = zeros(u)
      print("p = ",p,": ")
      @time differentiate!(sbp, di, u, res)
    end
    println()
  end

  context("Timing SummationByParts.differentiate! (TriSBP, vector field method)") do
    # warm-up 
    sbp = TriSBP{Float64}(degree=1)
    u = ones(Float64, (5,sbp.numnodes,2))
    di = 1
    res = zeros(u)
    print("\tIgnore this --->")
    @time differentiate!(sbp, di, u, res)
    numelem = 10000
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      u = ones(Float64, (5,sbp.numnodes,numelem))
      di = 1
      res = zeros(u)
      print("p = ",p,": ")
      @time differentiate!(sbp, di, u, res)
    end
    println()
  end
  
  context("Timing SummationByParts.differentiate! (TetSBP, vector field method)") do
    # warm-up 
    sbp = TetSBP{Float64}(degree=1)
    u = ones(Float64, (5,sbp.numnodes,2))
    di = 1
    res = zeros(u)
    print("\tIgnore this --->")
    @time differentiate!(sbp, di, u, res)
    numelem = 10000
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (5,sbp.numnodes,numelem))
      di = 1
      res = zeros(u)
      print("p = ",p,": ")
      @time differentiate!(sbp, di, u, res)
    end
    println()
  end

  context("Timing SummationByParts.directionaldifferentiate! (TriSBP, scalar field method)") do
    function elemloop{T}(sbp::SBPOperator{T}, dir::Array{T,1}, 
                         u::AbstractArray{T,2}, i::Int, numelem::Int)
      Ddir = zero(T)
      for k = 1:numelem 
        Ddir = directionaldifferentiate(sbp, dir, view(u,:,k), i)
      end
    end
    # warm-up 
    sbp = TriSBP{Float64}(degree=1)
    u = ones(Float64, (sbp.numnodes,2))
    i = rand(1:sbp.numnodes)
    dir = rand(2)
    print("\tIgnore this --->")
    @time elemloop(sbp, dir, u, i, 2)
    numelem = 10000
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,numelem))
      print("p = ",p,": ")
      @time elemloop(sbp, dir, u, i, numelem)
    end
    println()
  end

  context("Timing SummationByParts.directionaldifferentiate! (TriSBP, vector field method)") do
    function elemloop{T}(sbp::SBPOperator{T}, dir::Array{T,1}, 
                         u::AbstractArray{T,3}, i::Int, numelem::Int)
      Ddir = zeros(T, size(u,1))
      for k = 1:numelem 
        Ddir = directionaldifferentiate(sbp, dir, view(u,:,:,k), i)
      end
    end
    # warm-up 
    sbp = TriSBP{Float64}(degree=1)
    u = ones(Float64, (5,sbp.numnodes,2))
    i = rand(1:sbp.numnodes)
    dir = rand(2)
    print("\tIgnore this --->")
    @time elemloop(sbp, dir, u, i, 2)
    numelem = 10000
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      u = ones(Float64, (5,sbp.numnodes,numelem))
      print("p = ",p,": ")
      @time elemloop(sbp, dir, u, i, numelem)
    end
    println()
  end

  context("Timing SummationByParts.directionaldifferentiate! (TetSBP, scalar field method)") do
    function elemloop{T}(sbp::SBPOperator{T}, dir::Array{T,1}, 
                         u::AbstractArray{T,2}, i::Int, numelem::Int)
      Ddir = zero(T)
      for k = 1:numelem 
        Ddir = directionaldifferentiate(sbp, dir, view(u,:,k), i)
      end
    end
    # warm-up 
    sbp = TetSBP{Float64}(degree=1)
    u = ones(Float64, (sbp.numnodes,2))
    i = rand(1:sbp.numnodes)
    dir = rand(3)
    print("\tIgnore this --->")
    @time elemloop(sbp, dir, u, i, 2)
    numelem = 10000
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,numelem))
      print("p = ",p,": ")
      @time elemloop(sbp, dir, u, i, numelem)
    end
    println()
  end

  context("Timing SummationByParts.directionaldifferentiate! (TetSBP, vector field method)") do
    function elemloop{T}(sbp::SBPOperator{T}, dir::Array{T,1}, 
                         u::AbstractArray{T,3}, i::Int, numelem::Int)
      Ddir = zeros(T, size(u,1))
      for k = 1:numelem 
        Ddir = directionaldifferentiate(sbp, dir, view(u,:,:,k), i)
      end
    end
    # warm-up 
    sbp = TetSBP{Float64}(degree=1)
    u = ones(Float64, (5,sbp.numnodes,2))
    i = rand(1:sbp.numnodes)
    dir = rand(3)
    print("\tIgnore this --->")
    @time elemloop(sbp, dir, u, i, 2)
    numelem = 10000
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (5,sbp.numnodes,numelem))
      print("p = ",p,": ")
      @time elemloop(sbp, dir, u, i, numelem)
    end
    println()
  end

  context("Timing SummationByParts.volumeintegrate! (TriSBP, scalar field method)") do
    # warm-up 
    sbp = TriSBP{Float64}(degree=1)
    u = ones(Float64, (sbp.numnodes,2))
    res = zeros(u)
    print("\tIgnore this --->")
    @time volumeintegrate!(sbp, u, res)
    numelem = 10000
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,numelem))
      res = zeros(u)
      print("p = ",p,": ")
      @time volumeintegrate!(sbp, u, res)
    end
    println()
  end

  context("Timing SummationByParts.volumeintegrate! (TetSBP, scalar field method)") do
    # warm-up 
    sbp = TetSBP{Float64}(degree=1)
    u = ones(Float64, (sbp.numnodes,2))
    res = zeros(u)
    print("\tIgnore this --->")
    @time volumeintegrate!(sbp, u, res)
    numelem = 10000
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,numelem))
      res = zeros(u)
      print("p = ",p,": ")
      @time volumeintegrate!(sbp, u, res)
    end
    println()
  end

  context("Timing SummationByParts.volumeintegrate! (TriSBP, vector field method)") do
    # warm-up 
    sbp = TriSBP{Float64}(degree=1)
    u = ones(Float64, (5,sbp.numnodes,2))
    res = zeros(u)
    print("\tIgnore this --->")
    @time volumeintegrate!(sbp, u, res)
    numelem = 10000
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      u = ones(Float64, (5,sbp.numnodes,numelem))
      res = zeros(u)
      print("p = ",p,": ")
      @time volumeintegrate!(sbp, u, res)
    end
    println()
  end

  context("Timing SummationByParts.volumeintegrate! (TetSBP, vector field method)") do
    # warm-up 
    sbp = TetSBP{Float64}(degree=1)
    u = ones(Float64, (5,sbp.numnodes,2))
    res = zeros(u)
    print("\tIgnore this --->")
    @time volumeintegrate!(sbp, u, res)
    numelem = 10000
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (5,sbp.numnodes,numelem))
      res = zeros(u)
      print("p = ",p,": ")
      @time volumeintegrate!(sbp, u, res)
    end
    println()
  end

  context("Timing SummationByParts.boundaryintegrate! (TriSBP, scalar field method)") do
    function bndryflux{T}(u::T, dξdx::AbstractArray{T,2}, nrm::AbstractArray{T,1})
      return u*(nrm[1]*(dξdx[1,1] + dξdx[1,2]) + nrm[2]*(dξdx[2,1] + dξdx[2,2]))
    end
    # warm-up
    sbp = TriSBP{Float64}(degree=1)
    x = zeros(Float64, (2,sbp.numnodes,1))
    vtx = [0. 0.; 1. 0.; 0. 1.]
    x[:,:,1] = calcnodes(sbp, vtx)
    dξdx = zeros(Float64, (2,2,sbp.numnodes,1))
    jac = zeros(Float64, (sbp.numnodes,1))
    mappingjacobian!(sbp, x, dξdx, jac)
    bndryfaces = Array(Boundary, 1)
    bndryfaces[1] = Boundary(1,1)
    u = ones(sbp.numnodes,1)
    res = zeros(u)
    print("\tIgnore this --->")
    @time boundaryintegrate!(sbp, bndryfaces, u, dξdx, bndryflux, res)

    numelem = 10000
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      x = zeros(Float64, (2,sbp.numnodes,numelem))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      for k = 1:numelem
        x[:,:,k] = calcnodes(sbp, vtx)
      end
      dξdx = zeros(Float64, (2,2,sbp.numnodes,numelem))
      jac = zeros(Float64, (sbp.numnodes,numelem))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, numelem)
      for k = 1:numelem
        bndryfaces[k] = Boundary(k,rand(1:3))
      end
      u = ones(sbp.numnodes, numelem)
      res = zeros(u)
      print("p = ",p,": ")
      @time boundaryintegrate!(sbp, bndryfaces, u, dξdx, bndryflux, res)
    end
  end

  context("Timing SummationByParts.boundaryintegrate! (TetSBP, scalar field method)") do
    function bndryflux{T}(u::T, dξdx::AbstractArray{T,2}, nrm::AbstractArray{T,1})
      return u*(nrm[1]*(dξdx[1,1] + dξdx[1,2] + dξdx[1,3]) + 
                nrm[2]*(dξdx[2,1] + dξdx[2,2] + dξdx[2,3]) +
                nrm[3]*(dξdx[3,1] + dξdx[3,2] + dξdx[3,3]))
    end
    # warm-up
    sbp = TetSBP{Float64}(degree=1)
    x = zeros(Float64, (3,sbp.numnodes,1))
    vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
    x[:,:,1] = calcnodes(sbp, vtx)
    dξdx = zeros(Float64, (3,3,sbp.numnodes,1))
    jac = zeros(Float64, (sbp.numnodes,1))
    mappingjacobian!(sbp, x, dξdx, jac)
    bndryfaces = Array(Boundary, 1)
    bndryfaces[1] = Boundary(1,1)
    u = ones(sbp.numnodes,1)
    res = zeros(u)
    print("\tIgnore this --->")
    @time boundaryintegrate!(sbp, bndryfaces, u, dξdx, bndryflux, res)

    numelem = 10000
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      x = zeros(Float64, (3,sbp.numnodes,numelem))
      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      for k = 1:numelem
        x[:,:,k] = calcnodes(sbp, vtx)
      end
      dξdx = zeros(Float64, (3,3,sbp.numnodes,numelem))
      jac = zeros(Float64, (sbp.numnodes,numelem))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, numelem)
      for k = 1:numelem
        bndryfaces[k] = Boundary(k,rand(1:3))
      end
      u = ones(sbp.numnodes, numelem)
      res = zeros(u)
      print("p = ",p,": ")
      @time boundaryintegrate!(sbp, bndryfaces, u, dξdx, bndryflux, res)
    end
  end

  context("Timing SummationByParts.boundaryintegrate! (TriSBP, vector field method)") do
    function bndryflux{T}(u::AbstractArray{T,1}, dξdx::AbstractArray{T,2}, 
                          nrm::AbstractArray{T,1}, flux::AbstractArray{T,1})
      #return u*(nrm[1]*(dξdx[1,1] + dξdx[1,2]) + nrm[2]*(dξdx[2,1] + dξdx[2,2]))
      tmp = (nrm[1]*(dξdx[1,1] + dξdx[1,2]) + nrm[2]*(dξdx[2,1] + dξdx[2,2]))
      for field = 1:size(u,1)
        flux[field] = u[field]*tmp
      end
    end
    # warm-up
    sbp = TriSBP{Float64}(degree=1)
    x = zeros(Float64, (2,sbp.numnodes,1))
    vtx = [0. 0.; 1. 0.; 0. 1.]
    x[:,:,1] = calcnodes(sbp, vtx)
    dξdx = zeros(Float64, (2,2,sbp.numnodes,1))
    jac = zeros(Float64, (sbp.numnodes,1))
    mappingjacobian!(sbp, x, dξdx, jac)
    bndryfaces = Array(Boundary, 1)
    bndryfaces[1] = Boundary(1,1)
    u = ones(5,sbp.numnodes,1)
    res = zeros(u)
    print("\tIgnore this --->")
    @time boundaryintegrate!(sbp, bndryfaces, u, dξdx, bndryflux, res)

    numelem = 10000
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      x = zeros(Float64, (2,sbp.numnodes,numelem))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      for k = 1:numelem
        x[:,:,k] = calcnodes(sbp, vtx)
      end
      dξdx = zeros(Float64, (2,2,sbp.numnodes,numelem))
      jac = zeros(Float64, (sbp.numnodes,numelem))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, numelem)
      for k = 1:numelem
        bndryfaces[k] = Boundary(k,rand(1:3))
      end
      u = ones(5, sbp.numnodes, numelem)
      res = zeros(u)
      print("p = ",p,": ")
      @time boundaryintegrate!(sbp, bndryfaces, u, dξdx, bndryflux, res)
    end
  end

  context("Timing SummationByParts.boundaryintegrate! (TetSBP, vector field method)") do
    function bndryflux{T}(u::AbstractArray{T,1}, dξdx::AbstractArray{T,2},
                          nrm::AbstractArray{T,1}, flux::AbstractArray{T,1})
      flux = u*(nrm[1]*(dξdx[1,1] + dξdx[1,2] + dξdx[1,3]) + 
                nrm[2]*(dξdx[2,1] + dξdx[2,2] + dξdx[2,3]) +
                nrm[3]*(dξdx[3,1] + dξdx[3,2] + dξdx[3,3]))
    end
    # warm-up
    sbp = TetSBP{Float64}(degree=1)
    x = zeros(Float64, (3,sbp.numnodes,1))
    vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
    x[:,:,1] = calcnodes(sbp, vtx)
    dξdx = zeros(Float64, (3,3,sbp.numnodes,1))
    jac = zeros(Float64, (sbp.numnodes,1))
    mappingjacobian!(sbp, x, dξdx, jac)
    bndryfaces = Array(Boundary, 1)
    bndryfaces[1] = Boundary(1,1)
    u = ones(5,sbp.numnodes,1)
    res = zeros(u)
    print("\tIgnore this --->")
    @time boundaryintegrate!(sbp, bndryfaces, u, dξdx, bndryflux, res)

    numelem = 10000
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      x = zeros(Float64, (3,sbp.numnodes,numelem))
      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      for k = 1:numelem
        x[:,:,k] = calcnodes(sbp, vtx)
      end
      dξdx = zeros(Float64, (3,3,sbp.numnodes,numelem))
      jac = zeros(Float64, (sbp.numnodes,numelem))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, numelem)
      for k = 1:numelem
        bndryfaces[k] = Boundary(k,rand(1:3))
      end
      u = ones(5, sbp.numnodes, numelem)
      res = zeros(u)
      print("p = ",p,": ")
      @time boundaryintegrate!(sbp, bndryfaces, u, dξdx, bndryflux, res)
    end
  end

  context("Timing SummationByParts.boundaryintegrate! (TriSBP, scalar field method)") do
    function bndryflux{T}(u::T, x::AbstractArray{T,1}, dξdx::AbstractArray{T,2}, 
                          nrm::AbstractArray{T,1})
      return u*(nrm[1]*(dξdx[1,1] + dξdx[1,2]) + nrm[2]*(dξdx[2,1] + dξdx[2,2]))
    end
    # warm-up
    sbp = TriSBP{Float64}(degree=1)
    x = zeros(Float64, (2,sbp.numnodes,1))
    vtx = [0. 0.; 1. 0.; 0. 1.]
    x[:,:,1] = calcnodes(sbp, vtx)
    dξdx = zeros(Float64, (2,2,sbp.numnodes,1))
    jac = zeros(Float64, (sbp.numnodes,1))
    mappingjacobian!(sbp, x, dξdx, jac)
    bndryfaces = Array(Boundary, 1)
    bndryfaces[1] = Boundary(1,1)
    u = ones(sbp.numnodes,1)
    res = zeros(u)
    print("\tIgnore this --->")
    @time boundaryintegrate!(sbp, bndryfaces, u, x, dξdx, bndryflux, res)

    numelem = 10000
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      x = zeros(Float64, (2,sbp.numnodes,numelem))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      for k = 1:numelem
        x[:,:,k] = calcnodes(sbp, vtx)
      end
      dξdx = zeros(Float64, (2,2,sbp.numnodes,numelem))
      jac = zeros(Float64, (sbp.numnodes,numelem))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, numelem)
      for k = 1:numelem
        bndryfaces[k] = Boundary(k,rand(1:3))
      end
      u = ones(sbp.numnodes, numelem)
      res = zeros(u)
      print("p = ",p,": ")
      @time boundaryintegrate!(sbp, bndryfaces, u, x, dξdx, bndryflux, res)
    end
  end

  context("Timing SummationByParts.boundaryintegrate! (TetSBP, scalar field method)") do
    function bndryflux{T}(u::T, x::AbstractArray{T,1}, dξdx::AbstractArray{T,2},
                          nrm::AbstractArray{T,1})
      return u*(nrm[1]*(dξdx[1,1] + dξdx[1,2] + dξdx[1,3]) + 
                nrm[2]*(dξdx[2,1] + dξdx[2,2] + dξdx[2,3]) +
                nrm[3]*(dξdx[3,1] + dξdx[3,2] + dξdx[3,3]))
    end
    # warm-up
    sbp = TetSBP{Float64}(degree=1)
    x = zeros(Float64, (3,sbp.numnodes,1))
    vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
    x[:,:,1] = calcnodes(sbp, vtx)
    dξdx = zeros(Float64, (3,3,sbp.numnodes,1))
    jac = zeros(Float64, (sbp.numnodes,1))
    mappingjacobian!(sbp, x, dξdx, jac)
    bndryfaces = Array(Boundary, 1)
    bndryfaces[1] = Boundary(1,1)
    u = ones(sbp.numnodes,1)
    res = zeros(u)
    print("\tIgnore this --->")
    @time boundaryintegrate!(sbp, bndryfaces, u, x, dξdx, bndryflux, res)

    numelem = 10000
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      x = zeros(Float64, (3,sbp.numnodes,numelem))
      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      for k = 1:numelem
        x[:,:,k] = calcnodes(sbp, vtx)
      end
      dξdx = zeros(Float64, (3,3,sbp.numnodes,numelem))
      jac = zeros(Float64, (sbp.numnodes,numelem))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, numelem)
      for k = 1:numelem
        bndryfaces[k] = Boundary(k,rand(1:3))
      end
      u = ones(sbp.numnodes, numelem)
      res = zeros(u)
      print("p = ",p,": ")
      @time boundaryintegrate!(sbp, bndryfaces, u, x, dξdx, bndryflux, res)
    end
  end

  context("Timing SummationByParts.boundaryintegrate! (TriSBP, vector field method)") do
    function bndryflux{T}(u::AbstractArray{T,1}, x::AbstractArray{T,1},
                          dξdx::AbstractArray{T,2}, nrm::AbstractArray{T,1},
                          flux::AbstractArray{T,1})
      tmp = (nrm[1]*(dξdx[1,1] + dξdx[1,2]) + nrm[2]*(dξdx[2,1] + dξdx[2,2]))
      for field = 1:size(u,1)
        flux[field] = u[field]*tmp
      end
    end
    # warm-up
    sbp = TriSBP{Float64}(degree=1)
    x = zeros(Float64, (2,sbp.numnodes,1))
    vtx = [0. 0.; 1. 0.; 0. 1.]
    x[:,:,1] = calcnodes(sbp, vtx)
    dξdx = zeros(Float64, (2,2,sbp.numnodes,1))
    jac = zeros(Float64, (sbp.numnodes,1))
    mappingjacobian!(sbp, x, dξdx, jac)
    bndryfaces = Array(Boundary, 1)
    bndryfaces[1] = Boundary(1,1)
    u = ones(5,sbp.numnodes,1)
    res = zeros(u)
    print("\tIgnore this --->")
    @time boundaryintegrate!(sbp, bndryfaces, u, x, dξdx, bndryflux, res)

    numelem = 10000
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      x = zeros(Float64, (2,sbp.numnodes,numelem))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      for k = 1:numelem
        x[:,:,k] = calcnodes(sbp, vtx)
      end
      dξdx = zeros(Float64, (2,2,sbp.numnodes,numelem))
      jac = zeros(Float64, (sbp.numnodes,numelem))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, numelem)
      for k = 1:numelem
        bndryfaces[k] = Boundary(k,rand(1:3))
      end
      u = ones(5, sbp.numnodes, numelem)
      res = zeros(u)
      print("p = ",p,": ")
      @time boundaryintegrate!(sbp, bndryfaces, u, x, dξdx, bndryflux, res)
    end
  end

  context("Timing SummationByParts.boundaryintegrate! (TetSBP, vector field method)") do
    function bndryflux{T}(u::AbstractArray{T,1}, x::AbstractArray{T,1},
                          dξdx::AbstractArray{T,2}, nrm::AbstractArray{T,1},
                          flux::AbstractArray{T,1})
      flux = u*(nrm[1]*(dξdx[1,1] + dξdx[1,2] + dξdx[1,3]) + 
                nrm[2]*(dξdx[2,1] + dξdx[2,2] + dξdx[2,3]) +
                nrm[3]*(dξdx[3,1] + dξdx[3,2] + dξdx[3,3]))
    end
    # warm-up
    sbp = TetSBP{Float64}(degree=1)
    x = zeros(Float64, (3,sbp.numnodes,1))
    vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
    x[:,:,1] = calcnodes(sbp, vtx)
    dξdx = zeros(Float64, (3,3,sbp.numnodes,1))
    jac = zeros(Float64, (sbp.numnodes,1))
    mappingjacobian!(sbp, x, dξdx, jac)
    bndryfaces = Array(Boundary, 1)
    bndryfaces[1] = Boundary(1,1)
    u = ones(5,sbp.numnodes,1)
    res = zeros(u)
    print("\tIgnore this --->")
    @time boundaryintegrate!(sbp, bndryfaces, u, x, dξdx, bndryflux, res)

    numelem = 10000
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      x = zeros(Float64, (3,sbp.numnodes,numelem))
      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      for k = 1:numelem
        x[:,:,k] = calcnodes(sbp, vtx)
      end
      dξdx = zeros(Float64, (3,3,sbp.numnodes,numelem))
      jac = zeros(Float64, (sbp.numnodes,numelem))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, numelem)
      for k = 1:numelem
        bndryfaces[k] = Boundary(k,rand(1:3))
      end
      u = ones(5, sbp.numnodes, numelem)
      res = zeros(u)
      print("p = ",p,": ")
      @time boundaryintegrate!(sbp, bndryfaces, u, x, dξdx, bndryflux, res)
    end
  end

end

  context("Timing SummationByParts.edgestabilize! (TriSBP, scalar field method)") do
    function stabscale{T}(u::T, dξdx::AbstractArray{T,2}, nrm::AbstractArray{T,1})
      return one(T)
    end
    # warm-up
    sbp = TriSBP{Float64}(degree=1)
    x = zeros(Float64, (2,sbp.numnodes,1))
    vtx = [0. 0.; 1. 0.; 0. 1.]
    x[:,:,1] = calcnodes(sbp, vtx)
    dξdx = zeros(Float64, (2,2,sbp.numnodes,1))
    jac = zeros(Float64, (sbp.numnodes,1))
    mappingjacobian!(sbp, x, dξdx, jac)
    α = zeros(dξdx)
    for i = 1:sbp.numnodes
      for di1 = 1:2
        for di2 = 1:2
          α[di1,di2,i,1] = (dξdx[di1,1,i,1].*dξdx[di2,1,i,1] + 
                            dξdx[di1,2,i,1].*dξdx[di2,2,i,1])*jac[i,1]
        end
      end
    end
    ifaces = Array(Interface, 1)
    ifaces[1] = Interface(1,1,1,2)
    u = ones(sbp.numnodes, 1)
    res = zeros(u)
    print("\tIgnore this --->")
    @time edgestabilize!(sbp, ifaces, u, x, dξdx, jac, α, stabscale, res)

    numelem = 10000
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      x = zeros(Float64, (2,sbp.numnodes,numelem))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      for k = 1:numelem
        x[:,:,k] = calcnodes(sbp, vtx)
      end
      dξdx = zeros(Float64, (2,2,sbp.numnodes,numelem))
      jac = zeros(Float64, (sbp.numnodes,numelem))
      mappingjacobian!(sbp, x, dξdx, jac)
      α = zeros(dξdx)
      for k = 1:2
        for i = 1:sbp.numnodes
          for di1 = 1:2
            for di2 = 1:2
              α[di1,di2,i,k] = (dξdx[di1,1,i,k].*dξdx[di2,1,i,k] + 
                                dξdx[di1,2,i,k].*dξdx[di2,2,i,k])*jac[i,k]
            end
          end
        end
      end
      ifaces = Array(Interface, numelem-1)
      for k = 1:numelem-1
        ifaces[k] = Interface(k,k+1,rand(1:3),rand(1:3))
      end
      u = ones(sbp.numnodes, numelem)
      res = zeros(u)
      print("p = ",p,": ")
      @time edgestabilize!(sbp, ifaces, u, x, dξdx, jac, α, stabscale, res)
    end
  end

  context("Timing SummationByParts.edgestabilize! (TriSBP, vector field method)") do
    function stabscale{T}(u::AbstractArray{T,1}, dξdx::AbstractArray{T,2},
                          nrm::AbstractArray{T,1})
      return one(T)
    end
    # warm-up
    sbp = TriSBP{Float64}(degree=1)
    x = zeros(Float64, (2,sbp.numnodes,1))
    vtx = [0. 0.; 1. 0.; 0. 1.]
    x[:,:,1] = calcnodes(sbp, vtx)
    dξdx = zeros(Float64, (2,2,sbp.numnodes,1))
    jac = zeros(Float64, (sbp.numnodes,1))
    mappingjacobian!(sbp, x, dξdx, jac)
    α = zeros(dξdx)
    for i = 1:sbp.numnodes
      for di1 = 1:2
        for di2 = 1:2
          α[di1,di2,i,1] = (dξdx[di1,1,i,1].*dξdx[di2,1,i,1] + 
                            dξdx[di1,2,i,1].*dξdx[di2,2,i,1])*jac[i,1]
        end
      end
    end
    ifaces = Array(Interface, 1)
    ifaces[1] = Interface(1,1,1,2)
    u = ones(5, sbp.numnodes, 1)
    res = zeros(u)
    print("\tIgnore this --->")
    @time edgestabilize!(sbp, ifaces, u, x, dξdx, jac, α, stabscale, res)

    numelem = 10000
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      x = zeros(Float64, (2,sbp.numnodes,numelem))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      for k = 1:numelem
        x[:,:,k] = calcnodes(sbp, vtx)
      end
      dξdx = zeros(Float64, (2,2,sbp.numnodes,numelem))
      jac = zeros(Float64, (sbp.numnodes,numelem))
      mappingjacobian!(sbp, x, dξdx, jac)
      α = zeros(dξdx)
      for k = 1:2
        for i = 1:sbp.numnodes
          for di1 = 1:2
            for di2 = 1:2
              α[di1,di2,i,k] = (dξdx[di1,1,i,k].*dξdx[di2,1,i,k] + 
                                dξdx[di1,2,i,k].*dξdx[di2,2,i,k])*jac[i,k]
            end
          end
        end
      end
      ifaces = Array(Interface, numelem-1)
      for k = 1:numelem-1
        ifaces[k] = Interface(k,k+1,rand(1:3),rand(1:3))
      end
      u = ones(5, sbp.numnodes, numelem)
      res = zeros(u)
      print("p = ",p,": ")
      @time edgestabilize!(sbp, ifaces, u, x, dξdx, jac, α, stabscale, res)
    end
  end

context("Timing SummationByParts.getnbrnodeindex") do
  sbp = TriSBP{Float64}(degree=1)
  iface = Interface(1,1,1,2)
  i = rand(1:sbp.numfacenodes)
  print("\tIgnore this --->")
  @time SummationByParts.getnbrnodeindex(sbp, iface, i)    
  for p = 1:4
    sbp = TriSBP{Float64}(degree=p)
    iface = Interface(1,1,1,2)
    i = rand(1:sbp.numfacenodes)
    print("p = ",p,": ")
    @time SummationByParts.getnbrnodeindex(sbp, iface, i)    
  end
end

end