using SummationByParts
using SummationByParts.OrthoPoly
using SummationByParts.Cubature
using SummationByParts.SymCubatures
using FactCheck
using ArrayViews

facts("Checking run time of various operations...") do

if true
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

  # context("Timing SummationByParts.weakdifferentiate! (TriSBP, scalar, op)") do
  #   # warm-up 
  #   sbp = TriSBP{Float64}(degree=1)
  #   u = ones(Float64, (sbp.numnodes,2))
  #   di = 1
  #   res = zeros(u)
  #   print("\tIgnore this --->")
  #   @time weakdifferentiate!(sbp, di, u, res, op=:-, trans=true)
  #   numelem = 10000
  #   for p = 1:4
  #     sbp = TriSBP{Float64}(degree=p)
  #     u = ones(Float64, (sbp.numnodes,numelem))
  #     di = 1
  #     res = zeros(u)
  #     print("p = ",p,": ")
  #     @time weakdifferentiate!(sbp, di, u, res, op=:-, trans=true)
  #   end
  #   println()
  # end

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
        Ddir = directionaldifferentiate!(sbp, dir, view(u,:,k), i)
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
        directionaldifferentiate!(sbp, dir, view(u,:,:,k), i, Ddir)
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
        Ddir = directionaldifferentiate!(sbp, dir, view(u,:,k), i)
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
        directionaldifferentiate!(sbp, dir, view(u,:,:,k), i, Ddir)
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

end

end


