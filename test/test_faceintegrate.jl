@testset "Testing SummationByParts Module (face-data integration methods)..." begin

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing integratefunctional! ($(string($TSBP)) scalar field method)" begin
        # build single element grid and verify the accuracy of boundary integration
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          x = zeros(Float64, (1,sbp.numnodes,1))
          xf = zeros(Float64, (1,sbpface.numnodes,2))
          vtx = reshape([-1.0; 1.0], (2,1))
          x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[1;]],(1,1)))
          xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          bndryfaces = Array{Boundary}(undef, 2)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(1,2)
          uface = zeros(Float64, (sbpface.numnodes, 2))
          for i = 0:p
            uface[:,:] = xf[1,:,:].^i
            fun = integratefunctional!(sbpface, bndryfaces, uface)            
            funexact = (1.0^(i) + (-1.0)^(i)) # sign normal is not accounted for
            @test ≈(fun, funexact, atol=1e-12)
          end
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing integratefunctional! ($(string($TSBP)) vector field method)" begin
        # build single element grid and verify the accuracy of boundary integration
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          x = zeros(Float64, (1,sbp.numnodes,1))
          xf = zeros(Float64, (1,sbpface.numnodes,2))
          vtx = reshape([-1.0; 1.0], (2,1))
          x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[1;]],(1,1)))
          xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          bndryfaces = Array{Boundary}(undef, 2)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(1,2)
          uface = zeros(Float64, (2, sbpface.numnodes, 2))
          fun = zeros(2)
          for i = 0:p
            uface[1,:,:] = xf[1,:,:].^i
            uface[2,:,:] .= 1.0
            fill!(fun, 0.0)
            integratefunctional!(sbpface, bndryfaces, uface, fun)
            funexact = (1.0^(i) + (-1.0)^(i)) # sign normal is not accounted for
            @test ≈(fun[1], funexact, atol=1e-12)
            @test ≈(fun[2], 2.0, atol=1e-12)
          end
        end
      end
    end
  end
      
  @testset "Testing integratefunctional! (TriSBP, scalar field method)" begin
    # build a two element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      #x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = zeros(Float64, (sbpface.numnodes, 4))
      for d = 0:2*p
        for j = 0:d
          i = d-j
          # the function be integrated is (x+1)^i (y+1)^j
          uface[:,:] = ((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j)
          uface .*= 0.5 #scale!(uface, 0.5) # 0.5 factor accounts for tranformation to ref space
          fun = integratefunctional!(sbpface, bndryfaces, uface)
          funexact = (2^(i+1)-1)*(1 + 2^j)/(i+1) + (2^(j+1)-1)*(1 + 2^i)/(j+1)
          @test ≈(fun, funexact, atol=1e-12)
        end
      end
    end
  end

  @testset "Testing integratefunctional! (TriSBP, vector field method)" begin
    # build a two element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      #x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = zeros(Float64, (2, sbpface.numnodes, 4))
      fun = zeros(2)
      for d = 0:2*p
        for j = 0:d
          i = d-j
          # the function be integrated is (x+1)^i (y+1)^j; 0.5 factor accounts
          # for the transformation to ref space
          uface[1,:,:] = 0.5*((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j)
          uface[2,:,:] .= 0.5 # integrate constant, gives perimeter
          fill!(fun, 0.0)
          integratefunctional!(sbpface, bndryfaces, uface, fun)
          funexact = (2^(i+1)-1)*(1 + 2^j)/(i+1) + (2^(j+1)-1)*(1 + 2^i)/(j+1)
          @test ≈(fun[1], funexact, atol=1e-12)
          @test ≈(fun[2], 4.0, atol=1e-12)
        end
      end
    end
  end

  @testset "Testing integratefunctional! (TriSparseFace, scalar field method)" begin
    # build a two element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      #x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = zeros(Float64, (sbpface.numnodes, 4))
      for d = 0:2*p
        for j = 0:d
          i = d-j
          # the function be integrated is (x+1)^i (y+1)^j
          uface[:,:] = ((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j)
          uface .*= 0.5 #scale!(uface, 0.5) # 0.5 factor accounts for tranformation to ref space
          fun = integratefunctional!(sbpface, bndryfaces, uface)
          funexact = (2^(i+1)-1)*(1 + 2^j)/(i+1) + (2^(j+1)-1)*(1 + 2^i)/(j+1)
          @test ≈(fun, funexact, atol=1e-12)
        end
      end
    end
  end

  @testset "Testing integratefunctional! (TriSparseFace, vector field method)" begin
    # build a two element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      #x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = zeros(Float64, (2, sbpface.numnodes, 4))
      fun = zeros(2)
      for d = 0:2*p
        for j = 0:d
          i = d-j
          # the function be integrated is (x+1)^i (y+1)^j; 0.5 factor accounts
          # for the transformation to ref space
          uface[1,:,:] = 0.5*((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j)
          uface[2,:,:] .= 0.5 # integrate constant, gives perimeter
          fill!(fun, 0.0)
          integratefunctional!(sbpface, bndryfaces, uface, fun)
          funexact = (2^(i+1)-1)*(1 + 2^j)/(i+1) + (2^(j+1)-1)*(1 + 2^i)/(j+1)
          @test ≈(fun[1], funexact, atol=1e-12)
          @test ≈(fun[2], 4.0, atol=1e-12)
        end
      end
    end
  end
  
  @testset "Testing integratefunctional! (TetSBP, scalar field method)" begin
    # build a four element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      uface = zeros(Float64, (sbpface.numnodes, 12))
      for d = 0:2*p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            # the function be integrated is (x+1)^i (y+1)^j (z+1)^k
            uface[:,:] = ((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j).*((xf[3,:,:].+1).^k)
            uface .*= 0.25 #scale!(uface, 0.25) # 0.25 factor accounts for tranformation to ref space
            fun = integratefunctional!(sbpface, bndryfaces, uface)
            funexact = (2^(j+1)-1)*(2^(k+1)-1)*(1 + 2^i)/((j+1)*(k+1)) +
            (2^(i+1)-1)*(2^(k+1)-1)*(1 + 2^j)/((i+1)*(k+1)) + 
            (2^(i+1)-1)*(2^(j+1)-1)*(1 + 2^k)/((i+1)*(j+1))
            @test ≈(fun, funexact, atol=1e-12)
          end
        end
      end
    end
  end

  @testset "Testing integratefunctional! (TetSBP, vector field method)" begin
    # build a four element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      uface = zeros(Float64, (2, sbpface.numnodes, 12))
      fun = zeros(2)
      for d = 0:2*p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            # the function be integrated is (x+1)^i (y+1)^j (z+1)^k
            uface[1,:,:] = ((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j).*((xf[3,:,:].+1).^k)
            uface[2,:,:] .= 1.0 # integrate constant, gives surface area
            uface .*= 0.25 #scale!(uface, 0.25) # 0.25 factor accounts for tranformation to ref space
            fill!(fun, 0.0)
            integratefunctional!(sbpface, bndryfaces, uface, fun)
            funexact = (2^(j+1)-1)*(2^(k+1)-1)*(1 + 2^i)/((j+1)*(k+1)) +
            (2^(i+1)-1)*(2^(k+1)-1)*(1 + 2^j)/((i+1)*(k+1)) + 
            (2^(i+1)-1)*(2^(j+1)-1)*(1 + 2^k)/((i+1)*(j+1))
            @test ≈(fun[1], funexact, atol=1e-12)
            @test ≈(fun[2], 6.0, atol=1e-12)
          end
        end
      end
    end
  end

  @testset "Testing integratefunctional! (TetSparseFace, scalar field method)" begin
    # build a four element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      uface = zeros(Float64, (sbpface.numnodes, 12))
      for d = 0:2*p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            # the function be integrated is (x+1)^i (y+1)^j (z+1)^k
            uface[:,:] = ((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j).*((xf[3,:,:].+1).^k)
            uface .*= 0.25 #scale!(uface, 0.25) # 0.25 factor accounts for tranformation to ref space
            fun = integratefunctional!(sbpface, bndryfaces, uface)
            funexact = (2^(j+1)-1)*(2^(k+1)-1)*(1 + 2^i)/((j+1)*(k+1)) +
            (2^(i+1)-1)*(2^(k+1)-1)*(1 + 2^j)/((i+1)*(k+1)) + 
            (2^(i+1)-1)*(2^(j+1)-1)*(1 + 2^k)/((i+1)*(j+1))
            @test ≈(fun, funexact, atol=1e-12)
          end
        end
      end
    end
  end

  @testset "Testing integratefunctional! (TetSparseFace, vector field method)" begin
    # build a four element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      uface = zeros(Float64, (2, sbpface.numnodes, 12))
      fun = zeros(2)
      for d = 0:2*p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            # the function be integrated is (x+1)^i (y+1)^j (z+1)^k
            uface[1,:,:] = ((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j).*((xf[3,:,:].+1).^k)
            uface[2,:,:] .= 1.0 # integrate constant, gives surface area
            uface .*= 0.25 #scale!(uface, 0.25) # 0.25 factor accounts for tranformation to ref space
            fill!(fun, 0.0)
            integratefunctional!(sbpface, bndryfaces, uface, fun)
            funexact = (2^(j+1)-1)*(2^(k+1)-1)*(1 + 2^i)/((j+1)*(k+1)) +
            (2^(i+1)-1)*(2^(k+1)-1)*(1 + 2^j)/((i+1)*(k+1)) + 
            (2^(i+1)-1)*(2^(j+1)-1)*(1 + 2^k)/((i+1)*(j+1))
            @test ≈(fun[1], funexact, atol=1e-12)
            @test ≈(fun[2], 6.0, atol=1e-12)
          end
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing integrateBoundaryFunctional! ($(string($TSBP)) scalar field method)" begin
        # build single element grid and verify the accuracy of boundary integration
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          x = zeros(Float64, (1,sbp.numnodes,1))
          xf = zeros(Float64, (1,sbpface.numnodes,2))
          vtx = reshape([-1.0; 1.0], (2,1))
          x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[1;]],(1,1)))
          xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          bndryfaces = Array{Boundary}(undef, 2)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(1,2)
          uface = zeros(Float64, (sbpface.numnodes, 2))
          for i = 0:p
            uface[:,:] = xf[1,:,:].^i
            fun = 0.0
            for (bindex, bndry) in enumerate(bndryfaces)
              fun += integrateBoundaryFunctional!(sbpface, bndry.face,
                                                  view(uface,:,bindex))
            end
            funexact = (1.0^(i) + (-1.0)^(i)) # sign normal is not accounted for
            @test ≈(fun, funexact, atol=1e-12)
          end
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing integrateBoundaryFunctional! ($(string($TSBP)) vector field method)" begin
        # build single element grid and verify the accuracy of boundary integration
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          x = zeros(Float64, (1,sbp.numnodes,1))
          xf = zeros(Float64, (1,sbpface.numnodes,2))
          vtx = reshape([-1.0; 1.0], (2,1))
          x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[1;]],(1,1)))
          xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          bndryfaces = Array{Boundary}(undef, 2)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(1,2)
          uface = zeros(Float64, (2, sbpface.numnodes, 2))
          fun = zeros(2)
          for i = 0:p
            uface[1,:,:] = xf[1,:,:].^i
            uface[2,:,:] .= 1.0
            fill!(fun, 0.0)
            for (bindex, bndry) in enumerate(bndryfaces)          
              integrateBoundaryFunctional!(sbpface, bndry.face,
                                           view(uface,:,:,bindex), fun)
            end
            funexact = (1.0^(i) + (-1.0)^(i)) # sign normal is not accounted for
            @test ≈(fun[1], funexact, atol=1e-12)
            @test ≈(fun[2], 2.0, atol=1e-12)
          end
        end
      end
    end
  end

  @testset "Testing integrateBoundaryFunctional! (TriSBP, scalar field method)" begin
    # build a two element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      #x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = zeros(Float64, (sbpface.numnodes, 4))
      for d = 0:2*p
        for j = 0:d
          i = d-j
          # the function be integrated is (x+1)^i (y+1)^j
          uface[:,:] = ((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j)
          uface .*= 0.5 #scale!(uface, 0.5) # 0.5 factor accounts for tranformation to ref space
          fun = 0.0
          for (bindex, bndry) in enumerate(bndryfaces)
            fun += integrateBoundaryFunctional!(sbpface, bndry.face,
                                                view(uface,:,bindex))
          end
          funexact = (2^(i+1)-1)*(1 + 2^j)/(i+1) + (2^(j+1)-1)*(1 + 2^i)/(j+1)
          @test ≈(fun, funexact, atol=1e-12)
        end
      end
    end
  end

  @testset "Testing integrateBoundaryFunctional! (TriSBP, vector field method)" begin
    # build a two element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      #x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = zeros(Float64, (2, sbpface.numnodes, 4))
      fun = zeros(2)
      for d = 0:2*p
        for j = 0:d
          i = d-j
          # the function be integrated is (x+1)^i (y+1)^j; 0.5 factor accounts
          # for the transformation to ref space
          uface[1,:,:] = 0.5*((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j)
          uface[2,:,:] .= 0.5 # integrate constant, gives perimeter
          fill!(fun, 0.0)
          for (bindex, bndry) in enumerate(bndryfaces)          
            integrateBoundaryFunctional!(sbpface, bndry.face,
                                         view(uface,:,:,bindex), fun)
          end
          funexact = (2^(i+1)-1)*(1 + 2^j)/(i+1) + (2^(j+1)-1)*(1 + 2^i)/(j+1)
          @test ≈(fun[1], funexact, atol=1e-12)
          @test ≈(fun[2], 4.0, atol=1e-12)
        end
      end
    end
  end

  @testset "Testing integrateBoundaryFunctional! (TriSparseFace, scalar field method)" begin
    # build a two element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      #x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = zeros(Float64, (sbpface.numnodes, 4))
      for d = 0:2*p
        for j = 0:d
          i = d-j
          # the function be integrated is (x+1)^i (y+1)^j
          uface[:,:] = ((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j)
          uface .*= 0.5 #scale!(uface, 0.5) # 0.5 factor accounts for tranformation to ref space
          fun = 0.0
          for (bindex, bndry) in enumerate(bndryfaces)
            fun += integrateBoundaryFunctional!(sbpface, bndry.face,
                                                view(uface,:,bindex))
          end
          funexact = (2^(i+1)-1)*(1 + 2^j)/(i+1) + (2^(j+1)-1)*(1 + 2^i)/(j+1)
          @test ≈(fun, funexact, atol=1e-12)
        end
      end
    end
  end

  @testset "Testing integrateBoundaryFunctional! (TriSparseFace, vector field method)" begin
    # build a two element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      #x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = zeros(Float64, (2, sbpface.numnodes, 4))
      fun = zeros(2)
      for d = 0:2*p
        for j = 0:d
          i = d-j
          # the function be integrated is (x+1)^i (y+1)^j; 0.5 factor accounts
          # for the transformation to ref space
          uface[1,:,:] = 0.5*((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j)
          uface[2,:,:] .= 0.5 # integrate constant, gives perimeter
          fill!(fun, 0.0)
          for (bindex, bndry) in enumerate(bndryfaces)          
            integrateBoundaryFunctional!(sbpface, bndry.face,
                                         view(uface,:,:,bindex), fun)
          end
          funexact = (2^(i+1)-1)*(1 + 2^j)/(i+1) + (2^(j+1)-1)*(1 + 2^i)/(j+1)
          @test ≈(fun[1], funexact, atol=1e-12)
          @test ≈(fun[2], 4.0, atol=1e-12)
        end
      end
    end
  end

  @testset "Testing integrateBoundaryFunctional! (TetSBP, scalar field method)" begin
    # build a four element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      uface = zeros(Float64, (sbpface.numnodes, 12))
      for d = 0:2*p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            # the function be integrated is (x+1)^i (y+1)^j (z+1)^k
            uface[:,:] = ((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j).*((xf[3,:,:].+1).^k)
            uface .*= 0.25 #scale!(uface, 0.25) # 0.25 factor accounts for tranformation to ref space
            fun = 0.0
            for (bindex, bndry) in enumerate(bndryfaces)
              fun += integrateBoundaryFunctional!(sbpface, bndry.face,
                                                  view(uface,:,bindex))
            end
            funexact = (2^(j+1)-1)*(2^(k+1)-1)*(1 + 2^i)/((j+1)*(k+1)) +
            (2^(i+1)-1)*(2^(k+1)-1)*(1 + 2^j)/((i+1)*(k+1)) + 
            (2^(i+1)-1)*(2^(j+1)-1)*(1 + 2^k)/((i+1)*(j+1))
            @test ≈(fun, funexact, atol=1e-12)
          end
        end
      end
    end
  end

  @testset "Testing integrateBoundaryFunctional! (TetSBP, vector field method)" begin
    # build a four element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      uface = zeros(Float64, (2, sbpface.numnodes, 12))
      fun = zeros(2)
      for d = 0:2*p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            # the function be integrated is (x+1)^i (y+1)^j (z+1)^k
            uface[1,:,:] = ((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j).*((xf[3,:,:].+1).^k)
            uface[2,:,:] .= 1.0 # integrate constant, gives surface area
            uface .*= 0.25 #scale!(uface, 0.25) # 0.25 factor accounts for tranformation to ref space
            fill!(fun, 0.0)
            for (bindex, bndry) in enumerate(bndryfaces)
              integrateBoundaryFunctional!(sbpface, bndry.face,
                                           view(uface,:,:,bindex), fun)
            end
            funexact = (2^(j+1)-1)*(2^(k+1)-1)*(1 + 2^i)/((j+1)*(k+1)) +
            (2^(i+1)-1)*(2^(k+1)-1)*(1 + 2^j)/((i+1)*(k+1)) + 
            (2^(i+1)-1)*(2^(j+1)-1)*(1 + 2^k)/((i+1)*(j+1))
            @test ≈(fun[1], funexact, atol=1e-12)
            @test ≈(fun[2], 6.0, atol=1e-12)
          end
        end
      end
    end
  end

  @testset "Testing integrateBoundaryFunctional! (TetSparseFace, scalar field method)" begin
    # build a four element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      uface = zeros(Float64, (sbpface.numnodes, 12))
      for d = 0:2*p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            # the function be integrated is (x+1)^i (y+1)^j (z+1)^k
            uface[:,:] = ((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j).*((xf[3,:,:].+1).^k)
            uface .*= 0.25 #scale!(uface, 0.25) # 0.25 factor accounts for tranformation to ref space
            fun = 0.0
            for (bindex, bndry) in enumerate(bndryfaces)
              fun += integrateBoundaryFunctional!(sbpface, bndry.face,
                                                  view(uface,:,bindex))
            end
            funexact = (2^(j+1)-1)*(2^(k+1)-1)*(1 + 2^i)/((j+1)*(k+1)) +
            (2^(i+1)-1)*(2^(k+1)-1)*(1 + 2^j)/((i+1)*(k+1)) + 
            (2^(i+1)-1)*(2^(j+1)-1)*(1 + 2^k)/((i+1)*(j+1))
            @test ≈(fun, funexact, atol=1e-12)
          end
        end
      end
    end
  end

  @testset "Testing integrateBoundaryFunctional! (TetSparseFace, vector field method)" begin
    # build a four element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      uface = zeros(Float64, (2, sbpface.numnodes, 12))
      fun = zeros(2)
      for d = 0:2*p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            # the function be integrated is (x+1)^i (y+1)^j (z+1)^k
            uface[1,:,:] = ((xf[1,:,:].+1).^i).*((xf[2,:,:].+1).^j).*((xf[3,:,:].+1).^k)
            uface[2,:,:] .= 1.0 # integrate constant, gives surface area
            uface .*= 0.25  #scale!(uface, 0.25) # 0.25 factor accounts for tranformation to ref space
            fill!(fun, 0.0)
            for (bindex, bndry) in enumerate(bndryfaces)
              integrateBoundaryFunctional!(sbpface, bndry.face,
                                           view(uface,:,:,bindex), fun)
            end
            funexact = (2^(j+1)-1)*(2^(k+1)-1)*(1 + 2^i)/((j+1)*(k+1)) +
            (2^(i+1)-1)*(2^(k+1)-1)*(1 + 2^j)/((i+1)*(k+1)) + 
            (2^(i+1)-1)*(2^(j+1)-1)*(1 + 2^k)/((i+1)*(j+1))
            @test ≈(fun[1], funexact, atol=1e-12)
            @test ≈(fun[2], 6.0, atol=1e-12)
          end
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing boundaryintegrate! and interiorfaceintegrate! ($(string($TSBP)) scalar field method)" begin
        # build a two element grid and verify that a constant integrated over all
        # faces is zero
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          ifaces = Array{Interface}(undef, 1)
          ifaces[1] = Interface(1,2,2,2,1)
          bndryfaces = Array{Boundary}(undef, 2)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(2,1)
          uface = ones(Float64, (sbpface.numnodes, 1))
          ubndry = ones(Float64, (sbpface.numnodes, 2))
          ubndry[:,1] *= -1.0
          res = zeros(Float64, (sbp.numnodes, 2))
          boundaryintegrate!(sbpface, bndryfaces, ubndry, res)
          interiorfaceintegrate!(sbpface, ifaces, uface, res)
          @test ≈(sum(res), 0.0, atol=1e-13)
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing boundaryintegrate! and interiorfaceintegrate! ($(string($TSBP)) vector field method)" begin
        # build a two element grid and verify that a constant integrated over all
        # faces is zero
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          ifaces = Array{Interface}(undef, 1)
          ifaces[1] = Interface(1,2,2,2,1)
          bndryfaces = Array{Boundary}(undef, 2)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(2,1)
          uface = ones(Float64, (2, sbpface.numnodes, 1))
          ubndry = ones(Float64, (2, sbpface.numnodes, 2))
          ubndry[:,:,1] *= -1.0
          res = zeros(Float64, (2, sbp.numnodes, 2))
          boundaryintegrate!(sbpface, bndryfaces, ubndry, res)
          interiorfaceintegrate!(sbpface, ifaces, uface, res)
          @test ≈(sum(res), 0.0, atol=1e-13)
        end        
      end
    end
  end

  @testset "Testing boundaryintegrate! and interiorfaceintegrate! (TriSBP scalar field method)" begin
    # build a two element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, [-1. -1.; 1. -1.; -1. 1.])
      ifaces = Array{Interface}(undef, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = 0.5.*ones(Float64, (sbpface.numnodes, 1))
      ubndry = 0.5.*ones(Float64, (sbpface.numnodes, 4))
      ubndry[:,1] *= -1.0
      ubndry[:,2] *= -1.0
      res = zeros(Float64, (sbp.numnodes, 2))
      boundaryintegrate!(sbpface, bndryfaces, ubndry, res)
      interiorfaceintegrate!(sbpface, ifaces, uface, res)
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

  @testset "Testing boundaryintegrate! and interiorfaceintegrate! (TriSBP vector field method)" begin
    # build a two element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, [-1. -1.; 1. -1.; -1. 1.])
      ifaces = Array{Interface}(undef, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = 0.5.*ones(Float64, (2, sbpface.numnodes, 1))
      ubndry = 0.5.*ones(Float64, (2, sbpface.numnodes, 4))
      ubndry[:,:,1] *= -1.0
      ubndry[:,:,2] *= -1.0
      res = zeros(Float64, (2, sbp.numnodes, 2))
      boundaryintegrate!(sbpface, bndryfaces, ubndry, res)
      interiorfaceintegrate!(sbpface, ifaces, uface, res)
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

  @testset "Testing boundaryintegrate! and interiorfaceintegrate! (TriSparseFace scalar field method)" begin
    # build a two element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      ifaces = Array{Interface}(undef, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = 0.5.*ones(Float64, (sbpface.numnodes, 1))
      ubndry = 0.5.*ones(Float64, (sbpface.numnodes, 4))
      ubndry[:,1] *= -1.0
      ubndry[:,2] *= -1.0
      res = zeros(Float64, (sbp.numnodes, 2))
      boundaryintegrate!(sbpface, bndryfaces, ubndry, res)
      interiorfaceintegrate!(sbpface, ifaces, uface, res)
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

  @testset "Testing boundaryintegrate! and interiorfaceintegrate! (TriSparseFace vector field method)" begin
    # build a two element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      ifaces = Array{Interface}(undef, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = 0.5.*ones(Float64, (2, sbpface.numnodes, 1))
      ubndry = 0.5.*ones(Float64, (2, sbpface.numnodes, 4))
      ubndry[:,:,1] *= -1.0
      ubndry[:,:,2] *= -1.0
      res = zeros(Float64, (2, sbp.numnodes, 2))
      boundaryintegrate!(sbpface, bndryfaces, ubndry, res)
      interiorfaceintegrate!(sbpface, ifaces, uface, res)
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

  @testset "Testing boundaryintegrate! and interiorfaceintegrate! (TetSBP scalar field method)" begin
    # build a five element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4      
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      ifaces = Array{Interface}(undef, 4)
      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      bndryfaces = Array{Boundary}(undef, 12)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)
      uface = ones(Float64, (sbpface.numnodes, 4))
      ubndry = ones(Float64, (sbpface.numnodes, 12))
      uface[:,2] *= -1.0
      uface[:,3] *= -1.0
      uface[:,4] *= -1.0
      ubndry[:,1:3] *= -1.0
      ubndry[:,4] *= -1.0
      ubndry[:,8] *= -1.0
      ubndry[:,10] *= -1.0
      res = zeros(Float64, (sbp.numnodes, 5))
      boundaryintegrate!(sbpface, bndryfaces, ubndry, res)
      interiorfaceintegrate!(sbpface, ifaces, uface, res)
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

  @testset "Testing boundaryintegrate! and interiorfaceintegrate! (TetSBP vector field method)" begin
    # build a five element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4      
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      ifaces = Array{Interface}(undef, 4)
      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      bndryfaces = Array{Boundary}(undef, 12)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)
      uface = ones(Float64, (2, sbpface.numnodes, 4))
      ubndry = ones(Float64, (2, sbpface.numnodes, 12))
      uface[:,:,2] *= -1.0
      uface[:,:,3] *= -1.0
      uface[:,:,4] *= -1.0
      ubndry[:,:,1:3] *= -1.0
      ubndry[:,:,4] *= -1.0
      ubndry[:,:,8] *= -1.0
      ubndry[:,:,10] *= -1.0
      res = zeros(Float64, (2, sbp.numnodes, 5))
      boundaryintegrate!(sbpface, bndryfaces, ubndry, res)
      interiorfaceintegrate!(sbpface, ifaces, uface, res)
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

  @testset "Testing boundaryintegrate! and interiorfaceintegrate! (TetSparseFace scalar field method)" begin
    # build a five element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      ifaces = Array{Interface}(undef, 4)
      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      bndryfaces = Array{Boundary}(undef, 12)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)
      uface = ones(Float64, (sbpface.numnodes, 4))
      ubndry = ones(Float64, (sbpface.numnodes, 12))
      uface[:,2] *= -1.0
      uface[:,3] *= -1.0
      uface[:,4] *= -1.0
      ubndry[:,1:3] *= -1.0
      ubndry[:,4] *= -1.0
      ubndry[:,8] *= -1.0
      ubndry[:,10] *= -1.0
      res = zeros(Float64, (sbp.numnodes, 5))
      boundaryintegrate!(sbpface, bndryfaces, ubndry, res)
      interiorfaceintegrate!(sbpface, ifaces, uface, res)
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

  @testset "Testing boundaryintegrate! and interiorfaceintegrate! (TetSparseFace vector field method)" begin
    # build a five element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      ifaces = Array{Interface}(undef, 4)
      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      bndryfaces = Array{Boundary}(undef, 12)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)
      uface = ones(Float64, (2, sbpface.numnodes, 4))
      ubndry = ones(Float64, (2, sbpface.numnodes, 12))
      uface[:,:,2] *= -1.0
      uface[:,:,3] *= -1.0
      uface[:,:,4] *= -1.0
      ubndry[:,:,1:3] *= -1.0
      ubndry[:,:,4] *= -1.0
      ubndry[:,:,8] *= -1.0
      ubndry[:,:,10] *= -1.0
      res = zeros(Float64, (2, sbp.numnodes, 5))
      boundaryintegrate!(sbpface, bndryfaces, ubndry, res)
      interiorfaceintegrate!(sbpface, ifaces, uface, res)
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing boundaryFaceIntegrate! and interiorFaceIntegrate! ($(string($TSBP)) scalar field method)" begin
        # build a two element grid and verify that a constant integrated over all
        # faces is zero
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          ifaces = Array{Interface}(undef, 1)
          ifaces[1] = Interface(1,2,2,2,1)
          bndryfaces = Array{Boundary}(undef, 2)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(2,1)
          uface = ones(Float64, (sbpface.numnodes, 1))
          ubndry = ones(Float64, (sbpface.numnodes, 2))
          ubndry[:,1] *= -1.0
          res = zeros(Float64, (sbp.numnodes, 2))
          for (bindex, bndry) in enumerate(bndryfaces)      
            boundaryFaceIntegrate!(sbpface, bndry.face, view(ubndry,:,bindex),
                                   view(res,:,bndry.element))
          end      
          interiorFaceIntegrate!(sbpface, ifaces[1], view(uface,:,1),
                                 view(res,:,ifaces[1].elementL),
                                 view(res,:,ifaces[1].elementR))
          @test ≈(sum(res), 0.0, atol=1e-13)
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing boundaryFaceIntegrate! and interiorFaceIntegrate! ($(string($TSBP)) vector field method)" begin
        # build a two element grid and verify that a constant integrated over all
        # faces is zero
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          ifaces = Array{Interface}(undef, 1)
          ifaces[1] = Interface(1,2,2,2,1)
          bndryfaces = Array{Boundary}(undef, 2)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(2,1)
          uface = ones(Float64, (2, sbpface.numnodes, 1))
          ubndry = ones(Float64, (2, sbpface.numnodes, 2))
          ubndry[:,:,1] *= -1.0
          res = zeros(Float64, (2, sbp.numnodes, 2))
          for (bindex, bndry) in enumerate(bndryfaces)
            boundaryFaceIntegrate!(sbpface, bndry.face, view(ubndry,:,:,bindex),
                                   view(res,:,:,bndry.element))
          end
          interiorFaceIntegrate!(sbpface, ifaces[1], view(uface,:,:,1),
                                 view(res,:,:,ifaces[1].elementL),
                                 view(res,:,:,ifaces[1].elementR))
          @test ≈(sum(res), 0.0, atol=1e-13)
        end        
      end
    end
  end

  @testset "Testing boundaryFaceIntegrate! and interiorFaceIntegrate! (TriSBP scalar field method)" begin
    # build a two element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, [-1. -1.; 1. -1.; -1. 1.])
      ifaces = Array{Interface}(undef, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = 0.5.*ones(Float64, (sbpface.numnodes, 1))
      ubndry = 0.5.*ones(Float64, (sbpface.numnodes, 4))
      ubndry[:,1] *= -1.0
      ubndry[:,2] *= -1.0
      res = zeros(Float64, (sbp.numnodes, 2))
      for (bindex, bndry) in enumerate(bndryfaces)      
        boundaryFaceIntegrate!(sbpface, bndry.face, view(ubndry,:,bindex),
                               view(res,:,bndry.element))
      end      
      interiorFaceIntegrate!(sbpface, ifaces[1], view(uface,:,1),
                             view(res,:,ifaces[1].elementL),
                             view(res,:,ifaces[1].elementR))
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

  @testset "Testing boundaryFaceIntegrate! and interiorFaceIntegrate! (TriSBP vector field method)" begin
    # build a two element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, [-1. -1.; 1. -1.; -1. 1.])
      ifaces = Array{Interface}(undef, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = 0.5.*ones(Float64, (2, sbpface.numnodes, 1))
      ubndry = 0.5.*ones(Float64, (2, sbpface.numnodes, 4))
      ubndry[:,:,1] *= -1.0
      ubndry[:,:,2] *= -1.0
      res = zeros(Float64, (2, sbp.numnodes, 2))
      for (bindex, bndry) in enumerate(bndryfaces)
        boundaryFaceIntegrate!(sbpface, bndry.face, view(ubndry,:,:,bindex),
                               view(res,:,:,bndry.element))
      end
      interiorFaceIntegrate!(sbpface, ifaces[1], view(uface,:,:,1),
                             view(res,:,:,ifaces[1].elementL),
                             view(res,:,:,ifaces[1].elementR))
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

  @testset "Testing boundaryFaceIntegrate! and interiorFaceIntegrate! (TriSparseFace scalar field method)" begin
    # build a two element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      ifaces = Array{Interface}(undef, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = 0.5.*ones(Float64, (sbpface.numnodes, 1))
      ubndry = 0.5.*ones(Float64, (sbpface.numnodes, 4))
      ubndry[:,1] *= -1.0
      ubndry[:,2] *= -1.0
      res = zeros(Float64, (sbp.numnodes, 2))
      for (bindex, bndry) in enumerate(bndryfaces)      
        boundaryFaceIntegrate!(sbpface, bndry.face, view(ubndry,:,bindex),
                               view(res,:,bndry.element))
      end      
      interiorFaceIntegrate!(sbpface, ifaces[1], view(uface,:,1),
                             view(res,:,ifaces[1].elementL),
                             view(res,:,ifaces[1].elementR))
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

  @testset "Testing boundaryFaceIntegrate! and interiorFaceIntegrate! (TriSparseFace vector field method)" begin
    # build a two element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      ifaces = Array{Interface}(undef, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = 0.5.*ones(Float64, (2, sbpface.numnodes, 1))
      ubndry = 0.5.*ones(Float64, (2, sbpface.numnodes, 4))
      ubndry[:,:,1] *= -1.0
      ubndry[:,:,2] *= -1.0
      res = zeros(Float64, (2, sbp.numnodes, 2))
      for (bindex, bndry) in enumerate(bndryfaces)
        boundaryFaceIntegrate!(sbpface, bndry.face, view(ubndry,:,:,bindex),
                               view(res,:,:,bndry.element))
      end
      interiorFaceIntegrate!(sbpface, ifaces[1], view(uface,:,:,1),
                             view(res,:,:,ifaces[1].elementL),
                             view(res,:,:,ifaces[1].elementR))
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

  @testset "Testing boundaryFaceIntegrate! and interiorFaceIntegrate! (TetSBP scalar field method)" begin
    # build a five element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4      
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      ifaces = Array{Interface}(undef, 4)
      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      bndryfaces = Array{Boundary}(undef, 12)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)
      uface = ones(Float64, (sbpface.numnodes, 4))
      ubndry = ones(Float64, (sbpface.numnodes, 12))
      uface[:,2] *= -1.0
      uface[:,3] *= -1.0
      uface[:,4] *= -1.0
      ubndry[:,1:3] *= -1.0
      ubndry[:,4] *= -1.0
      ubndry[:,8] *= -1.0
      ubndry[:,10] *= -1.0
      res = zeros(Float64, (sbp.numnodes, 5))
      for (bindex, bndry) in enumerate(bndryfaces)      
        boundaryFaceIntegrate!(sbpface, bndry.face, view(ubndry,:,bindex),
                               view(res,:,bndry.element))
      end
      for (findex, face) in enumerate(ifaces)
        interiorFaceIntegrate!(sbpface, face, view(uface,:,findex),
                               view(res,:,face.elementL),
                               view(res,:,face.elementR))
      end
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

  @testset "Testing boundaryFaceIntegrate! and interiorFaceIntegrate! (TetSBP vector field method)" begin
    # build a five element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4      
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      ifaces = Array{Interface}(undef, 4)
      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      bndryfaces = Array{Boundary}(undef, 12)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)
      uface = ones(Float64, (2, sbpface.numnodes, 4))
      ubndry = ones(Float64, (2, sbpface.numnodes, 12))
      uface[:,:,2] *= -1.0
      uface[:,:,3] *= -1.0
      uface[:,:,4] *= -1.0
      ubndry[:,:,1:3] *= -1.0
      ubndry[:,:,4] *= -1.0
      ubndry[:,:,8] *= -1.0
      ubndry[:,:,10] *= -1.0
      res = zeros(Float64, (2, sbp.numnodes, 5))
      for (bindex, bndry) in enumerate(bndryfaces)      
        boundaryFaceIntegrate!(sbpface, bndry.face, view(ubndry,:,:,bindex),
                               view(res,:,:,bndry.element))
      end
      for (findex, face) in enumerate(ifaces)
        interiorFaceIntegrate!(sbpface, face, view(uface,:,:,findex),
                               view(res,:,:,face.elementL),
                               view(res,:,:,face.elementR))
      end
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

  @testset "Testing boundaryFaceIntegrate! and interiorFaceIntegrate! (TetSparseFace scalar field method)" begin
    # build a five element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      ifaces = Array{Interface}(undef, 4)
      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      bndryfaces = Array{Boundary}(undef, 12)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)
      uface = ones(Float64, (sbpface.numnodes, 4))
      ubndry = ones(Float64, (sbpface.numnodes, 12))
      uface[:,2] *= -1.0
      uface[:,3] *= -1.0
      uface[:,4] *= -1.0
      ubndry[:,1:3] *= -1.0
      ubndry[:,4] *= -1.0
      ubndry[:,8] *= -1.0
      ubndry[:,10] *= -1.0
      res = zeros(Float64, (sbp.numnodes, 5))
      for (bindex, bndry) in enumerate(bndryfaces)      
        boundaryFaceIntegrate!(sbpface, bndry.face, view(ubndry,:,bindex),
                               view(res,:,bndry.element))
      end
      for (findex, face) in enumerate(ifaces)
        interiorFaceIntegrate!(sbpface, face, view(uface,:,findex),
                               view(res,:,face.elementL),
                               view(res,:,face.elementR))
      end
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

  @testset "Testing boundaryFaceIntegrate! and interiorFaceIntegrate! (TetSparseFace vector field method)" begin
    # build a five element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      ifaces = Array{Interface}(undef, 4)
      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      bndryfaces = Array{Boundary}(undef, 12)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)
      uface = ones(Float64, (2, sbpface.numnodes, 4))
      ubndry = ones(Float64, (2, sbpface.numnodes, 12))
      uface[:,:,2] *= -1.0
      uface[:,:,3] *= -1.0
      uface[:,:,4] *= -1.0
      ubndry[:,:,1:3] *= -1.0
      ubndry[:,:,4] *= -1.0
      ubndry[:,:,8] *= -1.0
      ubndry[:,:,10] *= -1.0
      res = zeros(Float64, (2, sbp.numnodes, 5))
      for (bindex, bndry) in enumerate(bndryfaces)      
        boundaryFaceIntegrate!(sbpface, bndry.face, view(ubndry,:,:,bindex),
                               view(res,:,:,bndry.element))
      end
      for (findex, face) in enumerate(ifaces)
        interiorFaceIntegrate!(sbpface, face, view(uface,:,:,findex),
                               view(res,:,:,face.elementL),
                               view(res,:,:,face.elementR))
      end
      @test ≈(sum(res), 0.0, atol=1e-13)
    end
  end

end
