@testset "Testing SummationByParts Module (face-data interpolation methods)..." begin

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing boundaryinterpolate! ($(string($TSBP)) scalar field method)" begin
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          x = zeros(Float64, (1,sbp.numnodes,2))
          xf = zeros(Float64, (1,sbpface.numnodes,4))
          vtx = reshape([0.0; 1.0], (2,1))
          x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[1;]],(1,1)))
          xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          vtx = reshape([2.0; 1.0], (2,1))  # note the reversal !!!
          x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[1;]],(1,1)))
          xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          bndryfaces = Array{Boundary}(undef, 4)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(1,2)
          bndryfaces[3] = Boundary(2,1)
          bndryfaces[4] = Boundary(2,2)
          u = zeros(Float64, (sbp.numnodes, 2))
          uface = zeros(Float64, (sbpface.numnodes, 4))
          for i = 0:p
            u[:,:] = x[1,:,:].^i
            boundaryinterpolate!(sbpface, bndryfaces, u, uface)
            @test ≈(vec(uface[:,:]), vec(xf[1,:,:].^i), atol=5e-14)
          end
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing boundaryinterpolate! ($(string($TSBP)) vector field method)" begin
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          x = zeros(Float64, (1,sbp.numnodes,2))
          xf = zeros(Float64, (1,sbpface.numnodes,4))
          vtx = reshape([0.0; 1.0], (2,1))
          x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[1;]],(1,1)))
          xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          vtx = reshape([2.0; 1.0], (2,1))  # note the reversal !!!
          x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[1;]],(1,1)))
          xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          bndryfaces = Array{Boundary}(undef, 4)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(1,2)
          bndryfaces[3] = Boundary(2,1)
          bndryfaces[4] = Boundary(2,2)
          u = zeros(Float64, (2,sbp.numnodes, 2))
          uface = zeros(Float64, (2,sbpface.numnodes, 4))
          for i = 0:p
            u[1,:,:] = x[1,:,:].^i
            u[2,:,:] = 2.0.*u[1,:,:]
            boundaryinterpolate!(sbpface, bndryfaces, u, uface)
            @test ≈(vec(uface[1,:,:]), vec(xf[1,:,:].^i), atol=5e-14)
            @test ≈(vec(uface[2,:,:]), 2.0.*vec(xf[1,:,:].^i), atol=5e-14)
          end
        end
      end
    end
  end
        
  @testset "Testing boundaryinterpolate! (TriSBP, scalar field method)" begin
    # build a two element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      u = zeros(Float64, (sbp.numnodes, 2))
      uface = zeros(Float64, (sbpface.numnodes, 4))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          boundaryinterpolate!(sbpface, bndryfaces, u, uface)
          @test ≈(vec(uface[:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=5e-14)
        end
      end
    end
  end

  @testset "Testing boundaryinterpolate! (TriSBP, vector field method)" begin
    # build a two element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      u = zeros(Float64, (2, sbp.numnodes, 2))
      uface = zeros(Float64, (2, sbpface.numnodes, 4))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j)
          boundaryinterpolate!(sbpface, bndryfaces, u, uface)
          @test ≈(vec(uface[1,:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=5e-14)
          @test ≈(vec(uface[2,:,:]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=5e-14)
        end
      end
    end
  end

  @testset "Testing boundaryinterpolate! (TriSparseFace, scalar field method)" begin
    # build a two element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      u = zeros(Float64, (sbp.numnodes, 2))
      uface = zeros(Float64, (sbpface.numnodes, 4))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          boundaryinterpolate!(sbpface, bndryfaces, u, uface)
          @test ≈(vec(uface[:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=5e-14)
        end
      end
    end
  end

  @testset "Testing boundaryinterpolate! (TriSparseSBP, vector field method)" begin
    # build a two element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      u = zeros(Float64, (2, sbp.numnodes, 2))
      uface = zeros(Float64, (2, sbpface.numnodes, 4))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j)
          boundaryinterpolate!(sbpface, bndryfaces, u, uface)
          @test ≈(vec(uface[1,:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=5e-14)
          @test ≈(vec(uface[2,:,:]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=5e-14)
        end
      end
    end
  end
  
  @testset "Testing boundaryinterpolate! (TetSBP, scalar field method)" begin
    # build a four element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (3,sbp.numnodes,4))
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      u = zeros(Float64, (sbp.numnodes, 4))
      uface = zeros(Float64, (sbpface.numnodes, 12))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            boundaryinterpolate!(sbpface, bndryfaces, u, uface)
            @test ≈(vec(uface[:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=5e-14)
          end
        end
      end
    end
  end

  @testset "Testing boundaryinterpolate! (TetSBP, vector field method)" begin
    # build a four element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (3,sbp.numnodes,4))
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      u = zeros(Float64, (2, sbp.numnodes, 4))
      uface = zeros(Float64, (2, sbpface.numnodes, 12))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            boundaryinterpolate!(sbpface, bndryfaces, u, uface)
            @test ≈(vec(uface[1,:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=1e-13)
            @test ≈(vec(uface[2,:,:]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=1e-13)
          end
        end
      end
    end
  end

  @testset "Testing boundaryinterpolate! (TetSparseFace, scalar field method)" begin
    # build a four element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      x = zeros(Float64, (3,sbp.numnodes,4))
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      u = zeros(Float64, (sbp.numnodes, 4))
      uface = zeros(Float64, (sbpface.numnodes, 12))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            boundaryinterpolate!(sbpface, bndryfaces, u, uface)
            @test ≈(vec(uface[:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=5e-14)
          end
        end
      end
    end
  end

  @testset "Testing boundaryinterpolate! (TetSparseFace, vector field method)" begin
    # build a four element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      x = zeros(Float64, (3,sbp.numnodes,4))
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      u = zeros(Float64, (2, sbp.numnodes, 4))
      uface = zeros(Float64, (2, sbpface.numnodes, 12))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            boundaryinterpolate!(sbpface, bndryfaces, u, uface)
            @test ≈(vec(uface[1,:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=5e-14)
            @test ≈(vec(uface[2,:,:]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=5e-13)
          end
        end
      end
    end
  end


  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing boundaryFaceInterpolate! ($(string($TSBP)) scalar field method)" begin
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          x = zeros(Float64, (1,sbp.numnodes,2))
          xf = zeros(Float64, (1,sbpface.numnodes,4))
          vtx = reshape([0.0; 1.0], (2,1))
          x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[1;]],(1,1)))
          xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          vtx = reshape([2.0; 1.0], (2,1))  # note the reversal !!!
          x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[1;]],(1,1)))
          xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          bndryfaces = Array{Boundary}(undef, 4)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(1,2)
          bndryfaces[3] = Boundary(2,1)
          bndryfaces[4] = Boundary(2,2)
          u = zeros(Float64, (sbp.numnodes, 2))
          uface = zeros(Float64, (sbpface.numnodes, 4))
          for i = 0:p
            u[:,:] = x[1,:,:].^i
            for (bindex, bndry) in enumerate(bndryfaces)          
              boundaryFaceInterpolate!(sbpface, bndry.face,
                                       view(u,:,bndry.element),
                                       view(uface,:,bindex))
            end
            @test ≈(vec(uface[:,:]), vec(xf[1,:,:].^i), atol=5e-14)
          end
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing boundaryFaceInterpolate! ($(string($TSBP)) vector field method)" begin
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          x = zeros(Float64, (1,sbp.numnodes,2))
          xf = zeros(Float64, (1,sbpface.numnodes,4))
          vtx = reshape([0.0; 1.0], (2,1))
          x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[1;]],(1,1)))
          xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          vtx = reshape([2.0; 1.0], (2,1))  # note the reversal !!!
          x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[1;]],(1,1)))
          xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          bndryfaces = Array{Boundary}(undef, 4)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(1,2)
          bndryfaces[3] = Boundary(2,1)
          bndryfaces[4] = Boundary(2,2)
          u = zeros(Float64, (2,sbp.numnodes, 2))
          uface = zeros(Float64, (2,sbpface.numnodes, 4))
          for i = 0:p
            u[1,:,:] = x[1,:,:].^i
            u[2,:,:] = 2.0.*u[1,:,:]
            for (bindex, bndry) in enumerate(bndryfaces)
              boundaryFaceInterpolate!(sbpface, bndry.face,
                                       view(u,:,:,bndry.element),
                                       view(uface,:,:,bindex))
            end
            @test ≈(vec(uface[1,:,:]), vec(xf[1,:,:].^i), atol=5e-14)
            @test ≈(vec(uface[2,:,:]), 2.0.*vec(xf[1,:,:].^i), atol=5e-14)
          end
        end
      end
    end
  end

  @testset "Testing boundaryFaceInterpolate! (TriSBP, scalar field method)" begin
    # build a two element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      u = zeros(Float64, (sbp.numnodes, 2))
      uface = zeros(Float64, (sbpface.numnodes, 4))      
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          for (bindex, bndry) in enumerate(bndryfaces)          
            boundaryFaceInterpolate!(sbpface, bndry.face,
                                     view(u,:,bndry.element),
                                     view(uface,:,bindex))
          end
          @test ≈(vec(uface[:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=5e-14)
        end
      end
    end
  end

  @testset "Testing boundaryFaceInterpolate! (TriSBP, vector field method)" begin
    # build a two element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      u = zeros(Float64, (2, sbp.numnodes, 2))
      uface = zeros(Float64, (2, sbpface.numnodes, 4))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j)
          for (bindex, bndry) in enumerate(bndryfaces)
            boundaryFaceInterpolate!(sbpface, bndry.face,
                                     view(u,:,:,bndry.element),
                                     view(uface,:,:,bindex))
          end
          @test ≈(vec(uface[1,:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=5e-14)
          @test ≈(vec(uface[2,:,:]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=5e-14)
        end
      end
    end
  end

  @testset "Testing boundaryFaceInterpolate! (TriSparseFace, scalar field method)" begin
    # build a two element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      u = zeros(Float64, (sbp.numnodes, 2))
      uface = zeros(Float64, (sbpface.numnodes, 4))      
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          for (bindex, bndry) in enumerate(bndryfaces)          
            boundaryFaceInterpolate!(sbpface, bndry.face,
                                     view(u,:,bndry.element),
                                     view(uface,:,bindex))
          end
          @test ≈(vec(uface[:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=5e-14)
        end
      end
    end
  end

  @testset "Testing boundaryFaceInterpolate! (TriSparseFace, vector field method)" begin
    # build a two element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array{Boundary}(undef, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      u = zeros(Float64, (2, sbp.numnodes, 2))
      uface = zeros(Float64, (2, sbpface.numnodes, 4))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j)
          for (bindex, bndry) in enumerate(bndryfaces)
            boundaryFaceInterpolate!(sbpface, bndry.face,
                                     view(u,:,:,bndry.element),
                                     view(uface,:,:,bindex))
          end
          @test ≈(vec(uface[1,:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=5e-14)
          @test ≈(vec(uface[2,:,:]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=5e-14)
        end
      end
    end
  end

  @testset "Testing boundaryFaceInterpolate! (TetSBP, scalar field method)" begin
    # build a four element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (3,sbp.numnodes,4))
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      u = zeros(Float64, (sbp.numnodes, 4))
      uface = zeros(Float64, (sbpface.numnodes, 12))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            for (bindex, bndry) in enumerate(bndryfaces)          
              boundaryFaceInterpolate!(sbpface, bndry.face,
                                       view(u,:,bndry.element),
                                       view(uface,:,bindex))
            end
            @test ≈(vec(uface[:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=5e-14)
          end
        end
      end
    end
  end

  @testset "Testing boundaryFaceInterpolate! (TetSBP, vector field method)" begin
    # build a four element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (3,sbp.numnodes,4))
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      u = zeros(Float64, (2, sbp.numnodes, 4))
      uface = zeros(Float64, (2, sbpface.numnodes, 12))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            for (bindex, bndry) in enumerate(bndryfaces)
              boundaryFaceInterpolate!(sbpface, bndry.face,
                                       view(u,:,:,bndry.element),
                                       view(uface,:,:,bindex))
            end
            @test ≈(vec(uface[1,:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=5e-14)
            @test ≈(vec(uface[2,:,:]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=5e-13)
          end
        end
      end
    end
  end

  @testset "Testing boundaryFaceInterpolate! (TetSparseFace, scalar field method)" begin
    # build a four element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      x = zeros(Float64, (3,sbp.numnodes,4))
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      u = zeros(Float64, (sbp.numnodes, 4))
      uface = zeros(Float64, (sbpface.numnodes, 12))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            for (bindex, bndry) in enumerate(bndryfaces)          
              boundaryFaceInterpolate!(sbpface, bndry.face,
                                       view(u,:,bndry.element),
                                       view(uface,:,bindex))
            end
            @test ≈(vec(uface[:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=5e-14)
          end
        end
      end
    end
  end

  @testset "Testing boundaryFaceInterpolate! (TetSparseFace, vector field method)" begin
    # build a four element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      x = zeros(Float64, (3,sbp.numnodes,4))
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array{Boundary}(undef, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      u = zeros(Float64, (2, sbp.numnodes, 4))
      uface = zeros(Float64, (2, sbpface.numnodes, 12))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            for (bindex, bndry) in enumerate(bndryfaces)
              boundaryFaceInterpolate!(sbpface, bndry.face,
                                       view(u,:,:,bndry.element),
                                       view(uface,:,:,bindex))
            end
            @test ≈(vec(uface[1,:,:]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=5e-14)
            @test ≈(vec(uface[2,:,:]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=5e-13)
          end
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing interiorfaceinterpolate! ($(string($TSBP)) scalar field method)" begin
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          x = zeros(Float64, (1,sbp.numnodes,2))
          xf = zeros(Float64, (1,sbpface.numnodes,1))
          vtx = reshape([0.0; 1.0], (2,1))
          x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          vtx = reshape([2.0; 1.0], (2,1))  # note the reversal !!!
          x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
          ifaces = Array{Interface}(undef, 1)
          ifaces[1] = Interface(1,2,2,2,1)
          u = zeros(Float64, (sbp.numnodes, 2))
          uface = zeros(Float64, (2, sbpface.numnodes, 1))
          for i = 0:p
            u[:,:] = x[1,:,:].^i
            interiorfaceinterpolate!(sbpface, ifaces, u, uface)
            @test ≈(vec(uface[1,:,:]), vec(xf[1,:,:].^i), atol=5e-14)
            @test ≈(vec(uface[2,:,:]), vec(xf[1,:,:].^i), atol=5e-14)
          end
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing interiorfaceinterpolate! ($(string($TSBP)) vector field method)" begin
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          x = zeros(Float64, (1,sbp.numnodes,2))
          xf = zeros(Float64, (1,sbpface.numnodes,1))
          vtx = reshape([0.0; 1.0], (2,1))
          x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          vtx = reshape([2.0; 1.0], (2,1))  # note the reversal !!!
          x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
          ifaces = Array{Interface}(undef,1)
          ifaces[1] = Interface(1,2,2,2,1)          
          u = zeros(Float64, (2, sbp.numnodes, 2))
          uface = zeros(Float64, (2, 2, sbpface.numnodes, 1))
          for i = 0:p
            u[1,:,:] = x[1,:,:].^i
            u[2,:,:] = 2.0.*u[1,:,:]
            interiorfaceinterpolate!(sbpface, ifaces, u, uface)
            @test ≈(vec(uface[1,1,:,:]), vec(xf[1,:,:].^i), atol=5e-14)
            @test ≈(vec(uface[2,1,:,:]), 2.0.*vec(xf[1,:,:].^i), atol=5e-14)
            @test ≈(vec(uface[1,2,:,:]), vec(xf[1,:,:].^i), atol=5e-14)
            @test ≈(vec(uface[2,2,:,:]), 2.0.*vec(xf[1,:,:].^i), atol=5e-14)
          end
        end
      end
    end
  end

  @testset "Testing interiorfaceinterpolate! (TriSBP scalar field method)" begin
    # build a two element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, [-1. -1.; 1. -1.; -1. 1.])
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,1))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      ifaces = Array{Interface}(undef,1)
      ifaces[1] = Interface(1,2,2,3,1)
      u = zeros(Float64, (sbp.numnodes,2))
      uface = zeros(Float64, (2, sbpface.numnodes, 1))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          interiorfaceinterpolate!(sbpface, ifaces, u, uface)
          # check that interpolation from left and right elements is exact
          @test ≈(vec(uface[1,:,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)),
                                              atol=1e-13)
          #@test ≈(uface[2,sbpface.nbrperm[:,1],1] --> 
          #roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[2,:,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)),
                                              atol=1e-13)
        end
      end
    end
  end

  @testset "Testing interiorfaceinterpolate! (TriSBP vector field method)" begin
    # build a two element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, [-1. -1.; 1. -1.; -1. 1.])
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,1))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      ifaces = Array{Interface}(undef,1)
      ifaces[1] = Interface(1,2,2,3,1)
      u = zeros(Float64, (2,sbp.numnodes,2))
      uface = zeros(Float64, (2, 2, sbpface.numnodes, 1))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j)
          interiorfaceinterpolate!(sbpface, ifaces, u, uface)
          # check that interpolation from left and right elements is exact
          @test ≈(vec(uface[1,1,:,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[2,1,:,1]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[1,2,:,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[2,2,:,1]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
        end
      end
    end
  end

  @testset "Testing interiorfaceinterpolate! (TriSparseFace scalar field method)" begin
    # build a two element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,1))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      ifaces = Array{Interface}(undef,1)
      ifaces[1] = Interface(1,2,2,3,1)
      u = zeros(Float64, (sbp.numnodes,2))
      uface = zeros(Float64, (2, sbpface.numnodes, 1))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          interiorfaceinterpolate!(sbpface, ifaces, u, uface)
          # check that interpolation from left and right elements is exact
          @test ≈(vec(uface[1,:,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)),
                                              atol=1e-13)
          #@test ≈(uface[2,sbpface.nbrperm[:,1],1] --> 
          #roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[2,:,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)),
                                              atol=1e-13)
        end
      end
    end
  end

  @testset "Testing interiorfaceinterpolate! (TriSparseFace vector field method)" begin
    # build a two element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,1))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      ifaces = Array{Interface}(undef,1)
      ifaces[1] = Interface(1,2,2,3,1)
      u = zeros(Float64, (2,sbp.numnodes,2))
      uface = zeros(Float64, (2, 2, sbpface.numnodes, 1))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j)
          interiorfaceinterpolate!(sbpface, ifaces, u, uface)
          # check that interpolation from left and right elements is exact
          @test ≈(vec(uface[1,1,:,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[2,1,:,1]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[1,2,:,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[2,2,:,1]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
        end
      end
    end
  end

  @testset "Testing interiorfaceinterpolate! (TetSBP, scalar field method)" begin
    # build a five element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (3,sbp.numnodes,5))
      xf = zeros(Float64, (3,sbpface.numnodes,4))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      ifaces = Array{Interface}(undef,4)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])
            
      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 0.; 0. 0. 1.; 1. 1. 1; 0. 1. 0.]
      x[:,:,5] = SymCubatures.calcnodes(sbp.cub, vtx)

      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      u = zeros(Float64, (sbp.numnodes, 5))
      uface = zeros(Float64, (2, sbpface.numnodes, 4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            interiorfaceinterpolate!(sbpface, ifaces, u, uface)
            # check that interpolation from left and right elements is exact
            for f = 1:4
              @test ≈(vec(uface[1,:,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[2,:,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
            end
          end
        end
      end
    end
  end

  @testset "Testing interiorfaceinterpolate! (TetSBP, vector field method)" begin
    # build a five element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (3,sbp.numnodes,5))
      xf = zeros(Float64, (3,sbpface.numnodes,4))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      ifaces = Array{Interface}(undef,4)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])
            
      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 0.; 0. 0. 1.; 1. 1. 1; 0. 1. 0.]
      x[:,:,5] = SymCubatures.calcnodes(sbp.cub, vtx)

      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      u = zeros(Float64, (2, sbp.numnodes, 5))
      uface = zeros(Float64, (2, 2, sbpface.numnodes, 4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            interiorfaceinterpolate!(sbpface, ifaces, u, uface)
            # check that interpolation from left and right elements is exact
            for f = 1:4
              @test ≈(vec(uface[1,1,:,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[2,1,:,f]), 2.0.*vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[1,2,:,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[2,2,:,f]), 2.0.*vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
            end
          end
        end
      end
    end
  end

  @testset "Testing interiorfaceinterpolate! (TetSparseFace scalar field method)" begin
    # build a five element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      x = zeros(Float64, (3,sbp.numnodes,5))
      xf = zeros(Float64, (3,sbpface.numnodes,4))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      ifaces = Array{Interface}(undef,4)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])
            
      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 0.; 0. 0. 1.; 1. 1. 1; 0. 1. 0.]
      x[:,:,5] = SymCubatures.calcnodes(sbp.cub, vtx)

      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      u = zeros(Float64, (sbp.numnodes, 5))
      uface = zeros(Float64, (2, sbpface.numnodes, 4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            interiorfaceinterpolate!(sbpface, ifaces, u, uface)
            # check that interpolation from left and right elements is exact
            for f = 1:4
              @test ≈(vec(uface[1,:,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[2,:,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
            end
          end
        end
      end
    end
  end

  @testset "Testing interiorfaceinterpolate! (TetSparseFace vector field method)" begin
    # build a five element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      x = zeros(Float64, (3,sbp.numnodes,5))
      xf = zeros(Float64, (3,sbpface.numnodes,4))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      ifaces = Array{Interface}(undef,4)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])
            
      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 0.; 0. 0. 1.; 1. 1. 1; 0. 1. 0.]
      x[:,:,5] = SymCubatures.calcnodes(sbp.cub, vtx)

      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      u = zeros(Float64, (2, sbp.numnodes, 5))
      uface = zeros(Float64, (2, 2, sbpface.numnodes, 4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            interiorfaceinterpolate!(sbpface, ifaces, u, uface)
            # check that interpolation from left and right elements is exact
            for f = 1:4
              @test ≈(vec(uface[1,1,:,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[2,1,:,f]), 2.0.*vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[1,2,:,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[2,2,:,f]), 2.0.*vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
            end
          end
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing interiorFaceInterpolate! ($(string($TSBP)) scalar field method)" begin
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          x = zeros(Float64, (1,sbp.numnodes,2))
          xf = zeros(Float64, (1,sbpface.numnodes,1))
          vtx = reshape([0.0; 1.0], (2,1))
          x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          vtx = reshape([2.0; 1.0], (2,1))  # note the reversal !!!
          x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
          ifaces = Array{Interface}(undef,1)
          ifaces[1] = Interface(1,2,2,2,1)
          u = zeros(Float64, (sbp.numnodes, 2))
          uface = zeros(Float64, (sbpface.numnodes, 2, 1))
          for i = 0:p
            u[:,:] = x[1,:,:].^i
            for (findex, face) in enumerate(ifaces)          
              interiorFaceInterpolate!(sbpface, face, view(u,:,face.elementL),
                                       view(u,:,face.elementR),
                                       view(uface,:,1,findex),
                                       view(uface,:,2,findex))
            end
            @test ≈(vec(uface[:,1,1]), vec(xf[1,:,:].^i), atol=5e-14)
            @test ≈(vec(uface[:,2,1]), vec(xf[1,:,:].^i), atol=5e-14)
          end
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing interiorFaceInterpolate! ($(string($TSBP)) vector field method)" begin
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = getLineSegFace(p, sbp.cub, sbp.vtx)
          x = zeros(Float64, (1,sbp.numnodes,2))
          xf = zeros(Float64, (1,sbpface.numnodes,1))
          vtx = reshape([0.0; 1.0], (2,1))
          x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
          xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, reshape(vtx[[2;]],(1,1)))
          vtx = reshape([2.0; 1.0], (2,1))  # note the reversal !!!
          x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
          ifaces = Array{Interface}(undef,1)
          ifaces[1] = Interface(1,2,2,2,1)          
          u = zeros(Float64, (2, sbp.numnodes, 2))
          uface = zeros(Float64, (2, sbpface.numnodes, 2, 1))
          for i = 0:p
            u[1,:,:] = x[1,:,:].^i
            u[2,:,:] = 2.0.*u[1,:,:]
            for (findex, face) in enumerate(ifaces)
              interiorFaceInterpolate!(sbpface, face, view(u,:,:,face.elementL),
                                       view(u,:,:,face.elementR),
                                       view(uface,:,:,1,findex),
                                       view(uface,:,:,2,findex))
            end
            @test ≈(vec(uface[1,:,1,1]), vec(xf[1,:,:].^i), atol=5e-14)
            @test ≈(vec(uface[2,:,1,1]), 2.0.*vec(xf[1,:,:].^i), atol=5e-14)
            @test ≈(vec(uface[1,:,2,1]), vec(xf[1,:,:].^i), atol=5e-14)
            @test ≈(vec(uface[2,:,2,1]), 2.0.*vec(xf[1,:,:].^i), atol=5e-14)
          end
        end
      end
    end
  end

  @testset "Testing interiorFaceInterpolate! (TriSBP, scalar field method)" begin
    # build a two element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, [-1. -1.; 1. -1.; -1. 1.])
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,1))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      ifaces = Array{Interface}(undef,1)
      ifaces[1] = Interface(1,2,2,3,1)
      u = zeros(Float64, (sbp.numnodes,2))
      uface = zeros(Float64, (sbpface.numnodes, 2, 1))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          for (findex, face) in enumerate(ifaces)          
            interiorFaceInterpolate!(sbpface, face, view(u,:,face.elementL),
                                     view(u,:,face.elementR),
                                     view(uface,:,1,findex),
                                     view(uface,:,2,findex))
          end
          # check that interpolation from left and right elements is exact
          @test ≈(vec(uface[:,1,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)),
                                              atol=1e-13)
          #@test ≈(uface[2,sbpface.nbrperm[:,1],1] --> 
          #roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[:,2,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)),
                                              atol=1e-13)
        end
      end
    end
  end

  @testset "Testing interiorFaceInterpolate! (TriSBP, vector field method)" begin
    # build a two element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, [-1. -1.; 1. -1.; -1. 1.])
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,1))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      ifaces = Array{Interface}(undef,1)
      ifaces[1] = Interface(1,2,2,3,1)
      u = zeros(Float64, (2,sbp.numnodes,2))
      uface = zeros(Float64, (2, sbpface.numnodes, 2, 1))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j)
          for (findex, face) in enumerate(ifaces)          
            interiorFaceInterpolate!(sbpface, face, view(u,:,:,face.elementL),
                                     view(u,:,:,face.elementR),
                                     view(uface,:,:,1,findex),
                                     view(uface,:,:,2,findex))
          end
          # check that interpolation from left and right elements is exact
          @test ≈(vec(uface[1,:,1,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[2,:,1,1]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[1,:,2,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[2,:,2,1]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
        end
      end
    end
  end

  @testset "Testing interiorFaceInterpolate! (TriSparseFace scalar field method)" begin
    # build a two element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,1))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      ifaces = Array{Interface}(undef,1)
      ifaces[1] = Interface(1,2,2,3,1)
      u = zeros(Float64, (sbp.numnodes,2))
      uface = zeros(Float64, (sbpface.numnodes, 2, 1))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          for (findex, face) in enumerate(ifaces)          
            interiorFaceInterpolate!(sbpface, face, view(u,:,face.elementL),
                                     view(u,:,face.elementR),
                                     view(uface,:,1,findex),
                                     view(uface,:,2,findex))
          end
          # check that interpolation from left and right elements is exact
          @test ≈(vec(uface[:,1,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)),
                                              atol=1e-13)
          #@test ≈(uface[2,sbpface.nbrperm[:,1],1] --> 
          #roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[:,2,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)),
                                              atol=1e-13)
        end
      end
    end
  end

  @testset "Testing interiorFaceInterpolate! (TriSparseFace vector field method)" begin
    # build a two element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,1))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      ifaces = Array{Interface}(undef,1)
      ifaces[1] = Interface(1,2,2,3,1)
      u = zeros(Float64, (2,sbp.numnodes,2))
      uface = zeros(Float64, (2, sbpface.numnodes, 2, 1))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j)
          for (findex, face) in enumerate(ifaces)          
            interiorFaceInterpolate!(sbpface, face, view(u,:,:,face.elementL),
                                     view(u,:,:,face.elementR),
                                     view(uface,:,:,1,findex),
                                     view(uface,:,:,2,findex))
          end
          # check that interpolation from left and right elements is exact
          @test ≈(vec(uface[1,:,1,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[2,:,1,1]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[1,:,2,1]), vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @test ≈(vec(uface[2,:,2,1]), 2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
        end
      end
    end
  end

  @testset "Testing interiorFaceInterpolate! (TetSBP, scalar field method)" begin
    # build a five element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (3,sbp.numnodes,5))
      xf = zeros(Float64, (3,sbpface.numnodes,4))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      ifaces = Array{Interface}(undef,4)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])
            
      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 0.; 0. 0. 1.; 1. 1. 1; 0. 1. 0.]
      x[:,:,5] = SymCubatures.calcnodes(sbp.cub, vtx)

      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      u = zeros(Float64, (sbp.numnodes, 5))
      uface = zeros(Float64, (sbpface.numnodes, 2, 4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            for (findex, face) in enumerate(ifaces)          
              interiorFaceInterpolate!(sbpface, face, view(u,:,face.elementL),
                                       view(u,:,face.elementR),
                                       view(uface,:,1,findex),
                                       view(uface,:,2,findex))
            end
            # check that interpolation from left and right elements is exact
            for f = 1:4
              @test ≈(vec(uface[:,1,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[:,2,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
            end
          end
        end
      end
    end
  end

  @testset "Testing interiorFaceInterpolate! (TetSBP, vector field method)" begin
    # build a five element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (3,sbp.numnodes,5))
      xf = zeros(Float64, (3,sbpface.numnodes,4))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      ifaces = Array{Interface}(undef,4)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])
            
      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 0.; 0. 0. 1.; 1. 1. 1; 0. 1. 0.]
      x[:,:,5] = SymCubatures.calcnodes(sbp.cub, vtx)

      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      u = zeros(Float64, (2, sbp.numnodes, 5))
      uface = zeros(Float64, (2, sbpface.numnodes, 2, 4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            for (findex, face) in enumerate(ifaces)          
              interiorFaceInterpolate!(sbpface, face, view(u,:,:,face.elementL),
                                       view(u,:,:,face.elementR),
                                       view(uface,:,:,1,findex),
                                       view(uface,:,:,2,findex))
            end
            # check that interpolation from left and right elements is exact
            for f = 1:4
              @test ≈(vec(uface[1,:,1,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[2,:,1,f]), 2.0.*vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[1,:,2,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[2,:,2,f]), 2.0.*vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
            end
          end
        end
      end
    end
  end

  @testset "Testing interiorFaceInterpolate! (TetSparseFace scalar field method)" begin
    # build a five element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      x = zeros(Float64, (3,sbp.numnodes,5))
      xf = zeros(Float64, (3,sbpface.numnodes,4))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      ifaces = Array{Interface}(undef,4)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])
            
      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 0.; 0. 0. 1.; 1. 1. 1; 0. 1. 0.]
      x[:,:,5] = SymCubatures.calcnodes(sbp.cub, vtx)

      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      u = zeros(Float64, (sbp.numnodes, 5))
      uface = zeros(Float64, (sbpface.numnodes, 2, 4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            for (findex, face) in enumerate(ifaces)          
              interiorFaceInterpolate!(sbpface, face, view(u,:,face.elementL),
                                       view(u,:,face.elementR),
                                       view(uface,:,1,findex),
                                       view(uface,:,2,findex))
            end
            # check that interpolation from left and right elements is exact
            for f = 1:4
              @test ≈(vec(uface[:,1,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[:,2,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
            end
          end
        end
      end
    end
  end

  @testset "Testing interiorFaceInterpolate! (TetSparseFace vector field method)" begin
    # build a five element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:1
      sbp = getTetSBPDiagE(degree=p, faceopertype=:Omega)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
      x = zeros(Float64, (3,sbp.numnodes,5))
      xf = zeros(Float64, (3,sbpface.numnodes,4))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      ifaces = Array{Interface}(undef,4)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])
            
      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 0.; 0. 0. 1.; 1. 1. 1; 0. 1. 0.]
      x[:,:,5] = SymCubatures.calcnodes(sbp.cub, vtx)

      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      u = zeros(Float64, (2, sbp.numnodes, 5))
      uface = zeros(Float64, (2, sbpface.numnodes, 2, 4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            for (findex, face) in enumerate(ifaces)          
              interiorFaceInterpolate!(sbpface, face, view(u,:,:,face.elementL),
                                       view(u,:,:,face.elementR),
                                       view(uface,:,:,1,findex),
                                       view(uface,:,:,2,findex))
            end
            # check that interpolation from left and right elements is exact
            for f = 1:4
              @test ≈(vec(uface[1,:,1,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[2,:,1,f]), 2.0.*vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[1,:,2,f]), vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @test ≈(vec(uface[2,:,2,f]), 2.0.*vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
            end
          end
        end
      end
    end
  end

end
