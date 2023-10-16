@testset "Testing SummationByParts Module (weak differentiate methods" begin

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing weakdifferentiate! ($(string($TSBP)) scalar field method)" begin
        # build a two element grid, and verify that Qxi * 1 = 0
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes,2))
          di = 1
          res = zeros(size(u))
          weakdifferentiate!(sbp, di, u, res)
          @test ≈(res[:,1], zeros(sbp.numnodes); atol = 5e-13)
          @test ≈(res[:,2], zeros(sbp.numnodes); atol = 5e-13)
        end
      end 
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      @testset "Testing weakdifferentiate! ($(string($TSBP)) scalar field method)" begin
        # build a two element grid, and verify that Qxi * 1 = 0
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes,2))
          di = 1
          res = zeros(size(u))
          weakdifferentiate!(sbp, di, u, res)
          @test ≈(res[:,1], zeros(sbp.numnodes); atol = 5e-13)
          @test ≈(res[:,2], zeros(sbp.numnodes); atol = 5e-13)
        end
      end 
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      @testset "Testing weakdifferentiate! ($(string($TSBP)) scalar field method)" begin
        # build a single element grid, and verify that Qxi * 1 = 0
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes,1))
          di = 1
          res = zeros(size(u))
          weakdifferentiate!(sbp, di, u, res)
          @test ≈(res[:,1], zeros(sbp.numnodes); atol = 1e-12)
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing weakdifferentiate! ($(string($TSBP)) vector field method)" begin
        # build a two element grid, and verify that \int (Dxi * u) d\Omega = 0
        # or 1 if u = 1 or x, respectively
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          vtx = reshape([0.0; 1.0], (2,1))
          x = ones(Float64, (2,sbp.numnodes,2))
          x[2,:,1] = calcnodes(sbp, vtx)          
          vtx = reshape([1.0; 2.0], (2,1))
          x[2,:,2] = calcnodes(sbp, vtx)
          di = 1
          res = zeros(size(x))
          weakdifferentiate!(sbp, di, x, res)
          @test ≈(sum(res[1,:,1]), 0.0, atol=1e-14)
          @test ≈(sum(res[2,:,1]), 1.0, atol=1e-14)
          @test ≈(sum(res[1,:,2]), 0.0, atol=1e-14)
          @test ≈(sum(res[2,:,2]), 1.0, atol=1e-14)
        end
      end
    end
  end
  
  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      @testset "Testing weakdifferentiate! ($(string($TSBP)) vector field method)" begin
        # build a two element grid, and verify that \int (Dxi * x) d\Omega = 1 or 0,
        # depending on orientation of local coordinates
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          vtx = [0. 0.; 1. 0.; 0. 1.]
          x = zeros(Float64, (2,sbp.numnodes,2))
          x[:,:,1] = calcnodes(sbp, vtx)
          vtx = [1. 0.; 1. 1.; 0. 1.]
          x[:,:,2] = calcnodes(sbp, vtx)
          di = 1
          res = zeros(size(x))
          weakdifferentiate!(sbp, di, x, res)
          @test ≈(sum(res[1,:,1]), 1.0, atol=1e-14)
          @test ≈(sum(res[2,:,1]), 0.0, atol=1e-14)
          @test ≈(sum(res[1,:,2]), 0.0, atol=1e-14)
          @test ≈(sum(res[2,:,2]), 1.0, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      @testset "Testing weakdifferentiate! ($(string($TSBP)) vector field method)" begin
        # build a single element grid, and verify that \int (D * x) d\Omega = 2/3
        # or 0, depending on orientation of local coordinates
        Id = (2/3).* I(3)
        vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
        for di = 1:3
          for p = 1:4
            sbp = ($TSBP)(degree=p)
            x = zeros(Float64, (3,sbp.numnodes,1))
            x[:,:,1] = calcnodes(sbp, vtx)
            res = zeros(size(x))
            weakdifferentiate!(sbp, di, x, res)
            @test ≈(sum(res[1,:,1]), Id[di,1], atol=1e-13)
            @test ≈(sum(res[2,:,1]), Id[di,2], atol=1e-13)
            @test ≈(sum(res[3,:,1]), Id[di,3], atol=1e-13)
          end
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing weakdifferentiateElement! ($(string($TSBP)) scalar field method)" begin
        # verify that Qxi * 1 = 0
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes))
          di = 1
          res = zeros(size(u))
          weakDifferentiateElement!(sbp, di, u, res)
          @test ≈(res[:], zeros(sbp.numnodes), atol=5e-13)
        end
      end 
    end
  end
  
  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      @testset "Testing weakdifferentiateElement! ($(string($TSBP)) scalar field method)" begin
        # verify that Qxi * 1 = 0
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes))
          di = 1
          res = zeros(size(u))
          weakDifferentiateElement!(sbp, di, u, res)
          @test ≈(res[:], zeros(sbp.numnodes), atol=5e-13)
        end
      end 
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      @testset "Testing weakDifferentiateElement! ($(string($TSBP)) scalar field method)" begin
        # build a single element, and verify that Qxi * 1 = 0
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes))
          di = 1
          res = zeros(size(u))
          weakDifferentiateElement!(sbp, di, u, res)
          @test ≈(res[:], zeros(sbp.numnodes), atol=5e-13)
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing weakdifferentiateElement! ($(string($TSBP)) vector field method)" begin
        # build a two element grid, and verify that \int (Dxi * u) d\Omega = 0
        # or 1 if u = 1 or x, respectively
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          vtx = reshape([0.0; 1.0], (2,1))
          x = ones(Float64, (2,sbp.numnodes,2))
          x[2,:,1] = calcnodes(sbp, vtx)          
          vtx = reshape([1.0; 2.0], (2,1))
          x[2,:,2] = calcnodes(sbp, vtx)
          di = 1
          res = zeros(size(x))
          weakDifferentiateElement!(sbp, di, view(x,:,:,1), view(res,:,:,1))
          weakDifferentiateElement!(sbp, di, view(x,:,:,2), view(res,:,:,2))
          @test ≈(sum(res[1,:,1]), 0.0, atol=1e-14)
          @test ≈(sum(res[2,:,1]), 1.0, atol=1e-14)
          @test ≈(sum(res[1,:,2]), 0.0, atol=1e-14)
          @test ≈(sum(res[2,:,2]), 1.0, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      @testset "Testing weakdifferentiateElement! ($(string($TSBP)) vector field method)" begin  
        # build a two element grid, and verify that \int (Dxi * x) d\Omega = 1 or 0,
        # depending on orientation of local coordinates
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          vtx = [0. 0.; 1. 0.; 0. 1.]
          x = zeros(Float64, (2,sbp.numnodes,2))
          x[:,:,1] = calcnodes(sbp, vtx)
          vtx = [1. 0.; 1. 1.; 0. 1.]
          x[:,:,2] = calcnodes(sbp, vtx)
          di = 1
          res = zeros(size(x))
          weakDifferentiateElement!(sbp, di, view(x,:,:,1), view(res,:,:,1))
          weakDifferentiateElement!(sbp, di, view(x,:,:,2), view(res,:,:,2))
          @test ≈(sum(res[1,:,1]), 1.0, atol=1e-14)
          @test ≈(sum(res[2,:,1]), 0.0, atol=1e-14)
          @test ≈(sum(res[1,:,2]), 0.0, atol=1e-14)
          @test ≈(sum(res[2,:,2]), 1.0, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      @testset "Testing weakDifferentiateElement! ($(string($TSBP)) vector field method)" begin
        # build a single element grid, and verify that \int (D * x) d\Omega = 2/3
        # or 0, depending on orientation of local coordinates
        Id = (2/3) .* I(3)
        vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
        for di = 1:3
          for p = 1:4
            sbp = ($TSBP)(degree=p)
            x = zeros(Float64, (3,sbp.numnodes,1))
            x[:,:,1] = calcnodes(sbp, vtx)
            res = zeros(size(x))
            weakDifferentiateElement!(sbp, di, view(x,:,:,1), view(res,:,:,1))
            @test ≈(sum(res[1,:,1]), Id[di,1], atol=1e-14)
            @test ≈(sum(res[2,:,1]), Id[di,2], atol=1e-14)
            @test ≈(sum(res[3,:,1]), Id[di,3], atol=1e-14)
          end
        end
      end
    end
  end
  
end
