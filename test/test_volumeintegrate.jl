@testset "Testing SummationByParts Module (volume integrate methods)..." begin

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing volumeintegrate! ($(string($TSBP)) scalar field method)" begin
        # build a two element grid, and verify that ones*H*ones = vol
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes,2))
          res = zeros(size(u))
          volumeintegrate!(sbp, u, res)
          @test ≈(sum(res[:,1]), 2.0, atol=1e-14)
          @test ≈(sum(res[:,2]), 2.0, atol=1e-14)
        end
      end
    end
  end
  
  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      @testset "Testing volumeintegrate! ($(string($TSBP)) scalar field method)" begin
        # build a two element grid, and verify that ones*H*ones = vol
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes,2))
          res = zeros(size(u))
          volumeintegrate!(sbp, u, res)
          @test ≈(sum(res[:,1]), 2.0, atol=1e-14)
          @test ≈(sum(res[:,2]), 2.0, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      @testset "Testing volumeintegrate! ($(string($TSBP)) scalar field method)" begin
        # build a single element grid, and verify that ones*H*ones = vol
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes,1))
          res = zeros(size(u))
          volumeintegrate!(sbp, u, res)
          @test ≈(sum(res[:,1]), 4/3, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing volumeintegrate! ($(string($TSBP)) vector field method)" begin
        # build a two element grid, and verify that ones^T*H*ones = (1,2)
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (2,sbp.numnodes,2))
          u[1,:,:] *= 0.5
          res = zeros(size(u))
          volumeintegrate!(sbp, u, res)
          @test ≈(sum(res[1,:,1]), 1.0, atol=1e-14)
          @test ≈(sum(res[2,:,1]), 2.0, atol=1e-14)
          @test ≈(sum(res[1,:,2]), 1.0, atol=1e-14)
          @test ≈(sum(res[2,:,2]), 2.0, atol=1e-14)
        end
      end
    end
  end
  
  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      @testset "Testing volumeintegrate! ($(string($TSBP)) vector field method)" begin
        # build a two element grid, and verify that ones^T*H*ones = (1,2)
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (2,sbp.numnodes,2))
          u[1,:,:] *= 0.5
          res = zeros(size(u))
          volumeintegrate!(sbp, u, res)
          @test ≈(sum(res[1,:,1]), 1.0, atol=1e-14)
          @test ≈(sum(res[2,:,1]), 2.0, atol=1e-14)
          @test ≈(sum(res[1,:,2]), 1.0, atol=1e-14)
          @test ≈(sum(res[2,:,2]), 2.0, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      @testset "Testing volumeintegrate! ($(string($TSBP)) vector field method)" begin
        # build a single element grid, and verify that ones*H*ones = (1,2)
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (2,sbp.numnodes,1))
          u[1,:,:] *= 3/4
          u[2,:,:] *= 3/2
          res = zeros(size(u))
          volumeintegrate!(sbp, u, res)
          @test ≈(sum(res[1,:,1]), 1.0, atol=1e-14)
          @test ≈(sum(res[2,:,1]), 2.0, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing volumeIntegrateElement! ($(string($TSBP)) scalar field method)" begin
        # build a two element grid, and verify that ones*H*ones = vol
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes,2))
          res = zeros(size(u))
          volumeIntegrateElement!(sbp, view(u,:,1), view(res,:,1))
          volumeIntegrateElement!(sbp, view(u,:,2), view(res,:,2))
          @test ≈(sum(res[:,1]), 2.0, atol=1e-14)
          @test ≈(sum(res[:,2]), 2.0, atol=1e-14)
        end
      end
    end
  end
  
  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      @testset "Testing volumeIntegrateElement! ($(string($TSBP)) scalar field method)" begin
        # build a two element grid, and verify that ones*H*ones = vol
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes,2))
          res = zeros(size(u))
          volumeIntegrateElement!(sbp, view(u,:,1), view(res,:,1))
          volumeIntegrateElement!(sbp, view(u,:,2), view(res,:,2))
          @test ≈(sum(res[:,1]), 2.0, atol=1e-14)
          @test ≈(sum(res[:,2]), 2.0, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      @testset "Testing volumeIntegrateElement! ($(string($TSBP)) scalar field method)" begin
        # build a single element grid, and verify that ones*H*ones = vol
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes,1))
          res = zeros(size(u))
          volumeIntegrateElement!(sbp, view(u,:,1), view(res,:,1))
          @test ≈(sum(res[:,1]), 4/3, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      @testset "Testing volumeIntegrateElement! ($(string($TSBP)) vector field method)" begin
        # build a two element grid, and verify that ones^T*H*ones = (1,2)
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (2,sbp.numnodes,2))
          u[1,:,:] *= 0.5
          res = zeros(size(u))
          volumeIntegrateElement!(sbp, view(u,:,:,1), view(res,:,:,1))
          volumeIntegrateElement!(sbp, view(u,:,:,2), view(res,:,:,2))
          @test ≈(sum(res[1,:,1]), 1.0, atol=1e-14)
          @test ≈(sum(res[2,:,1]), 2.0, atol=1e-14)
          @test ≈(sum(res[1,:,2]), 1.0, atol=1e-14)
          @test ≈(sum(res[2,:,2]), 2.0, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      @testset "Testing volumeIntegrateElement! ($(string($TSBP)) vector field method)" begin
        # build a two element grid, and verify that ones^T*H*ones = (1,2)
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (2,sbp.numnodes,2))
          u[1,:,:] *= 0.5
          res = zeros(size(u))
          volumeIntegrateElement!(sbp, view(u,:,:,1), view(res,:,:,1))
          volumeIntegrateElement!(sbp, view(u,:,:,2), view(res,:,:,2))
          @test ≈(sum(res[1,:,1]), 1.0, atol=1e-14)
          @test ≈(sum(res[2,:,1]), 2.0, atol=1e-14)
          @test ≈(sum(res[1,:,2]), 1.0, atol=1e-14)
          @test ≈(sum(res[2,:,2]), 2.0, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      @testset "Testing volumeIntegrateElement! ($(string($TSBP)) vector field method)" begin
        # build a single element grid, and verify that ones*H*ones = (1,2)
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (2,sbp.numnodes,1))
          u[1,:,:] *= 3/4
          u[2,:,:] *= 3/2
          res = zeros(size(u))
          volumeIntegrateElement!(sbp, view(u,:,:,1), view(res,:,:,1))
          @test ≈(sum(res[1,:,1]), 1.0, atol=1e-14)
          @test ≈(sum(res[2,:,1]), 2.0, atol=1e-14)
        end
      end
    end
  end
  
end
