@testset "Testing SummationByParts Module (reverse-diff of volume integrate methods)..." begin

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre,
              getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE,
              getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      @testset "Testing volumeintegrate_rev! ($(string($TSBP)) scalar field method)" begin
        for p = 1:4
          # verify that H*u = (u^T*H)^T
          sbp = ($TSBP)(degree=p)
          u = rand(Float64, (sbp.numnodes,2))
          Hu = zeros(size(u))
          utH = zeros(size(u))
          volumeintegrate!(sbp, u, Hu)
          volumeintegrate_rev!(sbp, utH, u)
          @test ≈(Hu, utH, atol=1e-15)          
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre,
              getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE,
              getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      @testset "Testing volumeintegrate_rev! ($(string($TSBP)) vector field method)" begin
        for p = 1:4
          # verify that H*u = (u^T*H)^T
          sbp = ($TSBP)(degree=p)
          u = rand(Float64, (4,sbp.numnodes,2))
          Hu = zeros(size(u))
          utH = zeros(size(u))
          volumeintegrate!(sbp, u, Hu)
          volumeintegrate_rev!(sbp, utH, u)
          @test ≈(Hu, utH, atol=1e-15)          
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre,
              getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE,
              getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      @testset "Testing volumeIntegrateElement_rev! ($(string($TSBP)) scalar field method)" begin
        for p = 1:4
          # verify that H*u = (u^T*H)^T
          sbp = ($TSBP)(degree=p)
          u = rand(Float64, (sbp.numnodes))
          Hu = zeros(size(u))
          utH = zeros(size(u))
          volumeIntegrateElement!(sbp, u, Hu)
          volumeIntegrateElement_rev!(sbp, utH, u)
          @test ≈(Hu, utH, atol=1e-15)          
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre,
              getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE,
              getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      @testset "Testing volumeIntegrateElement_rev! ($(string($TSBP)) vector field method)" begin
        for p = 1:4
          # verify that H*u = (u^T*H)^T
          sbp = ($TSBP)(degree=p)
          u = rand(Float64, (4,sbp.numnodes))
          Hu = zeros(size(u))
          utH = zeros(size(u))
          volumeIntegrateElement!(sbp, u, Hu)
          volumeIntegrateElement_rev!(sbp, utH, u)
          @test ≈(Hu, utH, atol=1e-15)
        end
      end
    end
  end
  
end
