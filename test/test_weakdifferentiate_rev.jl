@testset "Testing SummationByParts Module (reverse-diff of weak differentiate methods)..." begin

  for TSBP = ((getLineSegSBPLobbato,1), (getLineSegSBPLegendre,1),
              (getTriSBPGamma,2), (getTriSBPOmega,2), (getTriSBPDiagE,2),
              (getTetSBPGamma,3), (getTetSBPOmega,3), (getTetSBPDiagE,3))
    @eval begin
      @testset "Testing weakdifferentiate_rev! ($(string($TSBP[1])) scalar field method)" begin
        for p = 1:4
          sbp = ($TSBP[1])(degree=p)
          v = rand(Float64, (sbp.numnodes,2))
          Qv = zeros(size(v))
          vQ = zeros(size(v))
          for di = 1:($TSBP[2])        
            # verify that v^T Q = (Q^T v)^T
            weakdifferentiate!(sbp, di, v, Qv, trans=true)
            weakdifferentiate_rev!(sbp, di, vQ, v, trans=false)
            @test ≈(Qv, vQ, atol=1e-15)
            # verify that v^T Q^T = (Qv)^T
            weakdifferentiate!(sbp, di, v, Qv, trans=false)
            weakdifferentiate_rev!(sbp, di, vQ, v, trans=true)
            @test ≈(Qv, vQ, atol=1e-15)
          end
        end
      end 
    end
  end

  for TSBP = ((getLineSegSBPLobbato,1), (getLineSegSBPLegendre,1),
              (getTriSBPGamma,2), (getTriSBPOmega,2), (getTriSBPDiagE,2),
              (getTetSBPGamma,3), (getTetSBPOmega,3), (getTetSBPDiagE,3))
    @eval begin
      @testset "Testing weakdifferentiate_rev! ($(string($TSBP[1])) vector field method)" begin
        for p = 1:4
          sbp = ($TSBP[1])(degree=p)
          v = rand(Float64, (4,sbp.numnodes,2))
          Qv = zeros(size(v))
          vQ = zeros(size(v))
          for di = 1:($TSBP[2])           
            # verify that v^T Q = (Q^T v)^T
            weakdifferentiate!(sbp, di, v, Qv, trans=true)
            weakdifferentiate_rev!(sbp, di, vQ, v, trans=false)
            @test ≈(Qv, vQ, atol=1e-15)
            # verify that v^T Q^T = (Qv)^T
            weakdifferentiate!(sbp, di, v, Qv, trans=false)
            weakdifferentiate_rev!(sbp, di, vQ, v, trans=true)
            @test ≈(Qv, vQ, atol=1e-15)
          end
        end
      end 
    end
  end

  for TSBP = ((getLineSegSBPLobbato,1), (getLineSegSBPLegendre,1),
              (getTriSBPGamma,2), (getTriSBPOmega,2), (getTriSBPDiagE,2),
              (getTetSBPGamma,3), (getTetSBPOmega,3), (getTetSBPDiagE,3))
    @eval begin
      @testset "Testing weakDifferentiateElement_rev! ($(string($TSBP[1])) scalar field method)" begin
        for p = 1:4
          sbp = ($TSBP[1])(degree=p)
          v = rand(Float64, (sbp.numnodes))
          Qv = zeros(size(v))
          vQ = zeros(size(v))
          for di = 1:($TSBP[2])         
            # verify that v^T Q = (Q^T v)^T
            weakDifferentiateElement!(sbp, di, v, Qv, SummationByParts.Add(), true)
            weakDifferentiateElement_rev!(sbp, di, vQ, v, SummationByParts.Add(), false)
            @test ≈(Qv, vQ, atol=1e-15)
            # verify that v^T Q^T = (Qv)^T
            weakDifferentiateElement!(sbp, di, v, Qv, SummationByParts.Add(), false)
            weakDifferentiateElement_rev!(sbp, di, vQ, v, SummationByParts.Add(), true)
            @test ≈(Qv, vQ, atol=1e-15)
          end
        end
      end 
    end
  end

  for TSBP = ((getLineSegSBPLobbato,1), (getLineSegSBPLegendre,1),
              (getTriSBPGamma,2), (getTriSBPOmega,2), (getTriSBPDiagE,2),
              (getTetSBPGamma,3), (getTetSBPOmega,3), (getTetSBPDiagE,3))
    @eval begin
      @testset "Testing weakDifferentiateElement_rev! ($(string($TSBP[1])) vector field method)" begin
        for p = 1:4
          sbp = ($TSBP[1])(degree=p)
          v = rand(Float64, (4,sbp.numnodes))
          Qv = zeros(size(v))
          vQ = zeros(size(v))
          for di = 1:($TSBP[2])           
            # verify that v^T Q = (Q^T v)^T
            weakDifferentiateElement!(sbp, di, v, Qv, SummationByParts.Add(), true)
            weakDifferentiateElement_rev!(sbp, di, vQ, v, SummationByParts.Add(), false)
            @test ≈(Qv, vQ, atol=1e-15)
            # verify that v^T Q^T = (Qv)^T
            weakDifferentiateElement!(sbp, di, v, Qv, SummationByParts.Add(), false)
            weakDifferentiateElement_rev!(sbp, di, vQ, v, SummationByParts.Add(), true)
            @test ≈(Qv, vQ, atol=1e-15)
          end
        end
      end 
    end
  end
  
end
