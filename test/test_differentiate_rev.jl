@testset "Testing SummationByParts Module (reverse-diff of differentiate methods)..." begin

  for TSBP = ((getLineSegSBPLobbato,1), (getLineSegSBPLegendre,1),
              (getTriSBPGamma,2), (getTriSBPOmega,2), (getTriSBPDiagE,2),
              (getTetSBPGamma,3), (getTetSBPOmega,3), (getTetSBPDiagE,3))
    @eval begin
      @testset "Testing differentiate_rev! ($(string($TSBP[1])) scalar field method)" begin
        for p = 1:4
          sbp = ($TSBP[1])(degree=p)
          u = rand(Float64, (sbp.numnodes,2))
          v = rand(Float64, (sbp.numnodes,2))
          res = zeros(size(u))
          vtDu = zero(Float64)
          utDtv = zero(Float64)
          for di = 1:($TSBP[2])
            # compute the vector-matrix-vector product in forward mode
            fill!(res, 0.0)
            differentiate!(sbp, di, u, res)
            vtDu = sum(v.*res)
            # compute the vector-matrix-vector
            fill!(res, 0.0)
            differentiate_rev!(sbp, di, res, v)
            utDtv = sum(u.*res)
            @test ≈(vtDu, utDtv, atol=1e-13)
          end
        end
      end
    end
  end

  for TSBP = ((getLineSegSBPLobbato,1), (getLineSegSBPLegendre,1),
              (getTriSBPGamma,2), (getTriSBPOmega,2), (getTriSBPDiagE,2),
              (getTetSBPGamma,3), (getTetSBPOmega,3), (getTetSBPDiagE,3))
    @eval begin
      @testset "Testing differentiate_rev! ($(string($TSBP[1])) vector field method)" begin
        for p = 1:4
          sbp = ($TSBP[1])(degree=p)
          u = rand(Float64, (4,sbp.numnodes,2))
          v = rand(Float64, (4,sbp.numnodes,2))
          res = zeros(size(u))
          vtDu = zero(Float64)
          utDtv = zero(Float64)
          for di = 1:($TSBP[2])
            # compute the vector-matrix-vector product in forward mode
            fill!(res, 0.0)
            differentiate!(sbp, di, u, res)
            vtDu = sum(v.*res)
            # compute the vector-matrix-vector
            fill!(res, 0.0)
            differentiate_rev!(sbp, di, res, v)
            utDtv = sum(u.*res)
            @test ≈(vtDu, utDtv, atol=1e-12)
          end
        end
      end
    end
  end

  for TSBP = ((getLineSegSBPLobbato,1), (getLineSegSBPLegendre,1),
              (getTriSBPGamma,2), (getTriSBPOmega,2), (getTriSBPDiagE,2),
              (getTetSBPGamma,3), (getTetSBPOmega,3), (getTetSBPDiagE,3))
    @eval begin
      @testset "Testing differentiateElement_rev! ($(string($TSBP[1])) scalar field method)" begin
        for p = 1:4
          sbp = ($TSBP[1])(degree=p)
          u = rand(Float64, (sbp.numnodes))
          v = rand(Float64, (sbp.numnodes))
          res = zeros(size(u))
          vtDu = zero(Float64)
          utDtv = zero(Float64)
          for di = 1:($TSBP[2])
            # compute the vector-matrix-vector product in forward mode
            fill!(res, 0.0)
            differentiateElement!(sbp, di, u, res)
            vtDu = sum(v.*res)
            # compute the vector-matrix-vector
            fill!(res, 0.0)
            differentiateElement_rev!(sbp, di, res, v)
            utDtv = sum(u.*res)
            @test ≈(vtDu, utDtv, atol=1e-13)
          end
        end
      end
    end
  end

  for TSBP = ((getLineSegSBPLobbato,1), (getLineSegSBPLegendre,1),
              (getTriSBPGamma,2), (getTriSBPOmega,2), (getTriSBPDiagE,2),
              (getTetSBPGamma,3), (getTetSBPOmega,3), (getTetSBPDiagE,3))
    @eval begin
      @testset "Testing differentiateElement_rev! ($(string($TSBP[1])) vector field method)" begin
        for p = 1:4
          sbp = ($TSBP[1])(degree=p)
          u = rand(Float64, (4,sbp.numnodes))
          v = rand(Float64, (4,sbp.numnodes))
          res = zeros(size(u))
          vtDu = zero(Float64)
          utDtv = zero(Float64)
          for di = 1:($TSBP[2])
            # compute the vector-matrix-vector product in forward mode
            fill!(res, 0.0)
            differentiateElement!(sbp, di, u, res)
            vtDu = sum(v.*res)
            # compute the vector-matrix-vector
            fill!(res, 0.0)
            differentiateElement_rev!(sbp, di, res, v)
            utDtv = sum(u.*res)
            @test ≈(vtDu, utDtv, atol=1e-12)
          end
        end
      end
    end
  end
  
end
