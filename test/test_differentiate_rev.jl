facts("Testing SummationByParts Module (reverse-diff of differentiate methods)...") do

  for TSBP = ((getLineSegSBPLobbato,1), (getLineSegSBPLegendre,1),
              (getTriSBPGamma,2), (getTriSBPOmega,2), (getTriSBPDiagE,2),
              (getTetSBPGamma,3), (getTetSBPOmega,3), (getTetSBPDiagE,3))
    @eval begin
      context("Testing differentiate_rev! ("string($TSBP[1])" scalar field method)") do
        for p = 1:4
          sbp = ($TSBP[1])(degree=p)
          u = rand(Float64, (sbp.numnodes,2))
          v = rand(Float64, (sbp.numnodes,2))
          res = zeros(u)
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
            @fact vtDu --> roughly(utDtv, atol=1e-15)
          end
        end
      end
    end
  end

  for TSBP = ((getLineSegSBPLobbato,1), (getLineSegSBPLegendre,1),
              (getTriSBPGamma,2), (getTriSBPOmega,2), (getTriSBPDiagE,2),
              (getTetSBPGamma,3), (getTetSBPOmega,3), (getTetSBPDiagE,3))
    @eval begin
      context("Testing differentiate_rev! ("string($TSBP[1])" vector field method)") do
        for p = 1:4
          sbp = ($TSBP[1])(degree=p)
          u = rand(Float64, (4,sbp.numnodes,2))
          v = rand(Float64, (4,sbp.numnodes,2))
          res = zeros(u)
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
            @fact vtDu --> roughly(utDtv, atol=1e-15)
          end
        end
      end
    end
  end

  for TSBP = ((getLineSegSBPLobbato,1), (getLineSegSBPLegendre,1),
              (getTriSBPGamma,2), (getTriSBPOmega,2), (getTriSBPDiagE,2),
              (getTetSBPGamma,3), (getTetSBPOmega,3), (getTetSBPDiagE,3))
    @eval begin
      context("Testing differentiateElement_rev! ("string($TSBP[1])" scalar field method)") do
        for p = 1:4
          sbp = ($TSBP[1])(degree=p)
          u = rand(Float64, (sbp.numnodes))
          v = rand(Float64, (sbp.numnodes))
          res = zeros(u)
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
            @fact vtDu --> roughly(utDtv, atol=1e-15)
          end
        end
      end
    end
  end

  for TSBP = ((getLineSegSBPLobbato,1), (getLineSegSBPLegendre,1),
              (getTriSBPGamma,2), (getTriSBPOmega,2), (getTriSBPDiagE,2),
              (getTetSBPGamma,3), (getTetSBPOmega,3), (getTetSBPDiagE,3))
    @eval begin
      context("Testing differentiateElement_rev! ("string($TSBP[1])" vector field method)") do
        for p = 1:4
          sbp = ($TSBP[1])(degree=p)
          u = rand(Float64, (4,sbp.numnodes))
          v = rand(Float64, (4,sbp.numnodes))
          res = zeros(u)
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
            @fact vtDu --> roughly(utDtv, atol=1e-15)
          end
        end
      end
    end
  end
  
end
