facts("Testing SummationByParts Module (reverse-diff of volume integrate methods)...") do

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE,
              getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      context("Testing volumeintegrate_rev! ("string($TSBP)" scalar field method)") do
        for p = 1:4
          # verify that H*u = (u^T*H)^T
          sbp = ($TSBP)(degree=p)
          u = rand(Float64, (sbp.numnodes,2))
          Hu = zeros(u)
          utH = zeros(u)
          volumeintegrate!(sbp, u, Hu)
          volumeintegrate_rev!(sbp, utH, u)
          @fact Hu --> roughly(utH, atol=1e-15)          
        end
      end
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE,
              getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      context("Testing volumeintegrate_rev! ("string($TSBP)" vector field method)") do
        for p = 1:4
          # verify that H*u = (u^T*H)^T
          sbp = ($TSBP)(degree=p)
          u = rand(Float64, (4,sbp.numnodes,2))
          Hu = zeros(u)
          utH = zeros(u)
          volumeintegrate!(sbp, u, Hu)
          volumeintegrate_rev!(sbp, utH, u)
          @fact Hu --> roughly(utH, atol=1e-15)          
        end
      end
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE,
              getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      context("Testing volumeIntegrateElement_rev! ("string($TSBP)" scalar field method)") do
        for p = 1:4
          # verify that H*u = (u^T*H)^T
          sbp = ($TSBP)(degree=p)
          u = rand(Float64, (sbp.numnodes))
          Hu = zeros(u)
          utH = zeros(u)
          volumeIntegrateElement!(sbp, u, Hu)
          volumeIntegrateElement_rev!(sbp, utH, u)
          @fact Hu --> roughly(utH, atol=1e-15)          
        end
      end
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE,
              getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      context("Testing volumeIntegrateElement_rev! ("string($TSBP)" vector field method)") do
        for p = 1:4
          # verify that H*u = (u^T*H)^T
          sbp = ($TSBP)(degree=p)
          u = rand(Float64, (4,sbp.numnodes))
          Hu = zeros(u)
          utH = zeros(u)
          volumeIntegrateElement!(sbp, u, Hu)
          volumeIntegrateElement_rev!(sbp, utH, u)
          @fact Hu --> roughly(utH, atol=1e-15)
        end
      end
    end
  end
  
end
