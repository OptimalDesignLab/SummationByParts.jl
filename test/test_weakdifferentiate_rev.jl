facts("Testing SummationByParts Module (reverse-diff of weak differentiate methods)...") do

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing weakdifferentiate_rev! ("string($TSBP)" scalar field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          v = rand(Float64, (sbp.numnodes,2))
          Qv = zeros(v)
          vQ = zeros(v)
          for di = 1:2            
            # verify that v^T Q = (Q^T v)^T
            weakdifferentiate!(sbp, di, v, Qv, trans=true)
            weakdifferentiate_rev!(sbp, di, vQ, v, trans=false)
            @fact Qv --> roughly(vQ, atol=1e-15)
            # verify that v^T Q^T = (Qv)^T
            weakdifferentiate!(sbp, di, v, Qv, trans=false)
            weakdifferentiate_rev!(sbp, di, vQ, v, trans=true)
            @fact Qv --> roughly(vQ, atol=1e-15)
          end
        end
      end 
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      context("Testing weakdifferentiate_rev! ("string($TSBP)" scalar field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          v = rand(Float64, (sbp.numnodes,2))
          Qv = zeros(v)
          vQ = zeros(v)
          for di = 1:3
            # verify that v^T Q = (Q^T v)^T
            weakdifferentiate!(sbp, di, v, Qv, trans=true)
            weakdifferentiate_rev!(sbp, di, vQ, v, trans=false)
            @fact Qv --> roughly(vQ, atol=1e-15)
            # verify that v^T Q^T = (Qv)^T
            weakdifferentiate!(sbp, di, v, Qv, trans=false)
            weakdifferentiate_rev!(sbp, di, vQ, v, trans=true)
            @fact Qv --> roughly(vQ, atol=1e-15)
          end
        end
      end 
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing weakdifferentiate_rev! ("string($TSBP)" vector field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          v = rand(Float64, (4,sbp.numnodes,2))
          Qv = zeros(v)
          vQ = zeros(v)
          for di = 1:2            
            # verify that v^T Q = (Q^T v)^T
            weakdifferentiate!(sbp, di, v, Qv, trans=true)
            weakdifferentiate_rev!(sbp, di, vQ, v, trans=false)
            @fact Qv --> roughly(vQ, atol=1e-15)
            # verify that v^T Q^T = (Qv)^T
            weakdifferentiate!(sbp, di, v, Qv, trans=false)
            weakdifferentiate_rev!(sbp, di, vQ, v, trans=true)
            @fact Qv --> roughly(vQ, atol=1e-15)
          end
        end
      end 
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      context("Testing weakdifferentiate_rev! ("string($TSBP)" vector field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          v = rand(Float64, (4,sbp.numnodes,2))
          Qv = zeros(v)
          vQ = zeros(v)
          for di = 1:3
            # verify that v^T Q = (Q^T v)^T
            weakdifferentiate!(sbp, di, v, Qv, trans=true)
            weakdifferentiate_rev!(sbp, di, vQ, v, trans=false)
            @fact Qv --> roughly(vQ, atol=1e-15)
            # verify that v^T Q^T = (Qv)^T
            weakdifferentiate!(sbp, di, v, Qv, trans=false)
            weakdifferentiate_rev!(sbp, di, vQ, v, trans=true)
            @fact Qv --> roughly(vQ, atol=1e-15)
          end
        end
      end 
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing weakDifferentiateElement_rev! ("string($TSBP)" scalar field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          v = rand(Float64, (sbp.numnodes))
          Qv = zeros(v)
          vQ = zeros(v)
          for di = 1:2            
            # verify that v^T Q = (Q^T v)^T
            weakDifferentiateElement!(sbp, di, v, Qv, SummationByParts.Add(), true)
            weakDifferentiateElement_rev!(sbp, di, vQ, v, SummationByParts.Add(), false)
            @fact Qv --> roughly(vQ, atol=1e-15)
            # verify that v^T Q^T = (Qv)^T
            weakDifferentiateElement!(sbp, di, v, Qv, SummationByParts.Add(), false)
            weakDifferentiateElement_rev!(sbp, di, vQ, v, SummationByParts.Add(), true)
            @fact Qv --> roughly(vQ, atol=1e-15)
          end
        end
      end 
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      context("Testing weakDifferentiateElement_rev! ("string($TSBP)" scalar field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          v = rand(Float64, (sbp.numnodes))
          Qv = zeros(v)
          vQ = zeros(v)
          for di = 1:3
            # verify that v^T Q = (Q^T v)^T
            weakDifferentiateElement!(sbp, di, v, Qv, SummationByParts.Add(), true)
            weakDifferentiateElement_rev!(sbp, di, vQ, v, SummationByParts.Add(), false)
            @fact Qv --> roughly(vQ, atol=1e-15)
            # verify that v^T Q^T = (Qv)^T
            weakDifferentiateElement!(sbp, di, v, Qv, SummationByParts.Add(), false)
            weakDifferentiateElement_rev!(sbp, di, vQ, v, SummationByParts.Add(), true)
            @fact Qv --> roughly(vQ, atol=1e-15)
          end
        end
      end 
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing weakDifferentiateElement_rev! ("string($TSBP)" vector field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          v = rand(Float64, (4,sbp.numnodes))
          Qv = zeros(v)
          vQ = zeros(v)
          for di = 1:2            
            # verify that v^T Q = (Q^T v)^T
            weakDifferentiateElement!(sbp, di, v, Qv, SummationByParts.Add(), true)
            weakDifferentiateElement_rev!(sbp, di, vQ, v, SummationByParts.Add(), false)
            @fact Qv --> roughly(vQ, atol=1e-15)
            # verify that v^T Q^T = (Qv)^T
            weakDifferentiateElement!(sbp, di, v, Qv, SummationByParts.Add(), false)
            weakDifferentiateElement_rev!(sbp, di, vQ, v, SummationByParts.Add(), true)
            @fact Qv --> roughly(vQ, atol=1e-15)
          end
        end
      end 
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      context("Testing weakDifferentiateElement_rev! ("string($TSBP)" vector field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          v = rand(Float64, (4,sbp.numnodes))
          Qv = zeros(v)
          vQ = zeros(v)
          for di = 1:2            
            # verify that v^T Q = (Q^T v)^T
            weakDifferentiateElement!(sbp, di, v, Qv, SummationByParts.Add(), true)
            weakDifferentiateElement_rev!(sbp, di, vQ, v, SummationByParts.Add(), false)
            @fact Qv --> roughly(vQ, atol=1e-15)
            # verify that v^T Q^T = (Qv)^T
            weakDifferentiateElement!(sbp, di, v, Qv, SummationByParts.Add(), false)
            weakDifferentiateElement_rev!(sbp, di, vQ, v, SummationByParts.Add(), true)
            @fact Qv --> roughly(vQ, atol=1e-15)
          end
        end
      end 
    end
  end
  
end
