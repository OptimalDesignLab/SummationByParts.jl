facts("Testing SummationByParts Module (volume integrate methods)...") do

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      context("Testing volumeintegrate! ("string($TSBP)" scalar field method)") do
        # build a two element grid, and verify that ones*H*ones = vol
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes,2))
          res = zeros(u)
          volumeintegrate!(sbp, u, res)
          @fact sum(res[:,1]) --> roughly(2.0, atol=1e-14)
          @fact sum(res[:,2]) --> roughly(2.0, atol=1e-14)
        end
      end
    end
  end
  
  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing volumeintegrate! ("string($TSBP)" scalar field method)") do
        # build a two element grid, and verify that ones*H*ones = vol
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes,2))
          res = zeros(u)
          volumeintegrate!(sbp, u, res)
          @fact sum(res[:,1]) --> roughly(2.0, atol=1e-14)
          @fact sum(res[:,2]) --> roughly(2.0, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      context("Testing volumeintegrate! ("string($TSBP)" scalar field method)") do
        # build a single element grid, and verify that ones*H*ones = vol
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes,1))
          res = zeros(u)
          volumeintegrate!(sbp, u, res)
          @fact sum(res[:,1]) --> roughly(4/3, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      context("Testing volumeintegrate! ("string($TSBP)" vector field method)") do
        # build a two element grid, and verify that ones^T*H*ones = (1,2)
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (2,sbp.numnodes,2))
          u[1,:,:] *= 0.5
          res = zeros(u)
          volumeintegrate!(sbp, u, res)
          @fact sum(res[1,:,1]) --> roughly(1.0, atol=1e-14)
          @fact sum(res[2,:,1]) --> roughly(2.0, atol=1e-14)
          @fact sum(res[1,:,2]) --> roughly(1.0, atol=1e-14)
          @fact sum(res[2,:,2]) --> roughly(2.0, atol=1e-14)
        end
      end
    end
  end
  
  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing volumeintegrate! ("string($TSBP)" vector field method)") do
        # build a two element grid, and verify that ones^T*H*ones = (1,2)
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (2,sbp.numnodes,2))
          u[1,:,:] *= 0.5
          res = zeros(u)
          volumeintegrate!(sbp, u, res)
          @fact sum(res[1,:,1]) --> roughly(1.0, atol=1e-14)
          @fact sum(res[2,:,1]) --> roughly(2.0, atol=1e-14)
          @fact sum(res[1,:,2]) --> roughly(1.0, atol=1e-14)
          @fact sum(res[2,:,2]) --> roughly(2.0, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      context("Testing volumeintegrate! ("string($TSBP)" vector field method)") do
        # build a single element grid, and verify that ones*H*ones = (1,2)
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (2,sbp.numnodes,1))
          u[1,:,:] *= 3/4
          u[2,:,:] *= 3/2
          res = zeros(u)
          volumeintegrate!(sbp, u, res)
          @fact sum(res[1,:,1]) --> roughly(1.0, atol=1e-14)
          @fact sum(res[2,:,1]) --> roughly(2.0, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      context("Testing volumeIntegrateElement! ("string($TSBP)" scalar field method)") do
        # build a two element grid, and verify that ones*H*ones = vol
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes,2))
          res = zeros(u)
          volumeIntegrateElement!(sbp, view(u,:,1), view(res,:,1))
          volumeIntegrateElement!(sbp, view(u,:,2), view(res,:,2))
          @fact sum(res[:,1]) --> roughly(2.0, atol=1e-14)
          @fact sum(res[:,2]) --> roughly(2.0, atol=1e-14)
        end
      end
    end
  end
  
  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing volumeIntegrateElement! ("string($TSBP)" scalar field method)") do
        # build a two element grid, and verify that ones*H*ones = vol
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes,2))
          res = zeros(u)
          volumeIntegrateElement!(sbp, view(u,:,1), view(res,:,1))
          volumeIntegrateElement!(sbp, view(u,:,2), view(res,:,2))
          @fact sum(res[:,1]) --> roughly(2.0, atol=1e-14)
          @fact sum(res[:,2]) --> roughly(2.0, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      context("Testing volumeIntegrateElement! ("string($TSBP)" scalar field method)") do
        # build a single element grid, and verify that ones*H*ones = vol
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (sbp.numnodes,1))
          res = zeros(u)
          volumeIntegrateElement!(sbp, view(u,:,1), view(res,:,1))
          @fact sum(res[:,1]) --> roughly(4/3, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getLineSegSBPLobbato, getLineSegSBPLegendre)
    @eval begin
      context("Testing volumeIntegrateElement! ("string($TSBP)" vector field method)") do
        # build a two element grid, and verify that ones^T*H*ones = (1,2)
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (2,sbp.numnodes,2))
          u[1,:,:] *= 0.5
          res = zeros(u)
          volumeIntegrateElement!(sbp, view(u,:,:,1), view(res,:,:,1))
          volumeIntegrateElement!(sbp, view(u,:,:,2), view(res,:,:,2))
          @fact sum(res[1,:,1]) --> roughly(1.0, atol=1e-14)
          @fact sum(res[2,:,1]) --> roughly(2.0, atol=1e-14)
          @fact sum(res[1,:,2]) --> roughly(1.0, atol=1e-14)
          @fact sum(res[2,:,2]) --> roughly(2.0, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing volumeIntegrateElement! ("string($TSBP)" vector field method)") do
        # build a two element grid, and verify that ones^T*H*ones = (1,2)
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (2,sbp.numnodes,2))
          u[1,:,:] *= 0.5
          res = zeros(u)
          volumeIntegrateElement!(sbp, view(u,:,:,1), view(res,:,:,1))
          volumeIntegrateElement!(sbp, view(u,:,:,2), view(res,:,:,2))
          @fact sum(res[1,:,1]) --> roughly(1.0, atol=1e-14)
          @fact sum(res[2,:,1]) --> roughly(2.0, atol=1e-14)
          @fact sum(res[1,:,2]) --> roughly(1.0, atol=1e-14)
          @fact sum(res[2,:,2]) --> roughly(2.0, atol=1e-14)
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      context("Testing volumeIntegrateElement! ("string($TSBP)" vector field method)") do
        # build a single element grid, and verify that ones*H*ones = (1,2)
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          u = ones(Float64, (2,sbp.numnodes,1))
          u[1,:,:] *= 3/4
          u[2,:,:] *= 3/2
          res = zeros(u)
          volumeIntegrateElement!(sbp, view(u,:,:,1), view(res,:,:,1))
          @fact sum(res[1,:,1]) --> roughly(1.0, atol=1e-14)
          @fact sum(res[2,:,1]) --> roughly(2.0, atol=1e-14)
        end
      end
    end
  end
  
end
