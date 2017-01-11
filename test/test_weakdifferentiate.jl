facts("Testing SummationByParts Module (weak differentiate methods)...") do

    for TSBP = (TriSBP, SparseTriSBP)
    @eval begin
      context("Testing weakdifferentiate! ("string($TSBP)" scalar field method)") do
        # build a two element grid, and verify that Qxi * 1 = 0
        for p = 1:4
          sbp = ($TSBP){Float64}(degree=p)
          u = ones(Float64, (sbp.numnodes,2))
          di = 1
          res = zeros(u)
          weakdifferentiate!(sbp, di, u, res)
          @fact res[:,1] --> roughly(zeros(sbp.numnodes), atol=5e-13)
          @fact res[:,2] --> roughly(zeros(sbp.numnodes), atol=5e-13)
        end
      end 
    end
  end

  context("Testing weakdifferentiate! (TetSBP, scalar field method)") do
    # build a single element grid, and verify that Qxi * 1 = 0
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,1))
      di = 1
      res = zeros(u)
      weakdifferentiate!(sbp, di, u, res)
      @fact res[:,1] --> roughly(zeros(sbp.numnodes), atol=1e-13)
    end
  end 

  context("Testing weakdifferentiate! (TriSBP, vector field method)") do
    # build a two element grid, and verify that \int (Dxi * x) d\Omega = 1 or 0,
    # depending on orientation of local coordinates
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x = zeros(Float64, (2,sbp.numnodes,2))
      x[:,:,1] = calcnodes(sbp, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = calcnodes(sbp, vtx)
      di = 1
      res = zeros(x)
      weakdifferentiate!(sbp, di, x, res)
      @fact sum(res[1,:,1]) --> roughly(1.0, atol=1e-14)
      @fact sum(res[2,:,1]) --> roughly(0.0, atol=1e-14)
      @fact sum(res[1,:,2]) --> roughly(0.0, atol=1e-14)
      @fact sum(res[2,:,2]) --> roughly(1.0, atol=1e-14)
    end
  end

  context("Testing weakdifferentiate! (TetSBP, vector field method)") do
    # build a single element grid, and verify that \int (D * x) d\Omega = 2/3
    # or 0, depending on orientation of local coordinates
    Id = (2/3).*eye(3,3)
    vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
    for di = 1:3
      for p = 1:4
        sbp = TetSBP{Float64}(degree=p)
        x = zeros(Float64, (3,sbp.numnodes,1))
        x[:,:,1] = calcnodes(sbp, vtx)
        res = zeros(x)
        weakdifferentiate!(sbp, di, x, res)
        @fact sum(res[1,:,1]) --> roughly(Id[di,1], atol=1e-14)
        @fact sum(res[2,:,1]) --> roughly(Id[di,2], atol=1e-14)
        @fact sum(res[3,:,1]) --> roughly(Id[di,3], atol=1e-14)
      end
    end
  end

  for TSBP = (TriSBP, SparseTriSBP)
    @eval begin
      context("Testing weakDifferentiateElement! ("string($TSBP)" scalar field method)") do
        # verify that Qxi * 1 = 0
        for p = 1:4
          sbp = ($TSBP){Float64}(degree=p)
          u = ones(Float64, (sbp.numnodes))
          di = 1
          res = zeros(u)
          weakDifferentiateElement!(sbp, di, u, res)
          @fact res[:] --> roughly(zeros(sbp.numnodes), atol=5e-13)
        end
      end 
    end
  end

  context("Testing weakDifferentiateElement! (TetSBP, scalar field method)") do
    # build a single element, and verify that Qxi * 1 = 0
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes))
      di = 1
      res = zeros(u)
      weakDifferentiateElement!(sbp, di, u, res)
      @fact res[:] --> roughly(zeros(sbp.numnodes), atol=1e-13)
    end
  end 

  context("Testing weakDifferentiateElement! (TriSBP, vector field method)") do
    # build a two element grid, and verify that \int (Dxi * x) d\Omega = 1 or 0,
    # depending on orientation of local coordinates
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x = zeros(Float64, (2,sbp.numnodes,2))
      x[:,:,1] = calcnodes(sbp, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = calcnodes(sbp, vtx)
      di = 1
      res = zeros(x)
      weakDifferentiateElement!(sbp, di, view(x,:,:,1), view(res,:,:,1))
      weakDifferentiateElement!(sbp, di, view(x,:,:,2), view(res,:,:,2))
      @fact sum(res[1,:,1]) --> roughly(1.0, atol=1e-14)
      @fact sum(res[2,:,1]) --> roughly(0.0, atol=1e-14)
      @fact sum(res[1,:,2]) --> roughly(0.0, atol=1e-14)
      @fact sum(res[2,:,2]) --> roughly(1.0, atol=1e-14)
    end
  end

  context("Testing weakDifferentiateElement! (TetSBP, vector field method)") do
    # build a single element grid, and verify that \int (D * x) d\Omega = 2/3
    # or 0, depending on orientation of local coordinates
    Id = (2/3).*eye(3,3)
    vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
    for di = 1:3
      for p = 1:4
        sbp = TetSBP{Float64}(degree=p)
        x = zeros(Float64, (3,sbp.numnodes,1))
        x[:,:,1] = calcnodes(sbp, vtx)
        res = zeros(x)
        weakDifferentiateElement!(sbp, di, view(x,:,:,1), view(res,:,:,1))
        @fact sum(res[1,:,1]) --> roughly(Id[di,1], atol=1e-14)
        @fact sum(res[2,:,1]) --> roughly(Id[di,2], atol=1e-14)
        @fact sum(res[3,:,1]) --> roughly(Id[di,3], atol=1e-14)
      end
    end
  end
  
end
