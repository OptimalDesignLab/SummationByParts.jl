facts("Testing SummationByParts Module (useoperators.jl file)...") do

  context("Testing SummationByParts.applyQ! (TriSymCub, scalar field method)") do
    # build a two element grid, and verify that Qxi * 1 = 0
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,2))
      di = 1
      res = zeros(u)
      applyQ!(sbp, di, u, res)
      @fact res[:,1] => roughly(zeros(sbp.numnodes), atol=1e-13)
      @fact res[:,2] => roughly(zeros(sbp.numnodes), atol=1e-13)
    end
  end 

  context("Testing SummationByParts.applyQ! (TetSymCub, scalar field method)") do
    # build a single element grid, and verify that Qxi * 1 = 0
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,1))
      di = 1
      res = zeros(u)
      applyQ!(sbp, di, u, res)
      @fact res[:,1] => roughly(zeros(sbp.numnodes), atol=1e-13)
    end
  end 

  context("Testing SummationByParts.applyQ! (TriSymCub, vector field method)") do
    # build a two element grid, and verify that \int (Dxi * x) d\Omega = 1 or 0,
    # depending on orientation of local coordinates
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x = zeros(Float64, (2,sbp.numnodes,2))
      x[1,:,1], x[2,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[1,:,2], x[2,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      di = 1
      res = zeros(x)
      applyQ!(sbp, di, x, res)
      @fact sum(res[1,:,1]) => roughly(1.0, atol=1e-15)
      @fact sum(res[2,:,1]) => roughly(0.0, atol=1e-15)
      @fact sum(res[1,:,2]) => roughly(0.0, atol=1e-15)
      @fact sum(res[2,:,2]) => roughly(1.0, atol=1e-15)
    end
  end

  context("Testing SummationByParts.applyQ! (TetSymCub, vector field method)") do
    # build a single element grid, and verify that \int (Dxi * x) d\Omega = 2/3
    # or 0, depending on orientation of local coordinates
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x = zeros(Float64, (3,sbp.numnodes,1))
      x[1,:,1], x[2,:,1], x[3,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      di = 1
      res = zeros(x)
      applyQ!(sbp, di, x, res)
      @fact sum(res[1,:,1]) => roughly(2/3, atol=1e-15)
      @fact sum(res[2,:,1]) => roughly(0.0, atol=1e-15)
      @fact sum(res[3,:,1]) => roughly(0.0, atol=1e-15)
    end
  end

  context("Testing SummationByParts.applyD! (TriSymCub, scalar field method)") do
    # build a two element grid, and verify that Dxi * 1 = 0
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,2))
      di = 1
      res = zeros(u)
      applyD!(sbp, di, u, res)
      @fact res[:,1] => roughly(zeros(sbp.numnodes), atol=1e-13)
      @fact res[:,2] => roughly(zeros(sbp.numnodes), atol=1e-13)
    end
  end 

  context("Testing SummationByParts.applyD! (TetSymCub, scalar field method)") do
    # build a single element grid, and verify that Dxi * 1 = 0
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,1))
      di = 1
      res = zeros(u)
      applyD!(sbp, di, u, res)
      @fact res[:,1] => roughly(zeros(sbp.numnodes), atol=5e-13)
    end
  end 

  context("Testing SummationByParts.applyD! (TriSymCub, vector field method)") do
    # build a two element grid, and verify that Dxi*x = 0.5 or 0, depending on
    # orientation of local coordinates
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x = zeros(Float64, (2,sbp.numnodes,2))
      x[1,:,1], x[2,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[1,:,2], x[2,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      di = 1
      res = zeros(x)
      applyD!(sbp, di, x, res)
      @fact res[1,:,1] => roughly(0.5.*ones(1,sbp.numnodes), atol=1e-13)
      @fact res[2,:,1] => roughly(zeros(1,sbp.numnodes), atol=1e-13)
      @fact res[1,:,2] => roughly(zeros(1,sbp.numnodes), atol=1e-13)
      @fact res[2,:,2] => roughly(0.5.*ones(1,sbp.numnodes), atol=1e-13)
    end
  end

  context("Testing SummationByParts.applyD! (TetSymCub, vector field method)") do
    # build a single element grid, and verify that Dxi * x = 0.5 or 0, depending
    # on orientation of local coordinates
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x = zeros(Float64, (3,sbp.numnodes,1))
      x[1,:,1], x[2,:,1], x[3,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      di = 1
      res = zeros(x)
      applyD!(sbp, di, x, res)
      @fact res[1,:,1] => roughly(0.5.*ones(1,sbp.numnodes), atol=1e-13)
      @fact res[2,:,1] => roughly(zeros(1,sbp.numnodes), atol=1e-13)
      @fact res[3,:,1] => roughly(zeros(1,sbp.numnodes), atol=1e-13)
    end
  end

  context("Testing SummationByParts.applyH! (TriSymCub, scalar field method)") do
    # build a two element grid, and verify that ones*H*ones = vol
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,2))
      res = zeros(u)
      applyH!(sbp, u, res)
      @fact sum(res[:,1]) => roughly(2.0, atol=1e-14)
      @fact sum(res[:,2]) => roughly(2.0, atol=1e-14)
    end
  end

  context("Testing SummationByParts.applyH! (TetSymCub, scalar field method)") do
    # build a single element grid, and verify that ones*H*ones = vol
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,1))
      res = zeros(u)
      applyH!(sbp, u, res)
      @fact sum(res[:,1]) => roughly(4/3, atol=1e-14)
    end
  end
  
  context("Testing SummationByParts.applyH! (TriSymCub, vector field method)") do
    # build a two element grid, and verify that ones^T*H*ones = (1,2)
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      u = ones(Float64, (2,sbp.numnodes,2))
      u[1,:,:] *= 0.5
      res = zeros(u)
      applyH!(sbp, u, res)
      @fact sum(res[1,:,1]) => roughly(1.0, atol=1e-14)
      @fact sum(res[2,:,1]) => roughly(2.0, atol=1e-14)
      @fact sum(res[1,:,2]) => roughly(1.0, atol=1e-14)
      @fact sum(res[2,:,2]) => roughly(2.0, atol=1e-14)
    end
  end

  context("Testing SummationByParts.applyH! (TetSymCub, vector field method)") do
    # build a single element grid, and verify that ones*H*ones = (1,2)
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (2,sbp.numnodes,1))
      u[1,:,:] *= 3/4
      u[2,:,:] *= 3/2
      res = zeros(u)
      applyH!(sbp, u, res)
      @fact sum(res[1,:,1]) => roughly(1.0, atol=1e-14)
      @fact sum(res[2,:,1]) => roughly(2.0, atol=1e-14)
    end
  end

  context("Testing SummationByParts.mappingjacobian! (TriSymCub method)") do
    # build a two element grid, and verify components of the Jacobian and its
    # determinant
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x = zeros(Float64, (2,sbp.numnodes,2))
      x[1,:,1], x[2,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[1,:,2], x[2,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      dxidx = zeros(Float64, (2,2,sbp.numnodes,2))
      jac = zeros(Float64, (sbp.numnodes,2))
      mappingjacobian!(sbp, x, dxidx, jac)
      # verify on element 1
      @fact dxidx[1,1,:,1] => roughly(0.5*ones(1,1,sbp.numnodes), atol=1e-13)
      @fact dxidx[1,2,:,1] => roughly(zeros(1,1,sbp.numnodes), atol=1e-13)
      @fact dxidx[2,2,:,1] => roughly(0.5*ones(1,1,sbp.numnodes), atol=1e-13)
      @fact dxidx[2,1,:,1] => roughly(zeros(1,1,sbp.numnodes), atol=1e-13)
      @fact jac[:,1] => roughly(0.25*ones(sbp.numnodes), atol=1e-13)
      # verify on element 2
      @fact dxidx[1,1,:,2] => roughly(0.5*ones(1,1,sbp.numnodes), atol=1e-13)
      @fact dxidx[1,2,:,2] => roughly(0.5*ones(1,1,sbp.numnodes), atol=1e-13)
      @fact dxidx[2,2,:,2] => roughly(zeros(1,1,sbp.numnodes), atol=1e-13)
      @fact dxidx[2,1,:,2] => roughly(-0.5*ones(1,1,sbp.numnodes), atol=1e-13)
      @fact jac[:,2] => roughly(0.25*ones(sbp.numnodes), atol=1e-13)
    end
  end

end