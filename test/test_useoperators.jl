facts("Testing SummationByParts Module (useoperators.jl file)...") do

  context("Testing SummationByParts.calcminnodedistance (TriSBP method)") do
    mindist = [1.0; 0.2357022603955159; 0.1487006728783353; 0.09492895652255572]
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      vtx = [0. 0.; 1. 0.; 0. 1.]
      @fact calcminnodedistance(sbp, vtx) --> roughly(mindist[p], atol=1e-13)
    end
  end

  context("Testing SummationByParts.calcminnodedistance (TetSBP method)") do
    mindist = [1.0; 0.4330127018922193; 0.2639696512367827; 0.1366241982649621]
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      @fact calcminnodedistance(sbp, vtx) --> roughly(mindist[p], atol=1e-13)
    end
  end

  context("Testing SummationByParts.weakdifferentiate! (TriSBP, scalar field method)") do
    # build a two element grid, and verify that Qxi * 1 = 0
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,2))
      di = 1
      res = zeros(u)
      weakdifferentiate!(sbp, di, u, res)
      @fact res[:,1] --> roughly(zeros(sbp.numnodes), atol=5e-13)
      @fact res[:,2] --> roughly(zeros(sbp.numnodes), atol=5e-13)
    end
  end 

  context("Testing SummationByParts.weakdifferentiate! (TetSBP, scalar field method)") do
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

  context("Testing SummationByParts.weakdifferentiate! (TriSBP, vector field method)") do
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

  context("Testing SummationByParts.weakdifferentiate! (TetSBP, vector field method)") do
    # build a single element grid, and verify that \int (Dxi * x) d\Omega = 2/3
    # or 0, depending on orientation of local coordinates
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x = zeros(Float64, (3,sbp.numnodes,1))
      x[:,:,1] = calcnodes(sbp, vtx)
      di = 1
      res = zeros(x)
      weakdifferentiate!(sbp, di, x, res)
      @fact sum(res[1,:,1]) --> roughly(2/3, atol=1e-15)
      @fact sum(res[2,:,1]) --> roughly(0.0, atol=1e-15)
      @fact sum(res[3,:,1]) --> roughly(0.0, atol=1e-15)
    end
  end

  context("Testing SummationByParts.differentiate! (TriSBP, scalar field method)") do
    # verify the accuracy of the differentiation operators
    for d = 1:4
      sbp = TriSBP{Float64}(degree=d)
      cub, vtx = tricubature(2*d-1, Float64)
      xy = calcnodes(sbp, vtx)
      x = zeros(Float64, (sbp.numnodes,1))
      y = zeros(x)
      x[:,1] = vec(xy[1,:]); y[:,1] = vec(xy[2,:])
      for r = 0:d
        for j = 0:r
          i = r-j
          u = (x.^i).*(y.^j)
          dudx = (i.*x.^max(0,i-1)).*(y.^j)          
          dudy = (x.^i).*(j.*y.^max(0,j-1))
          res = zeros(u)
          differentiate!(sbp, 1, u, res)
          @fact res --> roughly(dudx, atol=5e-13)
          res = zeros(u)
          differentiate!(sbp, 2, u, res)
          @fact res --> roughly(dudy, atol=5e-13)
        end
      end
    end
  end 

  context("Testing SummationByParts.differentiate! (TetSBP, scalar field method)") do
    # verify the accuracy of the differentiation operators
    for d = 1:4
      sbp = TetSBP{Float64}(degree=d)
      cub, vtx = tetcubature(2*d-1, Float64)
      xyz = calcnodes(sbp, vtx)
      x = zeros(Float64, (sbp.numnodes,1))
      y = zeros(x)
      z = zeros(x)
      x[:,1] = vec(xyz[1,:]); y[:,1] = vec(xyz[2,:]); z[:,1] = vec(xyz[3,:])
      for r = 0:d
        for k = 0:r
          for j = 0:r-k
            i = r-j-k
            u = (x.^i).*(y.^j).*(z.^k)
            dudx = (i.*x.^max(0,i-1)).*(y.^j).*(z.^k)
            dudy = (x.^i).*(j.*y.^max(0,j-1)).*(z.^k)
            dudz = (x.^i).*(y.^j).*(k.*z.^max(0,k-1))
            res = zeros(u)
            differentiate!(sbp, 1, u, res)
            @fact res --> roughly(dudx, atol=5e-13)
            res = zeros(u)
            differentiate!(sbp, 2, u, res)
            @fact res --> roughly(dudy, atol=5e-13)
            res = zeros(u)
            differentiate!(sbp, 3, u, res)
            @fact res --> roughly(dudz, atol=5e-13)
          end
        end
      end
    end
  end 

  context("Testing SummationByParts.differentiate! (TriSBP, vector field method)") do
    # build a two element grid, and verify that Dxi*x = 0.5 or 0, depending on
    # orientation of local coordinates
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x = zeros(Float64, (2,sbp.numnodes,2))
      x[:,:,1] = calcnodes(sbp, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = calcnodes(sbp, vtx)
      di = 1
      res = zeros(x)
      differentiate!(sbp, di, x, res)
      @fact vec(res[1,:,1]) --> roughly(0.5.*ones(sbp.numnodes), atol=5e-13)
      @fact vec(res[2,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-13)
      @fact vec(res[1,:,2]) --> roughly(zeros(sbp.numnodes), atol=5e-13)
      @fact vec(res[2,:,2]) --> roughly(0.5.*ones(sbp.numnodes), atol=5e-13)
    end
  end

  context("Testing SummationByParts.differentiate! (TetSBP, vector field method)") do
    # build a single element grid, and verify that Dxi * x = 0.5 or 0, depending
    # on orientation of local coordinates
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x = zeros(Float64, (3,sbp.numnodes,1))
      x[:,:,1] = calcnodes(sbp, vtx)
      di = 1
      res = zeros(x)
      differentiate!(sbp, di, x, res)
      @fact vec(res[1,:,1]) --> roughly(0.5.*ones(sbp.numnodes), atol=1e-13)
      @fact vec(res[2,:,1]) --> roughly(zeros(sbp.numnodes), atol=1e-13)
      @fact vec(res[3,:,1]) --> roughly(zeros(sbp.numnodes), atol=1e-13)
    end
  end

  context("Testing SummationByParts.directionaldifferentiate! (TriSBP, scalar field method)") do
    # build a single element grid, define u = x+y, and verify that Ddir = 2
    for p = 1:4 
      sbp = TriSBP{Float64}(degree=p)
      x = zeros(Float64, (2,sbp.numnodes))
      vtx = [-1. -1.; -1. 1.; 1. -1.]
      x[:,:] = calcnodes(sbp, vtx)
      u = ones(Float64, (sbp.numnodes))
      u = vec(x[1,:] + x[2,:])
      dir = [1.;1.]
      for i = 1:sbp.numnodes
        Ddir = directionaldifferentiate!(sbp, dir, u, i)
        @fact Ddir --> roughly(2.0, atol=1e-13)
      end
    end
  end

  context("Testing SummationByParts.directionaldifferentiate! (TetSBP, scalar field method)") do
    # build a single element grid, define u = x+y+z, and verify that Ddir = 3.0
    for p = 1:4 
      sbp = TetSBP{Float64}(degree=p)
      vtx = [-1. -1. -1.; 1. -1. -1.; -1. 1. -1.; -1. -1. 1.]
      x = zeros(Float64, (3,sbp.numnodes))
      x[:,:] = calcnodes(sbp, vtx)
      u = ones(Float64, (sbp.numnodes))
      u = vec(x[1,:] + x[2,:] + x[3,:])
      dir = [1.;1.;1.]
      for i = 1:sbp.numnodes
        Ddir = directionaldifferentiate!(sbp, dir, u, i)
        @fact Ddir --> roughly(3.0, atol=1e-13)
      end
    end
  end

  context("Testing SummationByParts.directionaldifferentiate! (TriSBP, vector field method)") do
    # build a single element grid, define u = x+y, and verify that Ddir = 2
    for p = 1:4 
      sbp = TriSBP{Float64}(degree=p)
      x = zeros(Float64, (2,sbp.numnodes))
      vtx = [-1. -1.; -1. 1.; 1. -1.]
      x[:,:] = calcnodes(sbp, vtx)
      u = ones(Float64, (2, sbp.numnodes))
      u[1,:] = vec(x[1,:] + x[2,:])
      u[2,:] = -u[1,:]
      dir = [1.;1.]
      Ddir = zeros(Float64, size(u,1))
      for i = 1:sbp.numnodes
        fill!(Ddir, 0.0)
        directionaldifferentiate!(sbp, dir, u, i, Ddir)
        @fact Ddir[1] --> roughly(2.0, atol=1e-13)
        @fact Ddir[2] --> roughly(-2.0, atol=1e-13)
      end
    end
  end

  context("Testing SummationByParts.directionaldifferentiate! (TetSBP, vector field method)") do
    # build a single element grid, define u = x+y+z, and verify that Ddir = 3.0
    for p = 1:4 
      sbp = TetSBP{Float64}(degree=p)
      vtx = [-1. -1. -1.; 1. -1. -1.; -1. 1. -1.; -1. -1. 1.]
      x = zeros(Float64, (3,sbp.numnodes))
      x[:,:] = calcnodes(sbp, vtx)
      u = ones(Float64, (2, sbp.numnodes))
      u[1,:] = vec(x[1,:] + x[2,:] + x[3,:])
      u[2,:] = -u[1,:]
      dir = [1.;1.;1.]
      Ddir = zeros(Float64, size(u,1))
      for i = 1:sbp.numnodes
        fill!(Ddir, 0.0)
        directionaldifferentiate!(sbp, dir, u, i, Ddir)
        @fact Ddir[1] --> roughly(3.0, atol=1e-13)
        @fact Ddir[2] --> roughly(-3.0, atol=1e-13)
      end
    end
  end

  context("Testing SummationByParts.volumeintegrate! (TriSBP, scalar field method)") do
    # build a two element grid, and verify that ones*H*ones = vol
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,2))
      res = zeros(u)
      volumeintegrate!(sbp, u, res)
      @fact sum(res[:,1]) --> roughly(2.0, atol=1e-14)
      @fact sum(res[:,2]) --> roughly(2.0, atol=1e-14)
    end
  end
  
  context("Testing SummationByParts.volumeintegrate! (TetSBP, scalar field method)") do
    # build a single element grid, and verify that ones*H*ones = vol
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (sbp.numnodes,1))
      res = zeros(u)
      volumeintegrate!(sbp, u, res)
      @fact sum(res[:,1]) --> roughly(4/3, atol=1e-14)
    end
  end
  
  context("Testing SummationByParts.volumeintegrate! (TriSBP, vector field method)") do
    # build a two element grid, and verify that ones^T*H*ones = (1,2)
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
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
  
  context("Testing SummationByParts.volumeintegrate! (TetSBP, vector field method)") do
    # build a single element grid, and verify that ones*H*ones = (1,2)
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      u = ones(Float64, (2,sbp.numnodes,1))
      u[1,:,:] *= 3/4
      u[2,:,:] *= 3/2
      res = zeros(u)
      volumeintegrate!(sbp, u, res)
      @fact sum(res[1,:,1]) --> roughly(1.0, atol=1e-14)
      @fact sum(res[2,:,1]) --> roughly(2.0, atol=1e-14)
    end
  end
    
  context("Testing SummationByParts.mappingjacobian! (TriSBP method)") do
    # build a two element grid, and verify components of the Jacobian and its
    # determinant
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x = zeros(Float64, (2,sbp.numnodes,2))
      x[:,:,1] = calcnodes(sbp, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = calcnodes(sbp, vtx)
      dξdx = zeros(Float64, (2,2,sbp.numnodes,2))
      jac = zeros(Float64, (sbp.numnodes,2))
      mappingjacobian!(sbp, x, dξdx, jac)
      # verify on element 1
      @fact vec(dξdx[1,1,:,1]) --> roughly(0.5*ones(sbp.numnodes), atol=5e-13)
      @fact vec(dξdx[1,2,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-13)
      @fact vec(dξdx[2,2,:,1]) --> roughly(0.5*ones(sbp.numnodes), atol=5e-13)
      @fact vec(dξdx[2,1,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-13)
      @fact vec(jac[:,1]) --> roughly(4.0*ones(sbp.numnodes), atol=5e-13)
      # verify on element 2
      @fact vec(dξdx[1,1,:,2]) --> roughly(0.5*ones(sbp.numnodes), atol=5e-13)
      @fact vec(dξdx[1,2,:,2]) --> roughly(0.5*ones(sbp.numnodes), atol=5e-13)
      @fact vec(dξdx[2,2,:,2]) --> roughly(zeros(sbp.numnodes), atol=5e-13)
      @fact vec(dξdx[2,1,:,2]) --> roughly(-0.5*ones(sbp.numnodes), atol=5e-13)
      @fact vec(jac[:,2]) --> roughly(4.0*ones(sbp.numnodes), atol=5e-13)
    end
  end
  
  context("Testing SummationByParts.mappingjacobian! (TetSBP method)") do
    # build one element grid, and verify components of the Jacobian and its
    # determinant
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      vtx = [0. 0. 0.; 2. 0. 0.; 0. 2. 0.; 0. 0. 2.]
      x = zeros(Float64, (3,sbp.numnodes,1))
      x[:,:,1] = calcnodes(sbp, vtx)
      dξdx = zeros(Float64, (3,3,sbp.numnodes,1))
      jac = zeros(Float64, (sbp.numnodes,1))
      mappingjacobian!(sbp, x, dξdx, jac)
      # dxi/dx = (1,0,0)
      @fact vec(dξdx[1,1,:,1]) --> roughly(ones(sbp.numnodes), atol=5e-12)
      @fact vec(dξdx[1,2,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-12)
      @fact vec(dξdx[1,3,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-12)
      # deta/dx = (0,1,0)
      @fact vec(dξdx[2,1,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-12)
      @fact vec(dξdx[2,2,:,1]) --> roughly(ones(sbp.numnodes), atol=5e-12)
      @fact vec(dξdx[2,3,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-12)
      # dzeta/dx = (0,1,0)
      @fact vec(dξdx[3,1,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-12)
      @fact vec(dξdx[3,2,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-12)
      @fact vec(dξdx[3,3,:,1]) --> roughly(ones(sbp.numnodes), atol=5e-12)
      # jac = 1
      @fact vec(jac[:,1]) --> roughly(ones(sbp.numnodes), atol=5e-12)
    end
  end
  
  # context("Testing metric invariants (TetSBP type)") do
  #   # build one element grid, and verify components of the Jacobian and its
  #   # determinant
  #   mapdegree = [1;2] 
  #   for p = 1:2 # metric invariants only satisfied for p=1 and p=2
  #     sbp = TetSBP{Float64}(degree=p)
  #     vtx = [0. 0. 0.; 2. 1. 0.; -1. 2. 0.; 0.5 0.5 2.]
  #     xi = zeros(Float64, (3,sbp.numnodes,1))
  #     xi[:,:,1] = calcnodes(sbp, vtx)
  #     for r = 0:p
  #       for k = 0:r
  #         for j = 0:r-k
  #           i = r-j-k
  #           x = deepcopy(xi)
  #           x += (xi[1,:,:].^i).*(xi[2,:,:].^j).*(xi[3,:,:].^k)
  #           dξdx = zeros(Float64, (3,3,sbp.numnodes,1))
  #           jac = zeros(Float64, (sbp.numnodes,1))
  #           mappingjacobian!(sbp, x, dξdx, jac)
  
  #           invariant = zeros(Float64, (3,sbp.numnodes,1))
  #           differentiate!(sbp, 1, slice(dξdx,1,:,:,:), invariant)
  #           differentiate!(sbp, 2, slice(dξdx,2,:,:,:), invariant)
  #           differentiate!(sbp, 3, slice(dξdx,3,:,:,:), invariant)
  #           @fact invariant --> roughly(zeros(Float64, (3,sbp.numnodes,1)), atol=1e-13)
  #         end
  #       end
  #     end
  #   end
  # end

end
  