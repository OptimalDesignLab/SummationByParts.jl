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
      @fact res[:,1] --> roughly(zeros(sbp.numnodes), atol=1e-13)
      @fact res[:,2] --> roughly(zeros(sbp.numnodes), atol=1e-13)
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
      @fact sum(res[1,:,1]) --> roughly(1.0, atol=1e-15)
      @fact sum(res[2,:,1]) --> roughly(0.0, atol=1e-15)
      @fact sum(res[1,:,2]) --> roughly(0.0, atol=1e-15)
      @fact sum(res[2,:,2]) --> roughly(1.0, atol=1e-15)
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
          @fact res --> roughly(dudx, atol=1e-13)
          res = zeros(u)
          differentiate!(sbp, 2, u, res)
          @fact res --> roughly(dudy, atol=1e-13)
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
      @fact vec(res[1,:,1]) --> roughly(0.5.*ones(sbp.numnodes), atol=1e-13)
      @fact vec(res[2,:,1]) --> roughly(zeros(sbp.numnodes), atol=1e-13)
      @fact vec(res[1,:,2]) --> roughly(zeros(sbp.numnodes), atol=1e-13)
      @fact vec(res[2,:,2]) --> roughly(0.5.*ones(sbp.numnodes), atol=1e-13)
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
      @fact vec(dξdx[1,1,:,1]) --> roughly(0.5*ones(sbp.numnodes), atol=1e-13)
      @fact vec(dξdx[1,2,:,1]) --> roughly(zeros(sbp.numnodes), atol=1e-13)
      @fact vec(dξdx[2,2,:,1]) --> roughly(0.5*ones(sbp.numnodes), atol=1e-13)
      @fact vec(dξdx[2,1,:,1]) --> roughly(zeros(sbp.numnodes), atol=1e-13)
      @fact vec(jac[:,1]) --> roughly(4.0*ones(sbp.numnodes), atol=1e-13)
      # verify on element 2
      @fact vec(dξdx[1,1,:,2]) --> roughly(0.5*ones(sbp.numnodes), atol=1e-13)
      @fact vec(dξdx[1,2,:,2]) --> roughly(0.5*ones(sbp.numnodes), atol=1e-13)
      @fact vec(dξdx[2,2,:,2]) --> roughly(zeros(sbp.numnodes), atol=1e-13)
      @fact vec(dξdx[2,1,:,2]) --> roughly(-0.5*ones(sbp.numnodes), atol=1e-13)
      @fact vec(jac[:,2]) --> roughly(4.0*ones(sbp.numnodes), atol=1e-13)
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

  context("Testing SummationByParts.boundaryintegrate! (TriSBP, scalar field method)") do
    # build a two element grid, and verify the accuracy of the boundary integration
    function bndryflux(u, dξdx, nrm)
      return u*sum(nrm.'*dξdx)
    end
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      x = zeros(Float64, (2,sbp.numnodes,2))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = calcnodes(sbp, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = calcnodes(sbp, vtx)
      dξdx = zeros(Float64, (2,2,sbp.numnodes,2))
      jac = zeros(Float64, (sbp.numnodes,2))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)

      u = zeros(Float64, (sbp.numnodes,2))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          res = zeros(u)
          flux = zeros(sbp.numfacenodes, size(bndryfaces,1))
          for k = 1:size(bndryfaces,1)
            kB = bndryfaces[k].element
            for l = 1:sbp.numfacenodes
              lB = sbp.facenodes[l, bndryfaces[k].face]
              flux[l,k] = bndryflux(u[lB,kB], dξdx[:,:,lB,kB],
                                    sbp.facenormal[:,bndryfaces[k].face])
            end
          end
          boundaryintegrate!(sbp, bndryfaces, flux, res)
          exact = 0.0
          if i == 0 && j != 0
            exact = 1.
          elseif j == 0 && i != 0
            exact = 1.
          elseif i != 0 && j != 0
            exact = 1/(j+1) + 1/(i+1)
          end
          #println("i,j = ",i,",",j,": i+j = ",i+j)
          @fact sum(res) --> roughly(exact, atol=1e-13)
        end
      end
    end
  end

  context("Testing SummationByParts.boundaryintegrate! (TetSBP, scalar field method)") do
    # build a two element grid, and verify the accuracy of the boundary integration
    function bndryflux(u, dξdx, nrm)
      return u*sum(nrm.'*dξdx)
    end
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      x = zeros(Float64, (3,sbp.numnodes,4))
      vtx = Float64[0 0 0; 1 0 0; 0 1 0; 0 0 1]
      x[:,:,1] = calcnodes(sbp, vtx)
      vtx = Float64[1 0 1; 0 0 1; 1 1 1; 1 0 0]
      x[:,:,2] = calcnodes(sbp, vtx)
      vtx = Float64[1 1 0; 0 1 0; 1 0 0; 1 1 1]
      x[:,:,3] = calcnodes(sbp, vtx)
      vtx = Float64[0 1 1; 1 1 1; 0 0 1; 0 1 0]
      x[:,:,4] = calcnodes(sbp, vtx)

      dξdx = zeros(Float64, (3,3,sbp.numnodes,4))
      jac = zeros(Float64, (sbp.numnodes,4))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, 12)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(1,4)
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,3)
      bndryfaces[6] = Boundary(2,4)
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,3)
      bndryfaces[9] = Boundary(3,4)
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,3)
      bndryfaces[12] = Boundary(4,4)

      u = zeros(Float64, (sbp.numnodes,4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            res = zeros(u)
            flux = zeros(sbp.numfacenodes, size(bndryfaces,1))
            for bindex = 1:size(bndryfaces,1)
              kB = bndryfaces[bindex].element
              for l = 1:sbp.numfacenodes
                lB = sbp.facenodes[l, bndryfaces[bindex].face]
                flux[l,bindex] = bndryflux(u[lB,kB], dξdx[:,:,lB,kB],
                                           sbp.facenormal[:,bndryfaces[bindex].face])
              end
            end
            boundaryintegrate!(sbp, bndryfaces, flux, res)
            exact = 0.0
            if !(i == j == k == 0)
              i != 0 ? exact += 1/((j+1)*(k+1)) : nothing
              j != 0 ? exact += 1/((i+1)*(k+1)) : nothing
              k != 0 ? exact += 1/((i+1)*(j+1)) : nothing
            end
            #println("i,j,k = ",i,",",j,",",k," : i+j+k = ",i+j+k)
            @fact sum(res) --> roughly(exact, atol=1e-13)
          end
        end
      end
    end
  end

  context("Testing SummationByParts.boundaryintegrate! (TriSBP, vector field method)") do
    # build a two element grid, and verify the accuracy of the boundary integration
    function bndryflux{T}(u::AbstractArray{T,1}, dξdx::AbstractArray{T,2}, 
                          nrm::AbstractArray{T,1}, flux::AbstractArray{T,1})
      tmp = sum(nrm.'*dξdx)
      for field = 1:size(u,1)
        flux[field] = u[field]*tmp
      end
    end

    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      x = zeros(Float64, (2,sbp.numnodes,2))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = calcnodes(sbp, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = calcnodes(sbp, vtx)
      dξdx = zeros(Float64, (2,2,sbp.numnodes,2))
      jac = zeros(Float64, (sbp.numnodes,2))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)

      u = zeros(Float64, (1,sbp.numnodes,2))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          res = zeros(u)
          flux = zeros(1, sbp.numfacenodes, size(bndryfaces,1))
          for k = 1:size(bndryfaces,1)
            kB = bndryfaces[k].element
            for l = 1:sbp.numfacenodes
              lB = sbp.facenodes[l, bndryfaces[k].face]
              bndryflux(u[:,lB,kB], dξdx[:,:,lB,kB],
                        sbp.facenormal[:,bndryfaces[k].face],
                        view(flux,:,l,k))
            end
          end
          boundaryintegrate!(sbp, bndryfaces, flux, res)
          exact = 0.0
          if i == 0 && j != 0
            exact = 1.
          elseif j == 0 && i != 0
            exact = 1.
          elseif i != 0 && j != 0
            exact = 1/(j+1) + 1/(i+1)
          end
          #println("i,j = ",i,",",j,": i+j = ",i+j)
          @fact sum(res) --> roughly(exact, atol=1e-13)
        end
      end
    end
  end

  context("Testing SummationByParts.boundaryintegrate! (TetSBP, vector field method)") do
    # build a two element grid, and verify the accuracy of the boundary integration
    function bndryflux{T}(u::AbstractArray{T,1}, dξdx::AbstractArray{T,2}, 
                          nrm::AbstractArray{T,1}, flux::AbstractArray{T,1})
      tmp = sum(nrm.'*dξdx)
      for field = 1:size(u,1)
        flux[field] = u[field]*tmp
      end
    end
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      x = zeros(Float64, (3,sbp.numnodes,4))
      vtx = Float64[0 0 0; 1 0 0; 0 1 0; 0 0 1]
      x[:,:,1] = calcnodes(sbp, vtx)
      vtx = Float64[1 0 1; 0 0 1; 1 1 1; 1 0 0]
      x[:,:,2] = calcnodes(sbp, vtx)
      vtx = Float64[1 1 0; 0 1 0; 1 0 0; 1 1 1]
      x[:,:,3] = calcnodes(sbp, vtx)
      vtx = Float64[0 1 1; 1 1 1; 0 0 1; 0 1 0]
      x[:,:,4] = calcnodes(sbp, vtx)

      dξdx = zeros(Float64, (3,3,sbp.numnodes,4))
      jac = zeros(Float64, (sbp.numnodes,4))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, 12)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(1,4)
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,3)
      bndryfaces[6] = Boundary(2,4)
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,3)
      bndryfaces[9] = Boundary(3,4)
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,3)
      bndryfaces[12] = Boundary(4,4)

      u = zeros(Float64, (1,sbp.numnodes,4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            res = zeros(u)
            flux = zeros(1,sbp.numfacenodes, size(bndryfaces,1))
            for bindex = 1:size(bndryfaces,1)
              kB = bndryfaces[bindex].element
              for l = 1:sbp.numfacenodes
                lB = sbp.facenodes[l, bndryfaces[bindex].face]
                bndryflux(u[:,lB,kB], dξdx[:,:,lB,kB],
                          sbp.facenormal[:,bndryfaces[bindex].face],
                          view(flux,:,l,bindex))
              end
            end
            boundaryintegrate!(sbp, bndryfaces, flux, res)
            exact = 0.0
            if !(i == j == k == 0)
              i != 0 ? exact += 1/((j+1)*(k+1)) : nothing
              j != 0 ? exact += 1/((i+1)*(k+1)) : nothing
              k != 0 ? exact += 1/((i+1)*(j+1)) : nothing
            end
            #println("i,j,k = ",i,",",j,",",k," : i+j+k = ",i+j+k)
            @fact sum(res) --> roughly(exact, atol=1e-13)
          end
        end
      end
    end
  end

  context("Testing SummationByParts.interiorfaceintegrate! (TriSBP, scalar field method)") do
    function fluxfunc(uL, uR, dξdxL, dξdxR, jacL, jacR, αL, αR, nrmL, nrmR)
      flux = sum(nrmL.'*dξdxL)
      flux >= zero(flux) ? flux *= uL : flux *= uR
      return flux
    end
    function bndryflux(u, dξdx, nrm)
      return u*sum(nrm.'*dξdx)
    end
    # build a two element grid and verify that interiorfaceintegrate does
    # nothing when given a continuous linear field
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      x = zeros(Float64, (2,sbp.numnodes,2))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = calcnodes(sbp, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = calcnodes(sbp, vtx)
      dξdx = zeros(Float64, (2,2,sbp.numnodes,2))
      jac = zeros(Float64, (sbp.numnodes,2))
      mappingjacobian!(sbp, x, dξdx, jac)
      α = zeros(dξdx)
      for k = 1:2
        for i = 1:sbp.numnodes
          for di1 = 1:2
            for di2 = 1:2
              α[di1,di2,i,k] = (dξdx[di1,1,i,k].*dξdx[di2,1,i,k] + 
                                dξdx[di1,2,i,k].*dξdx[di2,2,i,k])*jac[i,k]
            end
          end
        end
      end
      ifaces = Array(Interface, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      bndryfaces = Array(Boundary, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      u = zeros(Float64, (sbp.numnodes,2))
      u[:,:] = x[1,:,:] + x[2,:,:]
      res = zeros(u)
      Fξ = zeros(u)
      Fη = zeros(u)
      for k = 1:2
        for i = 1:sbp.numnodes
          Fξ[i,k] = u[i,k]*(dξdx[1,1,i,k] + dξdx[1,2,i,k])
          Fη[i,k] = u[i,k]*(dξdx[2,1,i,k] + dξdx[2,2,i,k])
        end
      end
      weakdifferentiate!(sbp, 1, Fξ, res, trans=true)
      weakdifferentiate!(sbp, 2, Fη, res, trans=true)
      res *= -1.0

      # compute boundary fluxes and integrate over boundary
      bflux = zeros(sbp.numfacenodes,size(bndryfaces,1))
      for findex = 1:size(bndryfaces,1)
        for i = 1:sbp.numfacenodes
          iB = sbp.facenodes[i, bndryfaces[findex].face]
          bflux[i,findex] = 
          bndryflux(u[iB,bndryfaces[findex].element], 
                    view(dξdx,:,:,iB,bndryfaces[findex].element),
                    view(sbp.facenormal,:,bndryfaces[findex].face))
        end
      end
      boundaryintegrate!(sbp, bndryfaces, bflux, res)

      # compute interior fluxes and integrate over interior faces
      flux = zeros(sbp.numfacenodes,size(ifaces,1))
      for findex = 1:size(ifaces,1)
        for i = 1:sbp.numfacenodes
          iL = sbp.facenodes[i, ifaces[findex].faceL]
          iR = sbp.facenodes[sbp.numfacenodes-i+1, ifaces[findex].faceR]
          flux[i,findex] = 
          fluxfunc(u[iL,ifaces[findex].elementL], u[iR,ifaces[findex].elementR],
                   view(dξdx,:,:,iL,ifaces[findex].elementL),
                   view(dξdx,:,:,iR,ifaces[findex].elementR),
                   jac[iL,ifaces[findex].elementL], jac[iR,ifaces[findex].elementR],
                   view(α,:,:,iL,ifaces[findex].elementL),
                   view(α,:,:,iR,ifaces[findex].elementR),
                   view(sbp.facenormal,:,ifaces[findex].faceL),
                   view(sbp.facenormal,:,ifaces[findex].faceR))
        end
      end
      interiorfaceintegrate!(sbp, ifaces, flux, res)
      for k = 1:2
        for i = 1:sbp.numnodes
          res[i,k] /= (sbp.w[i]/jac[i,k])
        end
      end
      @fact res --> roughly(2.0*ones(res), atol=1e-11)
    end
  end 

end
  