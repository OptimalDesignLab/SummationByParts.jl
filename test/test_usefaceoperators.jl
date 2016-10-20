facts("Testing SummationByParts Module (usefaceoperators.jl file)...") do

  context("Testing SummationByParts.boundaryinterpolate! (TriSBP, scalar field method)") do
    # build a two element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p, reorder=false)
      sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array(Boundary, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      u = zeros(Float64, (sbp.numnodes, 2))
      uface = zeros(Float64, (sbpface.numnodes, 4))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          boundaryinterpolate!(sbpface, bndryfaces, u, uface)
          @fact vec(uface[:,:]) -->
          roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-14)
        end
      end
    end
  end

  context("Testing SummationByParts.boundaryinterpolate! (TriSBP, vector field method)") do
    # build a two element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p, reorder=false)
      sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array(Boundary, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      u = zeros(Float64, (2, sbp.numnodes, 2))
      uface = zeros(Float64, (2, sbpface.numnodes, 4))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j)
          boundaryinterpolate!(sbpface, bndryfaces, u, uface)
          @fact vec(uface[1,:,:]) -->
          roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-14)
          @fact vec(uface[2,:,:]) -->
          roughly(2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-14)
        end
      end
    end
  end

  context("Testing SummationByParts.boundaryinterpolate! (TetSBP, scalar field method)") do
    # build a four element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p, reorder=false, internal=false)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (3,sbp.numnodes,4))
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array(Boundary, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      u = zeros(Float64, (sbp.numnodes, 4))
      uface = zeros(Float64, (sbpface.numnodes, 12))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            boundaryinterpolate!(sbpface, bndryfaces, u, uface)
            @fact vec(uface[:,:]) -->
            roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=1e-14)
          end
        end
      end
    end
  end

  context("Testing SummationByParts.boundaryinterpolate! (TetSBP, vector field method)") do
    # build a four element grid and verify that interpolation is exact for degree p
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p, reorder=false, internal=false)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (3,sbp.numnodes,4))
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array(Boundary, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      u = zeros(Float64, (2, sbp.numnodes, 4))
      uface = zeros(Float64, (2, sbpface.numnodes, 12))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            boundaryinterpolate!(sbpface, bndryfaces, u, uface)
            @fact vec(uface[1,:,:]) -->
            roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=1e-14)
            @fact vec(uface[2,:,:]) -->
            roughly(2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=1e-14)
          end
        end
      end
    end
  end

  context("Testing SummationByParts.integratefunctional! (TriSBP, scalar field method)") do
    # build a two element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p, reorder=false)
      sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      #x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array(Boundary, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = zeros(Float64, (sbpface.numnodes, 4))
      for d = 0:2*p
        for j = 0:d
          i = d-j
          # the function be integrated is (x+1)^i (y+1)^j
          uface[:,:] = ((xf[1,:,:]+1).^i).*((xf[2,:,:]+1).^j)
          scale!(uface, 0.5) # 0.5 factor accounts for tranformation to ref space
          fun = integratefunctional!(sbpface, bndryfaces, uface)
          funexact = (2^(i+1)-1)*(1 + 2^j)/(i+1) + (2^(j+1)-1)*(1 + 2^i)/(j+1)
          @fact fun --> roughly(funexact, atol=1e-14)
        end
      end
    end
  end

  context("Testing SummationByParts.integratefunctional! (TriSBP, vector field method)") do
    # build a two element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p, reorder=false)
      sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,4))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      #x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[[3;1],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[[1;2],:])
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      bndryfaces = Array(Boundary, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = zeros(Float64, (2, sbpface.numnodes, 4))
      fun = zeros(2)
      for d = 0:2*p
        for j = 0:d
          i = d-j
          # the function be integrated is (x+1)^i (y+1)^j; 0.5 factor accounts
          # for the transformation to ref space
          uface[1,:,:] = 0.5*((xf[1,:,:]+1).^i).*((xf[2,:,:]+1).^j)
          uface[2,:,:] = 0.5 # integrate constant, gives perimeter
          fill!(fun, 0.0)
          integratefunctional!(sbpface, bndryfaces, uface, fun)
          funexact = (2^(i+1)-1)*(1 + 2^j)/(i+1) + (2^(j+1)-1)*(1 + 2^i)/(j+1)
          @fact fun[1] --> roughly(funexact, atol=1e-14)
          @fact fun[2] --> roughly(4.0, atol=1e-14)
        end
      end
    end
  end

  context("Testing SummationByParts.integratefunctional! (TetSBP, scalar field method)") do
    # build a four element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p, reorder=false, internal=false)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array(Boundary, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      uface = zeros(Float64, (sbpface.numnodes, 12))
      for d = 0:2*p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            # the function be integrated is (x+1)^i (y+1)^j (z+1)^k
            uface[:,:] = ((xf[1,:,:]+1).^i).*((xf[2,:,:]+1).^j).*((xf[3,:,:]+1).^k)
            scale!(uface, 0.25) # 0.25 factor accounts for tranformation to ref space
            fun = integratefunctional!(sbpface, bndryfaces, uface)
            funexact = (2^(j+1)-1)*(2^(k+1)-1)*(1 + 2^i)/((j+1)*(k+1)) +
            (2^(i+1)-1)*(2^(k+1)-1)*(1 + 2^j)/((i+1)*(k+1)) + 
            (2^(i+1)-1)*(2^(j+1)-1)*(1 + 2^k)/((i+1)*(j+1))
            @fact fun --> roughly(funexact, atol=1e-14)
          end
        end
      end
    end
  end

  context("Testing SummationByParts.integratefunctional! (TetSBP, vector field method)") do
    # build a four element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p, reorder=false, internal=false)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      xf = zeros(Float64, (3,sbpface.numnodes,12))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      bndryfaces = Array(Boundary, 12)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)

      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,5] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,6] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      xf[:,:,7] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,8] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,9] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      xf[:,:,10] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,1],:])
      xf[:,:,11] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,2],:])
      xf[:,:,12] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,4],:])
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)

      uface = zeros(Float64, (2, sbpface.numnodes, 12))
      fun = zeros(2)
      for d = 0:2*p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            # the function be integrated is (x+1)^i (y+1)^j (z+1)^k
            uface[1,:,:] = ((xf[1,:,:]+1).^i).*((xf[2,:,:]+1).^j).*((xf[3,:,:]+1).^k)
            uface[2,:,:] = 1.0 # integrate constant, gives surface area
            scale!(uface, 0.25) # 0.25 factor accounts for tranformation to ref space
            fill!(fun, 0.0)
            integratefunctional!(sbpface, bndryfaces, uface, fun)
            funexact = (2^(j+1)-1)*(2^(k+1)-1)*(1 + 2^i)/((j+1)*(k+1)) +
            (2^(i+1)-1)*(2^(k+1)-1)*(1 + 2^j)/((i+1)*(k+1)) + 
            (2^(i+1)-1)*(2^(j+1)-1)*(1 + 2^k)/((i+1)*(j+1))
            @fact fun[1] --> roughly(funexact, atol=1e-14)
            @fact fun[2] --> roughly(6.0, atol=1e-14)
          end
        end
      end
    end
  end

  context("Testing SummationByParts.interiorfaceinterpolate! (TriSBP, scalar field method)") do
    # build a two element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p, reorder=false)
      sbpface = TriFace{Float64}(p, sbp.cub, [-1. -1.; 1. -1.; -1. 1.])
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,1))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      ifaces = Array(Interface, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      u = zeros(Float64, (sbp.numnodes,2))
      uface = zeros(Float64, (2, sbpface.numnodes, 1))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          interiorfaceinterpolate!(sbpface, ifaces, u, uface)
          # check that interpolation from left and right elements is exact
          @fact vec(uface[1,:,1]) --> roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)),
                                              atol=1e-13)
          #@fact uface[2,sbpface.nbrperm[:,1],1] --> 
          #roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @fact vec(uface[2,:,1]) --> roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)),
                                              atol=1e-13)
        end
      end
    end
  end

  context("Testing SummationByParts.interiorfaceinterpolate! (TriSBP, vector field method)") do
    # build a two element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p, reorder=false)
      sbpface = TriFace{Float64}(p, sbp.cub, [-1. -1.; 1. -1.; -1. 1.])
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,1))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      ifaces = Array(Interface, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      u = zeros(Float64, (2,sbp.numnodes,2))
      uface = zeros(Float64, (2, 2, sbpface.numnodes, 1))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j)
          interiorfaceinterpolate!(sbpface, ifaces, u, uface)
          # check that interpolation from left and right elements is exact
          @fact vec(uface[1,1,:,1]) --> roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @fact vec(uface[2,1,:,1]) --> roughly(2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @fact vec(uface[1,2,:,1]) --> roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @fact vec(uface[2,2,:,1]) --> roughly(2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
        end
      end
    end
  end

  context("Testing SummationByParts.interiorfaceinterpolate! (TetSBP, scalar field method)") do
    # build a five element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p, reorder=false, internal=false)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (3,sbp.numnodes,5))
      xf = zeros(Float64, (3,sbpface.numnodes,4))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      ifaces = Array(Interface, 4)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])
            
      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 0.; 0. 0. 1.; 1. 1. 1; 0. 1. 0.]
      x[:,:,5] = SymCubatures.calcnodes(sbp.cub, vtx)

      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      u = zeros(Float64, (sbp.numnodes, 5))
      uface = zeros(Float64, (2, sbpface.numnodes, 4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            interiorfaceinterpolate!(sbpface, ifaces, u, uface)
            # check that interpolation from left and right elements is exact
            for f = 1:4
              @fact vec(uface[1,:,f]) --> 
              roughly(vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @fact vec(uface[2,:,f]) --> 
              roughly(vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
            end
          end
        end
      end
    end
  end

  context("Testing SummationByParts.interiorfaceinterpolate! (TetSBP, vector field method)") do
    # build a five element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p, reorder=false, internal=false)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      x = zeros(Float64, (3,sbp.numnodes,5))
      xf = zeros(Float64, (3,sbpface.numnodes,4))
      facevtx = SymCubatures.getfacevertexindices(sbp.cub)
      ifaces = Array(Interface, 4)

      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])
            
      vtx = [1. 1. 0.; 0. 1. 0.; 1. 0. 0.; 1. 1. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,2] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 1.; 0. 0. 1.; 1. 1. 1.; 1. 0. 0.]
      x[:,:,3] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,3] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [0. 1. 1.; 1. 1. 1.; 0. 0. 1.; 0. 1. 0.]
      x[:,:,4] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,4] = SymCubatures.calcnodes(sbpface.cub, vtx[facevtx[:,3],:])

      vtx = [1. 0. 0.; 0. 0. 1.; 1. 1. 1; 0. 1. 0.]
      x[:,:,5] = SymCubatures.calcnodes(sbp.cub, vtx)

      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      u = zeros(Float64, (2, sbp.numnodes, 5))
      uface = zeros(Float64, (2, 2, sbpface.numnodes, 4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            interiorfaceinterpolate!(sbpface, ifaces, u, uface)
            # check that interpolation from left and right elements is exact
            for f = 1:4
              @fact vec(uface[1,1,:,f]) --> 
              roughly(vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @fact vec(uface[2,1,:,f]) --> 
              roughly(2.0.*vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @fact vec(uface[1,2,:,f]) --> 
              roughly(vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @fact vec(uface[2,2,:,f]) --> 
              roughly(2.0.*vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
            end
          end
        end
      end
    end
  end

  context("Testing boundaryintegrate! and interiorfaceintegrate! (TriSBP scalar field method)") do
    # build a two element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p, reorder=false)
      sbpface = TriFace{Float64}(p, sbp.cub, [-1. -1.; 1. -1.; -1. 1.])
      ifaces = Array(Interface, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      bndryfaces = Array(Boundary, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = 0.5.*ones(Float64, (sbpface.numnodes, 1))
      ubndry = 0.5.*ones(Float64, (sbpface.numnodes, 4))
      ubndry[:,1] *= -1.0
      ubndry[:,2] *= -1.0
      res = zeros(Float64, (sbp.numnodes, 2))
      boundaryintegrate!(sbpface, bndryfaces, ubndry, res)
      interiorfaceintegrate!(sbpface, ifaces, uface, res)
      @fact sum(res) --> roughly(0.0, atol=1e-13)
    end
  end

  context("Testing boundaryintegrate! and interiorfaceintegrate! (TriSBP vector field method)") do
    # build a two element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p, reorder=false)
      sbpface = TriFace{Float64}(p, sbp.cub, [-1. -1.; 1. -1.; -1. 1.])
      ifaces = Array(Interface, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      bndryfaces = Array(Boundary, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      uface = 0.5.*ones(Float64, (2, sbpface.numnodes, 1))
      ubndry = 0.5.*ones(Float64, (2, sbpface.numnodes, 4))
      ubndry[:,:,1] *= -1.0
      ubndry[:,:,2] *= -1.0
      res = zeros(Float64, (2, sbp.numnodes, 2))
      boundaryintegrate!(sbpface, bndryfaces, ubndry, res)
      interiorfaceintegrate!(sbpface, ifaces, uface, res)
      @fact sum(res) --> roughly(0.0, atol=1e-13)
    end
  end

  context("Testing boundaryintegrate! and interiorfaceintegrate! (TetSBP scalar field method)") do
    # build a five element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4      
      sbp = TetSBP{Float64}(degree=p, reorder=false, internal=false)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      ifaces = Array(Interface, 4)
      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      bndryfaces = Array(Boundary, 12)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)
      uface = ones(Float64, (sbpface.numnodes, 4))
      ubndry = ones(Float64, (sbpface.numnodes, 12))
      uface[:,2] *= -1.0
      uface[:,3] *= -1.0
      uface[:,4] *= -1.0
      ubndry[:,1:3] *= -1.0
      ubndry[:,4] *= -1.0
      ubndry[:,8] *= -1.0
      ubndry[:,10] *= -1.0
      res = zeros(Float64, (sbp.numnodes, 5))
      boundaryintegrate!(sbpface, bndryfaces, ubndry, res)
      interiorfaceintegrate!(sbpface, ifaces, uface, res)
      @fact sum(res) --> roughly(0.0, atol=1e-13)
    end
  end

  context("Testing boundaryintegrate! and interiorfaceintegrate! (TetSBP scalar field method)") do
    # build a five element grid and verify that a constant integrated over all
    # faces is zero
    for p = 1:4      
      sbp = TetSBP{Float64}(degree=p, reorder=false, internal=false)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      ifaces = Array(Interface, 4)
      ifaces[1] = Interface(1,5,3,2,1)
      ifaces[2] = Interface(2,5,3,4,3)
      ifaces[3] = Interface(3,5,3,1,2)
      ifaces[4] = Interface(4,5,3,3,3)
      bndryfaces = Array(Boundary, 12)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,2)
      bndryfaces[3] = Boundary(1,4)
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,2)
      bndryfaces[6] = Boundary(2,4)
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,2)
      bndryfaces[9] = Boundary(3,4)
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,2)
      bndryfaces[12] = Boundary(4,4)
      uface = ones(Float64, (2, sbpface.numnodes, 4))
      ubndry = ones(Float64, (2, sbpface.numnodes, 12))
      uface[:,:,2] *= -1.0
      uface[:,:,3] *= -1.0
      uface[:,:,4] *= -1.0
      ubndry[:,:,1:3] *= -1.0
      ubndry[:,:,4] *= -1.0
      ubndry[:,:,8] *= -1.0
      ubndry[:,:,10] *= -1.0
      res = zeros(Float64, (2, sbp.numnodes, 5))
      boundaryintegrate!(sbpface, bndryfaces, ubndry, res)
      interiorfaceintegrate!(sbpface, ifaces, uface, res)
      @fact sum(res) --> roughly(0.0, atol=1e-13)
    end
  end

  context("Testing SummationByParts.mappingjacobian! (TriFace method)") do
    # build a two element grid, and verify components of the Jacobian and its
    # determinant
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p, reorder=false)
      sbpface = TriFace{Float64}(p, sbp.cub, [-1. -1.; 1. -1.; -1. 1.])
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x = zeros(Float64, (2,sbp.numnodes,2))
      x[:,:,1] = SummationByParts.SymCubatures.calcnodes(sbp.cub, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SummationByParts.SymCubatures.calcnodes(sbp.cub, vtx)
      dξdx = zeros(Float64, (2,2,sbpface.numnodes,2,1))
      jac = zeros(Float64, (sbpface.numnodes,2,1))
      ifaces = Array(Interface, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      mappingjacobian!(sbpface, ifaces, x, dξdx, jac)
      # verify on element 1
      @fact vec(dξdx[1,1,:,1,1]) --> roughly(zeros(sbpface.numnodes), atol=5e-13)
      @fact vec(dξdx[1,2,:,1,1]) --> roughly(2.0*ones(sbpface.numnodes), atol=5e-13)
      @fact vec(dξdx[2,1,:,1,1]) --> roughly(-2.0*ones(sbpface.numnodes), atol=5e-13)
      @fact vec(dξdx[2,2,:,1,1]) --> roughly(-2.0*ones(sbpface.numnodes), atol=5e-13)
      @fact vec(jac[:,1,1]) --> roughly(4.0*ones(sbpface.numnodes), atol=5e-13)
      # verify on element 2
      @fact vec(dξdx[1,1,:,2,1]) --> roughly(zeros(sbpface.numnodes), atol=5e-13)
      @fact vec(dξdx[1,2,:,2,1]) --> roughly(-2.0*ones(sbpface.numnodes), atol=5e-13)
      @fact vec(dξdx[2,1,:,2,1]) --> roughly(2.0*ones(sbpface.numnodes), atol=5e-13)
      @fact vec(dξdx[2,2,:,2,1]) --> roughly(2.0*ones(sbpface.numnodes), atol=5e-13)
      @fact vec(jac[:,2,1]) --> roughly(4.0*ones(sbpface.numnodes), atol=5e-13)
    end
  end

  context("Testing SummationByParts.facenormal! (TriFace method)") do
    # build a curvilinear element, and verify that the geometric conservation
    # law holds
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p, reorder=false)
      sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
      function mapping(ξ)
        x = 1 - (1 - 0.5*(ξ[1]+1))^(p+1)
        y = 1 - (1 - 0.5*(ξ[2]+1))^(p+1)
        fac = 1/sqrt(2)
        return [fac*x - fac*y; fac*x + fac*y]
      end
      # set the coordinates of the Lagrangian nodes in reference and physical
      # space
      xref = zeros(1,p+2)
      for i = 0:(p+1)
        xref[1,i+1] = 2*i/(p+1) - 1.0
      end
      xlag = zeros(2,p+2,3)
      for i = 0:(p+1)
        xlag[:,i+1,1] = mapping([xref[1,i+1]; -1.0])
        xlag[:,i+1,2] = mapping([xref[1,p-i+2]; xref[1,i+1]])
        xlag[:,i+1,3] = mapping([-1.0; xref[1,p-i+2]])
      end
      # get the SBP nodes and the normal vector
      xsbp = zeros(2,sbpface.numnodes,3)
      nrm = zeros(2,sbpface.numnodes,3)
      divfree = zeros(2)
      for f = 1:3
        facenormal!(sbpface, p+1, slice(xlag,:,:,f), xref,
                    slice(xsbp,:,:,f), slice(nrm,:,:,f))
        divfree[1] += dot(vec(nrm[1,:,f]),sbpface.wface)
        divfree[2] += dot(vec(nrm[2,:,f]),sbpface.wface)
      end
      @fact divfree[1] --> roughly(0.0, atol=1e-15)
      @fact divfree[2] --> roughly(0.0, atol=1e-15)
    end
  end

  context("Testing SummationByParts.facenormal! (TetFace method)") do
    # build a curvilinear element, and verify that the geometric conservation
    # law holds
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p, reorder=false)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      function mapping(ξ)
        x = 1 - (1 - 0.5*(ξ[1]+1))^(p+1)
        y = 1 - (1 - 0.5*(ξ[2]+1))^(p+1)
        z = 1 - (1 - 0.5*(ξ[3]+1))^(p+1)
        fac = 1/sqrt(2)
        return [fac*x - fac*y; fac*x + fac*y; z]
      end
      # set the coordinates of the Lagrangian nodes in reference and physical
      # space
      numdof = div((p+2)*(p+3),2)
      xref = zeros(2,numdof)
      ptr = 1
      for r = 0:p+1
        for j = 0:r
          i = r-j
          xref[1,ptr] = 2*i/(p+1) - 1.0
          xref[2,ptr] = 2*j/(p+1) - 1.0
          ptr += 1
        end
      end
      xlag = zeros(3,numdof,4)
      for i = 1:numdof
        xlag[:,i,1] = mapping([xref[1,i]; xref[2,i]; -1.0])
        xlag[:,i,2] = mapping([xref[2,i]; -1.0; xref[1,i]])
        xlag[:,i,3] = mapping([xref[2,i]; xref[1,i];
                               -1.0 - xref[1,i] - xref[2,i]])
        xlag[:,i,4] = mapping([-1.0; xref[1,i]; xref[2,i]]) 
      end
      # get the SBP nodes and the normal vector
      xsbp = zeros(3,sbpface.numnodes,4)
      nrm = zeros(3,sbpface.numnodes,4)
      divfree = zeros(3)
      for f = 1:4
        facenormal!(sbpface, p+1, slice(xlag,:,:,f), xref,
                    slice(xsbp,:,:,f), slice(nrm,:,:,f))
        divfree[1] += dot(vec(nrm[1,:,f]),sbpface.wface)
        divfree[2] += dot(vec(nrm[2,:,f]),sbpface.wface)
        divfree[3] += dot(vec(nrm[3,:,f]),sbpface.wface)
      end
      @fact divfree[1] --> roughly(0.0, atol=1e-15)
      @fact divfree[2] --> roughly(0.0, atol=1e-15)
      @fact divfree[3] --> roughly(0.0, atol=1e-15)
    end
  end

  context("Testing SummationByParts.edgestabilize! (TriFace, scalar field method)") do
    # build a two element grid, and verify that polynomials of degree p or less
    # vanish when edgestabilize is applied
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p, reorder=false)
      sbpface = TriFace{Float64}(p, sbp.cub, [-1. -1.; 1. -1.; -1. 1.])
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x = zeros(Float64, (2,sbp.numnodes,2))
      x[:,:,1] = SummationByParts.SymCubatures.calcnodes(sbp.cub, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SummationByParts.SymCubatures.calcnodes(sbp.cub, vtx)
      dξdx = zeros(Float64, (2,2,sbpface.numnodes,2,1))
      jac = zeros(Float64, (sbpface.numnodes,2,1))
      ifaces = Array(Interface, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      mappingjacobian!(sbpface, ifaces, x, dξdx, jac)
      # build the normal vector to the face; (<dξdx,dηdx> , <dηdx,dηdx>)
      dirvec = zeros(Float64, (2,sbpface.numnodes,2,1))
      for LR = 1:2
        for i = 1:sbpface.numnodes
          dirvec[1,i,LR,1] = dξdx[1,1,i,LR,1]*dξdx[2,1,i,LR,1] + 
          dξdx[1,2,i,LR,1]*dξdx[2,2,i,LR,1]
          dirvec[2,i,LR,1] = dξdx[2,1,i,LR,1]*dξdx[2,1,i,LR,1] +
          dξdx[2,2,i,LR,1]*dξdx[2,2,i,LR,1]
        end
      end
      # scaling factor of 1
      tau = ones(Float64, (sbpface.numnodes,1))
      u = zeros(Float64, (sbp.numnodes,2))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          res = zeros(u)
          edgestabilize!(sbpface, ifaces, dirvec, tau, u, res)
          @fact res --> roughly(zeros(res), atol=1e-10) # !!! note the size of atol
        end
      end
    end
  end

end