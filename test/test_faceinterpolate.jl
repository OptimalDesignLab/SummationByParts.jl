facts("Testing SummationByParts Module (face-data interpolation methods)...") do

    context("Testing boundaryinterpolate! (TriSBP, scalar field method)") do
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

  context("Testing boundaryinterpolate! (TriSBP, vector field method)") do
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

  context("Testing boundaryinterpolate! (TetSBP, scalar field method)") do
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

  context("Testing boundaryinterpolate! (TetSBP, vector field method)") do
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

  context("Testing boundaryFaceInterpolate! (TriSBP, scalar field method)") do
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
          for (bindex, bndry) in enumerate(bndryfaces)          
            boundaryFaceInterpolate!(sbpface, bndry.face,
                                     view(u,:,bndry.element),
                                     view(uface,:,bindex))
          end
          @fact vec(uface[:,:]) -->
          roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-14)
        end
      end
    end
  end

  context("Testing boundaryFaceInterpolate! (TriSBP, vector field method)") do
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
          for (bindex, bndry) in enumerate(bndryfaces)
            boundaryFaceInterpolate!(sbpface, bndry.face,
                                     view(u,:,:,bndry.element),
                                     view(uface,:,:,bindex))
          end
          @fact vec(uface[1,:,:]) -->
          roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-14)
          @fact vec(uface[2,:,:]) -->
          roughly(2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-14)
        end
      end
    end
  end

  context("Testing boundaryFaceInterpolate! (TetSBP, scalar field method)") do
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
            for (bindex, bndry) in enumerate(bndryfaces)          
              boundaryFaceInterpolate!(sbpface, bndry.face,
                                       view(u,:,bndry.element),
                                       view(uface,:,bindex))
            end
            @fact vec(uface[:,:]) -->
            roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j).*(xf[3,:,:].^k)),
                    atol=1e-14)
          end
        end
      end
    end
  end

  context("Testing boundaryFaceInterpolate! (TetSBP, vector field method)") do
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
            for (bindex, bndry) in enumerate(bndryfaces)
              boundaryFaceInterpolate!(sbpface, bndry.face,
                                       view(u,:,:,bndry.element),
                                       view(uface,:,:,bindex))
            end
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

  context("Testing interiorfaceinterpolate! (TriSBP, scalar field method)") do
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

  context("Testing interiorfaceinterpolate! (TriSBP, vector field method)") do
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

  context("Testing interiorfaceinterpolate! (TetSBP, scalar field method)") do
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

  context("Testing interiorfaceinterpolate! (TetSBP, vector field method)") do
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

  context("Testing interiorFaceInterpolate! (TriSBP, scalar field method)") do
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
      uface = zeros(Float64, (sbpface.numnodes, 2, 1))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          for (findex, face) in enumerate(ifaces)          
            interiorFaceInterpolate!(sbpface, face, view(u,:,face.elementL),
                                     view(u,:,face.elementR),
                                     view(uface,:,1,findex),
                                     view(uface,:,2,findex))
          end
          # check that interpolation from left and right elements is exact
          @fact vec(uface[:,1,1]) --> roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)),
                                              atol=1e-13)
          #@fact uface[2,sbpface.nbrperm[:,1],1] --> 
          #roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @fact vec(uface[:,2,1]) --> roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)),
                                              atol=1e-13)
        end
      end
    end
  end

  context("Testing interiorFaceInterpolate! (TriSBP, vector field method)") do
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
      uface = zeros(Float64, (2, sbpface.numnodes, 2, 1))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j)
          for (findex, face) in enumerate(ifaces)          
            interiorFaceInterpolate!(sbpface, face, view(u,:,:,face.elementL),
                                     view(u,:,:,face.elementR),
                                     view(uface,:,:,1,findex),
                                     view(uface,:,:,2,findex))
          end
          # check that interpolation from left and right elements is exact
          @fact vec(uface[1,:,1,1]) --> roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @fact vec(uface[2,:,1,1]) --> roughly(2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @fact vec(uface[1,:,2,1]) --> roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @fact vec(uface[2,:,2,1]) --> roughly(2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
        end
      end
    end
  end

  context("Testing interiorFaceInterpolate! (TetSBP, scalar field method)") do
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
      uface = zeros(Float64, (sbpface.numnodes, 2, 4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            for (findex, face) in enumerate(ifaces)          
              interiorFaceInterpolate!(sbpface, face, view(u,:,face.elementL),
                                       view(u,:,face.elementR),
                                       view(uface,:,1,findex),
                                       view(uface,:,2,findex))
            end
            # check that interpolation from left and right elements is exact
            for f = 1:4
              @fact vec(uface[:,1,f]) --> 
              roughly(vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @fact vec(uface[:,2,f]) --> 
              roughly(vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
            end
          end
        end
      end
    end
  end

  context("Testing interiorFaceInterpolate! (TetSBP, vector field method)") do
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
      uface = zeros(Float64, (2, sbpface.numnodes, 2, 4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            for (findex, face) in enumerate(ifaces)          
              interiorFaceInterpolate!(sbpface, face, view(u,:,:,face.elementL),
                                       view(u,:,:,face.elementR),
                                       view(uface,:,:,1,findex),
                                       view(uface,:,:,2,findex))
            end
            # check that interpolation from left and right elements is exact
            for f = 1:4
              @fact vec(uface[1,:,1,f]) --> 
              roughly(vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @fact vec(uface[2,:,1,f]) --> 
              roughly(2.0.*vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @fact vec(uface[1,:,2,f]) --> 
              roughly(vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
              @fact vec(uface[2,:,2,f]) --> 
              roughly(2.0.*vec((xf[1,:,f].^i).*(xf[2,:,f].^j).*(xf[3,:,f].^k)),
                      atol=1e-13)
            end
          end
        end
      end
    end
  end

end
