facts("Testing SummationByParts Module (reverse-diff of face-data integration methods)...") do
  
  context("Testing integratefunctional_rev! (TriSBP, vector field method)") do
    # build a two element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
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
      eps_cmplx = 1e-60
      uface_bar = zeros(2, sbpface.numnodes, 4)
      uface_cmplx = zeros(Complex128, (2, sbpface.numnodes, 4))
      fun_cmplx = zeros(Complex128, (2))
      dfdu_cmplx = zeros(2, sbpface.numnodes, 4)
      for d = 0:2*p
        for j = 0:d
          i = d-j
          # the function to be integrated is (x+1)^i (y+1)^j; 0.5 factor
          # accounts for the transformation to ref space
          # Note: output vector only depend on the corresponding input indices
          for f = 1:4
            for n = 1:sbpface.numnodes
              uface_cmplx[1,:,:] = 0.5*((xf[1,:,:]+1).^i).*((xf[2,:,:]+1).^j)
              uface_cmplx[2,:,:] = 0.5 # integrate constant, gives perimeter
              uface_cmplx[1,n,f] += complex(0.0, eps_cmplx)
              uface_cmplx[2,n,f] += complex(0.0, eps_cmplx)
              fill!(fun_cmplx, complex(0.0,0.0))
              integratefunctional!(sbpface, bndryfaces, uface_cmplx, fun_cmplx)
              dfdu_cmplx[1,n,f] = imag(fun_cmplx[1])/eps_cmplx
              dfdu_cmplx[2,n,f] = imag(fun_cmplx[2])/eps_cmplx
            end
          end
          fun_bar = [1.0; 0.0]
          fill!(uface_bar,0.0)
          integratefunctional_rev!(sbpface, bndryfaces, uface_bar, fun_bar)
          @fact uface_bar[1,:,:] --> roughly(dfdu_cmplx[1,:,:], atol=1e-14)
          @fact sum(uface_bar[2,:,:]) --> roughly(0.0, atol=1e-14)
          fun_bar = [0.0; 1.0]
          fill!(uface_bar,0.0)
          integratefunctional_rev!(sbpface, bndryfaces, uface_bar, fun_bar)
          @fact sum(uface_bar[1,:,:]) --> roughly(0.0, atol=1e-14)
          @fact uface_bar[2,:,:] --> roughly(dfdu_cmplx[2,:,:], atol=1e-14)          
        end
      end
    end
  end

  context("Testing integratefunctional_rev! (TetSBP, vector field method)") do
    # build a four element grid and verify the accuracy of boundary integration
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p, internal=false)
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

      eps_cmplx = 1e-60
      uface_bar = zeros(2, sbpface.numnodes, 12)
      uface_cmplx = zeros(Complex128, (2, sbpface.numnodes, 12))
      fun_cmplx = zeros(Complex128, (2))
      dfdu_cmplx = zeros(2, sbpface.numnodes, 12)
      for d = 0:2*p
        for k = 0:d
          for j = 0:d-k
            i = d-k-j
            # the function to be integrated is (x+1)^i (y+1)^j (z+1)^k
            # Note: output vector only depend on the corresponding input indices
            for f = 1:12
              for n = 1:sbpface.numnodes            
                uface_cmplx[1,:,:] = (((xf[1,:,:]+1).^i).*((xf[2,:,:]+1).^j).*
                                      ((xf[3,:,:]+1).^k))
                uface_cmplx[2,:,:] = 1.0 # integrate constant, gives surface area
                scale!(uface_cmplx, 0.25) # 0.25 factor accounts for tranformation
                uface_cmplx[1,n,f] += complex(0.0, eps_cmplx)
                uface_cmplx[2,n,f] += complex(0.0, eps_cmplx)
                fill!(fun_cmplx, 0.0)
                integratefunctional!(sbpface, bndryfaces, uface_cmplx, fun_cmplx)
                dfdu_cmplx[1,n,f] = imag(fun_cmplx[1])/eps_cmplx
                dfdu_cmplx[2,n,f] = imag(fun_cmplx[2])/eps_cmplx
              end
            end
            fun_bar = [1.0; 0.0]
            fill!(uface_bar,0.0)
            integratefunctional_rev!(sbpface, bndryfaces, uface_bar, fun_bar)
            @fact uface_bar[1,:,:] --> roughly(dfdu_cmplx[1,:,:], atol=1e-14)
            @fact sum(uface_bar[2,:,:]) --> roughly(0.0, atol=1e-14)
            fun_bar = [0.0; 1.0]
            fill!(uface_bar,0.0)
            integratefunctional_rev!(sbpface, bndryfaces, uface_bar, fun_bar)
            @fact sum(uface_bar[1,:,:]) --> roughly(0.0, atol=1e-14)
            @fact uface_bar[2,:,:] --> roughly(dfdu_cmplx[2,:,:], atol=1e-14)
          end
        end
      end
    end
  end

  for TSBP = (TriSBP, SparseTriSBP, TetSBP)
    @eval begin
      context("Testing boundaryFunctional_rev! ("string($TSBP)" vector field method)") do
        for p = 1:4
          sbp = TriSBP{Float64}(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          uface = rand(Float64, (4,sbpface.numnodes))
          vface = zeros(uface)
          vfun = rand(Float64, (4))
          ufun = zeros(vfun)
          for face = 1:size(sbp.Q,3)+1
            fill!(ufun, 0.0)
            integrateBoundaryFunctional!(sbpface, face, uface, ufun)
            vtBu = sum(vfun.*ufun)
            fill!(vface, 0.0)
            integrateBoundaryFunctional_rev!(sbpface, face, vface, vfun)
            utBv = sum(uface.*vface)
            @fact vtBu --> roughly(utBv, atol=1e-15)
          end
        end
      end
    end
  end

  for TSBP = (TriSBP, SparseTriSBP, TetSBP)
    @eval begin
      context("Testing boundaryintegrate_rev! ("string($TSBP)" scalar field method)") do
        for p = 1:4
          sbp = TriSBP{Float64}(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          uface = rand(Float64, (sbpface.numnodes,4))
          vface = zeros(uface)
          vvol = rand(Float64, (sbp.numnodes,2))
          uvol = zeros(vvol)
          bndryfaces = Array(Boundary, 4)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(1,3)
          bndryfaces[3] = Boundary(2,1)
          bndryfaces[4] = Boundary(2,2)
          boundaryintegrate!(sbpface, bndryfaces, uface, uvol)
          vtRtBu = sum(vvol.*uvol)
          boundaryintegrate_rev!(sbpface, bndryfaces, vface, vvol)
          utBRv = sum(uface.*vface)
          @fact vtRtBu --> roughly(utBRv, atol=1e-15)
        end
      end
    end
  end

  for TSBP = (TriSBP, SparseTriSBP, TetSBP)
    @eval begin
      context("Testing boundaryintegrate_rev! ("string($TSBP)" vector field method)") do
        for p = 1:4
          sbp = TriSBP{Float64}(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          uface = rand(Float64, (4,sbpface.numnodes,4))
          vface = zeros(uface)
          vvol = rand(Float64, (4,sbp.numnodes,2))
          uvol = zeros(vvol)
          bndryfaces = Array(Boundary, 4)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(1,3)
          bndryfaces[3] = Boundary(2,1)
          bndryfaces[4] = Boundary(2,2)
          boundaryintegrate!(sbpface, bndryfaces, uface, uvol)
          vtRtBu = sum(vvol.*uvol)
          boundaryintegrate_rev!(sbpface, bndryfaces, vface, vvol)
          utBRv = sum(uface.*vface)
          @fact vtRtBu --> roughly(utBRv, atol=1e-15)
        end
      end
    end
  end

  for TSBP = (TriSBP, SparseTriSBP, TetSBP)
    @eval begin
      context("Testing boundaryFaceIntegrate_rev! ("string($TSBP)" scalar field method)") do
        for p = 1:4
          sbp = TriSBP{Float64}(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          uface = rand(Float64, (sbpface.numnodes))
          vface = zeros(uface)
          vvol = rand(Float64, (sbp.numnodes))
          uvol = zeros(vvol)
          for face = 1:size(sbp.Q,3)+1
            fill!(uvol, 0.0)
            boundaryFaceIntegrate!(sbpface, face, uface, uvol)
            vtRtBu = sum(vvol.*uvol)
            fill!(vface, 0.0)
            boundaryFaceIntegrate_rev!(sbpface, face, vface, vvol)
            utBRv = sum(uface.*vface)
            @fact vtRtBu --> roughly(utBRv, atol=1e-15)
          end
        end
      end
    end
  end

  for TSBP = (TriSBP, SparseTriSBP, TetSBP)
    @eval begin
      context("Testing boundaryFaceIntegrate_rev! ("string($TSBP)" vector field method)") do
        for p = 1:4
          sbp = TriSBP{Float64}(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          uface = rand(Float64, (4,sbpface.numnodes))
          vface = zeros(uface)
          vvol = rand(Float64, (4,sbp.numnodes))
          uvol = zeros(vvol)
          for face = 1:size(sbp.Q,3)+1
            fill!(uvol, 0.0)
            boundaryFaceIntegrate!(sbpface, face, uface, uvol)
            vtRtBu = sum(vvol.*uvol)
            fill!(vface, 0.0)
            boundaryFaceIntegrate_rev!(sbpface, face, vface, vvol)
            utBRv = sum(uface.*vface)
            @fact vtRtBu --> roughly(utBRv, atol=1e-15)
          end
        end
      end
    end
  end

  for TSBP = (TriSBP, SparseTriSBP, TetSBP)
    @eval begin
      context("Testing interiorfaceintegrate_rev! ("string($TSBP)" scalar field method)") do
        for p = 1:4
          sbp = TriSBP{Float64}(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          ifaces = Array(Interface, 1)
          ifaces[1] = Interface(1,2,2,3,1)
          uface = rand(sbpface.numnodes, 1)
          vface = zeros(uface)
          vvol = rand(sbp.numnodes,2)
          uvol = zeros(vvol)
          interiorfaceintegrate!(sbpface, ifaces, uface, uvol)
          vtRtBu = sum(vvol.*uvol)
          interiorfaceintegrate_rev!(sbpface, ifaces, vface, vvol)
          utBRv = sum(uface.*vface)
          @fact vtRtBu --> roughly(utBRv, atol=1e-15)
        end
      end
    end
  end

  for TSBP = (TriSBP, SparseTriSBP, TetSBP)
    @eval begin
      context("Testing interiorfaceintegrate_rev! ("string($TSBP)" vector field method)") do
        for p = 1:4
          sbp = TriSBP{Float64}(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          ifaces = Array(Interface, 1)
          ifaces[1] = Interface(1,2,2,3,1)
          uface = rand(4,sbpface.numnodes, 1)
          vface = zeros(uface)
          vvol = rand(4,sbp.numnodes,2)
          uvol = zeros(vvol)
          interiorfaceintegrate!(sbpface, ifaces, uface, uvol)
          vtRtBu = sum(vvol.*uvol)
          interiorfaceintegrate_rev!(sbpface, ifaces, vface, vvol)
          utBRv = sum(uface.*vface)
          @fact vtRtBu --> roughly(utBRv, atol=1e-15)
        end
      end
    end
  end

  for TSBP = (TriSBP, SparseTriSBP, TetSBP)
    @eval begin
      context("Testing interiorFaceIntegrate_rev! ("string($TSBP)" scalar field method)") do
        for p = 1:4
          sbp = TriSBP{Float64}(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          uface = rand(Float64, (sbpface.numnodes))
          vface = zeros(uface)
          vL = rand(Float64, (sbp.numnodes))
          vR = rand(Float64, (sbp.numnodes))
          uL = zeros(vL)
          uR = zeros(vR)
          face = Interface(1,2,2,3,1)
          interiorFaceIntegrate!(sbpface, face, uface, uL, uR)
          vtRtBu = sum(vL.*uL) + sum(vR.*uR)
          interiorFaceIntegrate_rev!(sbpface, face, vface, vL, vR)
          utBRv = sum(uface.*vface)
          @fact vtRtBu --> roughly(utBRv, atol=1e-15)
        end
      end
    end
  end

  for TSBP = (TriSBP, SparseTriSBP, TetSBP)
    @eval begin
      context("Testing interiorFaceIntegrate_rev! ("string($TSBP)" vector field method)") do
        for p = 1:4
          sbp = TriSBP{Float64}(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          uface = rand(Float64, (4,sbpface.numnodes))
          vface = zeros(uface)
          vL = rand(Float64, (4,sbp.numnodes))
          vR = rand(Float64, (4,sbp.numnodes))
          uL = zeros(vL)
          uR = zeros(vR)
          face = Interface(1,2,2,3,1)
          interiorFaceIntegrate!(sbpface, face, uface, uL, uR)
          vtRtBu = sum(vL.*uL) + sum(vR.*uR)
          interiorFaceIntegrate_rev!(sbpface, face, vface, vL, vR)
          utBRv = sum(uface.*vface)
          @fact vtRtBu --> roughly(utBRv, atol=1e-15)
        end
      end
    end
  end

end
