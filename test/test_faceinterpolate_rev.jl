facts("Testing SummationByParts Module (reverse-diff of face-data interpolation methods)...") do

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing boundaryinterpolate_rev! ("string($TSBP)" scalar field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          uvol = rand(Float64, (sbp.numnodes,2))
          vvol = zeros(uvol)
          vface = rand(Float64, (sbpface.numnodes,4))
          uface = zeros(vface)
          bndryfaces = Array(Boundary, 4)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(1,3)
          bndryfaces[3] = Boundary(2,1)
          bndryfaces[4] = Boundary(2,2)
          boundaryinterpolate!(sbpface, bndryfaces, uvol, uface)
          vtRu = sum(vface.*uface)
          boundaryinterpolate_rev!(sbpface, bndryfaces, vvol, vface)
          utRtv = sum(vvol.*uvol)
          @fact vtRu --> roughly(utRtv, atol=1e-15)
        end
      end
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing boundaryinterpolate_rev! ("string($TSBP)" vector field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          uvol = rand(Float64, (4,sbp.numnodes,2))
          vvol = zeros(uvol)
          vface = rand(Float64, (4,sbpface.numnodes,4))
          uface = zeros(vface)
          bndryfaces = Array(Boundary, 4)
          bndryfaces[1] = Boundary(1,1)
          bndryfaces[2] = Boundary(1,3)
          bndryfaces[3] = Boundary(2,1)
          bndryfaces[4] = Boundary(2,2)
          boundaryinterpolate!(sbpface, bndryfaces, uvol, uface)
          vtRu = sum(vface.*uface)
          boundaryinterpolate_rev!(sbpface, bndryfaces, vvol, vface)
          utRtv = sum(vvol.*uvol)
          @fact vtRu --> roughly(utRtv, atol=1e-15)
        end
      end
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing boundaryFaceInterpolate_rev! ("string($TSBP)" scalar field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          uvol = rand(Float64, (sbp.numnodes))
          vvol = zeros(uvol)
          vface = rand(Float64, (sbpface.numnodes))
          uface = zeros(vface)
          for face = 1:size(sbp.Q,3)+1
            fill!(uface, 0.0)
            boundaryFaceInterpolate!(sbpface, face, uvol, uface)
            vtRu = sum(vface.*uface)
            fill!(vvol, 0.0)
            boundaryFaceInterpolate_rev!(sbpface, face, vvol, vface)
            utRtv = sum(vvol.*uvol)
            @fact vtRu --> roughly(utRtv, atol=1e-15)
          end
        end
      end
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing boundaryFaceInterpolate_rev! ("string($TSBP)" vector field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          uvol = rand(Float64, (4,sbp.numnodes))
          vvol = zeros(uvol)
          vface = rand(Float64, (4,sbpface.numnodes))
          uface = zeros(vface)
          for face = 1:size(sbp.Q,3)+1
            fill!(uface, 0.0)
            boundaryFaceInterpolate!(sbpface, face, uvol, uface)
            vtRu = sum(vface.*uface)
            fill!(vvol, 0.0)
            boundaryFaceInterpolate_rev!(sbpface, face, vvol, vface)
            utRtv = sum(vvol.*uvol)
            @fact vtRu --> roughly(utRtv, atol=1e-15)
          end
        end
      end
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing interiorfaceinterpolate_rev! ("string($TSBP)" scalar field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          ifaces = Array(Interface, 1)
          ifaces[1] = Interface(1,2,2,3,1)
          uvol = rand(sbp.numnodes,2)
          vvol = zeros(uvol)
          vface = rand(2, sbpface.numnodes, 1)
          uface = zeros(vface)
          interiorfaceinterpolate!(sbpface, ifaces, uvol, uface)
          vtRu = sum(vface.*uface)
          interiorfaceinterpolate_rev!(sbpface, ifaces, vvol, vface)
          utRtv = sum(vvol.*uvol)
          @fact vtRu --> roughly(utRtv, atol=1e-15)          
        end
      end
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing interiorfaceinterpolate_rev! ("string($TSBP)" vector field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          ifaces = Array(Interface, 1)
          ifaces[1] = Interface(1,2,2,3,1)
          uvol = rand(4,sbp.numnodes,2)
          vvol = zeros(uvol)
          vface = rand(4,2, sbpface.numnodes, 1)
          uface = zeros(vface)
          interiorfaceinterpolate!(sbpface, ifaces, uvol, uface)
          vtRu = sum(vface.*uface)
          interiorfaceinterpolate_rev!(sbpface, ifaces, vvol, vface)
          utRtv = sum(vvol.*uvol)
          @fact vtRu --> roughly(utRtv, atol=1e-15)          
        end
      end
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing interiorFaceInterpolate_rev! ("string($TSBP)" scalar field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          uL = rand(Float64, (sbp.numnodes))
          uR = rand(Float64, (sbp.numnodes))
          vL = zeros(uL)
          vR = zeros(uR)
          vfaceL = rand(Float64, (sbpface.numnodes))
          vfaceR = rand(Float64, (sbpface.numnodes))
          ufaceL = zeros(vfaceL)
          ufaceR = zeros(vfaceR)
          face = Interface(1,2,2,3,1)
          interiorFaceInterpolate!(sbpface, face, uL, uR, ufaceL, ufaceR)
          vtRu_left = sum(vfaceL.*ufaceL)
          vtRu_right = sum(vfaceR.*ufaceR)
          interiorFaceInterpolate_rev!(sbpface, face, vL, vR, vfaceL, vfaceR)
          utRtv_left = sum(vL.*uL)
          utRtv_right = sum(vR.*uR)
          @fact vtRu_left --> roughly(utRtv_left, atol=1e-15)
          @fact vtRu_right --> roughly(utRtv_right, atol=1e-15)
        end
      end
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing interiorFaceInterpolate_rev! ("string($TSBP)" vector field method)") do
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
          uL = rand(Float64, (4,sbp.numnodes))
          uR = rand(Float64, (4,sbp.numnodes))
          vL = zeros(uL)
          vR = zeros(uR)
          vfaceL = rand(Float64, (4,sbpface.numnodes))
          vfaceR = rand(Float64, (4,sbpface.numnodes))
          ufaceL = zeros(vfaceL)
          ufaceR = zeros(vfaceR)
          face = Interface(1,2,2,3,1)
          interiorFaceInterpolate!(sbpface, face, uL, uR, ufaceL, ufaceR)
          vtRu_left = sum(vfaceL.*ufaceL)
          vtRu_right = sum(vfaceR.*ufaceR)
          interiorFaceInterpolate_rev!(sbpface, face, vL, vR, vfaceL, vfaceR)
          utRtv_left = sum(vL.*uL)
          utRtv_right = sum(vR.*uR)
          @fact vtRu_left --> roughly(utRtv_left, atol=1e-15)
          @fact vtRu_right --> roughly(utRtv_right, atol=1e-15)
        end
      end
    end
  end

end
