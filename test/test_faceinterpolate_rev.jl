@testset "Testing SummationByParts Module (reverse-diff of face-data interpolation methods)..." begin

  for TSBP = ((getLineSegSBPLobbato,getLineSegFace,4),
              (getLineSegSBPLegendre,getLineSegFace,4),
              (getTriSBPGamma,TriFace{Float64},4),
              (getTriSBPOmega,TriFace{Float64},4),
              (getTriSBPDiagE,getTriFaceForDiagE,4),
              (getTetSBPGamma,TetFace{Float64},4),
              (getTetSBPOmega,TetFace{Float64},4),
              (getTetSBPDiagE,getTetFaceForDiagE,4))
    @eval begin
      @testset "Testing boundaryinterpolate_rev! ($(string($TSBP[1])) scalar field method)" begin
        for p = 1:($TSBP[3])
          # sbp = ($TSBP[1])(degree=p)
          # sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          if $TSBP == (getTetSBPDiagE,getTetFaceForDiagE,4)
            sbp = ($TSBP[1])(degree=p, faceopertype=:Omega)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
          else
            sbp = ($TSBP[1])(degree=p)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          end
          uvol = rand(Float64, (sbp.numnodes,2))
          vvol = zeros(size(uvol))
          vface = rand(Float64, (sbpface.numnodes,4))
          uface = zeros(size(vface))
          bndryfaces = Array{Boundary}(undef,4)
          bndryfaces[1] = Boundary(1,1)
          if size(sbp.Q,3) < 2
            bndryfaces[2] = Boundary(1,2)
          else
            bndryfaces[2] = Boundary(1,3)
          end
          bndryfaces[3] = Boundary(2,1)
          bndryfaces[4] = Boundary(2,2)
          boundaryinterpolate!(sbpface, bndryfaces, uvol, uface)
          vtRu = sum(vface.*uface)
          boundaryinterpolate_rev!(sbpface, bndryfaces, vvol, vface)
          utRtv = sum(vvol.*uvol)
          @test ≈(vtRu, utRtv, atol=5e-14)
        end
      end
    end
  end

  for TSBP = ((getLineSegSBPLobbato,getLineSegFace,4),
              (getLineSegSBPLegendre,getLineSegFace,4),
              (getTriSBPGamma,TriFace{Float64},4),
              (getTriSBPOmega,TriFace{Float64},4),
              (getTriSBPDiagE,getTriFaceForDiagE,4),
              (getTetSBPGamma,TetFace{Float64},4),
              (getTetSBPOmega,TetFace{Float64},4),
              (getTetSBPDiagE,getTetFaceForDiagE,4))
    @eval begin
      @testset "Testing boundaryinterpolate_rev! ($(string($TSBP[1])) vector field method)" begin
        for p = 1:($TSBP[3])
          # sbp = ($TSBP[1])(degree=p)
          # sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          if $TSBP == (getTetSBPDiagE,getTetFaceForDiagE,4)
            sbp = ($TSBP[1])(degree=p, faceopertype=:Omega)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
          else
            sbp = ($TSBP[1])(degree=p)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          end
          uvol = rand(Float64, (4,sbp.numnodes,2))
          vvol = zeros(size(uvol))
          vface = rand(Float64, (4,sbpface.numnodes,4))
          uface = zeros(size(vface))
          bndryfaces = Array{Boundary}(undef,4)
          bndryfaces[1] = Boundary(1,1)
          if size(sbp.Q,3) < 2
            bndryfaces[2] = Boundary(1,2)
          else
            bndryfaces[2] = Boundary(1,3)
          end
          bndryfaces[3] = Boundary(2,1)
          bndryfaces[4] = Boundary(2,2)
          boundaryinterpolate!(sbpface, bndryfaces, uvol, uface)
          vtRu = sum(vface.*uface)
          boundaryinterpolate_rev!(sbpface, bndryfaces, vvol, vface)
          utRtv = sum(vvol.*uvol)
          @test ≈(vtRu, utRtv, atol=1e-12)
        end
      end
    end
  end

  for TSBP = ((getLineSegSBPLobbato,getLineSegFace,4),
              (getLineSegSBPLegendre,getLineSegFace,4),
              (getTriSBPGamma,TriFace{Float64},4),
              (getTriSBPOmega,TriFace{Float64},4),
              (getTriSBPDiagE,getTriFaceForDiagE,4),
              (getTetSBPGamma,TetFace{Float64},4),
              (getTetSBPOmega,TetFace{Float64},4),
              (getTetSBPDiagE,getTetFaceForDiagE,4))
    @eval begin
      @testset "Testing boundaryFaceInterpolate_rev! ($(string($TSBP[1])) scalar field method)" begin
        for p = 1:($TSBP[3])
          # sbp = ($TSBP[1])(degree=p)
          # sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          if $TSBP == (getTetSBPDiagE,getTetFaceForDiagE,4)
            sbp = ($TSBP[1])(degree=p, faceopertype=:Omega)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
          else
            sbp = ($TSBP[1])(degree=p)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          end
          uvol = rand(Float64, (sbp.numnodes))
          vvol = zeros(size(uvol))
          vface = rand(Float64, (sbpface.numnodes))
          uface = zeros(size(vface))
          for face = 1:size(sbp.Q,3)+1
            fill!(uface, 0.0)
            boundaryFaceInterpolate!(sbpface, face, uvol, uface)
            vtRu = sum(vface.*uface)
            fill!(vvol, 0.0)
            boundaryFaceInterpolate_rev!(sbpface, face, vvol, vface)
            utRtv = sum(vvol.*uvol)
            @test ≈(vtRu, utRtv, atol=1e-14)
          end
        end
      end
    end
  end

  for TSBP = ((getLineSegSBPLobbato,getLineSegFace,4),
              (getLineSegSBPLegendre,getLineSegFace,4),
              (getTriSBPGamma,TriFace{Float64},4),
              (getTriSBPOmega,TriFace{Float64},4),
              (getTriSBPDiagE,getTriFaceForDiagE,4),
              (getTetSBPGamma,TetFace{Float64},4),
              (getTetSBPOmega,TetFace{Float64},4),
              (getTetSBPDiagE,getTetFaceForDiagE,4))
    @eval begin
      @testset "Testing boundaryFaceInterpolate_rev! ($(string($TSBP[1])) vector field method)" begin
        for p = 1:($TSBP[3])
          # sbp = ($TSBP[1])(degree=p)
          # sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          if $TSBP == (getTetSBPDiagE,getTetFaceForDiagE,4)
            sbp = ($TSBP[1])(degree=p, faceopertype=:Omega)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
          else
            sbp = ($TSBP[1])(degree=p)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          end
          uvol = rand(Float64, (4,sbp.numnodes))
          vvol = zeros(size(uvol))
          vface = rand(Float64, (4,sbpface.numnodes))
          uface = zeros(size(vface))
          for face = 1:size(sbp.Q,3)+1
            fill!(uface, 0.0)
            boundaryFaceInterpolate!(sbpface, face, uvol, uface)
            vtRu = sum(vface.*uface)
            fill!(vvol, 0.0)
            boundaryFaceInterpolate_rev!(sbpface, face, vvol, vface)
            utRtv = sum(vvol.*uvol)
            @test ≈(vtRu, utRtv, atol=5e-13)
          end
        end
      end
    end
  end

  for TSBP = ((getLineSegSBPLobbato,getLineSegFace,4),
              (getLineSegSBPLegendre,getLineSegFace,4),
              (getTriSBPGamma,TriFace{Float64},4),
              (getTriSBPOmega,TriFace{Float64},4),
              (getTriSBPDiagE,getTriFaceForDiagE,4),
              (getTetSBPGamma,TetFace{Float64},4),
              (getTetSBPOmega,TetFace{Float64},4),
              (getTetSBPDiagE,getTetFaceForDiagE,4))
    @eval begin
      @testset "Testing interiorfaceinterpolate_rev! ($(string($TSBP[1])) scalar field method)" begin
        for p = 1:($TSBP[3])
          # sbp = ($TSBP[1])(degree=p)
          # sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          if $TSBP == (getTetSBPDiagE,getTetFaceForDiagE,4)
            sbp = ($TSBP[1])(degree=p, faceopertype=:Omega)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
          else
            sbp = ($TSBP[1])(degree=p)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          end
          ifaces = Array{Interface}(undef,1)
          if size(sbp.Q,3) < 2
            ifaces[1] = Interface(1,2,2,1,1)
          else
            ifaces[1] = Interface(1,2,2,3,1)            
          end
          uvol = rand(sbp.numnodes,2)
          vvol = zeros(size(uvol))
          vface = rand(2, sbpface.numnodes, 1)
          uface = zeros(size(vface))
          interiorfaceinterpolate!(sbpface, ifaces, uvol, uface)
          vtRu = sum(vface.*uface)
          interiorfaceinterpolate_rev!(sbpface, ifaces, vvol, vface)
          utRtv = sum(vvol.*uvol)
          @test ≈(vtRu, utRtv, atol=1e-14)          
        end
      end
    end
  end

  for TSBP = ((getLineSegSBPLobbato,getLineSegFace,4),
              (getLineSegSBPLegendre,getLineSegFace,4),
              (getTriSBPGamma,TriFace{Float64},4),
              (getTriSBPOmega,TriFace{Float64},4),
              (getTriSBPDiagE,getTriFaceForDiagE,4),
              (getTetSBPGamma,TetFace{Float64},4),
              (getTetSBPOmega,TetFace{Float64},4),
              (getTetSBPDiagE,getTetFaceForDiagE,4))
    @eval begin
      @testset "Testing interiorfaceinterpolate_rev! ($(string($TSBP[1])) vector field method)" begin
        for p = 1:($TSBP[3])
          # sbp = ($TSBP[1])(degree=p)
          # sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          if $TSBP == (getTetSBPDiagE,getTetFaceForDiagE,4)
            sbp = ($TSBP[1])(degree=p, faceopertype=:Omega)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
          else
            sbp = ($TSBP[1])(degree=p)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          end
          ifaces = Array{Interface}(undef,1)
          if size(sbp.Q,3) < 2
            ifaces[1] = Interface(1,2,2,1,1)
          else
            ifaces[1] = Interface(1,2,2,3,1)            
          end
          uvol = rand(4,sbp.numnodes,2)
          vvol = zeros(size(uvol))
          vface = rand(4,2, sbpface.numnodes, 1)
          uface = zeros(size(vface))
          interiorfaceinterpolate!(sbpface, ifaces, uvol, uface)
          vtRu = sum(vface.*uface)
          interiorfaceinterpolate_rev!(sbpface, ifaces, vvol, vface)
          utRtv = sum(vvol.*uvol)
          @test ≈(vtRu, utRtv, atol=1e-13)          
        end
      end
    end
  end

  for TSBP = ((getLineSegSBPLobbato,getLineSegFace,4),
              (getLineSegSBPLegendre,getLineSegFace,4),
              (getTriSBPGamma,TriFace{Float64},4),
              (getTriSBPOmega,TriFace{Float64},4),
              (getTriSBPDiagE,getTriFaceForDiagE,4),
              (getTetSBPGamma,TetFace{Float64},4),
              (getTetSBPOmega,TetFace{Float64},4),
              (getTetSBPDiagE,getTetFaceForDiagE,4))
    @eval begin
      @testset "Testing interiorFaceInterpolate_rev! ($(string($TSBP[1])) scalar field method)" begin
        for p = 1:($TSBP[3])
          # sbp = ($TSBP[1])(degree=p)
          # sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          if $TSBP == (getTetSBPDiagE,getTetFaceForDiagE,4)
            sbp = ($TSBP[1])(degree=p, faceopertype=:Omega)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
          else
            sbp = ($TSBP[1])(degree=p)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          end
          uL = rand(Float64, (sbp.numnodes))
          uR = rand(Float64, (sbp.numnodes))
          vL = zeros(size(uL))
          vR = zeros(size(uR))
          vfaceL = rand(Float64, (sbpface.numnodes))
          vfaceR = rand(Float64, (sbpface.numnodes))
          ufaceL = zeros(size(vfaceL))
          ufaceR = zeros(size(vfaceR))
          if size(sbp.Q,3) < 2
            face = Interface(1,2,2,1,1)
          else
            face = Interface(1,2,2,3,1)            
          end
          interiorFaceInterpolate!(sbpface, face, uL, uR, ufaceL, ufaceR)
          vtRu_left = sum(vfaceL.*ufaceL)
          vtRu_right = sum(vfaceR.*ufaceR)
          interiorFaceInterpolate_rev!(sbpface, face, vL, vR, vfaceL, vfaceR)
          utRtv_left = sum(vL.*uL)
          utRtv_right = sum(vR.*uR)
          @test ≈(vtRu_left, utRtv_left, atol=1e-14)
          @test ≈(vtRu_right, utRtv_right, atol=1e-14)
        end
      end
    end
  end

  for TSBP = ((getLineSegSBPLobbato,getLineSegFace,4),
              (getLineSegSBPLegendre,getLineSegFace,4),
              (getTriSBPGamma,TriFace{Float64},4),
              (getTriSBPOmega,TriFace{Float64},4),
              (getTriSBPDiagE,getTriFaceForDiagE,4),
              (getTetSBPGamma,TetFace{Float64},4),
              (getTetSBPOmega,TetFace{Float64},4),
              (getTetSBPDiagE,getTetFaceForDiagE,4))
    @eval begin
      @testset "Testing interiorFaceInterpolate_rev! ($(string($TSBP[1])) vector field method)" begin
        for p = 1:($TSBP[3])
          # sbp = ($TSBP[1])(degree=p)
          # sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          if $TSBP == (getTetSBPDiagE,getTetFaceForDiagE,4)
            sbp = ($TSBP[1])(degree=p, faceopertype=:Omega)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
          else
            sbp = ($TSBP[1])(degree=p)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          end
          uL = rand(Float64, (4,sbp.numnodes))
          uR = rand(Float64, (4,sbp.numnodes))
          vL = zeros(size(uL))
          vR = zeros(size(uR))
          vfaceL = rand(Float64, (4,sbpface.numnodes))
          vfaceR = rand(Float64, (4,sbpface.numnodes))
          ufaceL = zeros(size(vfaceL))
          ufaceR = zeros(size(vfaceR))
          if size(sbp.Q,3) < 2
            face = Interface(1,2,2,1,1)
          else
            face = Interface(1,2,2,3,1)            
          end
          interiorFaceInterpolate!(sbpface, face, uL, uR, ufaceL, ufaceR)
          vtRu_left = sum(vfaceL.*ufaceL)
          vtRu_right = sum(vfaceR.*ufaceR)
          interiorFaceInterpolate_rev!(sbpface, face, vL, vR, vfaceL, vfaceR)
          utRtv_left = sum(vL.*uL)
          utRtv_right = sum(vR.*uR)
          @test ≈(vtRu_left, utRtv_left, atol=5e-14)
          @test ≈(vtRu_right, utRtv_right, atol=5e-14)
        end
      end
    end
  end

end
