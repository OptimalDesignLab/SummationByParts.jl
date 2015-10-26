facts("Testing SummationByParts Module (buildfaceoperators.jl file)...") do

  context("Testing SummationByParts.buildfacereconstruction (TriSymCub method, faceonly=true)") do
    # this checks that polynomials of total degree d are reconstructed accurately
    for d = 1:4
      cub, vtx = tricubature(2*d-1, Float64, internal=false)
      facecub, tmp = quadrature(2*d, Float64, internal=true)
      R, perm = SummationByParts.buildfacereconstruction(facecub, cub, vtx, d,
                                                         faceonly=true)
      xy = SymCubatures.calcnodes(cub, vtx)
      # loop over all monomials
      for r = 0:d
        for j = 0:r
          i = r-j
          u = vec((xy[1,:].^i).*(xy[2,:].^j))
          # loop over each face
          for f = 1:3
            xyface = SymCubatures.calcnodes(facecub, vtx[[f;mod(f,3)+1],:])
            uface = vec((xyface[1,:].^i).*(xyface[2,:].^j))
            @fact R*u[perm[:,f]] --> roughly(uface, atol=1e-15)
          end
        end
      end
    end
  end

  context("Testing SummationByParts.buildfacereconstruction (TriSymCub method, faceonly=false)") do
    # this checks that polynomials of total degree d are reconstructed accurately
    for d = 1:4
      cub, vtx = tricubature(2*d-1, Float64, internal=true)
      facecub, tmp = quadrature(2*d, Float64, internal=true)
      R, perm = SummationByParts.buildfacereconstruction(facecub, cub, vtx, d,
                                                         faceonly=false)
      xy = SymCubatures.calcnodes(cub, vtx)
      # loop over all monomials
      for r = 0:d
        for j = 0:r
          i = r-j
          u = vec((xy[1,:].^i).*(xy[2,:].^j))
          # loop over each face
          for f = 1:3
            xyface = SymCubatures.calcnodes(facecub, vtx[[f;mod(f,3)+1],:])
            uface = vec((xyface[1,:].^i).*(xyface[2,:].^j))
            @fact R*u[perm[:,f]] --> roughly(uface, atol=1e-15)
          end
        end
      end
    end
  end

  context("Testing SummationByParts.TriFace constructor (faceonly=true)") do
    for d = 1:4
      cub, vtx = tricubature(2*d-1, Float64, internal=false)
      xy = SymCubatures.calcnodes(cub, vtx)
      face = TriFace{Float64}(degree=d, faceonly=true)
      # loop over monomials of degree <= d
      for r = 0:d
        for j = 0:r
          i = r-j
          u = vec((xy[1,:].^i).*(xy[2,:].^j))
          for q = 0:d
            for l = 0:q
              k = q-l
              v = vec((xy[1,:].^k).*(xy[2,:].^l))
              # compute the boundary integral of u*v*(nx+ny)
              bndryintegral = 0.0
              for f = 1:3
                bndryintegral += sum(face.normal[:,f])* 
                dot(face.interp.'*v[face.perm[:,f]], 
                    diagm(face.wface)*face.interp.'*u[face.perm[:,f]])
              end
              @fact bndryintegral -->
              roughly(((-1)^(j+l+1))*(1^(i+k+1) - (-1)^(i+k+1))/(i+k+1) +
                      2.*((-1)^(j+l))*(1^(r+q+1) - (-1)^(r+q+1))/(r+q+1) +
                      ((-1)^(i+k+1))*(1^(j+l+1) - (-1)^(j+l+1))/(j+l+1),
                      atol=1e-15)
            end
          end
        end
      end
    end
  end

  context("Testing SummationByParts.TriFace constructor (faceonly=false)") do
    for d = 1:4
      cub, vtx = tricubature(2*d-1, Float64, internal=true)
      xy = SymCubatures.calcnodes(cub, vtx)
      face = TriFace{Float64}(degree=d, faceonly=false)
      # loop over monomials of degree <= d
      for r = 0:d
        for j = 0:r
          i = r-j
          u = vec((xy[1,:].^i).*(xy[2,:].^j))
          for q = 0:d
            for l = 0:q
              k = q-l
              v = vec((xy[1,:].^k).*(xy[2,:].^l))
              # compute the boundary integral of u*v*(nx+ny)
              bndryintegral = 0.0
              for f = 1:3
                bndryintegral += sum(face.normal[:,f])* 
                dot(face.interp.'*v[face.perm[:,f]], 
                    diagm(face.wface)*face.interp.'*u[face.perm[:,f]])
              end
              @fact bndryintegral -->
              roughly(((-1)^(j+l+1))*(1^(i+k+1) - (-1)^(i+k+1))/(i+k+1) +
                      2.*((-1)^(j+l))*(1^(r+q+1) - (-1)^(r+q+1))/(r+q+1) +
                      ((-1)^(i+k+1))*(1^(j+l+1) - (-1)^(j+l+1))/(j+l+1),
                      atol=1e-15)
            end
          end
        end
      end
    end
  end

end