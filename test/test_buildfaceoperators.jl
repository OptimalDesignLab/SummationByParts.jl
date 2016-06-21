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

  context("Testing SummationByParts.buildfacereconstruction (TriSymCub method, faceonly=false, internal=false)") do
    # this checks that polynomials of total degree d are reconstructed accurately
    for d = 1:4
      cub, vtx = tricubature(2*d-1, Float64, internal=false)
      facecub, tmp = quadrature(2*d, Float64, internal=false)
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

  context("Testing SummationByParts.buildfacereconstruction (TetSymCub method, faceonly=true)") do
    # this checks that polynomials of total degree d are reconstructed accurately
    for d = 1:4
      cub, vtx = tetcubature(2*d-1, Float64, internal=false)
      facecub, tmp = tricubature(2*d, Float64, internal=true)
      R, perm = SummationByParts.buildfacereconstruction(facecub, cub, vtx, d,
                                                         faceonly=true)
      vtxface = [1 2 3; 1 4 2; 2 4 3; 1 3 4].'
      xyz = SymCubatures.calcnodes(cub, vtx)
      # loop over all monomials
      for r = 0:d
        for k = 0:r
          for j = 0:r-k
            i = r-j-k
            u = vec((xyz[1,:].^i).*(xyz[2,:].^j).*(xyz[3,:].^k))
            # loop over each face
            for f = 1:4
              xyzface = SymCubatures.calcnodes(facecub, vtx[vtxface[:,f],:])
              uface = vec((xyzface[1,:].^i).*(xyzface[2,:].^j).*(xyzface[3,:].^k))
              @fact R*u[perm[:,f]] --> roughly(uface, atol=1e-15)
            end
          end
        end
      end
    end
  end

  context("Testing SummationByParts.buildfacereconstruction (TetSymCub method, faceonly=false, internal=true)") do
    # this checks that polynomials of total degree d are reconstructed accurately
    for d = 1:2
      cub, vtx = tetcubature(2*d-1, Float64, internal=true)
      facecub, tmp = tricubature(2*d, Float64, internal=true)
      R, perm = SummationByParts.buildfacereconstruction(facecub, cub, vtx, d,
                                                         faceonly=false)
      vtxface = [1 2 3; 1 4 2; 2 4 3; 1 3 4].'
      xyz = SymCubatures.calcnodes(cub, vtx)
      # loop over all monomials
      for r = 0:d
        for k = 0:r
          for j = 0:r-k
            i = r-j-k
            u = vec((xyz[1,:].^i).*(xyz[2,:].^j).*(xyz[3,:].^k))
            # loop over each face
            for f = 1:4
              xyzface = SymCubatures.calcnodes(facecub, vtx[vtxface[:,f],:])
              uface = vec((xyzface[1,:].^i).*(xyzface[2,:].^j).*(xyzface[3,:].^k))
              @fact R*u[perm[:,f]] --> roughly(uface, atol=1e-15)
            end
          end
        end
      end
    end
  end

  context("Testing SummationByParts.buildfacederivative (TriSymCub method)") do
    # this checks that polynomials of total degree d are differentiated
    for d = 1:4
      cub, vtx = tricubature(2*d-1, Float64, internal=false)
      facecub, tmp = quadrature(2*d, Float64, internal=true)
      D, perm = SummationByParts.buildfacederivatives(facecub, cub, vtx, d)
      xy = SymCubatures.calcnodes(cub, vtx)
      # loop over all monomials
      for r = 0:d
        for j = 0:r
          i = r-j
          u = vec((xy[1,:].^i).*(xy[2,:].^j))
          # consider face 1
          xyface = SymCubatures.calcnodes(facecub, vtx[[1;2],:])
          dudt = vec(i.*(xyface[1,:].^max(i-1,0)).*xyface[2,:].^j)
          @fact D[:,:,1].'*u[perm[:,1]] --> roughly(dudt, atol=1e-13)
          dudn = vec(j.*(xyface[1,:].^i).*(xyface[2,:].^max(j-1,0)))
          @fact D[:,:,2].'*u[perm[:,1]] --> roughly(dudn, atol=1e-13)
          # consider face 2
          xyface = SymCubatures.calcnodes(facecub, vtx[[2;3],:])
          dudt = vec(j.*(xyface[1,:].^i).*(xyface[2,:].^max(j-1,0))) -
          vec(i.*(xyface[1,:].^max(i-1,0)).*xyface[2,:].^j)
          @fact D[:,:,1].'*u[perm[:,2]] --> roughly(dudt, atol=1e-13)
          dudn = -vec(i.*(xyface[1,:].^max(i-1,0)).*xyface[2,:].^j)
          @fact D[:,:,2].'*u[perm[:,2]] --> roughly(dudn, atol=1e-13)
          # consider face 3
          xyface = SymCubatures.calcnodes(facecub, vtx[[3;1],:])
          dudt = -vec(j.*(xyface[1,:].^i).*(xyface[2,:].^max(j-1,0)))
          @fact D[:,:,1].'*u[perm[:,3]] --> roughly(dudt, atol=1e-13)
          dudn = vec(i.*(xyface[1,:].^max(i-1,0)).*xyface[2,:].^j) -
          vec(j.*(xyface[1,:].^i).*(xyface[2,:].^max(j-1,0)))
          @fact D[:,:,2].'*u[perm[:,3]] --> roughly(dudn, atol=1e-13)
        end
      end
    end
  end

  context("Testing SummationByParts.TriFace constructor (internal=false)") do
    for d = 1:4
      cub, vtx = tricubature(2*d-1, Float64, internal=false)
      xy = SymCubatures.calcnodes(cub, vtx)
      face = TriFace{Float64}(d, cub, vtx)
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

  context("Testing SummationByParts.TriFace constructor (internal=false, vertices=true)") do
    for d = 1:4
      cub, vtx = tricubature(2*d-1, Float64, internal=false, vertices=true)
      xy = SymCubatures.calcnodes(cub, vtx)
      face = TriFace{Float64}(d, cub, vtx)
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

  context("Testing SummationByParts.TriFace constructor (internal=true)") do
    for d = 1:4
      cub, vtx = tricubature(2*d-1, Float64, internal=true)
      xy = SymCubatures.calcnodes(cub, vtx)
      face = TriFace{Float64}(d, cub, vtx)
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

  context("Testing SummationByParts.TetFace constructor (internal=false)") do
    function integral(a, b, i, j)
      if a < 0 || b < 0
        return 0.0
      elseif i == j
        return ((1 + (-1)^(a + b))/(1 + a + b) + (2*(-1 + (-1)^(a + b)))/(2 +
        a + b) +(1 + (-1)^(a + b))/(3 + a + b))/2.
      else
        return (2*((-1)^a + (-1)^b + 6*(-1)^(a + b)) + (-1)^b*(3 + 13*(-1)^a)*b
                + (-1)^b*(1 + 3*(-1)^a)*b^2 + (-1)^a*a^2*(1 + 3*(-1)^b + 2*(-1)^b*b) +
                (-1)^a*a*(3 + 13*(-1)^b + 12*(-1)^b*b + 2*(-1)^b*b^2))/ 
        ((1 + a)*(2 + a)*(1 + b)*(2 + b)*(3 + a + b))
      end
    end
    for d = 1:4
      cub, vtx = tetcubature(2*d-1, Float64, internal=false)
      xyz = SymCubatures.calcnodes(cub, vtx)
      face = TetFace{Float64}(d, cub, vtx)
      # i and j are coordinate indices, and a and b are powers that determine
      # the polynomial degree
      for i = 1:3
        for a = 0:d
          u = vec(xyz[i,:].^a)
          for j = 1:3
            for b = 0:d
              v = vec(xyz[j,:].^b)
              bndryintegral = 0.0
              for f = 1:4
                bndryintegral += sum(face.normal[:,f])*
                dot(face.interp.'*v[face.perm[:,f]],
                    diagm(face.wface)*face.interp.'*u[face.perm[:,f]])
              end
              @fact bndryintegral --> 
              roughly( a*integral(a-1,b,i,j) + b*integral(a,b-1,i,j) , atol=1e-15)
            end
          end
        end
      end
    end
  end

end