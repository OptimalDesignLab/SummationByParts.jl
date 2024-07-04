@testset "Testing SummationByParts Module (buildfaceoperators.jl file)..." begin
  
  @testset "Testing SummationByParts.buildfacereconstruction (LineSymCub method, faceonly=true)" begin
    # this checks that polynomials of total degree d are reconstructed accurately
    for d = 1:4
      cub, vtx = Cubature.quadrature(2*d-1, Float64, internal=false)
      facecub, tmp = Cubature.pointCubature()
      R, perm = SummationByParts.buildfacereconstruction(facecub, cub, vtx, d)
      x = SymCubatures.calcnodes(cub, vtx)
      # loop over all monomials
      for i = 0:d
        u = vec(x[1,:].^i)
        # loop over each face
        for f = 1:2
          xface = SymCubatures.calcnodes(facecub, reshape(vtx[f,:],(1,1)))
          uface = vec(xface[1,:].^i)
          @test ≈(R*u[perm[:,f]], uface, atol=1e-15)
        end
      end
    end
  end

  @testset "Testing SummationByParts.buildfacereconstruction (LineSymCub method, faceonly=false)" begin
    # this checks that polynomials of total degree d are reconstructed accurately
    for d = 1:4
      cub, vtx = Cubature.quadrature(2*d, Float64, internal=false)
      facecub, tmp = Cubature.pointCubature()
      R, perm = SummationByParts.buildfacereconstruction(facecub, cub, vtx, d)
      x = SymCubatures.calcnodes(cub, vtx)
      # loop over all monomials
      for i = 0:d
        u = vec(x[1,:].^i)
        # loop over each face
        for f = 1:2
          xface = SymCubatures.calcnodes(facecub, reshape(vtx[f,:],(1,1)))
          uface = vec(xface[1,:].^i)
          @test ≈(R*u[perm[:,f]], uface, atol=1e-15)
        end
      end
    end
  end
  
  @testset "Testing SummationByParts.buildfacereconstruction (TriSymCub method, faceonly=true)" begin
    # this checks that polynomials of total degree d are reconstructed accurately
    for d = 1:4
      cub, vtx = Cubature.getTriCubatureGamma(2*d-1, Float64)
      facecub, tmp = Cubature.quadrature(2*d, Float64, internal=true)
      R, perm = SummationByParts.buildfacereconstruction(facecub, cub, vtx, d)
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
            @test ≈(R*u[perm[:,f]], uface, atol=1e-14)
          end
        end
      end
    end
  end

  @testset "Testing SummationByParts.buildfacereconstruction (TriSymCub method, faceonly=false, internal=false)" begin
    # this checks that polynomials of total degree d are reconstructed accurately
    for d = 1:4
      cub, vtx = Cubature.getTriCubatureGamma(2*d-1, Float64)
      facecub, tmp = Cubature.quadrature(2*d, Float64, internal=false)
      R, perm = SummationByParts.buildfacereconstruction(facecub, cub, vtx, d)
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
            @test ≈(R*u[perm[:,f]], uface, atol=1e-14)
          end
        end
      end
    end
  end

  @testset "Testing SummationByParts.buildfacereconstruction (TriSymCub method, faceonly=false)" begin
    # this checks that polynomials of total degree d are reconstructed accurately
    for d = 1:4
      cub, vtx = Cubature.getTriCubatureOmega(2*d, Float64)
      facecub, tmp = Cubature.quadrature(2*d, Float64, internal=true)
      R, perm = SummationByParts.buildfacereconstruction(facecub, cub, vtx, d)
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
            @test ≈(R*u[perm[:,f]], uface, atol=1e-14)
          end
        end
      end
    end
  end

  @testset "Testing SummationByParts.buildfacereconstruction (TetSymCub method, faceonly=true)" begin
    # this checks that polynomials of total degree d are reconstructed accurately
    for d = 1:4
      cub, vtx = Cubature.getTetCubatureGamma(2*d-1, Float64)
      facecub, tmp = Cubature.getTriCubatureOmega(2*d, Float64)
      R, perm = SummationByParts.buildfacereconstruction(facecub, cub, vtx, d)
      vtxface = [1 2 3; 1 4 2; 2 4 3; 1 3 4]'
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
              @test ≈(R*u[perm[:,f]], uface, atol=1e-14)
            end
          end
        end
      end
    end
  end

  @testset "Testing SummationByParts.buildfacereconstruction (TetSymCub method, faceonly=false, internal=true)" begin
    # this checks that polynomials of total degree d are reconstructed accurately
    for d = 1:2
      cub, vtx = Cubature.getTetCubatureOmega(2*d-1, Float64)
      facecub, tmp = Cubature.getTriCubatureOmega(2*d, Float64)
      R, perm = SummationByParts.buildfacereconstruction(facecub, cub, vtx, d)
      vtxface = [1 2 3; 1 4 2; 2 4 3; 1 3 4]'
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
              @test ≈(R*u[perm[:,f]], uface, atol=1e-14)
            end
          end
        end
      end
    end
  end

  @testset "Testing SummationByParts.buildfacederivative (LineSymCub method)" begin
    # this checks that polynomials of total degree d are differentiated
    for d = 1:4
      cub, vtx = Cubature.quadrature(2*d-1, Float64, internal=false)
      facecub, tmp = Cubature.pointCubature()
      D, perm = SummationByParts.buildfacederivatives(facecub, cub, vtx, d)
      x = SymCubatures.calcnodes(cub, vtx)
      # loop over all monomials
      for i = 0:d
        u = vec(x[1,:].^i)
        # consider face 1
        xface = SymCubatures.calcnodes(facecub, reshape(vtx[1,:],(1,1)))
        dudn = vec(i.*xface[1,:].^max(i-1,0))
        @test ≈(D[:,:,1]'*u[perm[:,1]], dudn, atol=1e-13)
        # consider face 2
        xface = SymCubatures.calcnodes(facecub, reshape(vtx[2,:],(1,1)))
        dudn = -vec(i.*xface[1,:].^max(i-1,0))
        @test ≈(D[:,:,1]'*u[perm[:,2]], dudn, atol=1e-13)
      end
    end
  end

  @testset "Testing SummationByParts.buildfacederivative (TriSymCub method)" begin
    # this checks that polynomials of total degree d are differentiated
    for d = 1:4
      cub, vtx = Cubature.getTriCubatureGamma(2*d-1, Float64)
      facecub, tmp = Cubature.quadrature(2*d, Float64, internal=true)
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
          @test ≈(D[:,:,1]'*u[perm[:,1]], dudt, atol=1e-13)
          dudn = vec(j.*(xyface[1,:].^i).*(xyface[2,:].^max(j-1,0)))
          @test ≈(D[:,:,2]'*u[perm[:,1]], dudn, atol=1e-13)
          # consider face 2
          xyface = SymCubatures.calcnodes(facecub, vtx[[2;3],:])
          dudt = vec(j.*(xyface[1,:].^i).*(xyface[2,:].^max(j-1,0))) -
          vec(i.*(xyface[1,:].^max(i-1,0)).*xyface[2,:].^j)
          @test ≈(D[:,:,1]'*u[perm[:,2]], dudt, atol=1e-13)
          dudn = -vec(i.*(xyface[1,:].^max(i-1,0)).*xyface[2,:].^j)
          @test ≈(D[:,:,2]'*u[perm[:,2]], dudn, atol=1e-13)
          # consider face 3
          xyface = SymCubatures.calcnodes(facecub, vtx[[3;1],:])
          dudt = -vec(j.*(xyface[1,:].^i).*(xyface[2,:].^max(j-1,0)))
          @test ≈(D[:,:,1]'*u[perm[:,3]], dudt, atol=1e-13)
          dudn = vec(i.*(xyface[1,:].^max(i-1,0)).*xyface[2,:].^j) -
          vec(j.*(xyface[1,:].^i).*(xyface[2,:].^max(j-1,0)))
          @test ≈(D[:,:,2]'*u[perm[:,3]], dudn, atol=1e-13)
        end
      end
    end
  end

  @testset "Testing SummationByParts.getLineSegFace constructor (internal=false)" begin
    for d = 1:4
      cub, vtx = Cubature.quadrature(2*d-1, Float64, internal=false)
      x = SymCubatures.calcnodes(cub, vtx)
      face = getLineSegFace(d, cub, vtx)
      # loop over monomials of degree <= d
      for i = 0:d
        u = vec(x[1,:].^i)
        for j = 0:d
          v = vec(x[1,:].^j)
          # compute the boundary integral of u*v*nx
          bndryintegral = 0.0
          for f = 1:2
            bndryintegral += face.normal[1,f]*
            dot(face.interp'*v[face.perm[:,f]],
                diagm(face.wface)*face.interp'u[face.perm[:,f]])
          end
          @test ≈(bndryintegral, 1.0 - (-1)^(i+j), atol=1e-14)
        end
      end
    end
  end

  @testset "Testing SummationByParts.getLineSegFace constructor (internal=true)" begin
    for d = 1:4
      cub, vtx = Cubature.quadrature(2*d, Float64, internal=true)
      x = SymCubatures.calcnodes(cub, vtx)
      face = getLineSegFace(d, cub, vtx)
      # loop over monomials of degree <= d
      for i = 0:d
        u = vec(x[1,:].^i)
        for j = 0:d
          v = vec(x[1,:].^j)
          # compute the boundary integral of u*v*nx
          bndryintegral = 0.0
          for f = 1:2
            bndryintegral += face.normal[1,f]*
            dot(face.interp'*v[face.perm[:,f]],
                diagm(face.wface)*face.interp'u[face.perm[:,f]])
          end
          @test ≈(bndryintegral, 1.0 - (-1)^(i+j), atol=1e-14)
        end
      end
    end
  end

  @testset "Testing SummationByParts.TriFace constructor (internal=false)" begin
    for d = 1:4
      cub, vtx = Cubature.getTriCubatureGamma(2*d-1, Float64)
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
                dot(face.interp'*v[face.perm[:,f]], 
                    diagm(face.wface)*face.interp'*u[face.perm[:,f]])
              end
              @test ≈(bndryintegral, ((-1)^(j+l+1))*(1^(i+k+1) - (-1)^(i+k+1))/(i+k+1) +
                      2.0*((-1)^(j+l))*(1^(r+q+1) - (-1)^(r+q+1))/(r+q+1) +
                      ((-1)^(i+k+1))*(1^(j+l+1) - (-1)^(j+l+1))/(j+l+1),
                      atol=1e-14)
            end
          end
        end
      end
    end
  end

  @testset "Testing SummationByParts.TriFace constructor (internal=false, vertices=true)" begin
    for d = 1:4
      cub, vtx = Cubature.getTriCubatureGamma(2*d-1, Float64)
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
                dot(face.interp'*v[face.perm[:,f]], 
                    diagm(face.wface)*face.interp'*u[face.perm[:,f]])
              end
              @test ≈(bndryintegral, ((-1)^(j+l+1))*(1^(i+k+1) - (-1)^(i+k+1))/(i+k+1) +
                      2.0*((-1)^(j+l))*(1^(r+q+1) - (-1)^(r+q+1))/(r+q+1) +
                      ((-1)^(i+k+1))*(1^(j+l+1) - (-1)^(j+l+1))/(j+l+1),
                      atol=1e-14)
            end
          end
        end
      end
    end
  end

  @testset "Testing SummationByParts.TriFace constructor (internal=true)" begin
    for d = 1:4
      cub, vtx = Cubature.getTriCubatureOmega(2*d, Float64)
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
                dot(face.interp'*v[face.perm[:,f]], 
                    diagm(face.wface)*face.interp'*u[face.perm[:,f]])
              end
              @test ≈(bndryintegral, ((-1)^(j+l+1))*(1^(i+k+1) - (-1)^(i+k+1))/(i+k+1) +
                      2.0*((-1)^(j+l))*(1^(r+q+1) - (-1)^(r+q+1))/(r+q+1) +
                      ((-1)^(i+k+1))*(1^(j+l+1) - (-1)^(j+l+1))/(j+l+1),
                      atol=1e-14)
            end
          end
        end
      end
    end
  end

  @testset "Testing SummationByParts.getTriFaceForDiagE constructor (vertices=true)" begin
    for d = 1:4
      cub, vtx = Cubature.getTriCubatureDiagE(2*d, Float64, vertices=true)
      xy = SymCubatures.calcnodes(cub, vtx)
      face = getTriFaceForDiagE(d, cub, vtx, vertices=true)
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
                dot(v[face.perm[:,f]], diagm(face.wface)*u[face.perm[:,f]])
              end
              @test ≈(bndryintegral, ((-1)^(j+l+1))*(1^(i+k+1) - (-1)^(i+k+1))/(i+k+1) +
                      2.0*((-1)^(j+l))*(1^(r+q+1) - (-1)^(r+q+1))/(r+q+1) +
                      ((-1)^(i+k+1))*(1^(j+l+1) - (-1)^(j+l+1))/(j+l+1),
                      atol=1e-14)
            end
          end
        end
      end
    end
  end

  @testset "Testing SummationByParts.getTriFaceForDiagE constructor (vertices=false)" begin
    for d = 1:4
      cub, vtx = Cubature.getTriCubatureDiagE(2*d, Float64, vertices=false)
      xy = SymCubatures.calcnodes(cub, vtx)
      face = getTriFaceForDiagE(d, cub, vtx, vertices=false)
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
                dot(v[face.perm[:,f]], diagm(face.wface)*u[face.perm[:,f]])
              end
              @test ≈(bndryintegral, ((-1)^(j+l+1))*(1^(i+k+1) - (-1)^(i+k+1))/(i+k+1) +
                      2.0*((-1)^(j+l))*(1^(r+q+1) - (-1)^(r+q+1))/(r+q+1) +
                      ((-1)^(i+k+1))*(1^(j+l+1) - (-1)^(j+l+1))/(j+l+1),
                      atol=1e-14)
            end
          end
        end
      end
    end
  end

  @testset "Testing SummationByParts.TetFace constructor (internal=false)" begin
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
      cub, vtx = Cubature.getTetCubatureGamma(2*d-1, Float64)
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
                dot(face.interp'*v[face.perm[:,f]],
                    diagm(face.wface)*face.interp'*u[face.perm[:,f]])
              end
              @test ≈(bndryintegral, a*integral(a-1,b,i,j) + b*integral(a,b-1,i,j) , atol=1e-14)
            end
          end
        end
      end
    end
  end

  @testset "Testing SummationByParts.getTetFaceForDiagE constructor" begin
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
      cub, vtx = Cubature.getTetCubatureDiagE(2*d, Float64, faceopertype=:Omega)
      xyz = SymCubatures.calcnodes(cub, vtx)
      face = getTetFaceForDiagE(d, cub, vtx, faceopertype=:Omega)
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
                dot(v[face.perm[:,f]], diagm(face.wface)*u[face.perm[:,f]])
              end
              @test ≈(bndryintegral, a*integral(a-1,b,i,j) + b*integral(a,b-1,i,j) , atol=1e-14)
            end
          end
        end
      end
    end
  end

  @testset "Testing SummationByParts.getfaceextrapolation function" begin
    opertypes = [:DiagE,:Omega,:Gamma]
    for ioptype=2:length(opertypes)-1
      for id = 2:3
        for ip =1:4
          q = 2*ip-1
          opertype=opertypes[ioptype]
          faceopertype=:Omega
          vertices=true
          if faceopertype==:Omega || faceopertype==:Gamma 
            vertices=false 
          end
          T = Float64
          if opertype==:Omega 
            if id==2
              volcub, volvtx = getTriCubatureOmega(q, T)
            elseif id==3 
              volcub, volvtx = getTetCubatureOmega(q, T)
            end
          elseif opertype==:Gamma 
            if id==2
              volcub, volvtx = getTriCubatureGamma(q, T)
            elseif id==3 
              volcub, volvtx = getTetCubatureGamma(q, T)
            end
          elseif opertype==:DiagE 
            if id==2
              volcub, volvtx = getTriCubatureDiagE(q, T, vertices=vertices)
            elseif id==3
              volcub, volvtx = getTetCubatureDiagE(q, T, faceopertype=faceopertype)
            end
          end
          if id==2
            facecub, facevtx = SummationByParts.Cubature.quadrature(2*ip, T, internal=!vertices)  
          elseif id==3
            if opertype==:DiagE 
              facecub, facevtx = getTriCubatureForTetFaceDiagE(2*ip, T, faceopertype=faceopertype)  
            else 
              facecub, facevtx = getTriCubatureOmega(2*ip, T)
            end
          end

          Rs = SummationByParts.getfaceextrapolation(ip, q, id, opertype=opertype, faceopertype=faceopertype)
          if id==2
            xf1 = SymCubatures.calcnodes(facecub, volvtx[[1;2],:])
            xf2 = SymCubatures.calcnodes(facecub, volvtx[[2;3],:])
            xf3 = SymCubatures.calcnodes(facecub, volvtx[[3;1],:])
            xfs = [xf1, xf2, xf3]
          elseif id==3
            xf1 = SymCubatures.calcnodes(facecub, volvtx[[1;2;3],:])
            xf2 = SymCubatures.calcnodes(facecub, volvtx[[1;4;2],:])
            xf3 = SymCubatures.calcnodes(facecub, volvtx[[2;3;4],:])
            xf4 = SymCubatures.calcnodes(facecub, volvtx[[1;3;4],:])
            xfs = [xf1, xf2, xf3, xf4]
          end

          xv = SymCubatures.calcnodes(volcub, volvtx)
          for ifacet=1:id+1
            @test ≈(Rs[:,:,ifacet]*xv', xfs[ifacet]', atol=1e-13)
          end
        end
      end
    end
  end

end
