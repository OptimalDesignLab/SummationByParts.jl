@testset "Testing SummationByParts Module (buildoperators.jl file)..." begin

  # @testset "Testing SummationByParts.bndrynodalexpansion (TriSymCub method)" begin
  #   # check that P*E produces the identity matrix on the boundary
  #   for d = 1:4
  #     cub, vtx = Cubature.getTriCubatureGamma(2*d-1, Float64)
  #     E = SummationByParts.bndrynodalexpansion(cub, vtx, d)
  #     xy = SymCubatures.calcnodes(cub, vtx)
  #     x = vec(xy[1,:]); y = vec(xy[2,:])
  #     N = convert(Int, (d+1)*(d+2)/2 )
  #     P = zeros(cub.numnodes, N)
  #     ptr = 1
  #     for r = 0:d
  #       for j = 0:r
  #         i = r-j
  #         P[:,ptr] = OrthoPoly.proriolpoly(x, y, i, j)
  #         ptr += 1
  #       end
  #     end
  #     A = P*E
  #     numbndry = SymCubatures.getnumboundarynodes(cub)
  #     bndryindices = SymCubatures.getbndrynodeindices(cub)
  #     @test ≈(A[bndryindices,1:numbndry], I(numbndry), atol=1e-14)
  #   end
  # end
  
  # @testset "Testing SummationByParts.bndrynodalexpansion (TetSymCub method)" begin
  #   # check that P*E produces the identity matrix on the boundary
  #   for d = 1:4
  #     cub, vtx = Cubature.getTetCubatureGamma(2*d-1, Float64)
  #     E = SummationByParts.bndrynodalexpansion(cub, vtx, d)
  #     xyz = SymCubatures.calcnodes(cub, vtx)
  #     x = vec(xyz[1,:]); y = vec(xyz[2,:]); z = vec(xyz[3,:])
  #     N = convert(Int, (d+1)*(d+2)*(d+3)/6 )
  #     P = zeros(cub.numnodes, N)
  #     ptr = 1
  #     for r = 0:d
  #       for k = 0:r
  #         for j = 0:r-k
  #           i = r-j-k
  #           P[:,ptr] = OrthoPoly.proriolpoly(x, y, z, i, j, k)
  #           ptr += 1
  #         end
  #       end
  #     end
  #     A = P*E
  #     numbndry = SymCubatures.getnumboundarynodes(cub)
  #     bndryindices = SymCubatures.getbndrynodeindices(cub)
  #     @test ≈(A[bndryindices,1:numbndry], I(numbndry), atol=1e-14)
  #   end
  # end

  # @testset "Testing SummationByParts.nodalexpansion (TriSymCub method)" begin
  #   # check that P*E produces the identity matrix at the nodes
  #   e = [1;3;4;5]
  #   for d = 1:4
  #     cub, vtx = Cubature.getTriCubatureGamma(2*d-1, Float64)
  #     #@fact_throws SummationByParts.nodalexpansion(cub, vtx, d, e[d])
  #     C = SummationByParts.nodalexpansion(cub, vtx, d, e[d])
  #     xy = SymCubatures.calcnodes(cub, vtx)
  #     x = vec(xy[1,:]); y = vec(xy[2,:])
  #     N = convert(Int, (e[d]+1)*(e[d]+2)/2 )
  #     P = zeros(cub.numnodes, N)
  #     ptr = 1
  #     for r = 0:e[d]
  #       for j = 0:r
  #         i = r-j
  #         P[:,ptr] = OrthoPoly.proriolpoly(x, y, i, j)
  #         ptr += 1
  #       end
  #     end
  #     @test ≈(P*C, I(cub.numnodes), atol=1e-14)
  #   end
  # end

  # @testset "Testing SummationByParts.boundaryoperators (TriSymCub method)" begin
  #   # check Ex, Ey by comparing with appropriate integral of divergence
  #   for d = 1:4
  #     cub, vtx = Cubature.getTetCubatureGamma(2*d-1, Float64)
  #     w = SymCubatures.calcweights(cub)
  #     Ex, Ey = SummationByParts.boundaryoperators(cub, vtx, d)
  #     H = diagm(w)
  #     xy = SymCubatures.calcnodes(cub, vtx)     
  #     x = vec(xy[1,:]); y = vec(xy[2,:])
  #     for r = 0:2*d-1
  #       for j = 0:r
  #         i = r-j
  #         if i > d
  #           u = x.^d
  #           dudx = d.*x.^(d-1)
  #           dudy = zeros(size(x))
  #           v = (x.^(i-d)).*(y.^j)
  #           dvdx = (i-d).*(x.^(i-d-1)).*(y.^j)
  #           dvdy = (x.^(i-d)).*(j.*y.^max(0,j-1))
  #         elseif j > d
  #           u = y.^d
  #           dudx = zeros(size(y))
  #           dudy = d.*y.^(d-1)
  #           v = (x.^i).*(y.^(j-d))
  #           dvdx = (i.*x.^max(0,i-1)).*(y.^(j-d))
  #           dvdy = (x.^i).*((j-d).*y.^(j-d-1))
  #         else # i <= d, j <= d
  #           u = x.^i
  #           dudx = i.*x.^max(0,i-1)
  #           dudy = zeros(size(x))
  #           v = y.^j
  #           dvdx = zeros(size(y))
  #           dvdy = j.*y.^max(0,j-1)
  #         end
  #         @test ≈(u'*Ex*v, dudx'*H*v + u'*H*dvdx, atol=1e-14)
  #         @test ≈(u'*Ey*v, dudy'*H*v + u'*H*dvdy, atol=1e-14)
  #       end
  #     end
  #   end   
  # end

  # @testset "Testing SummationByParts.boundaryoperators (TetSymCub method)" begin
  #   # check Ex, Ey, Ez by comparing with appropriate integral of divergence
  #   for d = 1:4
  #     cub, vtx = Cubature.getTetCubatureGamma(2*d-1, Float64)
  #     w = SymCubatures.calcweights(cub)
  #     Ex, Ey, Ez = SummationByParts.boundaryoperators(cub, vtx, d)
  #     H = diagm(w)
  #     xyz = SymCubatures.calcnodes(cub, vtx)      
  #     x = vec(xyz[1,:]); y = vec(xyz[2,:]); z = vec(xyz[3,:])
  #     for r = 0:2*d-1
  #       for k = 0:r
  #         for j = 0:r-k
  #           i = r-j-k
  #           if i > d
  #             u = x.^d
  #             dudx = d.*x.^(d-1)
  #             dudy = zeros(size(x))
  #             dudz = zeros(size(x))
  #             v = (x.^(i-d)).*(y.^j).*(z.^k)
  #             dvdx = (i-d).*(x.^(i-d-1)).*(y.^j).*(z.^k)
  #             dvdy = (x.^(i-d)).*(j.*y.^max(0,j-1)).*(z.^k)
  #             dvdz = (x.^(i-d)).*(y.^j).*(k.*z.^max(0,k-1))
  #           elseif j > d
  #             u = y.^d
  #             dudx = zeros(size(y))
  #             dudy = d.*y.^(d-1)
  #             dudz = zeros(size(y))
  #             v = (x.^i).*(y.^(j-d)).*(z.^k)
  #             dvdx = (i.*x.^max(0,i-1)).*(y.^(j-d)).*(z.^k)
  #             dvdy = (x.^i).*((j-d).*y.^(j-d-1)).*(z.^k)
  #             dvdz = (x.^i).*(y.^(j-d)).*(k.*z.^max(0,k-1))
  #           elseif k > d
  #             u = z.^d
  #             dudx = zeros(size(z))
  #             dudy = zeros(size(z))
  #             dudz = d.*z.^(d-1)
  #             v = (x.^i).*(y.^j).*(z.^(k-d))
  #             dvdx = (i.*x.^max(0,i-1)).*(y.^j).*(z.^(k-d))
  #             dvdy = (x.^i).*(j.*y.^max(0,j-1)).*(z.^(k-d))
  #             dvdz = (x.^i).*(y.^j).*((k-d).*z.^(k-d-1))
  #           elseif i+j > d
  #             ku, kv = (i > j ? [0,k] : [k,0])
  #             u = (x.^i).*(z.^ku)
  #             dudx = (i.*x.^max(0,i-1)).*(z.^ku)
  #             dudy = zeros(size(x))
  #             dudz = (x.^i).*(ku.*z.^max(0,ku-1))
  #             v = (y.^j).*(z.^kv)
  #             dvdx = zeros(size(y))
  #             dvdy = (j.*y.^max(0,j-1)).*(z.^kv)
  #             dvdz = (y.^j).*(kv.*z.^max(0,kv-1))
  #           elseif i+k > d
  #             ju, jv = (i > k ? [0,j] : [j,0])
  #             u = (x.^i).*(y.^ju)
  #             dudx = (i.*x.^max(0,i-1)).*(y.^ju)
  #             dudy = (x.^i).*(ju.*y.^max(0,ju-1))              
  #             dudz = zeros(size(x))
  #             v = (z.^k).*(y.^jv)
  #             dvdx = zeros(size(z))
  #             dvdy = (z.^k).*(jv.*y.^max(0,jv-1))
  #             dvdz = (k.*z.^max(0,k-1)).*(y.^jv)
  #           else
  #             iu, iv = (j > k ? [0,i] : [i,0])
  #             u = (x.^iu).*(y.^j)
  #             dudx = (iu.*x.^max(0,iu-1)).*(y.^j)
  #             dudy = (x.^iu).*(j.*y.^max(0,j-1))
  #             dudz = zeros(size(x))              
  #             v = (x.^iv).*(z.^k)
  #             dvdx = (iv.*x.^max(0,iv-1)).*(z.^k)
  #             dvdy = zeros(size(x))
  #             dvdz = (x.^iv).*(k.*z.^max(0,k-1))
  #           end
  #           @test ≈(u'*Ex*v, dudx'*H*v + u'*H*dvdx, atol=1e-14)
  #           @test ≈(u'*Ey*v, dudy'*H*v + u'*H*dvdy, atol=1e-14)
  #           @test ≈(u'*Ez*v, dudz'*H*v + u'*H*dvdz, atol=1e-14)
  #         end
  #       end
  #     end   
  #   end
  # end

  # @testset "Testing SummationByParts.boundaryoperator! (TriFace method)" begin
  #   # check by comparing with Ex, Ey produced by boundaryoperators
  #   for d = 1:4
  #     cub, vtx = Cubature.getTriCubatureGamma(2*d-1, Float64)
  #     w = SymCubatures.calcweights(cub)
  #     Ex, Ey = SummationByParts.boundaryoperators(cub, vtx, d)
  #     face = TriFace{Float64}(d, cub, vtx)
  #     E = zeros(cub.numnodes,cub.numnodes)
  #     SummationByParts.boundaryoperator!(face, 1, E)
  #     @test ≈(E, Ex, atol=1e-14)
  #     SummationByParts.boundaryoperator!(face, 2, E)
  #     @test ≈(E, Ey, atol=1e-14)
  #   end
  # end

  # @testset "Testing SummationByParts.boundaryoperator! (TetFace method)" begin
  #   # check by comparing with Ex, Ey produced by boundaryoperators
  #   for d = 1:4
  #     cub, vtx = Cubature.getTetCubatureGamma(2*d-1, Float64)
  #     w = SymCubatures.calcweights(cub)
  #     Ex, Ey, Ez = SummationByParts.boundaryoperators(cub, vtx, d)
  #     face = TetFace{Float64}(d, cub, vtx)
  #     E = zeros(cub.numnodes,cub.numnodes)
  #     SummationByParts.boundaryoperator!(face, 1, E)
  #     @test ≈(E, Ex, atol=1e-14)
  #     SummationByParts.boundaryoperator!(face, 2, E)
  #     @test ≈(E, Ey, atol=1e-14)
  #     SummationByParts.boundaryoperator!(face, 3, E)
  #     @test ≈(E, Ez, atol=1e-14)
  #   end
  # end

  # @testset "Testing SummationByParts.boundarymassmatrix (TriSymCub method)" begin
  #   # check that mass matrix can be assembled into Ex and Ey
  #   for d = 1:4
  #     cub, vtx = Cubature.getTriCubatureGamma(2*d-1, Float64)
  #     Ex, Ey = SummationByParts.boundaryoperators(cub, vtx, d)
  #     Hbndry, bndindx = SummationByParts.boundarymassmatrix(cub, vtx, d)
  #     Ex_inject = zeros(size(Ex))
  #     Ey_inject = zeros(size(Ey))
  #     Ex_inject[bndindx[:,2],bndindx[:,2]] = Hbndry
  #     Ex_inject[bndindx[:,3],bndindx[:,3]] -= Hbndry
  #     Ey_inject[bndindx[:,1],bndindx[:,1]] -= Hbndry
  #     Ey_inject[bndindx[:,2],bndindx[:,2]] += Hbndry
  #     @test ≈(Ex_inject, Ex, atol=1e-14)
  #     @test ≈(Ey_inject, Ey, atol=1e-14)
  #   end
  # end

  # @testset "Testing SummationByParts.boundarymassmatrix (TetSymCub method)" begin
  #   # check that mass matrix can be assembled into Ex and Ey
  #   for d = 1:3
  #     cub, vtx = Cubature.getTetCubatureGamma(2*d-1, Float64)
  #     Ex, Ey, Ez = SummationByParts.boundaryoperators(cub, vtx, d)
  #     Hbndry, bndindx = SummationByParts.boundarymassmatrix(cub, vtx, d)
  #     Ex_inject = zeros(size(Ex))
  #     Ey_inject = zeros(size(Ey))
  #     Ez_inject = zeros(size(Ez))
  #     Ex_inject[bndindx[:,3],bndindx[:,3]] += Hbndry
  #     Ex_inject[bndindx[:,4],bndindx[:,4]] -= Hbndry
  #     Ey_inject[bndindx[:,3],bndindx[:,3]] += Hbndry
  #     Ey_inject[bndindx[:,2],bndindx[:,2]] -= Hbndry
  #     Ez_inject[bndindx[:,3],bndindx[:,3]] += Hbndry
  #     Ez_inject[bndindx[:,1],bndindx[:,1]] -= Hbndry
  #     @test ≈(Ex_inject, Ex, atol=1e-14)
  #     @test ≈(Ey_inject, Ey, atol=1e-14)
  #     @test ≈(Ez_inject, Ez, atol=1e-14)
  #   end
  # end

  # @testset "Testing SummationByParts.accuracyconstraints (TriSymCub method)" begin
  #   # check that the null-space of the constraint Jacobian is the correct size
  #   # this is not an adequate unit test.
  #   sizenull = [0, 0, 1, 3]
  #   for d = 1:4
  #     cub, vtx = Cubature.getTriCubatureGamma(2*d-1, Float64)
  #     face = TriFace{Float64}(d, cub, vtx)
  #     Q = zeros(cub.numnodes,cub.numnodes,2)
  #     SummationByParts.boundaryoperator!(face, 1, view(Q,:,:,1))
  #     SummationByParts.boundaryoperator!(face, 2, view(Q,:,:,2))
  #     Q .*= 0.5 #scale!(Q, 0.5)
  #     A, bx, by = SummationByParts.accuracyconstraints(cub, vtx, d, Q)
  #     @test size(nullspace(A),2) == sizenull[d]
  #   end
  # end

  # @testset "Testing SummationByParts.accuracyconstraints (TetSymCub method)" begin
  #   # check that the null-space of the constraint Jacobian is the correct size
  #   # this is not an adequate unit test.
  #   sizenull = [0, 0, 6, 45]
  #   for d = 1:4
  #     cub, vtx = Cubature.getTetCubatureGamma(2*d-1, Float64)
  #     face = TetFace{Float64}(d, cub, vtx)
  #     Q = zeros(cub.numnodes,cub.numnodes,3)
  #     SummationByParts.boundaryoperator!(face, 1, view(Q,:,:,1))
  #     SummationByParts.boundaryoperator!(face, 2, view(Q,:,:,2))
  #     SummationByParts.boundaryoperator!(face, 3, view(Q,:,:,3))
  #     Q .*= 0.5 #scale!(Q, 0.5)
  #     A, bx, by, bz = SummationByParts.accuracyconstraints(cub, vtx, d, Q)
  #     @test size(nullspace(A),2) == sizenull[d]
  #   end
  # end
  
  # @testset "Testing SummationByParts.commuteerror (TriSymCub method)" begin
  #   reducedsol = (Float64[], Float64[], Float64[0, 0])
  #   error = [0, 0, 0.5*2.324812265031167] # error based on particular Q
  #   for d = 1:3
  #     cub, vtx = Cubature.getTriCubatureGamma(2*d-1, Float64)
  #     face = TriFace{Float64}(d, cub, vtx)
  #     Q = zeros(cub.numnodes,cub.numnodes,2)
  #     SummationByParts.boundaryoperator!(face, 1, view(Q,:,:,1))
  #     SummationByParts.boundaryoperator!(face, 2, view(Q,:,:,2))
  #     Q .*= 0.5 #scale!(Q, 0.5)
  #     w = SymCubatures.calcweights(cub)
  #     A, bx, by = SummationByParts.accuracyconstraints(cub, vtx, d, Q)
  #     # build Q that satisfies the accuracy constraints
  #     x = A\bx; y = A\by
  #     for row = 2:cub.numnodes
  #       offset = convert(Int, (row-1)*(row-2)/2)
  #       for col = 1:row-1
  #         Q[row,col,1] += x[offset+col]
  #         Q[col,row,1] -= x[offset+col]
  #         Q[row,col,2] += y[offset+col]
  #         Q[col,row,2] -= y[offset+col]
  #       end
  #     end
  #     Z = nullspace(A)
  #     f, dfdx = SummationByParts.commuteerror(w, view(Q,:,:,1), view(Q,:,:,2),
  #                                             Z, reducedsol[d])
  #     @test ≈(f, error[d], atol=1e-14)
  #   end
  # end

  # @testset "Testing SummationByParts.buildoperators (TriSymCub method)" begin
  #   for d = 1:4
  #     cub, vtx = Cubature.getTriCubatureGamma(2*d-1, Float64)
  #     w, Q = SummationByParts.buildoperators(cub, vtx, d)
  #     Dx = diagm(1.0 ./w)*Q[:,:,1]
  #     Dy = diagm(1.0 ./w)*Q[:,:,2]
  #     xy = SymCubatures.calcnodes(cub, vtx)
  #     x = vec(xy[1,:]); y = vec(xy[2,:])
  #     for r = 0:d
  #       for j = 0:r
  #         i = r-j
  #         u = (x.^i).*(y.^j)
  #         dudx = (i.*x.^max(0,i-1)).*(y.^j)
  #         dudy = (x.^i).*(j.*y.^max(0,j-1))
  #         @test ≈(Dx*u, dudx; atol=5e-13)
  #         @test ≈(Dy*u, dudy; atol=5e-13)
  #       end
  #     end
  #   end
  # end

  # @testset "Testing SummationByParts.buildoperators (TriSymCub method, internal nodes)" begin
  #   for d = 1:4
  #     cub, vtx = Cubature.getTriCubatureOmega(2*d, Float64)
  #     w, Q = SummationByParts.buildoperators(cub, vtx, d)
  #     Dx = diagm(1.0 ./w)*Q[:,:,1]
  #     Dy = diagm(1.0 ./w)*Q[:,:,2]
  #     xy = SymCubatures.calcnodes(cub, vtx)  
  #     x = vec(xy[1,:]); y = vec(xy[2,:])
  #     for r = 0:d
  #       for j = 0:r
  #         i = r-j
  #         u = (x.^i).*(y.^j)
  #         dudx = (i.*x.^max(0,i-1)).*(y.^j)
  #         dudy = (x.^i).*(j.*y.^max(0,j-1))
  #         @test ≈(Dx*u, dudx, atol=5e-13)
  #         @test ≈(Dy*u, dudy, atol=5e-13)
  #       end
  #     end
  #   end
  # end

  # @testset "Testing SummationByParts.buildoperators (TetSymCub method)" begin
  #   tol = [1e-12; 1e-12; 1e-12; 5e-8]
  #   for d = 1:4
  #     cub, vtx = Cubature.getTetCubatureGamma(2*d-1, Float64)
  #     w, Q = SummationByParts.buildoperators(cub, vtx, d)
  #     Dx = diagm(1.0 ./w)*Q[:,:,1]
  #     Dy = diagm(1.0 ./w)*Q[:,:,2]
  #     Dz = diagm(1.0 ./w)*Q[:,:,3]
  #     xyz = SymCubatures.calcnodes(cub, vtx)      
  #     x = vec(xyz[1,:]); y = vec(xyz[2,:]); z = vec(xyz[3,:])
  #     for r = 0:d
  #       for k = 0:r
  #         for j = 0:r-k
  #           i = r-j-k
  #           u = (x.^i).*(y.^j).*(z.^k)
  #           dudx = (i.*x.^max(0,i-1)).*(y.^j).*(z.^k)
  #           dudy = (x.^i).*(j.*y.^max(0,j-1)).*(z.^k)
  #           dudz = (x.^i).*(y.^j).*(k.*z.^max(0,k-1))
  #           @test ≈(Dx*u, dudx, atol=tol[d])
  #           @test ≈(Dy*u, dudy, atol=tol[d])
  #           @test ≈(Dz*u, dudz, atol=tol[d])
  #         end
  #       end
  #     end
  #   end
  # end

  # @testset "Testing SummationByParts.buildoperators (TetSymCub method, internal nodes)" begin
  #   tol = [1e-12; 1e-12; 1e-12; 1e-12] #5e-8]
  #   for d = 1:4
  #     cub, vtx = Cubature.getTetCubatureOmega(2*d-1, Float64)
  #     w, Q = SummationByParts.buildoperators(cub, vtx, d)
  #     Dx = diagm(1.0 ./w)*Q[:,:,1]
  #     Dy = diagm(1.0 ./w)*Q[:,:,2]
  #     Dz = diagm(1.0 ./w)*Q[:,:,3]
  #     xyz = SymCubatures.calcnodes(cub, vtx)      
  #     x = vec(xyz[1,:]); y = vec(xyz[2,:]); z = vec(xyz[3,:])
  #     for r = 0:d
  #       for k = 0:r
  #         for j = 0:r-k
  #           i = r-j-k
  #           u = (x.^i).*(y.^j).*(z.^k)
  #           dudx = (i.*x.^max(0,i-1)).*(y.^j).*(z.^k)
  #           dudy = (x.^i).*(j.*y.^max(0,j-1)).*(z.^k)
  #           dudz = (x.^i).*(y.^j).*(k.*z.^max(0,k-1))
  #           @test ≈(Dx*u, dudx, atol=tol[d])
  #           @test ≈(Dy*u, dudy, atol=tol[d])
  #           @test ≈(Dz*u, dudz, atol=tol[d])
  #         end
  #       end
  #     end
  #   end
  # end

  # @testset "Testing SummationByParts.buildoperators (spectral-element method)" begin
  #   e = [1;3;4;5]
  #   for d = 1:4
  #     cub, vtx = Cubature.getTriCubatureGamma(2*d-1, Float64)
  #     w, Q = SummationByParts.buildoperators(cub, vtx, d, e[d])
  #     Dx = diagm(1.0 ./w)*Q[:,:,1]
  #     Dy = diagm(1.0 ./w)*Q[:,:,2]
  #     xy = SymCubatures.calcnodes(cub, vtx)
  #     x = vec(xy[1,:]); y = vec(xy[2,:])
  #     for r = 0:d
  #       for j = 0:r
  #         i = r-j
  #         u = (x.^i).*(y.^j)
  #         dudx = (i.*x.^max(0,i-1)).*(y.^j)
  #         dudy = (x.^i).*(j.*y.^max(0,j-1))
  #         @test ≈(Dx*u, dudx, atol=1e-12)
  #         @test ≈(Dy*u, dudy, atol=1e-12)
  #       end
  #     end
  #   end
  # end

  # @testset "Testing SummationByParts.buildsparseoperators (TriSymCub method)" begin
  #   for d = 1:4 
  #     cub, vtx = Cubature.getTriCubatureDiagE(2*d, Float64)
  #     w, Q = SummationByParts.buildsparseoperators(cub, vtx, d)
  #     Dx = diagm(1.0 ./w)*Q[:,:,1]
  #     Dy = diagm(1.0 ./w)*Q[:,:,2]
  #     xy = SymCubatures.calcnodes(cub, vtx)
  #     x = vec(xy[1,:]); y = vec(xy[2,:])
  #     for r = 0:d
  #       for j = 0:r
  #         i = r-j
  #         u = (x.^i).*(y.^j)
  #         dudx = (i.*x.^max(0,i-1)).*(y.^j)
  #         dudy = (x.^i).*(j.*y.^max(0,j-1))
  #         @test ≈(Dx*u, dudx; atol=5e-12)
  #         @test ≈(Dy*u, dudy; atol=5e-12)
  #       end
  #     end
  #   end
  # end

  # @testset "Testing SummationByParts.buildMinConditionOperators (TriSymCub method, vertices=true)" begin
  #   for d = 1:4 
  #     cub, vtx = Cubature.getTriCubatureDiagE(2*d, Float64)
  #     w, Q = SummationByParts.buildMinConditionOperators(cub, vtx, d)
  #     Dx = diagm(1.0 ./w)*Q[:,:,1]
  #     Dy = diagm(1.0 ./w)*Q[:,:,2]
  #     xy = SymCubatures.calcnodes(cub, vtx)  
  #     x = vec(xy[1,:]); y = vec(xy[2,:])
  #     for r = 0:d
  #       for j = 0:r
  #         i = r-j
  #         u = (x.^i).*(y.^j)
  #         dudx = (i.*x.^max(0,i-1)).*(y.^j)
  #         dudy = (x.^i).*(j.*y.^max(0,j-1))
  #         @test ≈(Dx*u, dudx, atol=5e-13)
  #         @test ≈(Dy*u, dudy, atol=5e-13)
  #       end
  #     end
  #   end
  # end

  # @testset "Testing SummationByParts.buildMinConditionOperators (TriSymCub method, vertices=false)" begin
  #   for d = 1:4 
  #     cub, vtx = Cubature.getTriCubatureDiagE(2*d, Float64, vertices=false)
  #     w, Q = SummationByParts.buildMinConditionOperators(cub, vtx, d,
  #                                                        vertices=false)
  #     Dx = diagm(1.0 ./w)*Q[:,:,1]
  #     Dy = diagm(1.0 ./w)*Q[:,:,2]
  #     xy = SymCubatures.calcnodes(cub, vtx)  
  #     x = vec(xy[1,:]); y = vec(xy[2,:])
  #     for r = 0:d
  #       for j = 0:r
  #         i = r-j
  #         u = (x.^i).*(y.^j)
  #         dudx = (i.*x.^max(0,i-1)).*(y.^j)
  #         dudy = (x.^i).*(j.*y.^max(0,j-1))
  #         @test ≈(Dx*u, dudx, atol=5e-13)
  #         @test ≈(Dy*u, dudy, atol=5e-13)
  #       end
  #     end
  #   end
  # end

  # @testset "Testing SummationByParts.buildMinConditionOperators (TetSymCub method)" begin
  #   tol = [1e-12; 1e-12; 1e-12; 1e-12]
  #   # tol = [1e-12; 1e-12] #; 1e-12; 1e-12]
  #   for d = 1:2 # d=3,4 are too slow for tests
  #     # cub, vtx = Cubature.getTetCubatureDiagE(2*d, Float64, vertices=false)
  #     cub, vtx = Cubature.getTetCubatureDiagE(2*d, Float64, faceopertype=:Omega)
  #     w, Q = SummationByParts.buildMinConditionOperators(cub, vtx, d, tol=1e-2)
  #     Dx = diagm(1.0 ./w)*Q[:,:,1]
  #     Dy = diagm(1.0 ./w)*Q[:,:,2]
  #     Dz = diagm(1.0 ./w)*Q[:,:,3]
  #     xyz = SymCubatures.calcnodes(cub, vtx)      
  #     x = vec(xyz[1,:]); y = vec(xyz[2,:]); z = vec(xyz[3,:])
  #     for r = 0:d
  #       for k = 0:r
  #         for j = 0:r-k
  #           i = r-j-k
  #           u = (x.^i).*(y.^j).*(z.^k)
  #           dudx = (i.*x.^max(0,i-1)).*(y.^j).*(z.^k)
  #           dudy = (x.^i).*(j.*y.^max(0,j-1)).*(z.^k)
  #           dudz = (x.^i).*(y.^j).*(k.*z.^max(0,k-1))
  #           @test ≈(Dx*u, dudx, atol=tol[d])
  #           @test ≈(Dy*u, dudy, atol=tol[d])
  #           @test ≈(Dz*u, dudz, atol=tol[d])
  #         end
  #       end
  #     end
  #   end
  # end

  # # @testset "Testing SummationByParts.buildsparseoperators (TetSymCub method, internal nodes)" begin
  # #   tol = [1e-12; 1e-12; 1e-12; 1e-12] #5e-8]
  # #   for d = 1:3
  # #     cub, vtx = tetcubature(2*d+1, Float64, internal=true)
  # #     w, Q = SummationByParts.buildsparseoperators(cub, vtx, d)
  # #     Dx = diagm(1./w)*Q[:,:,1]
  # #     Dy = diagm(1./w)*Q[:,:,2]
  # #     Dz = diagm(1./w)*Q[:,:,3]
  # #     xyz = SymCubatures.calcnodes(cub, vtx)      
  # #     x = vec(xyz[1,:]); y = vec(xyz[2,:]); z = vec(xyz[3,:])
  # #     for r = 0:d
  # #       for k = 0:r
  # #         for j = 0:r-k
  # #           i = r-j-k
  # #           u = (x.^i).*(y.^j).*(z.^k)
  # #           dudx = (i.*x.^max(0,i-1)).*(y.^j).*(z.^k)
  # #           dudy = (x.^i).*(j.*y.^max(0,j-1)).*(z.^k)
  # #           dudz = (x.^i).*(y.^j).*(k.*z.^max(0,k-1))
  # #           @fact Dx*u --> roughly(dudx, atol=tol[d])
  # #           @fact Dy*u --> roughly(dudy, atol=tol[d])
  # #           @fact Dz*u --> roughly(dudz, atol=tol[d])
  # #         end
  # #       end
  # #     end
  # #   end
  # # end

  # @testset "Testing SummationByParts.getnodepermutation (TriSymCub method)" begin
  #   # check that vertices are first and edge nodes are ordered correctly
  #   for d = 1:4
  #     cub, vtx = Cubature.getTriCubatureGamma(2*d-1, Float64)
  #     perm, faceperm = SummationByParts.getnodepermutation(cub, d)
  #     xy = SymCubatures.calcnodes(cub, vtx)
  #     x = vec(xy[1,perm])
  #     y = vec(xy[2,perm])
  #     # check vertices
  #     @test ≈(x[1:3], vtx[:,1], atol=1e-15)
  #     @test ≈(y[1:3], vtx[:,2], atol=1e-15)
  #     ptr = 3
  #     # check ordering of x nodes on first edge
  #     @test issorted(x[ptr+1:ptr+d-1]) == true
  #     @test ≈(y[ptr+1:ptr+d-1], -ones(d-1), atol=1e-15)
  #     ptr += (d-1)
  #     # check ordering of x and y nodes on second edge
  #     @test issorted(x[ptr+1:ptr+d-1], rev=true) == true
  #     @test issorted(y[ptr+1:ptr+d-1]) == true
  #     ptr += (d-1)
  #     # check ordering of y nodes on third edge
  #     @test ≈( x[ptr+1:ptr+d-1], -ones(d-1), atol=1e-15)
  #     @test issorted(y[ptr+1:ptr+d-1], rev=true) == true
  #   end
  # end

  @testset "Testing SummationByParts.getnodepermutation (TetSymCub method)" begin
    # check that vertices are first, edge nodes are ordered correctly, and face
    # nodes lie on the expected faces
    for d = 1:3 #4 not working
      cub, vtx = Cubature.getTetCubatureGamma(2*d-1, Float64)
      perm, faceperm = SummationByParts.getnodepermutation(cub, d)
      xyz = SymCubatures.calcnodes(cub, vtx)
      x = vec(xyz[1,perm])
      y = vec(xyz[2,perm])
      z = vec(xyz[3,perm])
      # check vertices
      @test ≈(x[1:4], vtx[:,1], atol=1e-15)
      @test ≈(y[1:4], vtx[:,2], atol=1e-15)
      @test ≈(z[1:4], vtx[:,3], atol=1e-15)
      ptr = 4
      # check ordering of x nodes on first edge
      @test issorted(x[ptr+1:ptr+d-1]) == true
      @test ≈(y[ptr+1:ptr+d-1], -ones(d-1), atol=1e-15)
      @test ≈(z[ptr+1:ptr+d-1], -ones(d-1), atol=1e-15)
      ptr += (d-1)
      # check ordering of x and y nodes on second edge
      @test issorted(x[ptr+1:ptr+d-1], rev=true) == true
      @test issorted(y[ptr+1:ptr+d-1]) == true
      @test ≈(z[ptr+1:ptr+d-1], -ones(d-1), atol=1e-15)
      ptr += (d-1)
      # check ordering of y nodes on third edge
      @test ≈(x[ptr+1:ptr+d-1], -ones(d-1), atol=1e-15)
      @test issorted(y[ptr+1:ptr+d-1], rev=true) == true
      @test ≈(z[ptr+1:ptr+d-1], -ones(d-1), atol=1e-15)
      ptr += (d-1)
      # check ordering of z nodes on fourth edge
      @test ≈(x[ptr+1:ptr+d-1], -ones(d-1), atol=1e-15)
      @test ≈(y[ptr+1:ptr+d-1], -ones(d-1), atol=1e-15)
      @test issorted(z[ptr+1:ptr+d-1]) == true
      ptr += (d-1)
      # check ordering of x and z on fifth edge
      @test issorted(x[ptr+1:ptr+d-1], rev=true) == true
      @test ≈(y[ptr+1:ptr+d-1], -ones(d-1), atol=1e-15)
      @test issorted(z[ptr+1:ptr+d-1]) == true
      ptr += (d-1)
      # check ordering of y and z on sixth edge
      @test ≈(x[ptr+1:ptr+d-1], -ones(d-1), atol=1e-15)
      @test issorted(y[ptr+1:ptr+d-1], rev=true) == true
      @test issorted(z[ptr+1:ptr+d-1]) == true
      ptr += (d-1)
      # check that face nodes lie on appropriate faces
      numface = div((d-1)*(d-2),2)
      # check that z = -1 on face 1
      @test ≈(z[ptr+1:ptr+numface], -ones(numface), atol=1e-15)
      ptr += numface
      # check that y = -1 on face 2
      @test ≈(y[ptr+1:ptr+numface], -ones(numface), atol=1e-15)
      ptr += numface
      # check that z = -1 - x - y on face 3
      @test ≈(z[ptr+1:ptr+numface], -ones(numface)
                                            -x[ptr+1:ptr+numface]
                                            -y[ptr+1:ptr+numface], atol=1e-15)
      ptr += numface
      # check that x = -1 on face 4
      @test ≈(x[ptr+1:ptr+numface], -ones(numface), atol=1e-15)
      ptr += numface
    end
  end

end
