@testset "Testing SummationByParts Module (utils.jl file)..." begin

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE,
              getTetSBPGamma, getTetSBPOmega)
    @eval begin
      @testset "Testing SummationByParts.getNumFaceNodes $(string($TSBP)) method" begin
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          @test getNumFaceNodes(sbp) == SymCubatures.getnumfacenodes(sbp.cub)
        end
      end
    end
  end
  
  @testset "Testing SummationByParts.calcminnodedistance (TriSBP method)" begin
    mindist = [1.0; 0.2357022603955159; 0.1487006728783353; 0.09492895652255572]
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      vtx = [0. 0.; 1. 0.; 0. 1.]
      @test ≈(calcminnodedistance(sbp, vtx), mindist[p], atol=1e-13)
    end
  end

  @testset "Testing SummationByParts.calcminnodedistance (TetSBP method)" begin
    mindist = [1.0; 0.4330127018922193; 0.2639696512367827; 0.1366241982649621]
    for p = 1:4
      sbp = getTetSBPGamma(degree=p)
      vtx = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
      @test ≈(calcminnodedistance(sbp, vtx), mindist[p], atol=1e-13)
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      @testset "Testing buildinterpolation $(string($TSBP)) method" begin
        # this checks that polynomials of total degree d are reconstructed accurately
        numpoints = 3
        for d = 1:4
          sbp = ($TSBP)(degree=d)
          x = 2.0*rand(2,numpoints) .- 1.0
          R = SummationByParts.buildinterpolation(sbp, x)
          xsbp = calcnodes(sbp)
          # loop over all monomials
          for r = 0:d
            for j = 0:r
              i = r-j
              u = vec((x[1,:].^i).*(x[2,:].^j))
              usbp = vec((xsbp[1,:].^i).*(xsbp[2,:].^j))
              uinterp = R*usbp
              @test ≈(uinterp, u, atol=1e-12)
            end
          end
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega)
    @eval begin
      @testset "Testing buildinterpolation $(string($TSBP)) method" begin
        # this checks that polynomials of total degree d are reconstructed accurately
        numpoints = 10
        for d = 1:4
          sbp = ($TSBP)(degree=d)
          x = 2.0*rand(3,numpoints) .- 1.0
          R = SummationByParts.buildinterpolation(sbp, x)
          xsbp = calcnodes(sbp)
          # loop over all monomials
          for r = 0:d
            for k = 0:r
              for j = 0:r-k
                i = r-j-k
                u = vec((x[1,:].^i).*(x[2,:].^j).*(x[3,:].^k))
                usbp = vec((xsbp[1,:].^i).*(xsbp[2,:].^j).*(xsbp[3,:].^k))
                uinterp = R*usbp
                @test ≈(uinterp, u, atol=1e-11)
              end
            end
          end
        end
      end
    end
  end

  @testset "Testing SummationByParts.permuteface!" begin
    # test 2D version
    permvec = [2, 1, 3, 4]
    vals = rand(5, 4)
    vals2 = copy(vals)  # permute a copy, leaving the original unchanged
    workarr = zeros(size(vals))

    SummationByParts.permuteface!(permvec, workarr, vals2)

    for i=1:4
      for j=1:5
        @test ≈(vals[j, i], vals2[j, permvec[i] ], atol=1e-13)
      end
    end

    # test 1D version
    vals = rand(4)
    vals2 = copy(vals)
    workarr = zeros(size(vals))
    SummationByParts.permuteface!(permvec, workarr, vals2)

    for i=1:4
      @test ≈(vals[i], vals2[permvec[i]], atol=1e-13)
    end

  end  # end @testset testing permuteface!)

  @testset "Testing SummationByParts.permuteinterface!" begin

    # test 3D version
    nfaces = 2
    ndofpernode = 5

    # create an SBP with non-trivial number of face nodes
    sbp = getTriSBPOmega(degree=3)
    ref_verts = [-1. 1 -1; -1 -1 1]
    sbpface = TriFace{Float64}(sbp.degree, sbp.cub, Matrix(ref_verts'))

    # create some data
    q = rand(ndofpernode, sbpface.numnodes, nfaces)
    q2 = copy(q)

    ifaces = Array{Interface}(undef, 2)
    for i=1:nfaces
      ifaces[i] = Interface(1,1,1,1,1)
    end

    SummationByParts.permuteinterface!(sbpface, ifaces, q2)

    for iface=1:nfaces
      for j=1:sbpface.numnodes
        for k=1:ndofpernode
          @test ≈(q[k, j, iface], q2[k, sbpface.nbrperm[j, 1], iface], atol=1e-13)
        end
      end
    end

    # test 2D version
    q = rand(sbpface.numnodes, nfaces)
    q2 = copy(q)
    SummationByParts.permuteinterface!(sbpface, ifaces, q2)
    for iface=1:nfaces
      for j=1:sbpface.numnodes
          @test ≈(q[j, iface], q2[sbpface.nbrperm[j, 1], iface], atol=1e-13)
      end
    end

  end  # end @testset Testing permuteinterface!)
        
  @testset "Testing SummationByParts.basispursuit!" begin
    # check that a sparse solution is produced for a given problem; this needs a
    # better test
    A = [0.3240559919365473 0.7670391351853123 0.2760840637306585 0.8293572108203571 0.7272847946349628 0.7517453486022216 0.7489674270729445 0.06473192617477141 0.41202501224251553 0.0459986733625628 0.7718001847233651 0.8622024510885822 0.21767987994746618 0.45714495101454133 0.6949713430771973;
         0.29803429306849205 0.4434206195098742 0.927298924512211 0.1919232898137111 0.1250762177473692 0.25670469205309976 0.5572691190255206 0.8340797045140584 0.578031022922054 0.5934658938673283 0.09352598461878059 0.5596421183290323 0.16277010184303786 0.6936829447663804 0.855561662977216;
    0.40301815537750696 0.17930610701088434 0.23351970607979333 0.5785142899752729 0.24218369036938636 0.9398474572781728 0.013075108326236595 0.7328553912656008 0.9956493310946866 0.9803588951389592 0.5555928536272843 0.7150648069092651 0.13955009358904968 0.32898589105322307 0.1877638137341433;
    0.7427577164373185 0.5656491407521753 0.20737863002402435 0.6747469757050888 0.5349678001619618 0.8435987394854569 0.9186696381949431 0.7570368254815638 0.8085807239457743 0.7015949207882941 0.20190996439947595 0.07196121725427318 0.1351874819308645 0.1322155409235981 0.05557516782948646;
    0.013605678113414399 0.11615532396259942 0.5612391398057421 0.2799012713985811 0.22205249009933548 0.37425224598561657 0.9852918434644029 0.7423720549697086 0.8815730477210602 0.5139365866782812 0.5838760058313315 0.6105184693096728 0.9544124347991654 0.938059683206228 0.9829633837233192]
    b = [0.12548517218415522; 0.8009107593937155; 0.28606421432116824; 
         0.7763203971528092; 0.2175208275793281]
    xexact = [0.7769106288562887; 0.0; 0.4792978185035988; 0.0; 0.0; 0.0; 
              0.09040334213110997; 0.13806679721708773; 0.0; 0.0; 
              -0.4343745397847615; 0.0; 0.0; 0.0; 0.0]
    x = zeros(size(A,2))

    SummationByParts.basispursuit!(A, b, x, rho=1.5, alpha=1.0, hist=false,
                                   abstol=1e-6, reltol=1e-6)
    # P = zeros(size(A,2),size(A,1))
    # idx = sortperm(abs(x), rev=true)
    # for i = 1:size(A,1)
    #   P[idx[i],i] = 1.0
    # end
    # AP = A*P
    # xexact = P*(AP\b)    
    # println("x = ",x)
    # println("xexact = ",xexact)
    @test ≈(x, xexact, rtol=1e-4)
    @test ≈(A*x, b, atol=1e-13)
  end

  @testset "Testing SummationByParts.calcSparseSolution!" begin
    # check that a sparse solution is produced for a given problem; this needs a
    # better test
    A = [0.3240559919365473 0.7670391351853123 0.2760840637306585 0.8293572108203571 0.7272847946349628 0.7517453486022216 0.7489674270729445 0.06473192617477141 0.41202501224251553 0.0459986733625628 0.7718001847233651 0.8622024510885822 0.21767987994746618 0.45714495101454133 0.6949713430771973;
         0.29803429306849205 0.4434206195098742 0.927298924512211 0.1919232898137111 0.1250762177473692 0.25670469205309976 0.5572691190255206 0.8340797045140584 0.578031022922054 0.5934658938673283 0.09352598461878059 0.5596421183290323 0.16277010184303786 0.6936829447663804 0.855561662977216;
    0.40301815537750696 0.17930610701088434 0.23351970607979333 0.5785142899752729 0.24218369036938636 0.9398474572781728 0.013075108326236595 0.7328553912656008 0.9956493310946866 0.9803588951389592 0.5555928536272843 0.7150648069092651 0.13955009358904968 0.32898589105322307 0.1877638137341433;
    0.7427577164373185 0.5656491407521753 0.20737863002402435 0.6747469757050888 0.5349678001619618 0.8435987394854569 0.9186696381949431 0.7570368254815638 0.8085807239457743 0.7015949207882941 0.20190996439947595 0.07196121725427318 0.1351874819308645 0.1322155409235981 0.05557516782948646;
    0.013605678113414399 0.11615532396259942 0.5612391398057421 0.2799012713985811 0.22205249009933548 0.37425224598561657 0.9852918434644029 0.7423720549697086 0.8815730477210602 0.5139365866782812 0.5838760058313315 0.6105184693096728 0.9544124347991654 0.938059683206228 0.9829633837233192]
    b = [0.12548517218415522; 0.8009107593937155; 0.28606421432116824; 
         0.7763203971528092; 0.2175208275793281]
    xexact = [0.7769106288562887; 0.0; 0.4792978185035988; 0.0; 0.0; 0.0; 
              0.09040334213110997; 0.13806679721708773; 0.0; 0.0; 
              -0.4343745397847615; 0.0; 0.0; 0.0; 0.0]
    x = zeros(size(A,2))

    SummationByParts.calcSparseSolution!(A, b, x)
    @test ≈(x, xexact, rtol=1e-13)
    @test ≈(A*x, b, atol=1e-13)
  end

  @testset "Testing absMatrix!" begin
    # construct a symmetric matrix
    n = 10
    Q, R = qr(rand(n,n))
    λ = randn(n)
    A = Q*diagm(λ)*Q'
    Acheck = Q*diagm(abs.(λ))*Q'
    Aabs = zeros(size(A))
    SummationByParts.absMatrix!(A, Aabs)
    @test ≈(Aabs, Acheck, rtol=1e-13)        
  end

  for T = (Float64, ComplexF64)
    @eval begin
      @testset "Testing calcMatrixEigs! and calcMatrixEigs_rev! for DataType $(string($T))" begin
        # construct a symmetric matrix whose eigenvalues are the parameters
        n = 10
        Q, R = qr(rand(n,n))
        x = rand(($T), n)
        idx = sortperm(x, lt=SummationByParts.compareEigs)
        x[:] = x[idx]
        A = Q*diagm(x)*Q'
        λ = zeros(($T), n)
        SummationByParts.calcMatrixEigs!(A, λ)
        @test ≈(λ, x, rtol=1e-14)            
        λ_bar = deepcopy(λ)
        A_bar = zeros(eltype(A), size(A))
        SummationByParts.calcMatrixEigs_rev!(A, λ, λ_bar, A_bar)
        dfdx = zeros(($T), n)
        for k = 1:n
          for i = 1:n
            for j = 1:n
              dfdx[k] += A_bar[i,j]*Q[i,k]*Q[j,k]
            end
          end
        end
        @test ≈(dfdx, x, rtol=1e-13)
      end
    end
  end

  @testset "Testing calcMatrixEigs! and calcMatrixEigs_rev! (skewsymmetric matrix" begin
    # construct a skew symmetric matrix whose eigenvalues are the parameters
    n = 10
    A = rand(n,n)
    A = 0.5*(A - A')
    eigvecval = eigen(A)
    x = eigvecval.values
    V = eigvecval.vectors
    idx = sortperm(x, lt=SummationByParts.compareEigs)
    x[:] = x[idx]
    Q = zeros(eltype(V), size(V))
    Q[:,:] = V[:,idx]
    @test ≈(real(Q*diagm(x)*Q'), A, rtol=1e-14)
    
    λ = zeros(ComplexF64, n)
    Ac = complex(A)
    SummationByParts.calcMatrixEigs!(Ac, λ)
    #println(λ)
    #println(x)
    #println(abs(λ - x))
    @test ≈(λ, x, rtol=1e-12)
    
    λ_bar = deepcopy(λ)
    A_bar = zeros(ComplexF64, size(Ac))
    SummationByParts.calcMatrixEigs_rev!(Ac, λ, λ_bar, A_bar)
    #println("norm(λ_bar) = ",norm(λ_bar))
    #println("norm(A_bar) = ",norm(A_bar))
    #println("norm(Q) = ",norm(Q))
    dfdx = zeros(ComplexF64, n)
    for k = 1:n
      for i = 1:n
        for j = 1:n
          dfdx[k] += A_bar[i,j]*conj(Q[i,k])*Q[j,k]
        end
      end
    end
    @test ≈(dfdx, x, rtol=1e-13)
  end

  @testset "Testing eigenvalueObj and eigenvalueObjGrad!" begin
    # The full matrix is the set of design variables here (Z = I, yperp = 0),
    # and we use a finite-difference approximation to test the gradient
    # (complex-step is not an option due to the complex arithmetic)

    numnodes = 10
    n = div(numnodes*(numnodes-1),2)
    p = 10
    x = rand(n)
    xperp = zeros(size(x))
    Znull = Matrix{Float64}(I, n, n)
    w = rand(numnodes)
    E = rand(numnodes,numnodes)
    E = 0.5*(E + E')
    obj = SummationByParts.eigenvalueObj(x, p, xperp, Znull, w, E)
    grad = zeros(n)
    SummationByParts.eigenvalueObjGrad!(x, p, xperp, Znull, w, E, grad)
    
    # find the finite-difference gradient and compare it to the reverse-mode
    # gradient

    epsfd = 1e-6
    for i = 1:n
      x[i] += epsfd
      grad_fd = SummationByParts.eigenvalueObj(x, p, xperp, Znull, w, E)
      grad_fd = (grad_fd - obj)/epsfd
      @test ≈(grad_fd, grad[i], rtol=1e-3)
      x[i] -= epsfd
    end
  end

  @testset "Testing conditionObj and conditionObjGrad!" begin
    # The full matrix is the set of design variables here (Z = I, yperp = 0),
    # and we use a finite-difference approximation to test the gradient
    # (complex-step is not an option due to the complex arithmetic)

    numnodes = 10
    n = div(numnodes*(numnodes-1),2)
    p = 1
    x = rand(n)
    xperp = zeros(size(x))
    Znull = Matrix{Float64}(I, n, n)
    w = rand(numnodes)
    E = rand(numnodes,numnodes)
    E = 0.5*(E + E')
    obj = SummationByParts.conditionObj(x, p, xperp, Znull, E)
    grad = zeros(n)
    SummationByParts.conditionObjGrad!(x, p, xperp, Znull, E, grad)
    
    # find the finite-difference gradient and compare it to the reverse-mode
    # gradient
    epsfd = 1e-6
    for i = 1:n
      x[i] += epsfd
      grad_fd = SummationByParts.conditionObj(x, p, xperp, Znull, E)
      grad_fd = (grad_fd - obj)/epsfd
      @test ≈(grad_fd, grad[i], rtol=1e-1)
      x[i] -= epsfd
    end
  end

end
