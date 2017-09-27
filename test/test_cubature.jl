facts("Testing Cubature Module...") do

  context("Testing Cubature.cubatureresidual (TriSymCub method)") do
    # Check that a P1 vertices-only rule produces zero residual
    vtx = Float64[-1 -1; 1 -1; -1 1]
    cub = SymCubatures.TriSymCub{Float64}()
    #x,y = SymCubatures.calcnodes(cub, vtx)
    SymCubatures.setweights!(cub, [2/3])
    F, dF = Cubature.cubatureresidual(cub, 1)
    @fact F --> roughly(zeros(3), atol=1e-15)
  end

  context("Testing Cubature.cubatureresidual (TetSymCub method)") do
    # Check that a P1 vertices-only rule produces zero residual
    vtx = Float64[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
    cub = SymCubatures.TetSymCub{Float64}()
    #x,y,z = SymCubatures.calcnodes(cub, vtx)
    SymCubatures.setweights!(cub, [1/3])
    F, dF = Cubature.cubatureresidual(cub, 1)
    @fact F --> roughly(zeros(4), atol=1e-15)
  end

  context("Testing Cubature.solvecubature (TriSymCub method)") do
    # recover cubature for centroid only rule
    cub = SymCubatures.TriSymCub{Float64}(vertices=false, centroid=true)
    SymCubatures.setweights!(cub, [0.5])
    mask = 1:(cub.numparams+cub.numweights)
    Cubature.solvecubature!(cub, 1, mask, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    @fact w --> roughly([2], atol=1e-15)

    # recover P1 vertices-only rule
    cub = SymCubatures.TriSymCub{Float64}()
    SymCubatures.setweights!(cub, [0.5])
    mask = 1:(cub.numparams+cub.numweights)
    Cubature.solvecubature!(cub, 1, mask, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    @fact w --> roughly([2/3, 2/3, 2/3], atol=1e-15)

    # recover cubature for centroid+S21 (exact to degree 3)
    cub = SymCubatures.TriSymCub{Float64}(vertices=false, centroid=true, numS21=1)
    SymCubatures.setweights!(cub, [0.5, 0.5])
    SymCubatures.setparams!(cub, [0.25])
    mask = 1:(cub.numparams+cub.numweights)
    Cubature.solvecubature!(cub, 3, mask, tol=1e-14)
    w  = SymCubatures.calcweights(cub)
    @fact w --> roughly([50/48, 50/48, 50/48, -18/16], atol=1e-15)

    # create P2 element cubature with 1 bubble node at the centroid
    cub = SymCubatures.TriSymCub{Float64}(midedges=true, centroid=true)
    SymCubatures.setweights!(cub, [0.5, 0.5, 0.5])
    mask = 1:(cub.numparams+cub.numweights)
    Cubature.solvecubature!(cub, 3, mask, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    @fact w --> roughly([1/10, 1/10, 1/10, 4/15, 4/15, 4/15, 18/20], atol=1e-15)

    # create P3 element cubature with 3 bubble nodes
    cub = SymCubatures.TriSymCub{Float64}(numedge=1, numS21=1)
    SymCubatures.setweights!(cub, [0.5, 0.5, 0.5])
    SymCubatures.setparams!(cub, [0.25, 0.25])
    mask = 1:(cub.numparams+cub.numweights)
    Cubature.solvecubature!(cub, 5, mask, tol=1e-14)
    w  = SymCubatures.calcweights(cub)
    @fact w --> roughly([0.029745826049641155,0.029745826049641155,0.029745826049641155,0.44155411568082154,0.44155411568082154,0.44155411568082154,0.097683362468102,0.097683362468102,0.097683362468102,0.097683362468102,0.097683362468102,0.097683362468102], atol=1e-15)
  end

  context("Testing Cubature.solvecubature (TetSymCub method)") do
    # recover P1 vertices-only rule
    cub = SymCubatures.TetSymCub{Float64}()
    SymCubatures.setweights!(cub, [0.5])
    mask = 1:(cub.numparams+cub.numweights)
    Cubature.solvecubature!(cub, 1, mask, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    @fact w --> roughly([1/3, 1/3, 1/3, 1/3], atol=1e-15)

    # recover cubature for vertex+centroid (exact to degree 2)
    cub = SymCubatures.TetSymCub{Float64}(vertices=true, centroid=true)
    SymCubatures.setweights!(cub, [0.1 0.1])
    mask = 1:(cub.numparams+cub.numweights)
    Cubature.solvecubature!(cub, 2, mask, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    @fact w --> roughly([4/60, 4/60, 4/60, 4/60, 16/15], atol=1e-14)

    # recover cubature for S31+centroid (exact to degree 3)
    cub = SymCubatures.TetSymCub{Float64}(vertices=false, centroid=true,
                                          numS31=1)
    SymCubatures.setweights!(cub, [0.1 0.1])
    SymCubatures.setparams!(cub, [0.1])
    mask = 1:(cub.numparams+cub.numweights)
    Cubature.solvecubature!(cub, 3, mask, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    @fact w --> roughly([3/5, 3/5, 3/5, 3/5, -16/15], atol=1e-14)

    # create a P2 element with 1 bubble node
    cub = SymCubatures.TetSymCub{Float64}(midedges=true, centroid=true) #numS31=1)
    wuni = 0.1 #(4.0/3.0)/orbits.numnodes
    SymCubatures.setweights!(cub, [wuni wuni wuni]) #[0.02 0.05 0.23])
    mask = 1:(cub.numparams+cub.numweights)
    Cubature.solvecubature!(cub, 3, mask, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    @fact w --> roughly([1/45.*ones(4); 4/45.*ones(6); 32/45], atol=1e-14)

    # create a P3 element with 4 bubble nodes
    cub = SymCubatures.TetSymCub{Float64}(facecentroid=true,
                                          numedge=1, numS31=1)
    SymCubatures.setweights!(cub, [0.1 0.1 0.1 0.1])
    SymCubatures.setparams!(cub, [1/6 1/5])
    mask = 1:(cub.numparams+cub.numweights)
    Cubature.solvecubature!(cub, 5, mask, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    @fact w -->
    roughly([0.004421633248304814,0.004421633248304814,0.004421633248304814,0.004421633248304814,0.06935370366814599,0.06935370366814599,0.06935370366814599,0.06935370366814599,0.2065316361160523,0.2065316361160523,0.2065316361160523,0.2065316361160523,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603], atol=1e-14)

    # create a P4 element with 11 bubble nodes
    cub = SymCubatures.TetSymCub{Float64}(midedges=true, centroid=true, numedge=1,
                                          numfaceS21=1, numS31=1, numS22=1)
    SymCubatures.setweights!(cub, [0.001 0.004 0.005 0.02 0.08 0.06 0.1])
    SymCubatures.setparams!(cub, [0.28 0.22 0.75 0.45])
    mask = 1:(cub.numparams+cub.numweights)
    Cubature.solvecubature!(cub, 7, mask, tol=1e-14)
    w = SymCubatures.calcweights(cub)
    @fact w -->
    roughly([0.0015106273303336273,0.0015106273303336273,0.0015106273303336273,0.0015106273303336273,0.060490542374353584,0.060490542374353584,0.060490542374353584,0.060490542374353584,0.004038881996228382,0.004038881996228382,0.004038881996228382,0.004038881996228382,0.004038881996228382,0.004038881996228382,0.10344930834722398,0.10344930834722398,0.10344930834722398,0.10344930834722398,0.10344930834722398,0.10344930834722398,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.08113091859465722], atol=1e-14)
  end

  context("Testing Cubature.quadrature (internal=false)") do
    # verify by checking that polynomials are integrated exactly
    for q = 1:10
      quad, vtx = quadrature(q, Float64)
      x = SymCubatures.calcnodes(quad, vtx)
      w = SymCubatures.calcweights(quad)
      f = 0.5*dot(w, (0.5.*vec(x) + 0.5).^q)
      @fact f --> roughly(1/(q+1), atol=1e-15)      
    end
  end 

  context("Testing Cubature.quadrature (internal=true)") do
    # verify by checking that polynomials are integrated exactly
    for q = 1:10
      quad, vtx = quadrature(q, Float64, internal=true)
      x = SymCubatures.calcnodes(quad, vtx)
      w = SymCubatures.calcweights(quad)
      f = 0.5*dot(w, (0.5.*vec(x) + 0.5).^q)
      @fact f --> roughly(1/(q+1), atol=1e-15)
    end
  end

  context("Testing Cubature.getTriCubatureGamma") do
    # test using Float32, because this has not been done above
    cub, vtx = getTriCubatureGamma(1, Float64)
    w = SymCubatures.calcweights(cub)
    @fact w --> roughly(Float32[2/3, 2/3, 2/3], atol=1e-7)
    @fact SymCubatures.getnumboundarynodes(cub) --> 3
    cub, vtx = getTriCubatureGamma(3, Float64)
    w = SymCubatures.calcweights(cub)
    @fact w --> roughly(Float32[1/10, 1/10, 1/10, 4/15, 4/15, 4/15, 18/20],
                       atol=1e-7)
    @fact SymCubatures.getnumboundarynodes(cub) --> 6
    cub, vtx = getTriCubatureGamma(5, Float64)
    w = SymCubatures.calcweights(cub)
    @fact w -->
    roughly(Float32[0.029745826049641155,0.029745826049641155,0.029745826049641155,0.44155411568082154,0.44155411568082154,0.44155411568082154,0.097683362468102,0.097683362468102,0.097683362468102,0.097683362468102,0.097683362468102,0.097683362468102],
    atol=1e-7)
    @fact SymCubatures.getnumboundarynodes(cub) --> 9
  end

  context("Testing Cubature.getTriCubatureOmega") do
    cub, vtx = getTriCubatureOmega(2, Float64)
    @fact cub.weights --> roughly([2/3], atol=1e-14)
    @fact cub.params --> roughly([1/3], atol=1e-14)

    cub, vtx = getTriCubatureOmega(4, Float64)
    @fact cub.weights --> roughly([0.44676317935602283;
                                   0.2199034873106437], atol=1e-14)
    @fact cub.params --> roughly([0.8918969818319298;
                                  0.18315242701954149], atol=1e-14)

    cub, vtx = getTriCubatureOmega(5, Float64)
    @fact cub.weights --> roughly([0.11550472674301035;
                                   0.20924480696331949;
                                   0.39801697799105223], atol=1e-14)
    @fact cub.params --> roughly([0.13862330627662678;
                                  0.14215944055500324;
                                  0.6226442585632832], atol=1e-14)

    cub, vtx = getTriCubatureOmega(7, Float64)
    @fact cub.weights --> roughly([0.045386157905236965;
                                   0.1458284149509071;
                                   0.2543369199180239;
                                   0.11055758694624956], atol=1e-14)
    @fact cub.params --> roughly([0.08433122881886425;
                                  0.9485893782350207;
                                  0.4841719475189572;
                                  0.4231241172761849;
                                  0.0959626827429292], atol=1e-14)
  end

  context("Testing Cubature.getTetCubatureGamma") do
    # test using Float32, because this has not been done above
    cub, vtx = getTetCubatureGamma(1, Float32)
    w = SymCubatures.calcweights(cub)
    @fact w --> roughly(Float32[1/3, 1/3, 1/3, 1/3], atol=1e-7)
    @fact SymCubatures.getnumboundarynodes(cub) --> 4
    cub, vtx = getTetCubatureGamma(3, Float32)
    w = SymCubatures.calcweights(cub)
    @fact w --> roughly([Float32(2/90).*ones(Float32, 4);
                         Float32(8/90).*ones(Float32, 6);
                         Float32[32/45]], atol=1e-7)
    @fact SymCubatures.getnumboundarynodes(cub) --> 10
    cub, vtx = getTetCubatureGamma(5, Float64)
    w = SymCubatures.calcweights(cub)
    @fact w -->
    roughly([0.004421633248304814,0.004421633248304814,0.004421633248304814,0.004421633248304814,0.06935370366814599,0.06935370366814599,0.06935370366814599,0.06935370366814599,0.2065316361160523,0.2065316361160523,0.2065316361160523,0.2065316361160523,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603,0.017675453433610603], atol=1e-14)
    @fact SymCubatures.getnumboundarynodes(cub) --> 20
    cub, vtx = getTetCubatureGamma(7, Float64)
    w = SymCubatures.calcweights(cub)
    @fact w -->
    roughly([0.0015106273303336273,0.0015106273303336273,0.0015106273303336273,0.0015106273303336273,0.060490542374353584,0.060490542374353584,0.060490542374353584,0.060490542374353584,0.004038881996228382,0.004038881996228382,0.004038881996228382,0.004038881996228382,0.004038881996228382,0.004038881996228382,0.10344930834722398,0.10344930834722398,0.10344930834722398,0.10344930834722398,0.10344930834722398,0.10344930834722398,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.005696088152131421,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.02424296133613638,0.08113091859465722], atol=1e-14)
    @fact SymCubatures.getnumboundarynodes(cub) --> 34    
  end

  context("Testing Cubature.getTetCubatureOmega") do
    cub, vtx = getTetCubatureOmega(2, Float64)
    @fact cub.weights --> roughly([1/3], atol=1e-14)
    @fact cub.params --> roughly([(1 - sqrt(5)/5)*3/4], atol=1e-14)

    cub, vtx = getTetCubatureOmega(3, Float64)
    @fact cub.weights --> roughly([0.06483158243276162;
                                   0.17900116726703835], atol=1e-14)
    @fact cub.params --> roughly([0.22511815489558668;
                                  0.18771315212883505], atol=1e-14)
  end
  
  context("Testing Cubature.equivalenceconstant (tricubature method)") do
    λ = [4.0; 6.165789254884331; 7.136363791526314; 8.644675072049024]
    for p = 1:4
      cub, vtx = getTriCubatureGamma(2*p-1, Float64)
      @fact Cubature.equivalenceconstant(cub, vtx, 2*p-1) -->
      roughly(λ[p], atol=1e-13)
    end
  end

  context("Testing Cubature.equivalenceconstant (tetcubature method)") do
    λ = [5.0; 10.5; 11.77599246416848; 16.61463996715285]
    for p = 1:4
      cub, vtx = getTetCubatureGamma(2*p-1, Float64)
      @fact Cubature.equivalenceconstant(cub, vtx, 2*p-1) -->
      roughly(λ[p], atol=1e-13)
    end
  end

end
