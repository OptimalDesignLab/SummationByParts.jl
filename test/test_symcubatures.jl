@testset "Testing SymCubatures Module..." begin

  for T = (Float32, Float64, ComplexF64 )
    @eval begin
      @testset "Testing PointSymCub Inner Constructor for DataType $(string($T))" begin
        quad = PointSymCub{($T)}()
        @test quad.numparams == 0
        @test quad.numweights == 1
        @test quad.numnodes == 1
      end
    end
  end
  
  for T = (Float32, Float64, ComplexF64 )
    @eval begin
      @testset "Testing LineSymCub Inner Constructor for DataType $(string($T))" begin
        quad = LineSymCub{($T)}(vertices=false)
        @test quad.numparams == 0
        @test quad.numweights == 0
        @test quad.numnodes  == 0
        quad = LineSymCub{($T)}(numedge=2)
        @test quad.numparams  == 2
        @test quad.numweights  == 3
        @test quad.numnodes  == 2+2*2
        quad = LineSymCub{($T)}(numedge=3, centroid=true)
        @test quad.numparams  == 3
        @test quad.numweights  == 5
        @test quad.numnodes  == 2+2*3+1
      end
    end
  end
  
  for T = (Float32, Float64, ComplexF64 )
    @eval begin
      @testset "Testing TriSymCub Inner Constructor for DataType $(string($T))" begin
        tricub = TriSymCub{($T)}(vertices=false)
        @test tricub.numparams == 0
        @test tricub.numweights == 0
        @test tricub.numnodes == 0
        tricub = TriSymCub{($T)}(numedge=2,midedges=true)
        @test tricub.numparams == 2
        @test tricub.numweights == 4
        @test tricub.numnodes == 3+3+2*6
        tricub = TriSymCub{($T)}(midedges=true, numS21=3)
        @test tricub.numparams == 3
        @test tricub.numweights == 5
        @test tricub.numnodes == 3+3+3*3
        tricub = TriSymCub{($T)}(numedge=2, centroid=true, numS111=2)
        @test tricub.numparams == 2+4
        @test tricub.numweights == 6
        @test tricub.numnodes == 3+2*6+1+2*6
      end
    end
  end
  
  for T = (Float32, Float64, ComplexF64 )
    @eval begin
      @testset "Testing TetSymCub Inner Constructor for DataType $(string($T))" begin
        tetcub = TetSymCub{($T)}(vertices=false)
        @test tetcub.numparams == 0
        @test tetcub.numweights == 0
        @test tetcub.numnodes == 0
        tetcub = TetSymCub{($T)}(numedge=2,midedges=true)
        @test tetcub.numparams == 2
        @test tetcub.numweights == 4
        @test tetcub.numnodes == 4+6+2*12
        tetcub = TetSymCub{($T)}(midedges=true, numS31=3)
        @test tetcub.numparams == 3
        @test tetcub.numweights == 5
        @test tetcub.numnodes == 4+6+4*3
        tetcub = TetSymCub{($T)}(midedges=true, numedge=1, numfaceS21=2)
        @test tetcub.numparams == 3
        @test tetcub.numweights == 5
        @test tetcub.numnodes == 4+6+1*12+2*12
        tetcub = TetSymCub{($T)}(midedges=true, centroid=true, numS22=2)
        @test tetcub.numparams == 2
        @test tetcub.numweights == 5
        @test tetcub.numnodes == 4+6+1+2*6
        tetcub = TetSymCub{($T)}(vertices=false, numfaceS111=2, numS211=1,
                                 numS1111=2)
        @test tetcub.numparams == 2*2 + 2 + 2*3
        @test tetcub.numweights == 5
        @test tetcub.numnodes == 2*6*4 + 12 + 2*24
      end
    end
  end

  @testset "Testing getnumboundarynodes (PointSymCub method)" begin
    point = PointSymCub{Float64}()
    @test SymCubatures.getnumboundarynodes(point) == 1
  end
  
  @testset "Testing getnumboundarynodes (LineSymCub method)" begin
    quad = LineSymCub{Float64}() # vertex only rule
    @test SymCubatures.getnumboundarynodes(quad) == 2
    quad = LineSymCub{Float64}(numedge = 2, centroid=true, vertices=false)
    @test SymCubatures.getnumboundarynodes(quad) == 0
  end

  @testset "Testing getnumboundarynodes (TriSymCub method)" begin
    tricub = TriSymCub{Float64}() # vertex only rule
    @test SymCubatures.getnumboundarynodes(tricub) == 3
    tricub = TriSymCub{Float64}(numedge = 2, midedges=true, numS21 = 4)
    @test SymCubatures.getnumboundarynodes(tricub) == 3+3+2*6
  end

  @testset "Testing getnumboundarynodes (TetSymCub method)" begin
    tetcub = TetSymCub{Float64}() # vertex only rule
    @test SymCubatures.getnumboundarynodes(tetcub) == 4
    tetcub = TetSymCub{Float64}(numedge=2, midedges=true, numfaceS21=3,
                                numS31 = 4)
    @test SymCubatures.getnumboundarynodes(tetcub) == 4+6+2*12+3*12
    tetcub = TetSymCub{Float64}(numfaceS111=2, numS1111=1)
    @test SymCubatures.getnumboundarynodes(tetcub) == 4+2*4*6
  end

  @testset "Testing getnumfacenodes (PointSymCub method)" begin
    point = PointSymCub{Float64}()
    @test SymCubatures.getnumfacenodes(point) == 1
  end
  
  @testset "Testing getnumfacenodes (LineSymCub method)" begin
    quad = LineSymCub{Float64}() # vertex only rule
    @test SymCubatures.getnumfacenodes(quad) == 1
    quad = LineSymCub{Float64}(numedge = 2, centroid=true, vertices=false)
    @test SymCubatures.getnumfacenodes(quad) == 0
  end

  @testset "Testing getnumfacenodes (TriSymCub method)" begin
    tricub = TriSymCub{Float64}() # vertex only rule
    @test SymCubatures.getnumfacenodes(tricub) == 2
    tricub = TriSymCub{Float64}(numedge = 2, midedges=true, numS21 = 4)
    @test SymCubatures.getnumfacenodes(tricub) == 2+2*2+1
  end

  @testset "Testing getnumfacenodes (TetSymCub method)" begin
    tetcub = TetSymCub{Float64}() # vertex only rule
    @test SymCubatures.getnumfacenodes(tetcub) == 3
    tetcub = TetSymCub{Float64}(numedge=2, midedges=true, numfaceS21=3,
                                numS31 = 4, numfaceS111=1)
    @test SymCubatures.getnumfacenodes(tetcub) == 3+3+2*6+3*3+6
  end

  @testset "Testing getbndrynodeindices (PointSymCub method)" begin
    point = PointSymCub{Float64}()
    bndryindices = SymCubatures.getbndrynodeindices(point)
    @test bndryindices == [1;]
  end

  @testset "Testing getbndrynodeindices (LineSymCub method)" begin
    quad = LineSymCub{Float64}(numedge=1, centroid=true)
    bndryindices = SymCubatures.getbndrynodeindices(quad)
    @test bndryindices == [1;2]
  end

  @testset "Testing getbndrynodeindices (TriSymCub method)" begin
    tricub = TriSymCub{Float64}(numedge=1, midedges=true, numS21 = 4)
    bndryindices = SymCubatures.getbndrynodeindices(tricub)
    @test bndryindices == [1; 2; 3; 4; 5; 6; 19; 20; 21; 22; 23; 24]
  end

  @testset "Testing getbndrynodeindices (TetSymCub method)" begin
    tetcub = TetSymCub{Float64}(numedge=2, facecentroid=true, numS31=1, 
                                midedges=true, numS22=1, numfaceS21=1,
                                numfaceS111=1)
    bndryindices = SymCubatures.getbndrynodeindices(tetcub)
    # @test bndryindices == [1; 2; 3; 4; 
    #                         5; 6; 7; 8;
    #                         13; 14; 15; 16; 17; 18;
    #                         25; 26; 27; 28; 29; 30; 31; 32; 33; 34; 35; 36;
    #                         37; 38; 39; 40; 41; 42; 43; 44; 45; 46; 47; 48;
    #                         49; 50; 51; 52; 53; 54; 55; 56; 57; 58; 59; 60;
    #                         61; 62; 63; 64; 65; 66; 67; 68; 69; 70; 71; 72;
    #                         73; 74; 75; 76; 77; 78; 79; 80; 81; 82; 83; 84]
    @test bndryindices == [1:4;
                          9:14;
                          21:32;
                          33:44;
                          45:56;
                          57:80;
                          81:84]
  end

  @testset "Testing getinteriornodeindices (PointSymCub method)" begin
    point = PointSymCub{Float64}()
    indices = SymCubatures.getinteriornodeindices(point)
    @test indices == Int[]
  end

  @testset "Testing getinteriornodeindices (LineSymCub method)" begin
    quad = LineSymCub{Float64}(numedge=1, centroid=true)
    indices = SymCubatures.getinteriornodeindices(quad)
    @test indices == [3;4;5]
  end

  @testset "Testing getinteriornodeindices (TriSymCub method)" begin
    tricub = TriSymCub{Float64}(numedge=1, midedges=true, numS21=4,
                                centroid=true)
    indices = SymCubatures.getinteriornodeindices(tricub)
    @test indices == [7:18; 25]
  end

  # @testset "Testing getinteriornodeindices (TetSymCub method)" begin
  #   tetcub = TetSymCub{Float64}(facecentroid=true, numS31=1, 
  #                               midedges=true, numS22=1,
  #                               numedge=2, numfaceS21=1, numS211=2,
  #                               numfaceS111=1, numS1111=1,
  #                               centroid=true)
  #   indices = SymCubatures.getinteriornodeindices(tetcub)
  #   @test indices == [9:12; 19:24; 61:84; 109:132; 133]
  # end

  @testset "Testing getfacevertexindices (PointSymCub method)" begin
    point = PointSymCub{Float64}()
    facevtx = SymCubatures.getfacevertexindices(point)
    @test facevtx[1,1] == 1
  end

  @testset "Testing getfacevertexindices (LineSymCub method)" begin
    quad = LineSymCub{Float64}(numedge=1, centroid=true)
    facevtx = SymCubatures.getfacevertexindices(quad)
    @test facevtx[1,1] == 1
    @test facevtx[1,2] == 2
  end

  @testset "Testing getfacevertexindices (TriSymCub method)" begin
    tricub = TriSymCub{Float64}(numedge=1, midedges=true, numS21 = 4)
    facevtx = SymCubatures.getfacevertexindices(tricub)
    @test facevtx[:,1] == [1; 2]
    @test facevtx[:,2] == [2; 3]
    @test facevtx[:,3] == [3; 1]
  end

  @testset "Testing getfacevertexindices (TetSymCub method)" begin
    tetcub = TetSymCub{Float64}(numedge=1, facecentroid=true, numS31=1, 
                                midedges=true, numS22=1, numfaceS21=1,
                                centroid=true)
    facevtx = SymCubatures.getfacevertexindices(tetcub)
    @test facevtx[:,1] == [1; 2; 3]
    @test facevtx[:,2] == [1; 4; 2]
    @test facevtx[:,3] == [2; 4; 3]
    @test facevtx[:,4] == [1; 3; 4]
  end

  @testset "Testing getfacenodeindices (PointSymCub method)" begin
    point = PointSymCub{Float64}()
    bndryindices = SymCubatures.getfacenodeindices(point)
    @test bndryindices == reshape([1;],(1,1))
  end

  @testset "Testing getfacenodeindices (LineSymCub method)" begin
    quad = LineSymCub{Float64}(numedge=1, centroid=true)
    bndryindices = SymCubatures.getfacenodeindices(quad)
    @test bndryindices == [1 2]
  end

  @testset "Testing getfacenodeindices (TriSymCub method)" begin
    tricub = TriSymCub{Float64}(numedge=1, midedges=true, numS21 = 4)
    bndryindices = SymCubatures.getfacenodeindices(tricub)
    @test bndryindices == [1 2 3;
                            2 3 1;
                            4 5 6;
                            19 21 23;
                            20 22 24]
  end

  @testset "Testing getfacenodeindices (TetSymCub method)" begin
    tetcub = TetSymCub{Float64}(facecentroid=true, numS31=1,
                                midedges=true, numS22=1,
                                numedge=1, numfaceS21=1, numS211=2,
                                numfaceS111=1, numS1111=2,
                                centroid=true)
    bndryindices = SymCubatures.getfacenodeindices(tetcub)
    # @test bndryindices == [1 1 2 1; # vertices
    #                         2 4 4 3; # ...
    #                         3 2 3 4; # ...
    #                         5 6 7 8; # face centroids
    #                         13 16 17 15; 14 17 18 18; 15 13 14 16; # midedges
    #                         25 31 33 30; 26 32 34 29; 27 34 36 35; # edge nodes
    #                         28 33 35 36; 29 26 28 32; 30 25 27 31; # ...
    #                         37 40 43 46; 38 41 44 47; 39 42 45 48; # face S21
    #                         73 79 85 91; 74 80 86 92; 75 81 87 93; # face S111
    #                         76 82 88 94; 77 83 89 95; 78 84 90 96] # ...
    @test bndryindices == [1 1 2 1; # vertices
                            2 4 4 3; # ...
                            3 2 3 4; # ...
                            9 12 13 11; 10 13 14 14; 11 9 10 12; #midedges
                            21 24 27 30; 22 25 28 31; 23 26 29 32; #face S21
                            33 39 41 38; 34 40 42 37; 35 42 44 43; #edge
                            36 41 43 44; 37 34 36 40; 38 33 35 39; #...
                            69 75 81 87; 70 76 82 88; 71 77 83 89; #faceS111
                            72 78 84 90; 73 79 85 91; 74 80 86 92; #...
                            93 94 95 96] # facecentroid
  end

  @testset "Testing findleftperm!" begin
    # This uses the matrix for computing the Tetrahedron's S22 orbit to test the
    # findleftperm! method
    alpha = pi
    A = [alpha alpha (0.5-alpha) (0.5-alpha);
         alpha (0.5-alpha) alpha (0.5-alpha);
         alpha (0.5-alpha) (0.5-alpha) alpha;
         (0.5-alpha) alpha alpha (0.5-alpha);
         (0.5-alpha) alpha (0.5-alpha) alpha;
         (0.5-alpha) (0.5-alpha) alpha alpha]
    n = size(A,1)
    m = size(A,2)
    permR = shuffle([1:m;])
    perm = zeros(Int, (n))

    # First, find a valid left permutation, and check it
    success = SymCubatures.findleftperm!(A, permR, perm)
    @test success == true
    @test A[perm,:] == A[:,permR]

    # Next, transpose the matrix, and check that we recover the right permutation
    permL = deepcopy(perm)
    resize!(perm, m)
    success = SymCubatures.findleftperm!(A', permL, perm)
    @test success == true
    @test A[:,perm] == A[permL,:]
    @test perm == permR

    # Finally, check for failure in a 3x2 case that has no valid left permutation
    A = [1.0 2.0; 3.0 4.0; 5.0 6.0]
    n = size(A,1)
    m = size(A,2)
    resize!(permR, m)
    resize!(perm, n)
    permR = shuffle([1:m;])
    while permR == [1:m;]
      permR = shuffle([1:m;])
    end
    success = SymCubatures.findleftperm!(A, permR, perm)
    @test success == false
  end

  @testset "Test getpermutation (PointSymCub method)" begin
    point = PointSymCub{Float64}()
    SymCubatures.setparams!(point, rand(point.numparams))
    vtxperm = [1;]
    perm = SymCubatures.getpermutation(point, vtxperm)
    vtx = Float64[1 1;]
    x = SymCubatures.calcnodes(point, vtx)
    vtx = Float64[1 1;]
    xr = SymCubatures.calcnodes(point, vtx)
    @test ≈(x[:,perm], xr; atol=eps())
  end

  @testset "Test getpermutation (LineSymCub method)" begin
    quad = LineSymCub{Float64}(vertices=true, numedge=2, centroid=true)
    SymCubatures.setparams!(quad, rand(quad.numparams))
    vtxperm = [2; 1]
    perm = SymCubatures.getpermutation(quad, vtxperm)
    vtx = Float64[0 0; 1 1]
    x = SymCubatures.calcnodes(quad, vtx)
    vtx = Float64[1 1; 0 0]
    xr = SymCubatures.calcnodes(quad, vtx)
    @test ≈(x[:,perm], xr, atol=eps())
  end

  @testset "Test getpermutation (TriSymCub method)" begin
    tricub = TriSymCub{Float64}(vertices=true, numedge=2, midedges=true,
                                numS21 = 3, numS111=4, centroid=true)
    SymCubatures.setparams!(tricub, rand(tricub.numparams))
    vtxperm = [3; 2; 1]
    perm = SymCubatures.getpermutation(tricub, vtxperm)
    vtx = Float64[0 0 0; 1 0 0; 0.5 0.5 0.5]
    x = SymCubatures.calcnodes(tricub, vtx)
    vtx = Float64[0.5 0.5 0.5; 1 0 0; 0 0 0]
    xr = SymCubatures.calcnodes(tricub, vtx)
    @test ≈(x[:,perm], xr, atol=eps())
  end

  @testset "Test getpermutation (TetSymCub method)" begin
    tetcub = TetSymCub{Float64}(vertices=true, midedges=true, centroid=true,
                                facecentroid=true, numedge=2, numfaceS21=2,
                                numS31=2, numS22=2, numfaceS111=1, numS211=2,
                                numS1111=1)
    SymCubatures.setparams!(tetcub, rand(tetcub.numparams))
    vtxperm = [3; 4; 1; 2]
    perm = SymCubatures.getpermutation(tetcub, vtxperm)
    vtx = Float64[0 0 0; 1 0 0; 0 1 0; 0 0 1]
    x = SymCubatures.calcnodes(tetcub, vtx)   
    vtx = Float64[0 1 0; 0 0 1; 0 0 0; 1 0 0]
    xr = SymCubatures.calcnodes(tetcub, vtx)
    @test ≈(x[:,perm], xr, atol=eps())
  end  

  @testset "Test getfacebasedpermutation (PointSymCub method)" begin
    point = PointSymCub{Float64}()
    perm = SymCubatures.getfacebasedpermutation(point)
    @test perm[:,1]==[1;]
  end

  @testset "Test getfacebasedpermutation (LineSymCub method)" begin
    quad = LineSymCub{Float64}(numedge=1, centroid=true)
    perm = SymCubatures.getfacebasedpermutation(quad)
    @test perm[:,1] == [1:5;]
    @test perm[:,2] == [2; 1; 4; 3; 5]
    perm = SymCubatures.getfacebasedpermutation(quad, faceonly=true)
    @test perm[:,1] == [1;]
    @test perm[:,2] == [2;]
  end

  @testset "Test getfacebasedpermutation (TriSymCub method)" begin
    tricub = TriSymCub{Float64}(numedge=1, midedges=true, numS21 = 1, numS111=1,
                                centroid=true)
    perm = SymCubatures.getfacebasedpermutation(tricub)
    @test perm[:,1] == [1:22;]
    @test perm[:,2] == [2; 3; 1; 5; 6; 4; 8; 9; 7;
                         12; 13; 14; 15; 10; 11; 
                         18; 19; 20; 21; 16; 17; 
                         22]
    @test perm[:,3] == [3; 1; 2; 6; 4; 5; 9; 7; 8;
                         14; 15; 10; 11; 12; 13; 
                         20; 21; 16; 17; 18; 19;
                         22]
    perm = SymCubatures.getfacebasedpermutation(tricub, faceonly=true)
    @test perm[:,1] == [1; 2; 4; 10; 11]
    @test perm[:,2] == [2; 3; 5; 12; 13]
    @test perm[:,3] == [3; 1; 6; 14; 15]
  end

  @testset "Test getfacebasedpermutation (TetSymCub method)" begin
    # This test works by providing different sets of vertices to calcnodes, and
    # then checking for equivalence with a base set of vertices reordered using
    # the face-based permutation
    tetcub = TetSymCub{Float64}(vertices=true, midedges=true, centroid=true,
                                facecentroid=true, numedge=2, numfaceS21=2,
                                numS31=2, numS22=2, numS211=2, numfaceS111=2,
                                numS1111=1)
    SymCubatures.setparams!(tetcub, rand(tetcub.numparams))
    perm = SymCubatures.getfacebasedpermutation(tetcub)

    # get the nodes for each vertex arrangement
    vtx = Float64[0 0 0; 1 0 0; 0 1 0; 0 0 1]
    xf1 = SymCubatures.calcnodes(tetcub, vtx)   
    vtx = Float64[0 0 0; 0 0 1; 1 0 0; 0 1 0]
    xf2 = SymCubatures.calcnodes(tetcub, vtx)
    vtx = Float64[1 0 0; 0 0 1; 0 1 0; 0 0 0]
    xf3 = SymCubatures.calcnodes(tetcub, vtx)
    vtx = Float64[0 0 0; 0 1 0; 0 0 1; 1 0 0]
    xf4 = SymCubatures.calcnodes(tetcub, vtx)

    # now, reorder xf1 and check against the other vertex arrangments
    @test ≈(xf1[:,perm[:,1]], xf1, atol=eps())
    @test ≈(xf1[:,perm[:,2]], xf2, atol=eps())
    @test ≈(xf1[:,perm[:,3]], xf3, atol=eps())
    @test ≈(xf1[:,perm[:,4]], xf4, atol=eps())
  end

  @testset "Test getneighbourpermutation (PointSymCub method)" begin
    point = PointSymCub{Float64}()
    perm = SymCubatures.getneighbourpermutation(point)
    @test perm[1,1] == 1
  end

  @testset "Test getneighbourpermutation (LineSymCub method)" begin
    quad = LineSymCub{Float64}(numedge=1, centroid=true)
    perm = SymCubatures.getneighbourpermutation(quad)
    @test perm[:,1] == [2; 1; 4; 3; 5]
    quad = LineSymCub{Float64}(vertices=false, numedge=1, centroid=true)
    perm = SymCubatures.getneighbourpermutation(quad)
    @test perm[:,1] == [2; 1; 3]
  end

  @testset "Test getneighbourpermutation (TriSymCub method)" begin
    # This test constructs a face in 3D based on the vertices [0,0,0], [1,1,1],
    # and [-1,1,1].  It then computes the cubature nodes for this face using a
    # reference orientation, and three different neigbhour orientations.  The
    # permuted reference orientation should match the neighbour orientations.
    tricub = TriSymCub{Float64}(vertices=true, numedge=2, midedges=true,
                                numS21=3, numS111=4, centroid=true)
    SymCubatures.setparams!(tricub, rand(tricub.numparams))
    perm = SymCubatures.getneighbourpermutation(tricub)   
    vtx = Float64[0 0 0; 1 1 1; -1 1 1]
    x = SymCubatures.calcnodes(tricub, vtx)
    vtx = Float64[0 0 0; -1 1 1; 1 1 1]
    x1 = SymCubatures.calcnodes(tricub, vtx)
    vtx = Float64[1 1 1; 0 0 0; -1 1 1]
    x2 = SymCubatures.calcnodes(tricub, vtx)
    vtx = Float64[-1 1 1; 1 1 1; 0 0 0]
    x3 = SymCubatures.calcnodes(tricub, vtx)
    @test ≈(x[:,perm[:,1]], x1, atol=eps())
    @test ≈(x[:,perm[:,2]], x2, atol=eps())
    @test ≈(x[:,perm[:,3]], x3, atol=eps())
  end

  for T = (Float32, Float64, ComplexF64 )
    @eval begin
      @testset "Testing calcnodes (PointSymCub method) for DataType $(string($T))" begin
        vtx = ($T)[-1 1 0]
        point = PointSymCub{($T)}()
        @test SymCubatures.calcnodes(point, vtx) == vtx'
      end
    end
  end

  for T = (Float32, Float64, ComplexF64 )
    @eval begin
      @testset "Testing calcnodes (LineSymCub method) for DataType $(string($T))" begin
        vtx = reshape(($T)[-1; 1], (2,1))
        quad = LineSymCub{($T)}()
        @test SymCubatures.calcnodes(quad, vtx) == vtx[:,1]'

        quad = LineSymCub{($T)}(numedge=1)
        alpha = ($T)(1/pi)
        A = ($T)[alpha (1-alpha);
                 (1-alpha) alpha]
        SymCubatures.setparams!(quad, [alpha])
        x = SymCubatures.calcnodes(quad, vtx)
        @test x == [vtx[:,1]' (A*vtx[:,1])']

        # uniformly spaced 7 nodes between (-1,1) (no vertices)
        quad = LineSymCub{($T)}(numedge=3, centroid=true, vertices=false)
        SymCubatures.setparams!(quad, ($T)[5/8 6/8 7/8])
        x = SymCubatures.calcnodes(quad, vtx)
        @test ≈(sort(real(vec(x)))', range(-0.75,0.75,7)',
                                               atol=eps(real(one($T))) )
      end
    end
  end           

  for T = (Float32, Float64, ComplexF64 )
    @eval begin
      @testset "Testing calcnodes (TriSymCub method) for DataType $(string($T))" begin
        vtx = ($T)[-1 -1; 1 -1; -1 1]
        tricub = TriSymCub{($T)}()
        @test SymCubatures.calcnodes(tricub, vtx) == vtx[:,:]'

        tricub = TriSymCub{($T)}(numedge=1)
        alpha = ($T)(1/pi)
        A = ($T)[alpha (1-alpha) 0;
                 (1-alpha) alpha 0;
                 0 alpha (1-alpha);
                 0 (1-alpha) alpha;
                 (1-alpha) 0 alpha;
                 alpha 0 (1-alpha)]
        SymCubatures.setparams!(tricub, [alpha])
        x = SymCubatures.calcnodes(tricub, vtx)
        @test vec(x[1,:]) == [vtx[:,1]; (A*vtx[:,1])]
        @test vec(x[2,:]) == [vtx[:,2]; (A*vtx[:,2])]

        tricub = TriSymCub{($T)}(numS21=1)
        SymCubatures.setparams!(tricub, ($T)[2/3])
        x = SymCubatures.calcnodes(tricub, vtx)
        @test ≈(vec(x[1,:]), [vtx[:,1]; fill(sum(vtx[:,1])/3, (3) )],
                                      atol=eps(real(one($T))) )
        @test ≈(vec(x[2,:]), [vtx[:,2]; fill(sum(vtx[:,2])/3, (3) )],
                                      atol=eps(real(one($T))) )

        tricub = TriSymCub{($T)}(numS111=1)
        alpha = ($T)(1/4)
        beta = ($T)(3/4)
        SymCubatures.setparams!(tricub, ($T)[1/4 3/4])
        alpha *= 0.5
        beta *= 0.5
        A = ($T)[alpha beta (1-alpha-beta);
                 beta alpha (1-alpha-beta);
                 (1-alpha-beta) alpha beta;
                 (1-alpha-beta) beta alpha;
                 beta (1-alpha-beta) alpha;
                 alpha (1-alpha-beta) beta]
        x = SymCubatures.calcnodes(tricub, vtx)
        @test vec(x[1,:]) == [vtx[:,1]; (A*vtx[:,1])]
        @test vec(x[2,:]) == [vtx[:,2]; (A*vtx[:,2])]
      end
    end
  end

  for T = (Float32, Float64, ComplexF64 )
    @eval begin
      @testset "Testing calcnodes (TetSymCub method) for DataType $(string($T))" begin
        vtx = ($T)[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
        tetcub = TetSymCub{($T)}()
        @test SymCubatures.calcnodes(tetcub, vtx) == vtx[:,:]'

        tetcub = TetSymCub{($T)}(numedge=1)
        alpha = ($T)(0.25)
        A = ($T)[alpha (1-alpha) 0 0;
                 (1-alpha) alpha 0 0;
                 0 alpha (1-alpha) 0;
                 0 (1-alpha) alpha 0;
                 (1-alpha) 0 alpha 0;
                 alpha 0 (1-alpha) 0;
                 alpha 0 0 (1-alpha);
                 (1-alpha) 0 0 alpha;
                 0 alpha 0 (1-alpha);
                 0 (1-alpha) 0 alpha;
                 0 0 alpha (1-alpha);
                 0 0 (1-alpha) alpha]

        SymCubatures.setparams!(tetcub, [alpha])
        x = SymCubatures.calcnodes(tetcub, vtx)
        @test vec(x[1,:]) == [vtx[:,1]; (A*vtx[:,1])]
        @test vec(x[2,:]) == [vtx[:,2]; (A*vtx[:,2])]
        @test vec(x[3,:]) == [vtx[:,3]; (A*vtx[:,3])]

        tetcub = TetSymCub{($T)}(numfaceS21=1)
        alpha = ($T)(0.1)
        SymCubatures.setparams!(tetcub, [alpha])
        alpha *= 0.5
        A = ($T)[alpha alpha (1-2*alpha);
                 (1-2*alpha) alpha alpha;
                 alpha (1-2*alpha) alpha]
        facevtx = [1 1 2 1;
                   2 4 4 3;
                   3 2 3 4]
        x = SymCubatures.calcnodes(tetcub, vtx)
        @test vec(x[1,:]) == [vtx[:,1]; A*vtx[facevtx[:,1],1]; A*vtx[facevtx[:,2],1]; A*vtx[facevtx[:,3],1]; A*vtx[facevtx[:,4],1]]
        @test vec(x[2,:]) == [vtx[:,2]; A*vtx[facevtx[:,1],2]; A*vtx[facevtx[:,2],2]; A*vtx[facevtx[:,3],2]; A*vtx[facevtx[:,4],2]]
        @test vec(x[3,:]) == [vtx[:,3]; A*vtx[facevtx[:,1],3]; A*vtx[facevtx[:,2],3]; A*vtx[facevtx[:,3],3]; A*vtx[facevtx[:,4],3]]

        tetcub = TetSymCub{($T)}(numS31=1)
        SymCubatures.setparams!(tetcub, ($T)[3/4])
        x = SymCubatures.calcnodes(tetcub, vtx)
        @test vec(x[1,:]) == [vtx[:,1]; fill(sum(vtx[:,1])/4, (4) )]
        @test vec(x[2,:]) == [vtx[:,2]; fill(sum(vtx[:,2])/4, (4) )]
        @test vec(x[3,:]) == [vtx[:,3]; fill(sum(vtx[:,3])/4, (4) )]

        tetcub = TetSymCub{($T)}(numS22=1)
        alpha = ($T)(6/7)
        SymCubatures.setparams!(tetcub, [alpha])
        alpha *= 0.5
        A = ($T)[alpha alpha (0.5-alpha) (0.5-alpha);
                 alpha (0.5-alpha) alpha (0.5-alpha);
                 alpha (0.5-alpha) (0.5-alpha) alpha;
                 (0.5-alpha) alpha alpha (0.5-alpha);
                 (0.5-alpha) alpha (0.5-alpha) alpha;
                 (0.5-alpha) (0.5-alpha) alpha alpha]
        x = SymCubatures.calcnodes(tetcub, vtx)
        @test ≈(vec(x[1,:]), [vtx[:,1]; A*vtx[:,1]]) 
        @test ≈(vec(x[2,:]), [vtx[:,2]; A*vtx[:,2]])
        @test ≈(vec(x[3,:]), [vtx[:,3]; A*vtx[:,3]])
      end
    end
  end

  for T = (Float32, Float64, ComplexF64 )
    @eval begin
      @testset "Testing calcweights (PointSymCub method) for DataType $(string($T))" begin
        point = PointSymCub{($T)}()
        wgt = ($T)(1)
        SymCubatures.setweights!(point, [wgt])
        @test ≈(SymCubatures.calcweights(point), ($T)[1])
      end
    end
  end

  for T = (Float32, Float64, ComplexF64 )
    @eval begin
      @testset "Testing calcweights (LineSymCub method) for DataType $(string($T))" begin
        quad = LineSymCub{($T)}()
        wgt = ($T)(1)
        SymCubatures.setweights!(quad, [wgt])
        @test ≈(SymCubatures.calcweights(quad), fill(wgt, (2)))

        quad = LineSymCub{($T)}(numedge=2, centroid=true, vertices=false)
        w = ($T)[1/3 1/4 1/5]
        SymCubatures.setweights!(quad, w)
        @test ≈(SymCubatures.calcweights(quad), 
                [w[1]*ones(($T), (2))
                 w[2]*ones(($T), (2))
                 w[3]], atol=1e-15)
      end
    end
  end

  for T = (Float32, Float64, ComplexF64 )
    @eval begin
      @testset "Testing calcweights (TriSymCub method) for DataType $(string($T))" begin
        tricub = TriSymCub{($T)}()
        wgt = ($T)(1/3)
        SymCubatures.setweights!(tricub, [wgt])
        @test ≈(SymCubatures.calcweights(tricub), fill(wgt, (3)))

        tricub = TriSymCub{($T)}(midedges=true, numedge=1, numS21=2, numS111=1)
        w = ($T)[1/3 1/4 1/5 1/6 1/7 1/8]
        SymCubatures.setweights!(tricub, w)
        @test ≈(SymCubatures.calcweights(tricub), 
                [w[1]*ones(($T), (3))
                 w[2]*ones(($T), (3))
                 w[3]*ones(($T), (3))
                 w[4]*ones(($T), (3))
                 w[5]*ones(($T), (6))
                 w[6]*ones(($T), (6))], atol=1e-15)
      end
    end
  end

  for T = (Float32, Float64, ComplexF64 )
    @eval begin
      @testset "Testing calcweights (TetSymCub method) for DataType $(string($T))" begin
        tetcub = TetSymCub{($T)}()
        wgt = ($T)(1/3)
        SymCubatures.setweights!(tetcub, [wgt])
        @test ≈(SymCubatures.calcweights(tetcub), fill(wgt, (4)))

        tetcub = TetSymCub{($T)}(numS31=2,
                                 midedges=true, numS22=1,
                                 numedge=1, numfaceS21=1, numS211=2,
                                 numfaceS111=1, numS1111=2)
        w = ($T)[1/3 1/4 1/5 1/6 1/7 1/8 1/9 1/10 1/11 1/12 1/13 1/14]
        SymCubatures.setweights!(tetcub, w)
        @test ≈(SymCubatures.calcweights(tetcub), 
                [w[1]*ones(($T), (4))
                 w[2]*ones(($T), (4))
                 w[3]*ones(($T), (4))
                 w[4]*ones(($T), (6))
                 w[5]*ones(($T), (6))
                 w[6]*ones(($T), (12))
                 w[7]*ones(($T), (12))
                 w[8]*ones(($T), (12))
                 w[9]*ones(($T), (12))
                 w[10]*ones(($T), (24))
                 w[11]*ones(($T), (24))
                 w[12]*ones(($T), (24))], atol=1e-15)
      end
    end
  end

  @testset "Testing calcjacobianofnodes (TriSymCub method)" begin
    # loop over parameters and check Jacobian using complex step
    vtx = Float64[-1 -1; 1 -1; -1 1]
    tricub = TriSymCub{Float64}(midedges=true, numedge=2, numS21=1, numS111=1)
    SymCubatures.setparams!(tricub, [1/3, 2/3, 1/4, 1/5, 4/5])
    Jac = SymCubatures.calcjacobianofnodes(tricub, vtx)
    Jac_cs = zeros(size(Jac))
    tricub_cmplx = TriSymCub{ComplexF64}(midedges=true, numedge=2, numS21=1,
                                         numS111=1)
    params_cmplx = tricub.params .+ 0im
    eps_step = 1e-60
    for i = 1:tricub.numparams
      params_cmplx[i] += eps_step*im
      SymCubatures.setparams!(tricub_cmplx, params_cmplx)
      xc = SymCubatures.calcnodes(tricub_cmplx,
                                  convert(Array{ComplexF64}, vtx))
      Jac_cs[1:tricub.numnodes,i] = imag(xc[1,:])/eps_step
      Jac_cs[tricub.numnodes+1:2*tricub.numnodes,i] = imag(xc[2,:])/eps_step
      params_cmplx[i] -= eps_step*im
    end
    @test ≈(Jac, -Jac_cs, atol=1e-15)
  end

  @testset "Testing calcjacobianofnodes (TetSymCub method)" begin
    # loop over parameters and check Jacobian using complex step
    vtx = Float64[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
    tetcub = TetSymCub{Float64}(midedges=true, numS31=1, numS22=1, numedge=2,
                                numfaceS21=1, numS211=2, numfaceS111=1,
                                numS1111=2)
    SymCubatures.setparams!(tetcub, rand(tetcub.numparams))
    Jac = SymCubatures.calcjacobianofnodes(tetcub, vtx)
    Jac_cs = zeros(size(Jac))
    tetcub_cmplx = TetSymCub{ComplexF64}(midedges=true, numS31=1, numS22=1,
                                         numedge=2, numfaceS21=1, numS211=2,
                                         numfaceS111=1, numS1111=2)
    params_cmplx = tetcub.params .+ 0im
    eps_step = 1e-60
    for i = 1:tetcub.numparams
      params_cmplx[i] += eps_step*im
      SymCubatures.setparams!(tetcub_cmplx, params_cmplx)
      xc = SymCubatures.calcnodes(tetcub_cmplx,
                                  convert(Array{ComplexF64}, vtx))
      Jac_cs[1:tetcub.numnodes,i] = imag(xc[1,:])/eps_step
      Jac_cs[tetcub.numnodes+1:2*tetcub.numnodes,i] = imag(xc[2,:])/eps_step
      Jac_cs[2*tetcub.numnodes+1:3*tetcub.numnodes,i] = imag(xc[3,:])/eps_step
      params_cmplx[i] -= eps_step*im
    end
    @test ≈(Jac, -Jac_cs, atol=1e-15)
  end

  @testset "Testing calcjacobianofweights (TriSymCub method)" begin
    # loop over weights and check Jacobian using complex step
    tricub = TriSymCub{Float64}(midedges=true, numedge=1, numS21=2, numS111=1)
    w = Float64[1/3 1/4 1/5 1/6 1/7 1/8]
    SymCubatures.setweights!(tricub, w)
    Jac = SymCubatures.calcjacobianofweights(tricub)
    Jac_cs = zeros(size(Jac))
    tricub_cmplx = TriSymCub{ComplexF64}(midedges=true, numedge=1, numS21=2,
                                         numS111=1)
    weights_cmplx = tricub.weights .+ 0im
    eps_step = 1e-60
    for i = 1:tricub.numweights
      weights_cmplx[i] += eps_step*im
      SymCubatures.setweights!(tricub_cmplx, weights_cmplx)
      wc = SymCubatures.calcweights(tricub_cmplx)
      Jac_cs[:,i] = imag(wc)/eps_step
      weights_cmplx[i] -= eps_step*im
    end
    @test ≈(Jac, Jac_cs, atol=1e-15)
  end

  @testset "Testing calcjacobianofweights (TetSymCub method)" begin
    # loop over weights and check Jacobian using complex step
    tetcub = TetSymCub{Float64}(numS31=2,
                                midedges=true, numS22=1,
                                numedge=1, numfaceS21=1, numS211=2,
                                numfaceS111=1, numS1111=2)
    w = Float64[1/3 1/4 1/5 1/6 1/7 1/8 1/9 1/10 1/11 1/12 1/13 1/14]
    SymCubatures.setweights!(tetcub, w)
    Jac = SymCubatures.calcjacobianofweights(tetcub)
    Jac_cs = zeros(size(Jac))
    tetcub_cmplx = TetSymCub{ComplexF64}(numS31=2,
                                         midedges=true, numS22=1,
                                         numedge=1, numfaceS21=1, numS211=2,
                                         numfaceS111=1, numS1111=2)
    weights_cmplx = tetcub.weights .+ 0im
    eps_step = 1e-60
    for i = 1:tetcub.numweights
      weights_cmplx[i] += eps_step*im
      SymCubatures.setweights!(tetcub_cmplx, weights_cmplx)
      wc = SymCubatures.calcweights(tetcub_cmplx)
      Jac_cs[:,i] = imag(wc)/eps_step
      weights_cmplx[i] -= eps_step*im
    end
    @test ≈(Jac, Jac_cs, atol=1e-15)
  end

  @testset "Testing calcjacobian (TriSymCub method)" begin
    # This is essentially just a copy of the code in calcjacobian, so not good test
    vtx = Float64[-1 -1; 1 -1; -1 1]
    tricub = TriSymCub{Float64}(midedges=true, numedge=1, numS21=1, numS111=1)
    SymCubatures.setparams!(tricub, [1/3, 2/3, 3/4, 4/5])
    SymCubatures.setweights!(tricub, [1/3 1/4 1/5 1/6 1/7])
    Jac_params = SymCubatures.calcjacobianofnodes(tricub, vtx)
    Jac_weights = SymCubatures.calcjacobianofweights(tricub)
    Jac = SymCubatures.calcjacobian(tricub, vtx)
    numnodes = tricub.numnodes
    numparams = tricub.numparams
    numweights = tricub.numweights
    @test ≈(Jac[1:2*numnodes,1:numparams], Jac_params, atol=1e-15)
    @test ≈(Jac[1:2*numnodes,numparams+1:end], zeros(Float64, (2*numnodes, numweights)), atol=1e-15)
    @test ≈(Jac[2*numnodes+1:end,1:numparams], zeros(Float64, (numnodes, numparams)), atol=1e-15)
    @test ≈(Jac[2*numnodes+1:end,numparams+1:end], Jac_weights, atol=1e-15)
  end

  @testset "Testing calcjacobian (TetSymCub method)" begin
    # This is essentially just a copy of the code in calcjacobian, so not good test
    vtx = Float64[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
    tetcub = TetSymCub{Float64}(midedges=true, numedge=1, numfaceS21=2, numS31=1,
                                numS22=2, numS211=1, numfaceS111=1, numS1111=1)
    SymCubatures.setparams!(tetcub, [1/3, 2/3, 3/4, 4/5, 5/6, 6/7, 7/8, 8/9,
                                     9/10, 10/11, 11/12, 12/13, 13/14])
    SymCubatures.setweights!(tetcub, [1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 1/9, 1/10,
                                      1/11, 1/12, 1/13])
    Jac_params = SymCubatures.calcjacobianofnodes(tetcub, vtx)
    Jac_weights = SymCubatures.calcjacobianofweights(tetcub)
    Jac = SymCubatures.calcjacobian(tetcub, vtx)
    numnodes = tetcub.numnodes
    numparams = tetcub.numparams
    numweights = tetcub.numweights
    @test ≈(Jac[1:3*numnodes,1:numparams], Jac_params, atol=1e-15)
    @test ≈(Jac[1:3*numnodes,numparams+1:end], zeros(Float64, (3*numnodes, numweights)), atol=1e-15)
    @test ≈(Jac[3*numnodes+1:end,1:numparams], zeros(Float64, (numnodes, numparams)), atol=1e-15)
    @test ≈(Jac[3*numnodes+1:end,numparams+1:end], Jac_weights, atol=1e-15)
  end

  @testset "Testing getInternalParamMask (TriSymCub method)" begin
    tricub = TriSymCub{Float64}(vertices=true, midedges=true, centroid=true, numedge=2,
                                numS21=3, numS111=2)
    mask = SymCubatures.getInternalParamMask(tricub)
    @test mask == [1; 2; 3; 6; 7; 8; 9]    
  end

  @testset "Testing getInternalParamMask (TetSymCub method)" begin    
    tetcub = TetSymCub{Float64}(vertices=true, midedges=true, centroid=true,
                                facecentroid=true, numedge=2, numfaceS21=3, numS31=1,
                                numS22=2, numS211=1, numfaceS111=1, numS1111=1)
    mask = SymCubatures.getInternalParamMask(tetcub)
    @test mask == [1; 2; 3; 9; 10; 13; 14; 15]
  end
    
end
