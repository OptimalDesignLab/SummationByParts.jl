facts("Testing SymCubatures Module...") do
  
  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing LineSymCub Inner Constructor for DataType "string($T)) do
        quad = LineSymCub{($T)}(vertices=false)
        @fact quad.numparams --> 0
        @fact quad.numweights --> 0
        @fact quad.numnodes --> 0
        quad = LineSymCub{($T)}(numedge=2)
        @fact quad.numparams --> 2
        @fact quad.numweights --> 3
        @fact quad.numnodes --> 2+2*2
        quad = LineSymCub{($T)}(numedge=3, centroid=true)
        @fact quad.numparams --> 3
        @fact quad.numweights --> 5
        @fact quad.numnodes --> 2+2*3+1
      end
    end
  end
  
  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing TriSymCub Inner Constructor for DataType "string($T)) do
        tricub = TriSymCub{($T)}(vertices=false)
        @fact tricub.numparams --> 0
        @fact tricub.numweights --> 0
        @fact tricub.numnodes --> 0
        tricub = TriSymCub{($T)}(numedge=2,midedges=true)
        @fact tricub.numparams --> 2
        @fact tricub.numweights --> 4
        @fact tricub.numnodes --> 3+3+2*6
        tricub = TriSymCub{($T)}(midedges=true, numS21=3)
        @fact tricub.numparams --> 3
        @fact tricub.numweights --> 5
        @fact tricub.numnodes --> 3+3+3*3
        tricub = TriSymCub{($T)}(numedge=2, centroid=true, numS111=2)
        @fact tricub.numparams --> 2+4
        @fact tricub.numweights --> 6
        @fact tricub.numnodes --> 3+2*6+1+2*6
      end
    end
  end
  
  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing TetSymCub Inner Constructor for DataType "string($T)) do
        tetcub = TetSymCub{($T)}(vertices=false)
        @fact tetcub.numparams --> 0
        @fact tetcub.numweights --> 0
        @fact tetcub.numnodes --> 0
        tetcub = TetSymCub{($T)}(numedge=2,midedges=true)
        @fact tetcub.numparams --> 2
        @fact tetcub.numweights --> 4
        @fact tetcub.numnodes --> 4+6+2*12
        tetcub = TetSymCub{($T)}(midedges=true, numS31=3)
        @fact tetcub.numparams --> 3
        @fact tetcub.numweights --> 5
        @fact tetcub.numnodes --> 4+6+4*3
        tetcub = TetSymCub{($T)}(midedges=true, numedge=1, numfaceS21=2)
        @fact tetcub.numparams --> 3
        @fact tetcub.numweights --> 5
        @fact tetcub.numnodes --> 4+6+1*12+2*12
        tetcub = TetSymCub{($T)}(midedges=true, centroid=true, numS22=2)
        @fact tetcub.numparams --> 2
        @fact tetcub.numweights --> 5
        @fact tetcub.numnodes --> 4+6+1+2*6
      end
    end
  end
  
  context("Testing getnumboundarynodes (LineSymCub method)") do
    quad = LineSymCub{Float64}() # vertex only rule
    @fact SymCubatures.getnumboundarynodes(quad) --> 2
    quad = LineSymCub{Float64}(numedge = 2, centroid=true, vertices=false)
    @fact SymCubatures.getnumboundarynodes(quad) --> 0
  end

  context("Testing getnumboundarynodes (TriSymCub method)") do
    tricub = TriSymCub{Float64}() # vertex only rule
    @fact SymCubatures.getnumboundarynodes(tricub) --> 3
    tricub = TriSymCub{Float64}(numedge = 2, midedges=true, numS21 = 4)
    @fact SymCubatures.getnumboundarynodes(tricub) --> 3+3+2*6
  end

  context("Testing getnumboundarynodes (TetSymCub method)") do
    tetcub = TetSymCub{Float64}() # vertex only rule
    @fact SymCubatures.getnumboundarynodes(tetcub) --> 4
    tetcub = TetSymCub{Float64}(numedge=2, midedges=true, numfaceS21=3,
                                numS31 = 4)
    @fact SymCubatures.getnumboundarynodes(tetcub) --> 4+6+2*12+3*12
  end

  context("Testing getnumfacenodes (LineSymCub method)") do
    quad = LineSymCub{Float64}() # vertex only rule
    @fact SymCubatures.getnumfacenodes(quad) --> 1
    quad = LineSymCub{Float64}(numedge = 2, centroid=true, vertices=false)
    @fact SymCubatures.getnumfacenodes(quad) --> 0
  end

  context("Testing getnumfacenodes (TriSymCub method)") do
    tricub = TriSymCub{Float64}() # vertex only rule
    @fact SymCubatures.getnumfacenodes(tricub) --> 2
    tricub = TriSymCub{Float64}(numedge = 2, midedges=true, numS21 = 4)
    @fact SymCubatures.getnumfacenodes(tricub) --> 2+2*2+1
  end

  context("Testing getnumfacenodes (TetSymCub method)") do
    tetcub = TetSymCub{Float64}() # vertex only rule
    @fact SymCubatures.getnumfacenodes(tetcub) --> 3
    tetcub = TetSymCub{Float64}(numedge=2, midedges=true, numfaceS21=3,
                                numS31 = 4)
    @fact SymCubatures.getnumfacenodes(tetcub) --> 3+3+2*6+3*3
  end

  context("Testing getbndrynodeindices (LineSymCub method)") do
    quad = LineSymCub{Float64}(numedge=1, centroid=true)
    bndryindices = SymCubatures.getbndrynodeindices(quad)
    @fact bndryindices --> [1;2]
  end

  context("Testing getbndrynodeindices (TriSymCub method)") do
    tricub = TriSymCub{Float64}(numedge=1, midedges=true, numS21 = 4)
    bndryindices = SymCubatures.getbndrynodeindices(tricub)
    @fact bndryindices --> [1; 2; 3; 4; 5; 6; 19; 20; 21; 22; 23; 24]
  end

  context("Testing getbndrynodeindices (TetSymCub method)") do
    tetcub = TetSymCub{Float64}(numedge=2, facecentroid=true, numS31=1, 
                                midedges=true, numS22=1, numfaceS21 = 1)
    bndryindices = SymCubatures.getbndrynodeindices(tetcub)
    @fact bndryindices --> [1; 2; 3; 4; 
                            5; 6; 7; 8;
                            13; 14; 15; 16; 17; 18;
                            25; 26; 27; 28; 29; 30; 31; 32; 33; 34; 35; 36;
                            37; 38; 39; 40; 41; 42; 43; 44; 45; 46; 47; 48;
                            49; 50; 51; 52; 53; 54; 55; 56; 57; 58; 59; 60]
  end

  context("Testing getinteriornodeindices (LineSymCub method)") do
    quad = LineSymCub{Float64}(numedge=1, centroid=true)
    indices = SymCubatures.getinteriornodeindices(quad)
    @fact indices --> [3;4;5]
  end

  context("Testing getinteriornodeindices (TriSymCub method)") do
    tricub = TriSymCub{Float64}(numedge=1, midedges=true, numS21=4,
                                centroid=true)
    indices = SymCubatures.getinteriornodeindices(tricub)
    @fact indices --> [7:18; 25]
  end

  context("Testing getinteriornodeindices (TetSymCub method)") do
    tetcub = TetSymCub{Float64}(numedge=2, facecentroid=true, numS31=1, 
                                midedges=true, numS22=1, numfaceS21=1,
                                centroid=true)
    indices = SymCubatures.getinteriornodeindices(tetcub)
    @fact indices --> [9:12; 19:24; 61]
  end

  context("Testing getfacenodeindices (LineSymCub method)") do
    quad = LineSymCub{Float64}(numedge=1, centroid=true)
    bndryindices = SymCubatures.getfacenodeindices(quad)
    @fact bndryindices --> [1 2]
  end

  context("Testing getfacenodeindices (TriSymCub method)") do
    tricub = TriSymCub{Float64}(numedge=1, midedges=true, numS21 = 4)
    bndryindices = SymCubatures.getfacenodeindices(tricub)
    @fact bndryindices --> [1 2 3;
                            2 3 1;
                            4 5 6;
                            19 21 23;
                            20 22 24]
  end

  context("Testing getfacenodeindices (TetSymCub method)") do
    tetcub = TetSymCub{Float64}(numedge=1, facecentroid=true, numS31=1, 
                                midedges=true, numS22=1, numfaceS21=1,
                                centroid=true)
    bndryindices = SymCubatures.getfacenodeindices(tetcub)
    @fact bndryindices --> [1 1 2 1; # vertices
                            2 4 4 3; #
                            3 2 3 4; # 
                            5 6 7 8; # face centroids
                            13 16 17 15; 14 17 18 18; 15 13 14 16; # midedges
                            25 31 33 30; 26 32 34 29; 27 34 36 35; 
                            28 33 35 36; 29 26 28 32; 30 25 27 31; # edge nodes
                            37 46 40 43; 38 47 41 44; 39 48 42 45] # face S21
  end

  context("Testing findleftperm!") do
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
    @fact success --> true
    @fact A[perm,:] --> A[:,permR]

    # Next, transpose the matrix, and check that we recover the right permutation
    permL = deepcopy(perm)
    resize!(perm, m)
    success = SymCubatures.findleftperm!(A.', permL, perm)
    @fact success --> true
    @fact A[:,perm] --> A[permL,:]
    @fact perm --> permR

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
    @fact success --> false
  end

  context("Test getfacebasedpermutation (LineSymCub method)") do
    quad = LineSymCub{Float64}(numedge=1, centroid=true)
    perm = SymCubatures.getfacebasedpermutation(quad)
    @fact perm[:,1] --> [1:5;]
    @fact perm[:,2] --> [2; 1; 4; 3; 5]
    perm = SymCubatures.getfacebasedpermutation(quad, faceonly=true)
    @fact perm[:,1] --> [1;]
    @fact perm[:,2] --> [2;]
  end

  context("Test getfacebasedpermutation (TriSymCub method)") do
    tricub = TriSymCub{Float64}(numedge=1, midedges=true, numS21 = 1, numS111=1,
                                centroid=true)
    perm = SymCubatures.getfacebasedpermutation(tricub)
    @fact perm[:,1] --> [1:22;]
    @fact perm[:,2] --> [2; 3; 1; 5; 6; 4; 8; 9; 7;
                         12; 13; 14; 15; 10; 11; 
                         18; 19; 20; 21; 16; 17; 
                         22]
    @fact perm[:,3] --> [3; 1; 2; 6; 4; 5; 9; 7; 8;
                         14; 15; 10; 11; 12; 13; 
                         20; 21; 16; 17; 18; 19;
                         22]
    perm = SymCubatures.getfacebasedpermutation(tricub, faceonly=true)
    @fact perm[:,1] --> [1; 2; 4; 10; 11]
    @fact perm[:,2] --> [2; 3; 5; 12; 13]
    @fact perm[:,3] --> [3; 1; 6; 14; 15]
  end

  context("Test getfacebasedpermutation (TetSymCub method)") do
    # This test works by providing different sets of vertices to calcnodes, and
    # then checking for equivalence with a base set of vertices reordered using
    # the face-based permutation
    tetcub = TetSymCub{Float64}(vertices=true, midedges=true, centroid=true,
                                facecentroid=true, numedge=2, numfaceS21=2,
                                numS31=2, numS22=2)
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
    @fact xf1[:,perm[:,1]] --> roughly(xf1, atol=eps())
    @fact xf1[:,perm[:,2]] --> roughly(xf2, atol=eps())
    @fact xf1[:,perm[:,3]] --> roughly(xf3, atol=eps())
    @fact xf1[:,perm[:,4]] --> roughly(xf4, atol=eps())
  end

  context("Test getneighbourpermutation (LineSymCub method)") do
    quad = LineSymCub{Float64}(numedge=1, centroid=true)
    perm = SymCubatures.getneighbourpermutation(quad)
    @fact perm[:,1] --> [2; 1; 4; 3; 5]
    quad = LineSymCub{Float64}(vertices=false, numedge=1, centroid=true)
    perm = SymCubatures.getneighbourpermutation(quad)
    @fact perm[:,1] --> [2; 1; 3]
  end

  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing calcnodes (LineSymCub method) for DataType "string($T)) do
        vtx = reshape(($T)[-1; 1], (2,1))
        quad = LineSymCub{($T)}()
        @fact SymCubatures.calcnodes(quad, vtx) --> vtx[:,1].'

        quad = LineSymCub{($T)}(numedge=1)
        alpha = ($T)(1/pi)
        A = ($T)[alpha (1-alpha);
                 (1-alpha) alpha]
        SymCubatures.setparams!(quad, [alpha])
        x = SymCubatures.calcnodes(quad, vtx)
        @fact x --> [vtx[:,1].' (A*vtx[:,1]).']

        # uniformly spaced 7 nodes between (-1,1) (no vertices)
        quad = LineSymCub{($T)}(numedge=3, centroid=true, vertices=false)
        SymCubatures.setparams!(quad, ($T)[5/8 6/8 7/8])
        x = SymCubatures.calcnodes(quad, vtx)
        @fact sort(real(vec(x))).' --> roughly(linspace(-0.75,0.75,7).',
                                               atol=eps(real(one($T))) )
      end
    end
  end           

  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing calcnodes (TriSymCub method) for DataType "string($T)) do
        vtx = ($T)[-1 -1; 1 -1; -1 1]
        tricub = TriSymCub{($T)}()
        @fact SymCubatures.calcnodes(tricub, vtx) --> vtx[:,:].'

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
        @fact vec(x[1,:]) --> [vtx[:,1]; (A*vtx[:,1])]
        @fact vec(x[2,:]) --> [vtx[:,2]; (A*vtx[:,2])]

        tricub = TriSymCub{($T)}(numS21=1)
        SymCubatures.setparams!(tricub, ($T)[2/3])
        x = SymCubatures.calcnodes(tricub, vtx)
        @fact vec(x[1,:]) --> roughly([vtx[:,1]; fill(sum(vtx[:,1])/3, (3) )],
                                      atol=eps(real(one($T))) )
        @fact vec(x[2,:]) --> roughly([vtx[:,2]; fill(sum(vtx[:,2])/3, (3) )],
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
        @fact vec(x[1,:]) --> [vtx[:,1]; (A*vtx[:,1])]
        @fact vec(x[2,:]) --> [vtx[:,2]; (A*vtx[:,2])]
      end
    end
  end

  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing calcnodes (TetSymCub method) for DataType "string($T)) do
        vtx = ($T)[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
        tetcub = TetSymCub{($T)}()
        @fact SymCubatures.calcnodes(tetcub, vtx) --> vtx[:,:].'

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
        @fact vec(x[1,:]) --> [vtx[:,1]; (A*vtx[:,1])]
        @fact vec(x[2,:]) --> [vtx[:,2]; (A*vtx[:,2])]
        @fact vec(x[3,:]) --> [vtx[:,3]; (A*vtx[:,3])]

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
        @fact vec(x[1,:]) --> [vtx[:,1]; A*vtx[facevtx[:,1],1]; A*vtx[facevtx[:,2],1]; A*vtx[facevtx[:,3],1]; A*vtx[facevtx[:,4],1]]
        @fact vec(x[2,:]) --> [vtx[:,2]; A*vtx[facevtx[:,1],2]; A*vtx[facevtx[:,2],2]; A*vtx[facevtx[:,3],2]; A*vtx[facevtx[:,4],2]]
        @fact vec(x[3,:]) --> [vtx[:,3]; A*vtx[facevtx[:,1],3]; A*vtx[facevtx[:,2],3]; A*vtx[facevtx[:,3],3]; A*vtx[facevtx[:,4],3]]

        tetcub = TetSymCub{($T)}(numS31=1)
        SymCubatures.setparams!(tetcub, ($T)[3/4])
        x = SymCubatures.calcnodes(tetcub, vtx)
        @fact vec(x[1,:]) --> [vtx[:,1]; fill(sum(vtx[:,1])/4, (4) )]
        @fact vec(x[2,:]) --> [vtx[:,2]; fill(sum(vtx[:,2])/4, (4) )]
        @fact vec(x[3,:]) --> [vtx[:,3]; fill(sum(vtx[:,3])/4, (4) )]

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
        @fact vec(x[1,:]) --> roughly([vtx[:,1]; A*vtx[:,1]]) 
        @fact vec(x[2,:]) --> roughly([vtx[:,2]; A*vtx[:,2]])
        @fact vec(x[3,:]) --> roughly([vtx[:,3]; A*vtx[:,3]])
      end
    end
  end

  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing calcweights (LineSymCub method) for DataType "string($T)) do
        quad = LineSymCub{($T)}()
        wgt = ($T)(1)
        SymCubatures.setweights!(quad, [wgt])
        @fact SymCubatures.calcweights(quad) --> roughly(fill(wgt, (2)))

        quad = LineSymCub{($T)}(numedge=2, centroid=true, vertices=false)
        w = ($T)[1/3 1/4 1/5]
        SymCubatures.setweights!(quad, w)
        @fact SymCubatures.calcweights(quad) -->
        roughly([w[1]*ones(($T), (2))
                 w[2]*ones(($T), (2))
                 w[3]], atol=1e-15)
      end
    end
  end

  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing calcweights (TriSymCub method) for DataType "string($T)) do
        tricub = TriSymCub{($T)}()
        wgt = ($T)(1/3)
        SymCubatures.setweights!(tricub, [wgt])
        @fact SymCubatures.calcweights(tricub) --> roughly(fill(wgt, (3)))

        tricub = TriSymCub{($T)}(midedges=true, numedge=1, numS21=2, numS111=1)
        w = ($T)[1/3 1/4 1/5 1/6 1/7 1/8]
        SymCubatures.setweights!(tricub, w)
        @fact SymCubatures.calcweights(tricub) -->
        roughly([w[1]*ones(($T), (3))
                 w[2]*ones(($T), (3))
                 w[3]*ones(($T), (3))
                 w[4]*ones(($T), (3))
                 w[5]*ones(($T), (6))
                 w[6]*ones(($T), (6))], atol=1e-15)
      end
    end
  end

  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing calcweights (TetSymCub method) for DataType "string($T)) do
        tetcub = TetSymCub{($T)}()
        wgt = ($T)(1/3)
        SymCubatures.setweights!(tetcub, [wgt])
        @fact SymCubatures.calcweights(tetcub) --> roughly(fill(wgt, (4)))

        tetcub = TetSymCub{($T)}(midedges=true, numedge=1, numfaceS21=1,
                                 numS31=2, numS22=1)
        w = ($T)[1/3 1/4 1/5 1/6 1/7 1/8 1/9]
        SymCubatures.setweights!(tetcub, w)
        @fact SymCubatures.calcweights(tetcub) -->
        roughly([w[1]*ones(($T), (4))
                 w[2]*ones(($T), (4))
                 w[3]*ones(($T), (4))
                 w[4]*ones(($T), (6))
                 w[5]*ones(($T), (6))
                 w[6]*ones(($T), (12))
                 w[7]*ones(($T), (12))], atol=1e-15)
      end
    end
  end

  context("Testing calcjacobianofnodes (TriSymCub method)") do
    # loop over parameters and check Jacobian using complex step
    vtx = Float64[-1 -1; 1 -1; -1 1]
    tricub = TriSymCub{Float64}(midedges=true, numedge=2, numS21=1, numS111=1)
    SymCubatures.setparams!(tricub, [1/3, 2/3, 1/4, 1/5, 4/5])
    Jac = SymCubatures.calcjacobianofnodes(tricub, vtx)
    Jac_cs = zeros(Jac)
    tricub_cmplx = TriSymCub{Complex128}(midedges=true, numedge=2, numS21=1,
                                         numS111=1)
    params_cmplx = tricub.params + 0im
    eps_step = 1e-60
    for i = 1:tricub.numparams
      params_cmplx[i] += eps_step*im
      SymCubatures.setparams!(tricub_cmplx, params_cmplx)
      xc = SymCubatures.calcnodes(tricub_cmplx,
                                  convert(Array{Complex128}, vtx))
      Jac_cs[1:tricub.numnodes,i] = imag(xc[1,:])/eps_step
      Jac_cs[tricub.numnodes+1:2*tricub.numnodes,i] = imag(xc[2,:])/eps_step
      params_cmplx[i] -= eps_step*im
    end
    @fact Jac --> roughly(Jac_cs, atol=1e-15)
  end

  context("Testing calcjacobianofnodes (TetSymCub method)") do
    # loop over parameters and check Jacobian using complex step
    vtx = Float64[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
    tetcub = TetSymCub{Float64}(midedges=true, numS31=1, numS22=1, numedge=2,
                                numfaceS21=1)
    SymCubatures.setparams!(tetcub, [1/4, 1/5, 1/3, 2/3, 1/10])
    Jac = SymCubatures.calcjacobianofnodes(tetcub, vtx)
    Jac_cs = zeros(Jac)
    tetcub_cmplx = TetSymCub{Complex128}(midedges=true, numS31=1, numS22=1,
                                         numedge=2, numfaceS21=1)
    params_cmplx = tetcub.params + 0im
    eps_step = 1e-60
    for i = 1:tetcub.numparams
      params_cmplx[i] += eps_step*im
      SymCubatures.setparams!(tetcub_cmplx, params_cmplx)
      xc = SymCubatures.calcnodes(tetcub_cmplx,
                                  convert(Array{Complex128}, vtx))
      Jac_cs[1:tetcub.numnodes,i] = imag(xc[1,:])/eps_step
      Jac_cs[tetcub.numnodes+1:2*tetcub.numnodes,i] = imag(xc[2,:])/eps_step
      Jac_cs[2*tetcub.numnodes+1:3*tetcub.numnodes,i] = imag(xc[3,:])/eps_step
      params_cmplx[i] -= eps_step*im
    end
    @fact Jac --> roughly(Jac_cs, atol=1e-15)
  end

  context("Testing calcjacobianofweights (TriSymCub method)") do
    # loop over weights and check Jacobian using complex step
    tricub = TriSymCub{Float64}(midedges=true, numedge=1, numS21=2, numS111=1)
    w = Float64[1/3 1/4 1/5 1/6 1/7 1/8]
    SymCubatures.setweights!(tricub, w)
    Jac = SymCubatures.calcjacobianofweights(tricub)
    Jac_cs = zeros(Jac)
    tricub_cmplx = TriSymCub{Complex128}(midedges=true, numedge=1, numS21=2,
                                         numS111=1)
    weights_cmplx = tricub.weights + 0im
    eps_step = 1e-60
    for i = 1:tricub.numweights
      weights_cmplx[i] += eps_step*im
      SymCubatures.setweights!(tricub_cmplx, weights_cmplx)
      wc = SymCubatures.calcweights(tricub_cmplx)
      Jac_cs[:,i] = imag(wc)/eps_step
      weights_cmplx[i] -= eps_step*im
    end
    @fact Jac --> roughly(Jac_cs, atol=1e-15)
  end

  context("Testing calcjacobianofweights (TetSymCub method)") do
    # loop over weights and check Jacobian using complex step
    tetcub = TetSymCub{Float64}(midedges=true, numedge=1, numfaceS21=2, numS31=2,
                                numS22=1)
    w = Float64[1/3 1/4 1/5 1/6 1/7 1/8 1/9 1/10]
    SymCubatures.setweights!(tetcub, w)
    Jac = SymCubatures.calcjacobianofweights(tetcub)
    Jac_cs = zeros(Jac)
    tetcub_cmplx = TetSymCub{Complex128}(midedges=true, numedge=1, numfaceS21=2, 
                                         numS31=2, numS22=1)
    weights_cmplx = tetcub.weights + 0im
    eps_step = 1e-60
    for i = 1:tetcub.numweights
      weights_cmplx[i] += eps_step*im
      SymCubatures.setweights!(tetcub_cmplx, weights_cmplx)
      wc = SymCubatures.calcweights(tetcub_cmplx)
      Jac_cs[:,i] = imag(wc)/eps_step
      weights_cmplx[i] -= eps_step*im
    end
    @fact Jac --> roughly(Jac_cs, atol=1e-15)
  end

  context("Testing calcjacobian (TriSymCub method)") do
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
    @fact Jac[1:2*numnodes,1:numparams] --> roughly(Jac_params, atol=1e-15)
    @fact Jac[1:2*numnodes,numparams+1:end] -->
    roughly(zeros(Float64, (2*numnodes, numweights)), atol=1e-15)
    @fact Jac[2*numnodes+1:end,1:numparams] -->
    roughly(zeros(Float64, (numnodes, numparams)), atol=1e-15)
    @fact Jac[2*numnodes+1:end,numparams+1:end] --> roughly(Jac_weights, atol=1e-15)
  end

  context("Testing calcjacobian (TetSymCub method)") do
    # This is essentially just a copy of the code in calcjacobian, so not good test
    vtx = Float64[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
    tetcub = TetSymCub{Float64}(midedges=true, numedge=1, numfaceS21=2, numS31=1,
                                numS22=2)
    SymCubatures.setparams!(tetcub, [1/3, 2/3, 3/4, 4/5, 5/6, 6/7])
    SymCubatures.setweights!(tetcub, [1/3 1/4 1/5 1/6 1/7 1/8 1/9 1/10])
    Jac_params = SymCubatures.calcjacobianofnodes(tetcub, vtx)
    Jac_weights = SymCubatures.calcjacobianofweights(tetcub)
    Jac = SymCubatures.calcjacobian(tetcub, vtx)
    numnodes = tetcub.numnodes
    numparams = tetcub.numparams
    numweights = tetcub.numweights
    @fact Jac[1:3*numnodes,1:numparams] --> roughly(Jac_params, atol=1e-15)
    @fact Jac[1:3*numnodes,numparams+1:end] -->
    roughly(zeros(Float64, (3*numnodes, numweights)), atol=1e-15)
    @fact Jac[3*numnodes+1:end,1:numparams] -->
    roughly(zeros(Float64, (numnodes, numparams)), atol=1e-15)
    @fact Jac[3*numnodes+1:end,numparams+1:end] --> roughly(Jac_weights, atol=1e-15)
  end

end
