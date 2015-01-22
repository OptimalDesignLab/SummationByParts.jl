# Tests for SymCubatures module
using SummationByParts.SymCubatures

facts("Testing SymCubatures Module...") do
  
  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing TriSymCub Inner Constructor for DataType "string($T)) do
        tricub = TriSymCub{($T)}(vertices=false)
        @fact tricub.numparams => 0
        @fact tricub.numweights => 0
        @fact tricub.numnodes => 0
        tricub = TriSymCub{($T)}(numedge=2,midedges=true)
        @fact tricub.numparams => 2
        @fact tricub.numweights => 4
        @fact tricub.numnodes => 3+3+2*6
        tricub = TriSymCub{($T)}(midedges=true, numS21=3)
        @fact tricub.numparams => 3
        @fact tricub.numweights => 5
        @fact tricub.numnodes => 3+3+3*3
      end
    end
  end

  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing TetSymCub Inner Constructor for DataType "string($T)) do
        tetcub = TetSymCub{($T)}(vertices=false)
        @fact tetcub.numparams => 0
        @fact tetcub.numweights => 0
        @fact tetcub.numnodes => 0
        tetcub = TetSymCub{($T)}(numedge=2,midedges=true)
        @fact tetcub.numparams => 2
        @fact tetcub.numweights => 4
        @fact tetcub.numnodes => 4+6+2*12
        tetcub = TetSymCub{($T)}(midedges=true, numS31=3)
        @fact tetcub.numparams => 3
        @fact tetcub.numweights => 5
        @fact tetcub.numnodes => 4+6+4*3
      end
    end
  end

  context("Testing getnumboundarynodes (TriSymCub method)") do
    tricub = TriSymCub{Float64}() # vertex only rule
    @fact SymCubatures.getnumboundarynodes(tricub) => 3
    tricub = TriSymCub{Float64}(numedge = 2, midedges=true, numS21 = 4)
    @fact SymCubatures.getnumboundarynodes(tricub) => 3+3+2*6
  end

  context("Testing getnumboundarynodes (TetSymCub method)") do
    tetcub = TetSymCub{Float64}() # vertex only rule
    @fact SymCubatures.getnumboundarynodes(tetcub) => 4
    tetcub = TetSymCub{Float64}(numedge = 2, midedges=true, numS31 = 4)
    @fact SymCubatures.getnumboundarynodes(tetcub) => 4+6+2*12
  end

  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing calcnodes (TriSymCub method) for DataType "string($T)) do
        vtx = ($T)[-1 -1; 1 -1; -1 1]
        tricub = TriSymCub{($T)}()
        @fact SymCubatures.calcnodes(tricub, vtx) => vtx[:,1], vtx[:,2]

        tricub = TriSymCub{($T)}(numedge=1)
        alpha = ($T)(1/pi)
        A = ($T)[alpha (1-alpha) 0;
                 (1-alpha) alpha 0;
                 0 alpha (1-alpha);
                 0 (1-alpha) alpha;
                 (1-alpha) 0 alpha;
                 alpha 0 (1-alpha)]
        SymCubatures.setparams!(tricub, [alpha])
        x, y = SymCubatures.calcnodes(tricub, vtx)
        @fact x.' => [vtx[:,1].' (A*vtx[:,1]).']
        @fact y.' => [vtx[:,2].' (A*vtx[:,2]).']

        tricub = TriSymCub{($T)}(numS21=1)
        SymCubatures.setparams!(tricub, ($T)[2/3])
        x, y = SymCubatures.calcnodes(tricub, vtx)
        @fact x.' => roughly([vtx[:,1].' fill(sum(vtx[:,1])/3, (1,3) )],
                             atol=eps(real(one($T))) )
        @fact y.' => roughly([vtx[:,2].' fill(sum(vtx[:,2])/3, (1,3) )],
                             atol=eps(real(one($T))) )
      end
    end
  end

  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing calcnodes (TetSymCub method) for DataType "string($T)) do
        vtx = ($T)[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
        tetcub = TetSymCub{($T)}()
        @fact SymCubatures.calcnodes(tetcub, vtx) => vtx[:,1], vtx[:,2], vtx[:,3]

        tetcub = TetSymCub{($T)}(numedge=1)
        alpha = ($T)(0.25)
        A = ($T)[alpha (1-alpha) 0 0;
                 (1-alpha) alpha 0 0;
                 0 alpha (1-alpha) 0;
                 0 (1-alpha) alpha 0;
                 0 0 alpha (1-alpha);
                 0 0 (1-alpha) alpha;
                 (1-alpha) 0 0 alpha;
                 alpha 0 0 (1-alpha);
                 alpha 0 (1-alpha) 0;
                 (1-alpha) 0 alpha 0;
                 0 alpha 0 (1-alpha)
                 0 (1-alpha) 0 alpha]
        SymCubatures.setparams!(tetcub, [alpha])
        x, y, z = SymCubatures.calcnodes(tetcub, vtx)
        @fact x.' => [vtx[:,1].' (A*vtx[:,1]).']
        @fact y.' => [vtx[:,2].' (A*vtx[:,2]).']
        @fact z.' => [vtx[:,3].' (A*vtx[:,3]).']

        tetcub = TetSymCub{($T)}(numS31=1)
        SymCubatures.setparams!(tetcub, ($T)[3/4])
        x, y, z = SymCubatures.calcnodes(tetcub, vtx)
        @fact x.' => [vtx[:,1].' fill(sum(vtx[:,1])/4, (1,4) )]
        @fact y.' => [vtx[:,2].' fill(sum(vtx[:,2])/4, (1,4) )]
        @fact z.' => [vtx[:,3].' fill(sum(vtx[:,3])/4, (1,4) )]
      end
    end
  end

  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing calcweights (TriSymCub method) for DataType "string($T)) do
        tricub = TriSymCub{($T)}()
        wgt = ($T)(1/3)
        SymCubatures.setweights!(tricub, [wgt])
        @fact SymCubatures.calcweights(tricub) => roughly(fill(wgt, (3)))

        tricub = TriSymCub{($T)}(midedges=true, numedge=1, numS21=2)
        w = ($T)[1/3 1/4 1/5 1/6 1/7]
        SymCubatures.setweights!(tricub, w)
        @fact SymCubatures.calcweights(tricub) =>
        roughly([w[1]*ones(($T), (3))
                 w[2]*ones(($T), (3))
                 w[3]*ones(($T), (6))
                 w[4]*ones(($T), (3))
                 w[5]*ones(($T), (3))], atol=1e-15)
      end
    end
  end

  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing calcweights (TetSymCub method) for DataType "string($T)) do
        tetcub = TetSymCub{($T)}()
        wgt = ($T)(1/3)
        SymCubatures.setweights!(tetcub, [wgt])
        @fact SymCubatures.calcweights(tetcub) => roughly(fill(wgt, (4)))

        tetcub = TetSymCub{($T)}(midedges=true, numedge=1, numS31=2)
        w = ($T)[1/3 1/4 1/5 1/6 1/7]
        SymCubatures.setweights!(tetcub, w)
        @fact SymCubatures.calcweights(tetcub) =>
        roughly([w[1]*ones(($T), (4))
                 w[2]*ones(($T), (6))
                 w[3]*ones(($T), (12))
                 w[4]*ones(($T), (4))
                 w[5]*ones(($T), (4))], atol=1e-15)
      end
    end
  end

  context("Testing calcjacobianofnodes (TriSymCub method)") do
    # loop over parameters and check Jacobian using complex step
    vtx = Float64[-1 -1; 1 -1; -1 1]
    tricub = TriSymCub{Float64}(midedges=true, numedge=2, numS21=1)
    SymCubatures.setparams!(tricub, [1/3, 2/3, 1/4])
    Jac = SymCubatures.calcjacobianofnodes(tricub, vtx)
    Jac_cs = zeros(Jac)
    tricub_cmplx = TriSymCub{Complex128}(midedges=true, numedge=2, numS21=1)
    params_cmplx = tricub.params + 0im
    eps_step = 1e-60
    for i = 1:tricub.numparams
      params_cmplx[i] += eps_step*im
      SymCubatures.setparams!(tricub_cmplx, params_cmplx)
      xc, yc = SymCubatures.calcnodes(tricub_cmplx,
                                      convert(Array{Complex128}, vtx))
      Jac_cs[1:tricub.numnodes,i] = imag(xc)/eps_step
      Jac_cs[tricub.numnodes+1:2*tricub.numnodes,i] = imag(yc)/eps_step
      params_cmplx[i] -= eps_step*im
    end
    @fact Jac => roughly(Jac_cs, atol=1e-15)
  end

  context("Testing calcjacobianofnodes (TetSymCub method)") do
    # loop over parameters and check Jacobian using complex step
    vtx = Float64[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
    tetcub = TetSymCub{Float64}(midedges=true, numedge=2, numS31=1)
    SymCubatures.setparams!(tetcub, [1/3, 2/3, 1/4])
    Jac = SymCubatures.calcjacobianofnodes(tetcub, vtx)
    Jac_cs = zeros(Jac)
    tetcub_cmplx = TetSymCub{Complex128}(midedges=true, numedge=2, numS31=1)
    params_cmplx = tetcub.params + 0im
    eps_step = 1e-60
    for i = 1:tetcub.numparams
      params_cmplx[i] += eps_step*im
      SymCubatures.setparams!(tetcub_cmplx, params_cmplx)
      xc, yc, zc = SymCubatures.calcnodes(tetcub_cmplx,
                                          convert(Array{Complex128}, vtx))
      Jac_cs[1:tetcub.numnodes,i] = imag(xc)/eps_step
      Jac_cs[tetcub.numnodes+1:2*tetcub.numnodes,i] = imag(yc)/eps_step
      Jac_cs[2*tetcub.numnodes+1:3*tetcub.numnodes,i] = imag(zc)/eps_step
      params_cmplx[i] -= eps_step*im
    end
    @fact Jac => roughly(Jac_cs, atol=1e-15)
  end

  context("Testing calcjacobianofweights (TriSymCub method)") do
    # loop over weights and check Jacobian using complex step
    tricub = TriSymCub{Float64}(midedges=true, numedge=1, numS21=2)
    w = Float64[1/3 1/4 1/5 1/6 1/7]
    SymCubatures.setweights!(tricub, w)
    Jac = SymCubatures.calcjacobianofweights(tricub)
    Jac_cs = zeros(Jac)
    tricub_cmplx = TriSymCub{Complex128}(midedges=true, numedge=1, numS21=2)
    weights_cmplx = tricub.weights + 0im
    eps_step = 1e-60
    for i = 1:tricub.numweights
      weights_cmplx[i] += eps_step*im
      SymCubatures.setweights!(tricub_cmplx, weights_cmplx)
      wc = SymCubatures.calcweights(tricub_cmplx)
      Jac_cs[:,i] = imag(wc)/eps_step
      weights_cmplx[i] -= eps_step*im
    end
    @fact Jac => roughly(Jac_cs, atol=1e-15)
  end

  context("Testing calcjacobianofweights (TetSymCub method)") do
    # loop over weights and check Jacobian using complex step
    tetcub = TetSymCub{Float64}(midedges=true, numedge=1, numS31=2)
    w = Float64[1/3 1/4 1/5 1/6 1/7]
    SymCubatures.setweights!(tetcub, w)
    Jac = SymCubatures.calcjacobianofweights(tetcub)
    Jac_cs = zeros(Jac)
    tetcub_cmplx = TetSymCub{Complex128}(midedges=true, numedge=1, numS31=2)
    weights_cmplx = tetcub.weights + 0im
    eps_step = 1e-60
    for i = 1:tetcub.numweights
      weights_cmplx[i] += eps_step*im
      SymCubatures.setweights!(tetcub_cmplx, weights_cmplx)
      wc = SymCubatures.calcweights(tetcub_cmplx)
      Jac_cs[:,i] = imag(wc)/eps_step
      weights_cmplx[i] -= eps_step*im
    end
    @fact Jac => roughly(Jac_cs, atol=1e-15)
  end

  context("Testing calcjacobian (TriSymCub method)") do
    # This is essentially just a copy of the code in calcjacobian, so not good test
    vtx = Float64[-1 -1; 1 -1; -1 1]
    tricub = TriSymCub{Float64}(midedges=true, numedge=1, numS21=1)
    SymCubatures.setparams!(tricub, [1/3, 2/3])
    SymCubatures.setweights!(tricub, [1/3 1/4 1/5 1/6])
    Jac_params = SymCubatures.calcjacobianofnodes(tricub, vtx)
    Jac_weights = SymCubatures.calcjacobianofweights(tricub)
    Jac = SymCubatures.calcjacobian(tricub, vtx)
    numnodes = tricub.numnodes
    numparams = tricub.numparams
    numweights = tricub.numweights
    @fact Jac[1:2*numnodes,1:numparams] => roughly(Jac_params, atol=1e-15)
    @fact Jac[1:2*numnodes,numparams+1:end] =>
    roughly(zeros(Float64, (2*numnodes, numweights)), atol=1e-15)
    @fact Jac[2*numnodes+1:end,1:numparams] =>
    roughly(zeros(Float64, (numnodes, numparams)), atol=1e-15)
    @fact Jac[2*numnodes+1:end,numparams+1:end] => roughly(Jac_weights, atol=1e-15)
  end

  context("Testing calcjacobian (TetSymCub method)") do
    # This is essentially just a copy of the code in calcjacobian, so not good test
    vtx = Float64[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
    tetcub = TetSymCub{Float64}(midedges=true, numedge=1, numS31=1)
    SymCubatures.setparams!(tetcub, [1/3, 2/3])
    SymCubatures.setweights!(tetcub, [1/3 1/4 1/5 1/6])
    Jac_params = SymCubatures.calcjacobianofnodes(tetcub, vtx)
    Jac_weights = SymCubatures.calcjacobianofweights(tetcub)
    Jac = SymCubatures.calcjacobian(tetcub, vtx)
    numnodes = tetcub.numnodes
    numparams = tetcub.numparams
    numweights = tetcub.numweights
    @fact Jac[1:3*numnodes,1:numparams] => roughly(Jac_params, atol=1e-15)
    @fact Jac[1:3*numnodes,numparams+1:end] =>
    roughly(zeros(Float64, (3*numnodes, numweights)), atol=1e-15)
    @fact Jac[3*numnodes+1:end,1:numparams] =>
    roughly(zeros(Float64, (numnodes, numparams)), atol=1e-15)
    @fact Jac[3*numnodes+1:end,numparams+1:end] => roughly(Jac_weights, atol=1e-15)
  end

end
