module Cubature
# routines for constructing cubatures

using ..OrthoPoly
using ..SymCubatures

export quadrature, tricubature, tetcubature

@doc """
### Cubature.cubatureresidual

This method computes the residuals, `F`, between a cubature, defined by `cub`,
and the true value of an integral.  Each residual corresponds with an orthogonal
polynomial on the simplex up to degree `q`.  The Jacobian, `dF`, of the
residual, with respect to the quadrature nodes and weights, is also returned.

**Inputs**

* `cub`: defines the nodes and weights of the cubature via symmetry orbits
* `q`: maximum degree of the othogonal polynomials used in the conditions

**Outputs**

* `F`: the accuracy conditions for orthogonal polynomials up to degree q
* `dF`: derivative of F with respect to x, y, (z,) w, in that order

"""->
function cubatureresidual{T}(cub::TriSymCub{T}, q::Int)
  # compute the nodes and weights defined by cub
  vtx = T[-1 -1; 1 -1; -1 1]
  x = SymCubatures.calcnodes(cub, vtx)
  w = SymCubatures.calcweights(cub)
  num_eq = convert(Int, (q+1)*(q+2)/2)
  # loop over orthogonal polynomials of degree r <= q and form conditions
  F = zeros(T, (num_eq) )
  F[1] = -2.0/sqrt(2.0)
  dF = zeros(T, (num_eq, 3*cub.numnodes) )
  ptr = 1
  for r = 0:q
    for j = 0:r
      i = r-j
      P = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      F[ptr] += (w.'*P)[1]
      dPdx, dPdy = OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      dF[ptr,:] = [w.*dPdx w.*dPdy P]
      ptr += 1
      #print("(i,j,k) = (",i,",",j,",",k,"): i+j+k = ",i+j+k,"\n")
    end
  end
  return F, dF 
end

function cubatureresidual{T}(cub::TetSymCub{T}, q::Int)
  # compute the nodes and weights defined by cub
  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
  x = SymCubatures.calcnodes(cub, vtx)
  w = SymCubatures.calcweights(cub)
  num_eq = convert(Int, (q+1)*(q+2)*(q+3)/6)
  # loop over orthogonal polynomials of degree r <= q and form conditions
  F = zeros(T, (num_eq) )
  F[1] = -2.0/sqrt(3.0)
  dF = zeros(T, (num_eq, 4*cub.numnodes) )
  ptr = 1
  for r = 0:q
    for k = 0:r
      for j = 0:r-k
        i = r-j-k
        P = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), vec(x[3,:]), i, j, k)
        F[ptr] += (w.'*P)[1]
        dPdx, dPdy, dPdz = OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]),
                                                     vec(x[3,:]), i, j, k)
        dF[ptr,:] = [w.*dPdx w.*dPdy w.*dPdz P]
        ptr += 1
        #print("(i,j,k) = (",i,",",j,",",k,"): i+j+k = ",i+j+k,"\n")
      end
    end
  end
  return F, dF 
end

@doc """
### Cubature.solvecubature!{SymCub{T}}
  
Attempts to solve for the nodes and weights of a cubature that is exact for
polynomials of degree r <= `q`.  The nodes and weights of the cubature are
defined by `cub`, which is a parametric abstract type (see symcubatures.jl).

**Inputs**

* `q`: maximum (desired) degree for which the cubature is exact
* `mask`: array of indicies of parameters and weights that are free
* `tol`: tolerance with which to solve the accuracy conditions
* `hist`: if true, print the residual-norm convergence history

**In/Outs**

* `cub`: on entry, defines the initial guess for the cubature nodes and weights.
  on exit, defines the nodes and weights that satisfy the desired accuracy.

"""->
function solvecubature!{T}(cub::SymCub{T}, q::Int, mask::AbstractArray{Int64,1};
                           tol=10*eps(typeof(real(one(T)))),
                           hist::Bool=false)
  @assert( length(mask) <= cub.numparams + cub.numweights )
  Jac = SymCubatures.calcjacobian(cub)

  # compute accuracy for initial guess 
  F, dF = Cubature.cubatureresidual(cub, q)
  res = norm(F)
  res0 = res
  res_old = res
  hist ? print("solvecubature!:\n") : nothing
  hist ? print("\titer ",0,": res norm = ",res,"\n") : nothing
  if (res < tol)
    return
  end

  # Levenberg–Marquardt loop
  maxiter = 200
  nu = 1000.0 #100.0
  v = zeros(T, (cub.numparams + cub.numweights) )
  v[1:cub.numparams] = cub.params
  v[cub.numparams+1:end] = cub.weights
  dv = zeros(v)
  for k = 1:maxiter
    JtJ = Jac.'*dF.'*dF*Jac
    H = JtJ + nu*diagm(diag(JtJ))
    g = -Jac.'*dF.'*F

    # solve only for those parameters and weights that are in mask
    fill!(dv, zero(T))
    Hred = H[mask,mask]
    dv[mask] = Hred\g[mask]

    # update cubature definition and check for convergence
    v += dv
    SymCubatures.setparams!(cub, v[1:cub.numparams])
    SymCubatures.setweights!(cub, v[cub.numparams+1:end])
    F, dF = Cubature.cubatureresidual(cub, q)
    res = norm(F)
    hist ? print("\titer ",k,": res norm = ",res,"\n") : nothing
    if res < tol
      #println("size(JtJ) = ",size(JtJ))
      #println("rank(JtJ) = ",rank(JtJ))
      return
    end

    # trust-region like update
    if res > res_old
      v -= dv
      SymCubatures.setparams!(cub, v[1:cub.numparams])
      SymCubatures.setweights!(cub, v[cub.numparams+1:end])
      F, dF = Cubature.cubatureresidual(cub, q)
      nu *= 4.0
    else
      nu /= 2.0
      res_old = res
    end

  end
  error("solvecubature failed to find solution in ",maxiter," iterations")
end

@doc """
### Cubature.solvecubatureweights!{SymCub{T}}
  
Attempts to solve for the weights of a cubature that is exact for
polynomials of degree r <= `q`.  The weights (and nodes) of the cubature are
defined by `cub`, which is a parametric abstract type (see symcubatures.jl).

**Inputs**

* `q`: maximum (desired) degree for which the cubature is exact
* `tol`: tolerance with which to solve the accuracy conditions
* `hist`: if true, print the residual-norm convergence history

**In/Outs**

* `cub`: on entry, defines the initial guess for the cubature nodes and weights.
  on exit, defines the nodes and weights that satisfy the desired accuracy.

"""->
function solvecubatureweights!{T}(cub::SymCub{T}, q::Int;
                                  tol=10*eps(typeof(real(one(T)))),
                                  hist::Bool=false)
  Jac = SymCubatures.calcjacobianofweights(cub)

  # compute accuracy for initial guess 
  F, dF = Cubature.cubatureresidual(cub, q)
  res = norm(F)
  res0 = res
  res_old = res
  hist ? print("solvecubatureweights!:\n") : nothing
  hist ? print("\titer ",0,": res norm = ",res,"\n") : nothing
  if (res < tol)
    return
  end
  
  # Levenberg–Marquardt loop
  maxiter = 200
  nu = 1000.0 #100.0
  v = zeros(T, (cub.numweights) )
  v = cub.weights
  for k = 1:maxiter
    JtJ = Jac.'*(dF[:,end-cub.numnodes+1:end].'*
      dF[:,end-cub.numnodes+1:end])*Jac
    H = JtJ + nu*diagm(diag(JtJ))
    g = -Jac.'*dF[:,end-cub.numnodes+1:end].'*F
    dv = H\g

    # update cubature definition and check for convergence
    v += dv
    SymCubatures.setweights!(cub, v)
    F, dF = Cubature.cubatureresidual(cub, q)
    res = norm(F)
    hist ? print("\titer ",k,": res norm = ",res,"\n") : nothing
    if res < tol
      #println("size(JtJ) = ",size(JtJ))
      #println("rank(JtJ) = ",rank(JtJ))
      return
    end

    # trust-region like update
    if res > res_old
      v -= dv
      SymCubatures.setweights!(cub, v)
      F, dF = Cubature.cubatureresidual(cub, q)
      nu *= 4.0
    else
      nu /= 2.0
      res_old = res
    end

  end
  error("solvecubatureweights failed to find solution in ",maxiter," iterations")
end

@doc """
### Cubature.quadrature{T}

This high-level function computes and returns a symmetric cubature of requested
accuracy on the interval [-1,1]

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `internal`: if true, all nodes are strictly internal (default false)
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the interval [-1,1]
* `vtx`: vertices, [-1,1]

"""->
function quadrature(q::Int, T=Float64; internal::Bool=false)
  if internal
    # all nodes are internal (LG quadrature)
    N = div(q+2,2)
    x, w = OrthoPoly.lgnodes(N, T)
    alpha = zeros(T, (div(N,2)))
    weight = zeros(T, (div(N+1,2)))
    for i = 1:div(N,2)
      alpha[i] = (1 + x[i])/2
      weight[i] = w[i]
    end
    if rem(N,2) == 0
      centroid = false
    else
      centroid = true
      weight[div(N+1,2)] = w[div(N+1,2)]
    end
    quad = LineSymCub{T}(vertices=false, centroid=centroid,
                         numedge=div(N,2))
    SymCubatures.setparams!(quad, alpha)
    SymCubatures.setweights!(quad, weight)
  else
    # vertices are included (LGL quadrature)
    N = div(q+4,2)
    x, w = OrthoPoly.lglnodes(N-1, T)
    alpha = zeros(T, (div(N-2,2)))
    weight = zeros(T, (div(N+1,2)))
    weight[1] = w[1]
    for i = 1:div(N-2,2)
      alpha[i] = (1 - x[i+1])/2
      weight[i+1] = w[i+1]
    end
    if rem(N,2) == 0 
      centroid=false
    else
      centroid=true
      weight[div(N+1,2)] = w[div(N+1,2)]
    end
    quad = LineSymCub{T}(vertices=true, centroid=centroid,
                         numedge=div(N-2,2))
    SymCubatures.setparams!(quad, alpha)
    SymCubatures.setweights!(quad, weight)
  end
  vtx = reshape(T[-1; 1], (2,1))
  return quad, vtx
end

@doc """
### Cubature.tricubature{T}

This high-level function computes and returns a symmetric cubature of requested
accuracy on the right triangle.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `internal`: if true, all nodes are strictly internal (default false)
* `vertices`: if true and `internal` is false, then vertices are not included
* `facequad`: if true, the cubatures' face nodes coincide with a quadrature
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""->
function tricubature(q::Int, T=Float64; internal::Bool=false,
                     vertices::Bool=true, facequad::Bool=false,
                     tol=10*eps(typeof(real(one(T)))))
  cub_degree = q
  mask = zeros(Int64, (0))
  if facequad
    # face nodes coincide with a quadrature rule
    if q <= 1
      cub = SymCubatures.TriSymCub{T}(vertices=true, midedges=true)
      #cub = SymCubatures.TriSymCub{T}(vertices=false, numedge=1)
      SymCubatures.setweights!(cub, T[1/3; 1/3])
      #SymCubatures.setweights!(cub, T[1/3])
      #SymCubatures.setparams!(cub, T[0.5*(1 + 1/sqrt(3))])
      cub_degree = 1
      #weights_only = true
      mask = SymCubatures.getInternalParamMask(cub)
      append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
    elseif q <= 3
      cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=1,
                                      centroid=true)
      SymCubatures.setweights!(cub, T[1/5; 1/5; 1/5])
      SymCubatures.setparams!(cub, T[0.5*(1 + 1/sqrt(5))])
      cub_degree = 3 
      mask = SymCubatures.getInternalParamMask(cub)
      append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
    elseif q <= 5
      cub = SymCubatures.TriSymCub{T}(vertices=true, midedges=true, numedge=1,
                                      numS21=1)
      SymCubatures.setweights!(cub, T[2/15; 2/15; 2/15; 2/15])
      SymCubatures.setparams!(cub, T[0.5; 0.5*(1 + sqrt(3/7))])
      cub_degree = 5
      mask = SymCubatures.getInternalParamMask(cub)
      append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
    elseif q <= 7
      # this does not work yet
      cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=2, numS21=2)
      SymCubatures.setweights!(cub, T[2/21; 2/21; 2/21; 2/21; 2/21])
      SymCubatures.setparams!(cub, T[0.25; 0.75;
                                     0.5*(1 + sqrt(1/3 + 2*sqrt(7)/21));
                                     0.5*(1 + sqrt(1/3 - 2*sqrt(7)/21))])
      cub_degree = 7
      mask = SymCubatures.getInternalParamMask(cub)
      append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))      
    end
  elseif internal
    # all nodes are internal
    if q <= 2
      # P1; 3rd order cubature
      cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=1)      
      SymCubatures.setweights!(cub, T[2/3])
      SymCubatures.setparams!(cub, T[1/3])
      cub_degree = 2
    elseif q <= 4
      # P2; 5th order cubature
      cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=2)
      SymCubatures.setweights!(cub, T[0.44676317935602283;
                                      0.2199034873106437])
      SymCubatures.setparams!(cub, T[0.8918969818319298;
                                     0.18315242701954149])
      cub_degree = 4
    elseif q <= 5
      # P3; 6th order cubature
      cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=true,
                                      numS21=1, numS111=1)
      SymCubatures.setweights!(cub, T[0.11550472674301035;
                                      0.20924480696331949;
                                      0.39801697799105223])
      SymCubatures.setparams!(cub, T[0.13862330627662678;
                                     0.14215944055500324;
                                     0.6226442585632832])
      cub_degree = 5
    elseif q <= 6
      # P3; 7th order cubature
      cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=false,
                                      numS21=2, numS111=1)
      SymCubatures.setweights!(cub, T[0.1016898127404136;
                                      0.23357255145275863;
                                      0.16570215123674722])
      SymCubatures.setparams!(cub, T[0.1261780289830045;
                                     0.49857349034182097;
                                     0.10629009968963399;
                                     0.6207049020675687])
      cub_degree = 6
    elseif q <= 7
      # P4; 8th order cubature
      cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=3, numS111=1)
      SymCubatures.setweights!(cub, T[0.045386157905236965;
                                      0.1458284149509071;
                                      0.2543369199180239;
                                      0.11055758694624956])
      SymCubatures.setparams!(cub, T[0.08433122881886425;
                                     0.9485893782350207;
                                     0.4841719475189572;
                                     0.4231241172761849;
                                     0.0959626827429292])
      cub_degree = 7
      # JEH The commented version below has much worse conditioning/CFL limit;
      # leaving it here for reference, and to avoid
      # SymCubatures.setweights!(cub, T[0.10482661091570668;
      #                                 0.2253930198733382;
      #                                 0.057547518977195254;
      #                                 0.13944975845021326])
      # SymCubatures.setparams!(cub, T[0.1290634461434249;
      #                                0.4731163893279408;
      #                                0.8413468069012109;
      #                                0.08815074437486997;
      #                                0.624003943088726])
    elseif q <= 8
      # P4; 9th order cubature
      cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=4, numS111=1)
      SymCubatures.setweights!(cub, T[0.05415768359243055;
                                      0.11606844251884181;
                                      0.16305903205856434;
                                      0.18042409969012593;
                                      0.07647870440335205])
      SymCubatures.setparams!(cub, T[0.09199464041197784;
                                     0.9547402313677645;
                                     0.353676230440559;
                                     0.8037465193903021;
                                     0.4612850345245523;
                                     0.06105167151116454])
      cub_degree = 8
    else
      error("polynomial degree must be <= 7 (presently)\n")
    end
    mask = 1:(cub.numparams+cub.numweights)
  else
    if vertices
      # at least (q+1)/2+1 nodes along each edge
      @assert( q >= 1 && q <= 9 && mod(q,2) == 1 )
      if q <= 1
        # P1 (vertices only); 2nd order cubature
        cub = SymCubatures.TriSymCub{T}(vertices=true) 
        SymCubatures.setweights!(cub, T[2/3])
        cub_degree = 1
      elseif q <= 3
        # P2 + 1 bubble node; 4th order cubature
        cub = SymCubatures.TriSymCub{T}(midedges=true, centroid=true)
        SymCubatures.setweights!(cub, T[9/10, 1/10, 4/15])
        cub_degree = 3
      elseif q <= 5
        # P3 + 3 bubble nodes; 6th order cubature
        cub = SymCubatures.TriSymCub{T}(numedge=1, numS21=1)
        SymCubatures.setweights!(cub, T[0.02974582604964118,0.4415541156808217,
                                        0.09768336246810204])
        SymCubatures.setparams!(cub, T[0.41469035132718185,0.29346955590904017])
        cub_degree = 5
      elseif q <= 7
        # P4 + 6 bubble nodes; 8th order cubature
        cub = SymCubatures.TriSymCub{T}(midedges=true, numedge=1, numS21=2)
        SymCubatures.setweights!(cub, T[0.012698412698412695,0.05079365079365077,
                                        0.2023354595827503,0.3151248578775673,
                                        0.04285714285714284])
        SymCubatures.setparams!(cub, T[0.2615831876594899,0.8495279234516212,
                                       0.2113248654051872])
        cub_degree = 7
      elseif q <= 9
        # P5 + 10 bubble nodes; 10th order cubature
        cub = SymCubatures.TriSymCub{T}(numedge=2, centroid=false, numS21=2,
                                        numS111=1)
        SymCubatures.setweights!(cub, T[0.005060870857201095,0.1656605522739661,
                                        0.11124762548151651,0.02601835402818015,
                                        0.022007571591738946,0.1443228834070724])
        SymCubatures.setparams!(cub, T[0.5264875797340474,0.2020312767621901,
                                       0.3647863788168577,0.12582399442561498,
                                       0.17313630713608186,0.6472196801547492])
        cub_degree = 9
      elseif q <= 13
        # P7; 14th order cubature
        cub = SymCubatures.TriSymCub{T}(numedge=3, centroid=false, numS21=3,
                                        numS111=3)
        SymCubatures.setweights!(cub, T[0.0013434826332230758,0.07754170837489897,
                                        0.0392937103109862,0.08132263825474927,
                                        0.009447719133224907,0.011139320053630379,
                                        0.007018823640551441,0.055270427044833675,
                                        0.06131670240670696,0.08938957126745721])
        SymCubatures.setparams!(cub, T[0.5749380303918797,0.11974220085793316,
                                       0.33343014155598055,0.20437477140696825,
                                       0.39432452613154606,0.06536359154212519,
                                       0.38875869465672286,0.10195028749023816,
                                       1.1491810826793598,0.09609353164480232,
                                       0.6657786329998556,1.0308822535578346])
        cub_degree = 13
      else
        error("polynomial degree must be <= 9 (presently)\n")
      end
    else 
      # do not include vertices in the cubature rule
      if q <= 2
        # P1 (faces only); 2nd order cubature
        cub = SymCubatures.TriSymCub{T}(vertices=false, midedges=true) 
        SymCubatures.setweights!(cub, T[2/3])
        cub_degree = 2
      elseif q <= 3
        # 2 edge nodes + 1 bubble node; 4th order cubature
        cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=true, numedge=1)
        SymCubatures.setweights!(cub, T[0.18333333333333333,9/10])
        SymCubatures.setparams!(cub, T[0.23888351606645322])
        cub_degree = 3
      elseif q <= 5
        # 3 edge nodes + 3 interior nodes; 6th order cubature
        cub = SymCubatures.TriSymCub{T}(vertices=false, midedges=true, numedge=1,
                                        numS21=1)
        SymCubatures.setweights!(cub, T[0.11281708967460279;
                                        0.44155411568082154;
                                        0.056147730655621154])
        SymCubatures.setparams!(cub, T[0.41469035132718185;
                                       0.12525844798726105])
        cub_degree = 5        
      elseif q <= 7
        # 4 edge nodes + 6 interior nodes; 8th order cubature
        cub = SymCubatures.TriSymCub{T}(vertices=false, numedge=2, numS21=2)
        SymCubatures.setweights!(cub, T[0.20233545958275004;
                                        0.31512485787756733;
                                        0.017098156200124687;
                                        0.05750501840304983])
        SymCubatures.setparams!(cub, T[0.2615831876594899;
                                       0.8495279234516212;
                                       0.06120727790001947;
                                       0.6801692177706541])
        cub_degree = 7                                 
      else
        error("polynomial degree must be <= 7 (presently)\n")
      end
    end  
    mask = 1:(cub.numparams+cub.numweights)
  end
  vtx = T[-1 -1; 1 -1; -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol)
  return cub, vtx
end

@doc """
### Cubature.tetcubature{T}

This high-level function computes and returns a symmetric cubature of requested
accuracy on the right tetrahedron.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `internal`: if true, all nodes are strictly internal (default false)
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right tetrahedron
* `vtx`: vertices for the right tetrahedron

"""->
function tetcubature(q::Int, T=Float64; internal::Bool=false,
                     tol=10*eps(typeof(real(one(T)))))
  mask = zeros(Int64, (0))
  if internal
    # all nodes are internal
    @assert( q >= 1 )
    if q <= 2
      # P1; 3rd order cubature
      cub = SymCubatures.TetSymCub{T}(vertices=false, numS31=1)
      SymCubatures.setweights!(cub, T[1/3])
      SymCubatures.setparams!(cub, T[(1 - sqrt(5)/5)*3/4])
    elseif q <= 3
      # P2; 4th order cubature
      cub = SymCubatures.TetSymCub{T}(vertices=false, numS31=1, numS22=1)
      SymCubatures.setweights!(cub, T[0.06483158243276162;
                                      0.17900116726703835])
      SymCubatures.setparams!(cub, T[0.22511815489558668;
                                     0.18771315212883505])
      #SymCubatures.setweights!(cub, T[0.1302091416313459;
      #                                0.13541612780132486])
      #SymCubatures.setparams!(cub, T[0.33398409622579817;
      #                               0.18658191164952043])
    elseif q <= 5
      # P3; 6th order cubature
      cub = SymCubatures.TetSymCub{T}(vertices=false, numS31=2, numS211=1)
      #SymCubatures.setweights!(cub, T[0.061630217648090097;
      #                                0.13793513058238085;
      #                                0.04458932836762084])
      #SymCubatures.setparams!(cub, T[0.24722530396402584;
      #                               0.9298082909679131;
      #                               0.11664936229736803;
      #                               0.6505900754758551])
      SymCubatures.setweights!(cub, T[0.1330522438256425;
                                      0.027128275352375174;
                                      0.05771760471843853])
      SymCubatures.setparams!(cub, T[0.9291909288071043;
                                     0.18740453102806648;
                                     0.12393623056047982;
                                     0.5651871191159269])
    elseif q <= 7
      # P4; 8th order cubature
      cub = SymCubatures.TetSymCub{T}(vertices=false, numS31=2, numS22=1,
                                      numS211=2)
      # minimum node distance for this one is 0.17819
      # SymCubatures.setweights!(cub, T[0.02832965568227839;
      #                                 0.048583147669377845;
      #                                 0.039177188602071006;
      #                                 0.03021820473246459;
      #                                 0.03566671096039236])
      # SymCubatures.setparams!(cub, T[0.18167711419341304;
      #                                0.5398647398205032;
      #                                0.7170540544966304;
      #                                0.0881323679975843;
      #                                0.5992257377201948;
      #                                0.4688384403943167;
      #                                1.0098301020743294])
      # minimum node distance for this one is 0.17937
      # SymCubatures.setweights!(cub, T[0.05138889021100641;
      #                                 0.028389644845730116;
      #                                 0.03630699901869689;
      #                                 0.03614773635849856;
      #                                 0.030217030224352005])
      # SymCubatures.setparams!(cub, T[0.5465722016176291;
      #                                0.18183572000191664;
      #                                0.7209450878473552;
      #                                0.46859987367504574;
      #                                0.0537103011531615;
      #                                0.08816163364683494;
      #                                0.5991524190716274])
      # minimum node distance for this one is 0.2000058
      SymCubatures.setweights!(cub, T[0.0680746674820208;
                                      0.028423617986180167;
                                      0.018310980956037618;
                                      0.025461800630426842;
                                      0.044327724846598505])
      SymCubatures.setparams!(cub, T[0.5821093910011628;
                                     0.18262906602002502;
                                     0.17262216834048155;
                                     0.08099388793344552;
                                     0.5908415591286749;
                                     0.4612405506278837;
                                     0.07057889678019784])
    end
    mask = 1:(cub.numparams+cub.numweights)
  else
    # at least (q+1)/2+1 nodes along each edge
    @assert( q >= 1 && q <= 7 && mod(q,2) == 1)
    if q <= 1
      # P1 (vertices only); 2nd order cubature
      cub = SymCubatures.TetSymCub{T}()
      SymCubatures.setweights!(cub, T[1/3])
    elseif q <= 3
      # P2 + 1 bubble node; 4th order cubature
      cub = SymCubatures.TetSymCub{T}(midedges=true, centroid=true)
      SymCubatures.setweights!(cub, T[1/45 4/45 32/45])
    elseif q <= 5
      # P3 + 4 bubble nodes; 6th order cubature
      cub = SymCubatures.TetSymCub{T}(facecentroid=true, numedge=1, numS31=1)
      SymCubatures.setweights!(cub, T[0.004421633248304776 0.20653163611605146
                                      0.06935370366814568 0.0176754534336105])
      SymCubatures.setparams!(cub, T[0.45720884759834435 0.30480589839889616])
    elseif q <= 7
      # P3 + 11 bubble nodes; 8th order cubature
      cub = SymCubatures.TetSymCub{T}(midedges=true, centroid=true, numedge=1,
                                      numfaceS21=1, numS31=1, numS22=1)
      SymCubatures.setweights!(cub, T[0.0015106273303336273,0.060490542374353584,
                                      0.004038881996228382, 0.10344930834722398,
                                      0.005696088152131421, 0.02424296133613638,
                                      0.08113091859465722])
      SymCubatures.setparams!(cub, T[0.28418700275470193,0.21742832019555544,
                                     0.25737274681480826,0.45008848310824695])
    else
      error("polynomial degree must be 1, 3, 5, or 7 (presently)\n")
    end
    mask = 1:(cub.numparams+cub.numweights)
  end
  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
  Cubature.solvecubature!(cub, q, mask, tol=tol)
  return cub, vtx
end

@doc """
### Cubature.equivalenceconstant{T}

Computes the equivalence constant for a given cubature; that is, it finds the
maximum eigenvalue for the matrix pk^T H pm, where H = diag(weights) and pk
denotes the orthogonal polynomial evaluated at the cubature points.

**Inputs**

* `cub`: symmetric cubature rule
* `vtx`: vertices of the right simplex
* `q`: maximum degree of polynomial for which the cubature is to be tested

**Outputs**

* `λmax`: maximum eigenvalue, which is the equivalence constant

"""->
function equivalenceconstant{T}(cub::TriSymCub{T}, vtx::Array{T,2}, q::Int)
  N = convert(Int, (q+1)*(q+2)/2) 
  P = zeros(T, (cub.numnodes,N) )
  x = SymCubatures.calcnodes(cub, vtx)
  ptr = 1
  for r = 0:q
    for j = 0:r
      i = r-j
      P[:,ptr] = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
      ptr += 1
    end
  end
  H = diagm(SymCubatures.calcweights(cub))
  A = P.'*H*P
  return eigmax(0.5.*(A + A.'))
end  

function equivalenceconstant{T}(cub::TetSymCub{T}, vtx::Array{T,2}, q::Int)
  N = convert(Int, (q+1)*(q+2)*(q+3)/6 )
  P = zeros(T, (cub.numnodes,N) )
  x = SymCubatures.calcnodes(cub, vtx)
  ptr = 1
  for r = 0:q
    for k = 0:r
      for j = 0:r-k
        i = r-j-k
        P[:,ptr] = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), vec(x[3,:]),
                                         i, j, k)
        ptr += 1
      end
    end
  end
  H = diagm(SymCubatures.calcweights(cub))
  A = P.'*H*P
  return eigmax(0.5.*(A + A.'))
end

end
