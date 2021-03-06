module Cubature
# routines for constructing cubatures

using ..OrthoPoly
using ..SymCubatures

export pointCubature
export quadrature, quadratureUniform
export getTriCubatureGamma, getTriCubatureOmega, getTriCubatureDiagE
export getTetCubatureGamma, getTetCubatureOmega, getTetCubatureDiagE

"""
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

"""
function cubatureresidual{T}(cub::LineSymCub{T}, q::Int)
  # compute the nodes and weights defined by cub
  vtx = reshape(T[-1; 1], (2,1))
  x = SymCubatures.calcnodes(cub, vtx)
  w = SymCubatures.calcweights(cub)
  num_eq = convert(Int, q+1)
  # loop over orthogonal polynomials of degree r <= q and form conditions
  F = zeros(T, (num_eq) )
  F[1] = -2.0/sqrt(2.0)
  dF = zeros(T, (num_eq, 2*cub.numnodes) )
  ptr = 1
  for r = 0:q
    P = OrthoPoly.jacobipoly(vec(x[1,:]), 0.0, 0.0, r)
    F[ptr] += (w.'*P)[1]
    dPdx = OrthoPoly.diffjacobipoly(vec(x[1,:]), 0.0, 0.0, r)
    dF[ptr,:] = [w.*dPdx P]
    ptr += 1
  end
  return F, dF 
end

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

"""
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

"""
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
    #dv[mask] = pinv(Hred,1e-6)*g[mask]
    
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
      #nu *= 0.25
      res_old = res
    end

  end
  error("solvecubature failed to find solution in ",maxiter," iterations")
end

"""
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

"""
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

"""
### Cubature.pointCubature

This returns a (trivial) point cubature and default vertex -1

**Inputs**

* `T`: the data type used to represent the cubature

**Outputs**

* `cub`: a symmetric cubature for point
* `vtx`: vertex, [-1]

"""
function pointCubature(T::Type=Float64)
  pt = PointSymCub{T}()
  SymCubatures.setweights!(pt, T[1;])
  vtx = reshape(T[-1;], (1,1))
  return pt, vtx
end

"""
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

"""
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

"""
### Cubature.quadratureUniform{T}

This high-level function computes and returns a uniform cubature of requested
accuracy on the interval [-1,1]

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `N`: number of nodes (N >= q+1)
* `T`: the data type used to represent the cubature

**Outputs**

* `cub`: a symmetric cubature for the interval [-1,1]
* `vtx`: vertices, [-1,1]

"""
function quadratureUniform(q::Int, N::Int, T=Float64; internal::Bool=false)
  @assert(N >= q+1)

  if rem(N,2) == 0 
    centroid=false
  else
    centroid=true
  end
  numedge = div(N-2,2)
  quad = LineSymCub{T}(vertices=true, centroid=centroid,
                       numedge=numedge)
  alpha = zeros(T, (numedge))
  dx = 1.0/(N-1)
  for i = 1:numedge
    alpha[i] = i*dx 
  end    
  SymCubatures.setparams!(quad, alpha)
  weight = (2.0/N)*ones(T, (numedge + 1 + centroid))
  SymCubatures.setweights!(quad, weight)
  solvecubatureweights!(quad, q)
  #mask = zeros(Int64, (0))
  #append!(mask, (quad.numparams+1):(quad.numparams+quad.numweights))
  #Cubature.solvecubature!(quad, q, mask, tol=1e-15)

  vtx = reshape(T[-1; 1], (2,1))
  return quad, vtx
end

"""
### Cubature.tricubature{T}

Deprecated; this function will be removed in the future

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function tricubature(q::Int, T=Float64)
  return getTriCubatureGamma(q, T)
end

"""
### Cubature.getTriCubatureGamma{T}

Returns a cubature rule and vertices for the SBP Gamma operators on triangles;
these are operators with p+1 nodes on each face, where, typically, p =
(`q`+1)/2.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function getTriCubatureGamma(q::Int, T=Float64;
                             tol=10*eps(typeof(real(one(T)))))
  @assert( q >= 1 && q <= 9 && mod(q,2) == 1 )
  cub_degree = q
  mask = zeros(Int64, (0))
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
  mask = 1:(cub.numparams+cub.numweights)
  vtx = T[-1 -1; 1 -1; -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol)
  return cub, vtx
end

"""
### Cubature.getTriCubatureOmega{T}

Returns a cubature rule and vertices for the SBP Omega operators on triangles;
these are cubatures that are analogous to Gauss-Legendre in 1D, and they are
strictly internal to the triangle. 

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function getTriCubatureOmega(q::Int, T=Float64;
                             tol=10*eps(typeof(real(one(T)))))
  cub_degree = q
  mask = zeros(Int64, (0))
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
    cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=true,
                                    numS21=3, numS111=1)
    SymCubatures.setweights!(cub, T[0.20643474106943666;
                                    0.19018326853457126;
                                    0.06491699524639659;
                                    0.05446062834886685;
                                    0.2886312153555818])
    SymCubatures.setparams!(cub, T[0.3411386155035168;
                                   0.9185851765854467;
                                   0.10109445663405793;
                                   0.5262256592692812;
                                   0.016789554819910586])      
    # cub = SymCubatures.TriSymCub{T}(vertices=false, numS21=4, numS111=1)
    # SymCubatures.setweights!(cub, T[0.05415768359243055;
    #                                 0.11606844251884181;
    #                                 0.16305903205856434;
    #                                 0.18042409969012593;
    #                                 0.07647870440335205])
    # SymCubatures.setparams!(cub, T[0.09199464041197784;
    #                                0.9547402313677645;
    #                                0.353676230440559;
    #                                0.8037465193903021;
    #                                0.4612850345245523;
    #                                0.06105167151116454])
    cub_degree = 8
  else
    error("polynomial degree must be <= 8 (presently)\n")
  end
  mask = 1:(cub.numparams+cub.numweights)
  vtx = T[-1 -1; 1 -1; -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol)
  return cub, vtx
end

"""
### Cubature.getTriCubatureDiagE{T}

Returns a cubature rule and vertices for the SBP DiagE operators on triangles;
these are cubatures that have nodes on the boundary that correspond with LG or
LGL quadrature rules, which then leads to diagonal E.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `vertices`: if true then vertices are included
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function getTriCubatureDiagE(q::Int, T=Float64; vertices::Bool=true,
                             tol=10*eps(typeof(real(one(T)))))
  cub_degree = q
  mask = zeros(Int64, (0))
  if vertices
    # include the vertices in the cubature
    if q <= 2
      cub = SymCubatures.TriSymCub{T}(midedges=true, centroid=true)
      SymCubatures.setweights!(cub, T[9/10, 1/10, 4/15])
      cub_degree = 3
    elseif q <= 4
      cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=1,
                                      midedges=false,
                                      numS21=1, numS111=0,
                                      centroid=false)
      SymCubatures.setparams!(cub, T[0.4257087142201423;
                                     0.5*(1 + sqrt(1/5))])
      SymCubatures.setweights!(cub, T[0.02504450602156441;
                                      0.4270317586471588;
                                      0.10729520099967835])
      cub_degree = 4
    elseif q <= 6
      cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=1,
                                      midedges=true,
                                      numS21=2, numS111=0, centroid=false)
      SymCubatures.setparams!(cub, T[0.8487720503426771;
                                     0.28401016818370567;
                                     0.5*(1 + sqrt(3/7))])
      SymCubatures.setweights!(cub, T[0.00913026457472031;
                                      0.06201621056804736;
                                      0.2987067788024998;
                                      0.20607267197855683;
                                      0.04537037036896131])
      cub_degree = 6
    elseif q <= 8
      cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=2,
                                      numS21=2, numS111=1)
      SymCubatures.setparams!(cub, T[0.20735501628574252;
                                     0.5306627609507977;
                                     0.5*(1 + sqrt(1/3 - 2*sqrt(7)/21));
                                     0.5*(1 + sqrt(1/3 + 2*sqrt(7)/21));
                                     0.6492809445444747;
                                     1.1741711342683223])
      SymCubatures.setweights!(cub, T[0.004361575620524937;
                                      0.11358876265867929;
                                      0.15984019213585915;
                                      0.027867771483031517;
                                      0.02207899054885172;
                                      0.14449130609379232])
      cub_degree = 8
    elseif q <= 12
      cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=3,
                                      numS21=3, numS111=3)
      SymCubatures.setparams!(cub, T[0.5747680454804257; 0.34168799534894295;
                                     0.1194662933444393; 0.6046496089512394;
                                     0.7958500907165711; 0.9358700742548033;
                                     0.38458694986468334; 0.10666597060509266;
                                     0.3025043145989276; 0.6718411052897879;
                                     0.09573484681006857;0.7519048745731671])
      SymCubatures.setweights!(cub, T[0.001271601169161372; 0.07780339049594198;
                                      0.08154407642794283; 0.03881033381664769;
                                      0.01101979920649565; 0.010101834223603934;
                                      0.006972694458714173; 0.056746881379720164;
                                      0.08690355889320943; 0.06187386421672949])  
      cub_degree = 12
    elseif q <= 16
      cub = SymCubatures.TriSymCub{Float64}(vertices=true, numedge=4,
                                            numS21=4, numS111=6)
      SymCubatures.setparams!
      (cub, T[0.0768946752469594; 0.6109907336234316; 0.4179369130153705;
              0.23221732669622028; 0.5826394788331936; 0.7389624749052223;
              0.8693869325527526; 0.9597669540832294; 0.930330150896981;
              0.6679157686119799; 0.8075697058065031; 1.1286554029515519;
              0.07116357646128006; 0.25084609708056116; 0.20705942876211147;
              0.46901560437791967; 0.06341212586405608; 0.504346239131436;
      0.7497163430497378; 1.0452430326898021])
      SymCubatures.setweights!
      (cub, T[0.0004934174763938973; 0.016280379328615396; 0.03461636692031435;
              0.0544553219641883; 0.04016179760262383; 0.005790378245871778;
              0.005099447551845056; 0.004407182316786121; 0.002824162676005338;
              0.05536224506473959; 0.03347658309313654; 0.025582552500936707;
              0.047069550256364515; 0.029816087172401348; 0.050901502809179516])
      cub_degree = 16
    else
      error("polynomial degree must be <= 16 (presently)\n")
    end

  else
    # do not include vertices in the cubature
    if q <= 2
      cub = SymCubatures.TriSymCub{T}(vertices=false, numedge=1,
                                      centroid=true)
      SymCubatures.setweights!(cub, T[1/6; 1.0])
      SymCubatures.setparams!(cub, T[0.5*(1 + 1/sqrt(3))])
      cub_degree = 2
    elseif q <= 4
      cub = SymCubatures.TriSymCub{T}(vertices=false, numedge=1,
                                      midedges=true, numS21=1, numS111=0,
                                      centroid=false)
      SymCubatures.setparams!(cub, T[0.4093612831422611;
                                     0.5*(1 + sqrt(3/5))])
      SymCubatures.setweights!(cub, T[0.11680963211922607;
                                      0.449036908703613;
                                      0.050410062902983145])
      cub_degree = 4
    elseif q <= 6
      cub = SymCubatures.TriSymCub{T}(vertices=false, numedge=2,
                                      midedges=false,
                                      numS21=1, numS111=1, centroid=false)
      SymCubatures.setparams!(cub, T[0.8478748148895112; 
                                     0.5*(1 + sqrt(3/7 - 2/7*sqrt(6/5)));
                                     0.5*(1 + sqrt(3/7 + 2/7*sqrt(6/5)));
                                     0.3154166831592224;
                                     0.2189659885706149])
      SymCubatures.setweights!(cub, T[0.30870082526604714;
                                      0.0539978996110132;
                                      0.01943211076834772;
                                      0.10555291032094156])
      cub_degree = 6
    elseif q <= 8
      cub = SymCubatures.TriSymCub{T}(vertices=false, numedge=2,
                                      midedges=true,
                                      numS21=0, numS111=2, centroid=true)
      SymCubatures.setparams!(cub, T[0.5*(1 +(1/3)*sqrt(5 - 2*sqrt(10/7)));
                                     0.5*(1 +(1/3)*sqrt(5 + 2*sqrt(10/7)));
                                     0.22099843842186342;
                                     1.1631912073287645;
                                     0.24591837943530243;
                                     0.12484120739275503])
      SymCubatures.setweights!(cub, T[0.03946097492219484;
                                      0.022750835455651795;
                                      0.008818681441998117;
                                      0.18379350807584505;
                                      0.055866234615190226;
                                      0.2542415177018296])
      cub_degree = 8
    else
      error("polynomial degree must be <= 8 (presently)\n")
    end
  end
  mask = SymCubatures.getInternalParamMask(cub)
  append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
  vtx = T[-1 -1; 1 -1; -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol)
  return cub, vtx
end

function getTriCubatureSparse(q::Int, T=Float64;
                              tol=10*eps(typeof(real(one(T)))))
  @assert( mod(q,4) == 0)
  function getGuess(ξ,η)
    return (ξ - ξ*η/3 -1), (η - ξ*η/3 - 1)
  end  
  cub_degree = q
  mask = zeros(Int64, (0))
  # get the underlying LGL nodes
  xlgl, wlgl = OrthoPoly.lglnodes(div(q+2,2), T)
  xlgl[:] += 1.0
  xlgl[:] *= 0.5
  # figure out how many edge, S21, and S111 nodes
  numedge = div(q,4)
  numS21 = numedge
  numS111 = div((numedge-1)*numedge,2)
  params = zeros(T, (numedge+numS21+2*numS111))
  weights = zeros(T, (1+numedge+numS21+numS111))
  ptr = 1
  wptr = 1
  weights[wptr] = wlgl[1]*wlgl[1]
  wptr += 1
  # S21 parameter guesses
  for i = 1:numS21
    x, y = getGuess(2*xlgl[i+1],2*xlgl[i+1])
    params[ptr] = (x + 1)
    weights[wptr] = wlgl[i+1]*wlgl[i+1]
    ptr += 1
    wptr += 1
  end
  # edge parameters
  for i = 1:numedge
    params[ptr] = xlgl[i+1]
    weights[wptr] = wlgl[i+1]*wlgl[1]
    ptr += 1
    wptr += 1
  end
  # S111 parameters
  for i = 2:numedge
    for j = i:numedge
      x, y = getGuess(2*xlgl[i+1],2*xlgl[j+1])
      params[ptr] = (x + 1)
      params[ptr+1] = (y + 1)
      weights[wptr] = wlgl[i+1]*wlgl[j+1]
      ptr += 2
      wptr += 1
    end
  end
  cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=numedge,
                                  numS21=numS21, numS111=numS111)
  #weights[:] = 2./cub.numnodes
  SymCubatures.setweights!(cub, weights)

  wts = SymCubatures.calcweights(cub)
  weights[:] *= 2./sum(wts)

  weights[:] = 0.1/cub.numnodes
  
  SymCubatures.setweights!(cub, weights)
  println("sum wts = ",sum(SymCubatures.calcweights(cub)))

  SymCubatures.setparams!(cub, params)
  mask = zeros(Int64, (0))
  mask = SymCubatures.getInternalParamMask(cub)
  param_mask = zeros(Int64, (size(mask)))
  param_mask[:] = mask[:]
  append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
  vtx = T[-1 -1; 1 -1; -1 1]
  success = false
  count = 0
  while success == false
    count += 1
    if count > 20
      break
    end
    SymCubatures.setweights!(cub, 0.05*rand(cub.numweights)/cub.numnodes)
    rand_params = params
    rand_params[param_mask] .*= (1.0 + 0.1*randn(size(param_mask)))
    SymCubatures.setparams!(cub, rand_params)
    try Cubature.solvecubature!(cub, cub_degree, mask, tol=tol, hist=true)
    catch
      success = false
    end
    if success
      if minimum(cub.weights) < 0.0
        success = false
        continue
      end
      break
    end
  end
    
  #Cubature.solvecubatureweights!(cub, cub_degree, tol=tol, hist=true)
  return cub, vtx
end

"""
### Cubature.tetcubature{T}

This high-level function computes and returns a symmetric cubature of requested
accuracy on the right tetrahedron.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `internal`: if true, all nodes are strictly internal (default false)
* `facequad`: if true, the cubatures' face nodes coincide with a quadrature
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right tetrahedron
* `vtx`: vertices for the right tetrahedron

"""
function tetcubature(q::Int, T=Float64; internal::Bool=false,
                     facequad::Bool=false,
                     tol=10*eps(typeof(real(one(T)))))

end

"""
### Cubature.getTetCubatureGamma{T}

Returns a cubature rule and vertices for the SBP Gamma operators on tetrahedra;
these are operators with (p+1)(p+2)/2 nodes on each face, where, typically, p =
(`q`+1)/2.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right tetrahedron
* `vtx`: vertices for the right tetrahedron

"""
function getTetCubatureGamma(q::Int, T=Float64;
                             tol=10*eps(typeof(real(one(T)))))
  cub_degree = -1
  if q <= 1
    # P1 (vertices only); 2nd order cubature
    cub = SymCubatures.TetSymCub{T}()
    SymCubatures.setweights!(cub, T[1/3])
    cub_degree = 1
  elseif q <= 3
    # P2 + 1 bubble node; 4th order cubature
    cub = SymCubatures.TetSymCub{T}(midedges=true, centroid=true)
    SymCubatures.setweights!(cub, T[1/45 4/45 32/45])
    cub_degree = 3
  elseif q <= 5
    # P3 + 4 bubble nodes; 6th order cubature
    cub = SymCubatures.TetSymCub{T}(facecentroid=true, numedge=1, numS31=1)
    SymCubatures.setweights!(cub, T[0.004421633248304776 0.20653163611605146
                                    0.06935370366814568 0.0176754534336105])
    SymCubatures.setparams!(cub, T[0.45720884759834435 0.30480589839889616])
    cub_degree = 5
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
    cub_degree = 7
  else
    error("polynomial degree must be 1, 3, 5, or 7 (presently)\n")
  end
  mask = 1:(cub.numparams+cub.numweights)
  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol)
  return cub, vtx
end

"""
### Cubature.getTetCubatureOmega{T}

Returns a cubature rule and vertices for the SBP Omega operators on tetrahedra;
these are cubatures that are analogous to Gauss-Legendre in 1D, and they are
strictly internal to the tet. 

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right tetrahedron
* `vtx`: vertices for the right tetrahedron

"""
function getTetCubatureOmega(q::Int, T=Float64;
                             tol=10*eps(typeof(real(one(T)))))
  cub_degree = -1
  if q <= 2
    # P1; 3rd order cubature
    cub = SymCubatures.TetSymCub{T}(vertices=false, numS31=1)
    SymCubatures.setweights!(cub, T[1/3])
    SymCubatures.setparams!(cub, T[(1 - sqrt(5)/5)*3/4])
    cub_degree = 2
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
    cub_degree = 3
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
    cub_degree = 5
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
    cub_degree = 7
  end
  mask = 1:(cub.numparams+cub.numweights)
  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol)
  return cub, vtx
end

"""
### Cubature.getTetCubatureDiagE{T}

Returns a cubature rule and vertices for the SBP DiagE operators on tetrahedra;
these are cubatures that have nodes on the boundary that are analogous to LG or
LGL quadrature rules, which then leads to diagonal E.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `vertices`: if true then vertices are included
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right tetrahedron
* `vtx`: vertices for the right tetrahedron

"""
function getTetCubatureDiagE(q::Int, T=Float64; vertices::Bool=false,
                             tol=10*eps(typeof(real(one(T)))))
  cub_degree = q
  mask = zeros(Int64, (0))
  if vertices
    # include the vertices in the cubature
    if q <= 2
      # 14 nodes
      cub = SymCubatures.TetSymCub{T}(vertices=true, midedges=true,
                                      facecentroid=true)      
      cub_degree = 2
    else
      error("polynomial degree must be <= 2 (presently)\n")
    end
  else
    if q <= 2
      # 13 nodes
      cub = SymCubatures.TetSymCub{T}(vertices=false, numfaceS21=1,
                                      centroid=true)
      SymCubatures.setparams!(cub, T[1/3])
      SymCubatures.setweights!(cub, T[2/30, 8/15])
      cub_degree = 2
    elseif q <= 4
      # 36 nodes
      # the boundary nodes are strictly internal to the face
      cub = SymCubatures.TetSymCub{T}(vertices=false, numfaceS21=2, numS211=1)
      SymCubatures.setparams!(cub, T[0.8918969818319298;
                                     0.18315242701954149;
                                     0.40398819659496876;
                                     0.18710499145686446])
      SymCubatures.setweights!(cub, T[0.02922867858673424;
                                      0.012914918864852366;
                                      0.06896751365952453])
      cub_degree = 4
      tol = 1e-14
    elseif q <= 6
      # 69 nodes
      cub = SymCubatures.TetSymCub{Float64}(vertices=false, centroid=true,
                                            numS31=2, numS22=0, 
                                            numfaceS21=2, numS211=1,
                                            numfaceS111=1,
                                            numS1111=0)
      SymCubatures.setparams!(cub, T[0.21759254267395756; 0.9034356836878749;
                                     0.1261780289830045; 0.49857349034182097;
                                     0.16668051173596254; 0.5518039588837356;
                                     0.10629009968963399; 0.6207049020675687])
      SymCubatures.setweights!(cub, T[0.020370182220389683; 0.08061671634468645;
                                      0.0027518135350883413; 0.01252025357308984;
                                      0.04988851015482234; 0.005364555718326347;
                                      0.01870947467719025])
      cub_degree = 6
      tol = 1e-13
    elseif q <= 8

      # 99 nodes
      numS31 = 1       # 4  symmetry
      numS22 = 1       # 6  symmetry
      numfaceS21 = 3   # 12 symmetry
      numS211 = 2      # 12 symmetry
      numfaceS111 = 1  # 24 symmetry
      numS1111 = 0     # 24 symmetry
      cub = SymCubatures.TetSymCub{Float64}(vertices=false, centroid=true,
                                            facecentroid=true,
                                            numS31=numS31, numS22=numS22, 
                                            numfaceS21=numfaceS21,
                                            numS211=numS211,
                                            numfaceS111=numfaceS111,
                                            numS1111=numS1111)
      SymCubatures.setparams!(cub, T[0.1836596600107025,0.8275896865477325,0.3411386155035204,0.9185851765854463,0.10109445663406191,0.10888330398707555,0.4087055731374049,0.4521142310502033,0.1419370312058845,0.5262256592692762,0.0167895548199152])
      SymCubatures.setweights!(cub, T[0.011702201915075788,0.007003762166809546,0.02826881609586945,0.005737494407876313,0.009215129850744767,0.001683701003412158,0.019901864184903844,0.0455776334197571,0.0008896293665706435,0.08215560123254996])
      
#       # 111 nodes
#       numS31 = 1       # 4  symmetry
#       numS22 = 1       # 6  symmetry
#       numfaceS21 = 3   # 12 symmetry
#       numS211 = 1      # 12 symmetry
#       numfaceS111 = 1  # 24 symmetry
#       numS1111 = 1     # 24 symmetry
#       cub = SymCubatures.TetSymCub{Float64}(vertices=false, centroid=true,
#                                             facecentroid=true,
#                                             numS31=numS31, numS22=numS22, 
#                                             numfaceS21=numfaceS21,
#                                             numS211=numS211,
#                                             numfaceS111=numfaceS111,
#                                             numS1111=numS1111)
# cub.params = [0.24307201749459587,0.8435468793138858,0.3411386155035204,0.9185851765854463,0.10109445663406191,0.09363694100517951,0.43677794723997304,0.5262256592692762,0.0167895548199152,0.42744456237056916,0.14167234604995277,0.47242481400604164]
# cub.weights = [0.011504416930286117,0.015525438231567767,0.027098221864205605,0.006185807464370513,0.007989877875501223,0.0019370641049635973,0.015428153593239246,0.0007527586394225552,0.02427788489726483,0.0833983001632923]
      
      # # 116 nodes
      # cub = SymCubatures.TetSymCub{Float64}(vertices=false, centroid=false,
      #                                       facecentroid=true,
      #                                       numS31=numS31, numS22=numS22, 
      #                                       numfaceS21=numfaceS21,
      #                                       numS211=numS211,
      #                                       numfaceS111=numfaceS111,
      #                                       numS1111=numS1111)
      # SymCubatures.setparams!(cub, T[0.1663986770026873,0.3411386155035204,0.9185851765854463,0.10109445663406191,0.883044331845347,0.04075024844212266,0.18257544241564874,0.4787226396482321,0.5262256592692762,0.0167895548199152,0.619218726084807,0.45939914058873205,0.7185455040125052])
# SymCubatures.setweights!(cub, T[0.014676083640931844,0.015741313792886773,0.010852498353539883,0.000534442381829274,0.0005966404611656181,0.016802665132381038,0.03586790372907627,0.0020038024751056808,0.016155111812484073])

      cub_degree = 8
      tol = 1e-14      
    else
      error("polynomial degree must be <= 4 (presently)\n")
    end
  end
  mask = SymCubatures.getInternalParamMask(cub)
  append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol)
  return cub, vtx
end
  
"""
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

"""
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
