module Cubature
# routines for constructing cubatures

using LinearAlgebra
using LeastSquaresOptim
using ..OrthoPoly
using ..SymCubatures
using ..Optimizer

export pointCubature
export quadrature, quadratureUniform
export getTriCubatureGamma, getTriCubatureOmega, getTriCubatureDiagE, getTriCubatureForTetFaceDiagE
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
* `calc_grad`: indicates whether to calculate the gradients the orthogonal polynomials

**Outputs**

* `F`: the accuracy conditions for orthogonal polynomials up to degree q
* `dF`: derivative of F with respect to x, y, (z,) w, in that order

"""
function cubatureresidual(cub::LineSymCub{T}, q::Int; compute_grad) where {T}
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
    F[ptr] += (w'*P)[1]
    if compute_grad
      dPdx = OrthoPoly.diffjacobipoly(vec(x[1,:]), 0.0, 0.0, r)
      dF[ptr,:] = [w.*dPdx P]
    end
    ptr += 1
  end
  return F, dF 
end

# function cubatureresidual(cub::TriSymCub{T}, q::Int; compute_grad=true) where {T}
#   # compute the nodes and weights defined by cub
#   vtx = T[-1 -1; 1 -1; -1 1]
#   x = SymCubatures.calcnodes(cub, vtx)
#   w = SymCubatures.calcweights(cub)
#   num_eq = convert(Int, (q+1)*(q+2)/2)

#   # loop over orthogonal polynomials of degree r <= q and form conditions
#   F = zeros(T, (num_eq) )
#   F[1] = -2.0/sqrt(2.0)
#   dF = zeros(T, (num_eq, 3*cub.numnodes) )
#   ptr = 1
#   for r = 0:q
#     for j = 0:r
#       i = r-j
#       P = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
#       F[ptr] += (w'*P)[1]
#       if compute_grad
#         dPdx, dPdy = OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
#         dF[ptr,:] = [w.*dPdx w.*dPdy P]
#       end
#       ptr += 1
#       #print("(i,j,k) = (",i,",",j,",",k,"): i+j+k = ",i+j+k,"\n")
#     end
#   end
#   return F, dF 
# end

function cubatureresidual(cub::TriSymCub{T}, q::Int; compute_grad=true) where {T}
  # compute the nodes and weights defined by cub
  vtx = T[-1 -1; 1 -1; -1 1]
  x = SymCubatures.calcnodes(cub, vtx)
  w = SymCubatures.calcweights(cub)
  num_eq = convert(Int, (q+1)*(q+2)/2)

  # loop over orthogonal polynomials of degree r <= q and form conditions
  F = zeros(T, (num_eq) )
  F[1] = -2.0/sqrt(2.0)
  dF = zeros(T, (num_eq, 3*cub.numnodes) )

  V, Vdx, Vdy = OrthoPoly.vandermonde(q, x[1,:], x[2,:], compute_grad=compute_grad)
  F = V'*w + F

  # P = Optimizer.preconditioner(Matrix(V'))
  # F = (P*V')*w + P*F

  # P = Optimizer.preconditioner_shannon(V)
  # F = (P'*V')*w + P'*F

  if compute_grad
    dF = [Vdx'*Diagonal(w) Vdy'*Diagonal(w) V']
    # dF = [P*Vdx'*Diagonal(w) P*Vdy'*Diagonal(w) P*V']
    # dF = [P'*Vdx'*Diagonal(w) P'*Vdy'*Diagonal(w) P'*V']
  end
  return F, dF 
end

# function cubatureresidual(cub::TriSymCub{T}, q::Int; compute_grad=true) where {T}
#   # compute the nodes and weights defined by cub
#   # vtx = T[-1 -1; 1 -1; -1 1]
#   vtx = T[0 0; 1 0; 0 1]
#   x = SymCubatures.calcnodes(cub, vtx)
#   w = SymCubatures.calcweights(cub)
#   num_eq = convert(Int, (q+1)*(q+2)/2)

#   # loop over orthogonal polynomials of degree r <= q and form conditions
#   F = zeros(T, (num_eq) )
#   # F[1] = -2.0/sqrt(2.0)
#   dF = zeros(T, (num_eq, 3*cub.numnodes) )

#   V, Vdx, Vdy, Vinteg = OrthoPoly.vandermonde_monomial(q, x[1,:], x[2,:], compute_grad=compute_grad, compute_integ=true)
#   F = V'*w - Vinteg
#   # P = Optimizer.preconditioner(Matrix(V'))
#   # F = (P*V')*w - P*Vinteg

#   if compute_grad
#     dF = [Vdx'*Diagonal(w) Vdy'*Diagonal(w) V']
#     # dF = [P*Vdx'*Diagonal(w) P*Vdy'*Diagonal(w) P*V']
#   end
#   return F, dF 
# end

function cubatureresidual(cub::TetSymCub{T}, q::Int; compute_grad=true) where {T}
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
        F[ptr] += (w'*P)[1]
        if compute_grad
          dPdx, dPdy, dPdz = OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]),vec(x[3,:]), i, j, k)
          dF[ptr,:] = [w.*dPdx w.*dPdy w.*dPdz P]
        end
        ptr += 1
        #print("(i,j,k) = (",i,",",j,",",k,"): i+j+k = ",i+j+k,"\n")
      end
    end
  end
  return F, dF 
end

# function cubatureresidual(cub::TetSymCub{T}, q::Int; compute_grad=true) where {T}
#   # compute the nodes and weights defined by cub
#   vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
#   x = SymCubatures.calcnodes(cub, vtx)
#   w = SymCubatures.calcweights(cub)
#   num_eq = convert(Int, (q+1)*(q+2)*(q+3)/6)
#   # loop over orthogonal polynomials of degree r <= q and form conditions
#   F = zeros(T, (num_eq) )
#   F[1] = -2.0/sqrt(3.0)
#   dF = zeros(T, (num_eq, 4*cub.numnodes) )

#   V, Vdx, Vdy, Vdz= OrthoPoly.vandermonde(q, x[1,:], x[2,:], x[3,:], compute_grad=compute_grad)
#   F = V'*w + F
#   # P = Optimizer.preconditioner(Matrix(V'))
#   # F = (P*V')*w + P*F

#   if compute_grad
#     dF = [Vdx'*Diagonal(w) Vdy'*Diagonal(w) Vdz'*Diagonal(w) V']
#     # dF = [P*Vdx'*Diagonal(w) P*Vdy'*Diagonal(w) P*V']
#   end
#   return F, dF 
# end

"""
### Cubature.solvecubature!{SymCub{T}}
  
Attempts to solve for the nodes and weights of a cubature that is exact for
polynomials of degree r <= `q`.  The nodes and weights of the cubature are
defined by `cub`, which is a parametric abstract type (see symcubatures.jl).

**Inputs**

* `cub`: symmetric cubature rule
* `q`: maximum (desired) degree for which the cubature is exact
* `mask`: array of indicies of parameters and weights that are free
* `tol`: tolerance with which to solve the accuracy conditions
* `hist`: if true, print the residual-norm convergence history
* `xinit`: initial parameter guess
* `delta1`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin <= 0.1
* `delta2`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin > 0.1

**In/Outs**

* `cub`: on entry, defines the initial guess for the cubature nodes and weights.
  on exit, defines the nodes and weights that satisfy the desired accuracy.

"""
function solvecubature!(cub::SymCub{T}, q::Int, mask::AbstractArray{Int64,1};
  tol=10*eps(typeof(real(one(T)))),hist::Bool=false, verbose::Bool=false, xinit=[],
  delta1::Float64=1e-2, delta2::Float64=1e-2) where {T}
  @assert( length(mask) <= cub.numparams + cub.numweights )

  # Particle Swarm Optimization 
  n = cub.numparams + cub.numweights
  v = zeros(T, n)
  v[1:cub.numparams] = cub.params
  v[cub.numparams+1:end] = cub.weights
  # xinit=[]
  if xinit==[]
    xinit = copy(v)
  end
  nperturb = 0
  k=0
  v1 = copy(v)
  v2 = copy(v)
  
  if xinit!=[]
    SymCubatures.setparams!(cub, xinit[1:cub.numparams])
    SymCubatures.setweights!(cub, xinit[cub.numparams+1:end])
    xinit[1:cub.numparams] = cub.params
    xinit[cub.numparams+1:end] = cub.weights
  end

  nperturb_all = 0
  iter_pso = 0
  iter_lma = 0
  fmin,xinit,k = Optimizer.levenberg_marquardt(Cubature.cubatureresidual,cub,q,mask, xinit=xinit, maxiter=200, tol=tol, nu=1000.0, verbose=0)
  iter_lma+=k
  # println(fmin)
  # println(xinit)

  fbest = 1.0
  xbest = copy(v)
  for i = 1:5000
    fmin1,v1,_,_,k_pso,nperturb = Optimizer.pso(Cubature.cubatureresidual, n, cub, q, mask,xinit=xinit, np=10,maxiter=800, 
                                            tol=tol, delta1=delta1, delta2=delta2, save_iter=false,verbose=verbose)
    fmin2,v2,k_lma = Optimizer.levenberg_marquardt(Cubature.cubatureresidual,cub,q,mask, xinit=v1, maxiter=200, tol=tol, nu=1000.0, verbose=verbose)

    iter_pso += k_pso
    iter_lma += k_lma
    nperturb_all += nperturb

    # if (fmin2 < fmin1)
    #   xinit = v2
    # else
    #   xinit = v1
    # end
    use_v2 = false
    fmin3,_,_,_,_,_ = Optimizer.pso(Cubature.cubatureresidual, n, cub, q, mask,xinit=v2, np=10, maxiter=100, 
                                            tol=tol, delta1=delta1, delta2=delta2, save_iter=false,verbose=false)
    if (fmin3 < fmin1)
      use_v2 = true
    end

    if use_v2
      xinit = v2
    else
      xinit = v1
    end

    if (fbest < 1.0 && fmin2-fbest > 0.10 && fmin1-fbest > 0.10)
      xinit = xbest
    elseif use_v2 #(fmin2 < fmin1)
      xbest = v2
      fbest = fmin2
    else
      xbest = v1
      fbest = fmin1
    end

    if (fmin2<1e-12||fmin1<1e-12) && (minimum(v2)>0||minimum(v1)>0)
      if verbose
        println(v2)
      end
      break
    end
    if verbose
      println(v2)
    end
  end
  v = v2

  SymCubatures.setparams!(cub, v[1:cub.numparams])
  SymCubatures.setweights!(cub, v[cub.numparams+1:end])
  F, _ = Cubature.cubatureresidual(cub, q)
  res = norm(F)
  # hist ? print("\titer ",k, ":  nrestart ",nrestart,":  res norm = ",res,"\n") : nothing
  hist ? println("----------------------------------------------------------------------------------------------") : nothing
  hist ? println("iter_pso = ",iter_pso, ":  iter_lma = ",iter_lma, ":  nperturb_pso = ",nperturb_all,":  res norm = ",res) : nothing
  hist ? println("----------------------------------------------------------------------------------------------") : nothing
  # hist ? print("res norm = ",res,"\n") : nothing
  if res < tol
    return
  end

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
function solvecubatureweights!(cub::SymCub{T}, q::Int;
                                  tol=10*eps(typeof(real(one(T)))),
                                  hist::Bool=false) where {T}
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
    JtJ = Jac'*(dF[:,end-cub.numnodes+1:end]'*
      dF[:,end-cub.numnodes+1:end])*Jac
    H = JtJ + nu*diagm(diag(JtJ))
    g = -Jac'*dF[:,end-cub.numnodes+1:end]'*F
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
    tol = 6e-14
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

Returns a cubature rule and vertices for facets of the SBP DiagE operators 
on tetrahedra; these should not be used for 2D problems as they do not 
satisfy the accruacy requirements along their edge.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `T`: the data type used to represent the cubature
* `faceopertype`: the operator type on the facets of the tetrahedron
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function getTriCubatureForTetFaceDiagE(q::Int, T=Float64; faceopertype::Symbol=:DiagE,
  tol=10*eps(typeof(real(one(T)))))
  cub_degree = q
  mask = zeros(Int64, (0))

  if faceopertype == :Omega
    cub, vtx = getTriCubatureOmega(q)
  else
    if q<=2 #3 nodes
      cub = SymCubatures.TriSymCub{T}(vertices=false, midedges=true, centroid=false)
      SymCubatures.setweights!(cub, T[0.6666666666666666])
      cub_degree = 2
      tol = 5e-15

      #4 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices=true, midedges=false, centroid=true)
      # SymCubatures.setweights!(cub, T[0.16666666666666624, 1.4999999999999973])
      # cub_degree = 2
      # tol = 5e-15
    elseif q<=4 #9 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = false,
                                      numedge = 0,
                                      numS21 = 1,
                                      numS111 = 0)
      SymCubatures.setparams!(cub, T[0.3771609693928902])
      SymCubatures.setweights!(cub, T[0.0410802706918665, 0.123950997736534, 0.5016353982382646])
      cub_degree = 4
      tol = 5e-15

      #7 nodes, but this leads to tet element with 26 nodes while the above gives 23 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = false,
      #                                 midedges = true,
      #                                 centroid = true,
      #                                 numedge = 0,
      #                                 numS21 = 1,
      #                                 numS111 = 0)
      # SymCubatures.setparams!(cub, T[0.22222222222222215])
      # SymCubatures.setweights!(cub, T[0.1523809523809521, 0.28928571428571387, 0.6749999999999992])
      # cub_degree = 4
      # tol = 5e-15
    elseif q<=6 #15 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      centroid = false,
                                      numedge = 1,
                                      numS21 = 2,
                                      numS111 = 0)
      SymCubatures.setparams!(cub, T[0.8506802519794945, 0.23722737279318576, 0.3077459416259917])
      SymCubatures.setweights!(cub, T[0.01426071861440897, 0.3303589772911334, 0.20376930605390392, 0.059138832353610636])
      cub_degree = 6
      tol = 5e-15
      
      # omega type facet nodes
      # cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=false,
      #                                 numS21=2, numS111=1)
      # SymCubatures.setweights!(cub, T[0.1016898127404136;
      #                                 0.23357255145275863;
      #                                 0.16570215123674722])
      # SymCubatures.setparams!(cub, T[0.1261780289830045;
      #                                0.49857349034182097;
      #                                0.10629009968963399;
      #                                0.6207049020675687])
      # cub_degree = 6
    elseif q<=8 #22 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = true,
                                      numedge = 1,
                                      numS21 = 1,
                                      numS111 = 1)
      SymCubatures.setparams!(cub, T[0.16099183834007516, 0.800367892880542, 0.6058255660767269, 0.21518364356973504])
      SymCubatures.setweights!(cub, T[0.006036623735435305, 0.040748592504715506, 0.09377442432581949, 0.02798160585550063, 0.19106553442197458, 0.2640382366372378])
      cub_degree = 8
      tol = 7e-15

      # omega type facet nodes
      # cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=true,
      # numS21=3, numS111=1)
      # SymCubatures.setweights!(cub, T[0.20643474106943666;
      #   0.19018326853457126;
      #   0.06491699524639659;
      #   0.05446062834886685;
      #   0.2886312153555818])
      # SymCubatures.setparams!(cub, T[0.3411386155035168;
      #   0.9185851765854467;
      #   0.10109445663405793;
      #   0.5262256592692812;
      #   0.016789554819910586])      
      # cub_degree = 8
      # tol = 6e-14
    elseif q<=10 #28 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = true,
                                      numedge = 1,
                                      numS21 = 3,
                                      numS111 = 1)
      SymCubatures.setparams!(cub, T[0.40849548761177584, 0.1693510205976773, 0.895541492424858, 0.9192858150590572, 0.5199457984883094, 0.07662224108118755])
      SymCubatures.setweights!(cub, T[0.0017341485435839296, 0.02676086117258247, 0.1584314816374832, 0.07114324960955178, 0.15198175500182232, 0.013615409089319942, 0.08154977731366873, 0.1988543936869974])
      cub_degree = 10
      tol = 5e-15
    else
      error("polynomial degree must be <= 10 (presently)\n")
    end
  end
  mask = SymCubatures.getInternalParamMask(cub)
  append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
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
    if q<=1 #6 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = false,
                                      numedge = 0,
                                      numS21 = 0,
                                      numS111 = 0)                
      SymCubatures.setweights!(cub, T[0.5525545450892143; 0.11411212157745271])
      cub_degree = 1
      tol = tol

    elseif q<=2 #7 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = true,
                                      numedge = 0,
                                      numS21 = 0,
                                      numS111 = 0)

      SymCubatures.setweights!(cub, T[0.04761904761904726, 0.47619047619047483, 0.4285714285714286])
      cub_degree = 2
      tol = tol

      # SymCubatures.setweights!(cub, T[9/10, 1/10, 4/15])
      # cub_degree = 3
      # tol = 1e-14
    elseif q <=3 #10 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      centroid = true,
                                      numedge = 1,
                                      numS21 = 0,
                                      numS111 = 0)

      SymCubatures.setparams!(cub, T[0.7236067977499789])
      SymCubatures.setweights!(cub, T[0.03333333333333329, 0.16666666666666605, 0.8999999999999969])
      cub_degree = 3
      tol=1e-14
    elseif q <= 4
      #12 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 1,
                                      numedge = 1,
                                      numS111 = 0,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.4257087142236166, 0.7236067977499789])
      SymCubatures.setweights!(cub, T[0.025044506019598876, 0.42703175864395354, 0.10729520100155711])
      # cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=1,
      #                                 midedges=false,
      #                                 numS21=1, numS111=0,
      #                                 centroid=false)
      # SymCubatures.setparams!(cub, T[0.4257087142201423;
      #                                0.5*(1 + sqrt(1/5))])
      # SymCubatures.setweights!(cub, T[0.02504450602156441;
      #                                 0.4270317586471588;
      #                                 0.10729520099967835])
      cub_degree = 4
      tol = 1e-14
    elseif q <=5 #15 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = false,
                                      numedge = 1,
                                      numS21 = 1,
                                      numS111 = 0)
      SymCubatures.setparams!(cub, T[0.41469035132718185; 0.8273268353539885])
      SymCubatures.setweights!(cub, T[0.014698618394803228, 0.09752600361864236, 0.44155411568082115, 0.056443964486199594])
      # SymCubatures.setweights!(cub, T[0.02323751046092639; 0.07171714179335946; 0.4415541156808216; 0.06507894936577964])
      cub_degree = 5
      tol = 1e-14
    elseif q <= 6 #18 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 2,
                                      numedge = 1,
                                      numS111 = 0,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.8487720503437628, 0.28401016819355557, 0.8273268353539885])
      SymCubatures.setweights!(cub, T[0.009130264572198617, 0.06201621057208993, 0.29870677879935964, 0.20607267198227744, 0.045370370370370325])
      cub_degree = 6
      tol = 1e-14
      # cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=1,
      #                                 midedges=true,
      #                                 numS21=2, numS111=0, centroid=false)
      # SymCubatures.setparams!(cub, T[0.8487720503426771;
      #                                0.28401016818370567;
      #                                0.5*(1 + sqrt(3/7))])
      # SymCubatures.setweights!(cub, T[0.00913026457472031;
      #                                 0.06201621056804736;
      #                                 0.2987067788024998;
      #                                 0.20607267197855683;
      #                                 0.04537037036896131])
      # tol = 1e-14
    elseif q <=7 #24 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      centroid = false,
                                      numedge = 2,
                                      numS21 = 1,
                                      numS111 = 1)
      SymCubatures.setparams!(cub, T[0.8461370386526059, 0.8825276619647324, 0.642615758240322, 0.33879488493771764, 0.19830545894321575])
      SymCubatures.setweights!(cub, T[0.00753161345765886, 0.30480530834098646, 0.02108627120396598, 0.045504525339300515, 0.1105740758907442])
      cub_degree = 7
      tol=5e-15
    elseif q <= 8 # 27 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 2,
                                      numedge = 2,
                                      numS111 = 1,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.5306627609684195, 0.20735501628561037, 0.8825276619647324, 0.6426157582403226, 0.17654792120316234, 0.6492809445301031])
      SymCubatures.setweights!(cub, T[0.004361575619973876, 0.15984019211140704, 0.1135887626556459, 0.022078990549428416, 0.02786777148585481, 0.14449130610453667])
      cub_degree = 8
      tol = 1e-14
      # cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=2,
      #                                 numS21=2, numS111=1)
      # SymCubatures.setparams!(cub, T[0.20735501628574252;
      #                                0.5306627609507977;
      #                                0.5*(1 + sqrt(1/3 - 2*sqrt(7)/21));
      #                                0.5*(1 + sqrt(1/3 + 2*sqrt(7)/21));
      #                                0.6492809445444747;
      #                                1.1741711342683223])
      # SymCubatures.setweights!(cub, T[0.004361575620524937;
      #                                 0.11358876265867929;
      #                                 0.15984019213585915;
      #                                 0.027867771483031517;
      #                                 0.02207899054885172;
      #                                 0.14449130609379232])  
      # cub_degree = 8

      # tol = 1e-14
    elseif q <=9 #34 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = false,
                                      numedge = 2,
                                      numS21 = 3,
                                      numS111 = 1)
      SymCubatures.setparams!(cub, T[0.5110998040469292, 0.9089231881407306, 0.19994428870176986, 0.9151119481392835, 0.7344243967353571, 0.1465093205286518, 0.5715830463482526])
      SymCubatures.setweights!(cub, T[0.0019241266200829663, 0.022987519027599365, 0.1866822178508963, 0.09270056884965712, 0.0985056912129899, 0.016540728394130698, 0.01779352639593344, 0.09759901676265635])
      cub_degree = 9
      tol=5e-15
    elseif q <= 10 #36 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 2,
                                      numedge = 2,
                                      numS111 = 2,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.1025947577289524, 0.514543490059506, 0.9151119481392835, 0.7344243967353571, 0.16739689557448953, 0.736084860943316, 0.37663587800802634, 0.15391188146518084])
      SymCubatures.setweights!(cub, T[0.001135617009180763, 0.022982971732153294, 0.03238151607768464, 0.1823790802107045, 0.009620012578487587, 0.019610669249206004, 0.09949978404801857, 0.08516327494275952])
      cub_degree = 10
      tol = 1e-14
      # 39 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = true,
      #                                 numS21 = 1,
      #                                 numedge = 2,
      #                                 numS111 = 3,
      #                                 centroid = false)
      # SymCubatures.setparams!(cub, T[0.1332594652103378, 0.9151119481392835, 0.7344243967353571, 0.8013438243116634, 0.14971621943108124, 0.16193803910480203, 0.4343229317011592, 0.4438585788661879, 0.5888276411796545])
      # SymCubatures.setweights!(cub, T[0.001911431546031246, 0.020226485585902412, 0.048604009344358975, 0.010871420814832881, 0.019753079058769735, 0.07138434601195054, 0.09622530312117668, 0.0997282210884573])
      # cub_degree = 10
      # tol = 1e-14
    elseif q <=11 #40 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      centroid = true,
                                      numedge = 3,
                                      numS21 = 4,
                                      numS111 = 1)
      SymCubatures.setparams!(cub, T[0.40821896354263243, 0.981013735325053, 0.16748106658466438, 0.8690323743735683, 0.9358700742548033, 0.7958500907165711, 0.6046496089512394, 0.10247525036525074, 0.5346096921202054])
      SymCubatures.setweights!(cub, T[0.0008860287032078429, 0.1399375485517579, 0.04956216973314194, 0.07404852901772899, 0.13638424140745345, 0.010480045148169446, 0.011623548319685177, 0.0012118355050860282, 0.08102106300279119, 0.1715254959057385])
      cub_degree = 11
      tol=5e-15
      # 43 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = false,
      #                                 centroid = true,
      #                                 numedge = 3,
      #                                 numS21 = 3,
      #                                 numS111 = 2)
      # SymCubatures.setparams!(cub, T[0.8571396116924792, 0.42137206563390533, 0.15416411447545023, 0.9358700742548033, 0.7958500907165711, 0.6046496089512394, 0.8193910152685325, 0.07542901871217762, 0.1301469442971666, 0.4770052452777402])
      # SymCubatures.setweights!(cub, T[0.0011717118884880214, 0.12473132707644294, 0.12730130899704278, 0.06274770142440542, 0.009242885505864128, 0.013103502174005579, 0.008618295447807097, 0.04030615822887277, 0.07902105209405838, 0.15039249113721462])
      # cub_degree = 11
      # tol=5e-15
      # 43 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = false,
      #                                 centroid = true,
      #                                 numedge = 3,
      #                                 numS21 = 5,
      #                                 numS111 = 1)
      # SymCubatures.setparams!(cub, T[0.8690570506162132, 0.4099529922719056, 0.1674893211310927, 0.9810521073173698, 0.40735472249904886, 0.9358700742548033, 0.7958500907165711, 0.6046496089512394, 0.10243650581409067, 1.3628948408784771])
      # SymCubatures.setweights!(cub, T[0.0008858048750959239, 0.13640002400783297, 0.04612209214143705, 0.07405649995653696, 0.049554179974936816, 0.0938448125463987, 0.010480865652606294, 0.011622127797623866, 0.0011870772410632817, 0.08102044324683055, 0.17154667586454087])
      # cub_degree = 11
      # tol=1e-14
    elseif q <= 12 # 48 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 3,
                                      numedge = 3,
                                      numS111 = 3,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.3416879952650103, 0.11946629335996972, 0.5747680456926054, 0.9358700742548033, 0.7958500907165711, 0.6046496089512394, 0.38458694989344594, 0.1066659705465824, 0.6718411050524655, 0.30250431498537617, 0.7519048745331144, 0.09573484698587167])
      SymCubatures.setweights!(cub, T[0.0012716011685488162, 0.08154407642608534, 0.03881033382161809, 0.07780339023363679, 0.006972694460371076, 0.010101834216810933, 0.011019799229776561, 0.05674688135939982, 0.08690355893881288, 0.06187386430321755])
      cub_degree = 12
      tol = 1e-14
      
      #48 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices=true, numedge=3,
      #                                 numS21=3, numS111=3)
      # SymCubatures.setparams!(cub, T[0.5747680454804257; 0.34168799534894295;
      #                                0.1194662933444393; 0.6046496089512394;
      #                                0.7958500907165711; 0.9358700742548033;
      #                                0.38458694986468334; 0.10666597060509266;
      #                                0.3025043145989276; 0.6718411052897879;
      #                                0.09573484681006857;0.7519048745731671])
      # SymCubatures.setweights!(cub, T[0.001271601169161372; 0.07780339049594198;
      #                                 0.08154407642794283; 0.03881033381664769;
      #                                 0.01101979920649565; 0.010101834223603934;
      #                                 0.006972694458714173; 0.056746881379720164;
      #                                 0.08690355889320943; 0.06187386421672949])  
      # cub_degree = 12
    elseif q <= 13 #55 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = true,
                                      numedge = 3,
                                      numS21 = 4,
                                      numS111 = 3)
      SymCubatures.setparams!(cub, T[0.3218395626601594, 0.4957183097055999, 0.9721805710603578, 0.10876332278702848, 0.94987899770573, 0.8385931397553689, 0.6815587319130891, 0.9943269168205083, 0.7196871268867262, 0.11795721403434152, 0.6250454566805667, 0.0943997083693319, 0.3344326521708524])
      SymCubatures.setweights!(cub, T[0.0007222803300461352, 0.0017482849407310453, 0.07989770555605481, 0.07525421818189708, 0.04552726413680475, 0.030615775227096503, 0.005199952283163629, 0.006832142397590008, 0.01150160780019001, 0.07355921265284922, 0.05745315943375207, 0.043511497615255, 0.11035798178531019])
      cub_degree = 13
      # # #61 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = true,
      #                                 centroid = true,
      #                                 numedge = 3,
      #                                 numS21 = 6,
      #                                 numS111 = 3)
      # SymCubatures.setparams!(cub, T[0.4757066927127499, 0.8640750851214627, 0.848071695696415, 0.09588932222312066, 0.9730573610986599, 0.2970876410785625, 0.94987899770573, 0.8385931397553689, 0.6815587319130891, 0.08772329380398057, 0.303910195608789, 0.583928734667523, 0.08412812390769837, 0.6214646072392227, 0.18555290129142354])
      # SymCubatures.setweights!(cub, T[0.0007495859766713502, 0.0004006124634786362, 0.09741666818077659, 0.018352298042798858, 0.0765136751062513, 0.024240908916246806, 0.050029973029268995, 0.06791617164641447, 0.004523675585682937, 0.00614216438323085, 0.01039041823865825, 0.038711997393405265, 0.031386407110491864, 0.05529527857696424, 0.11444067218367934])
      # cub_degree = 13
      tol=5e-15
    elseif q <= 14 # 57 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 5,
                                      numedge = 3,
                                      numS111 = 3,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.09123182474143035, 0.29679838270099523, 0.5601853746013142, 0.9547456706858821, 0.8556466743175067, 0.94987899770573, 0.8385931397553689, 0.6815587319130891, 0.09449105513358476, 0.29883210637720314, 0.28045625536995294, 0.5626405122741717, 0.08624428771540175, 0.6026786176032012])
      SymCubatures.setweights!(cub, T[0.0007538778164594649, 0.009414839808665405, 0.023018612111280527, 0.058120471428513544, 0.09454606798115053, 0.05373738699005786, 0.0727513309290953, 0.004178640246830719, 0.0072314035856924815, 0.008184710625980713, 0.040253997904260874, 0.06900626592215424, 0.04830702151580295])
      cub_degree = 14
      tol = 1e-14
    elseif q <= 15 #69 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      centroid = false,
                                      numedge = 4,
                                      numS21 = 6,
                                      numS111 = 4)
      SymCubatures.setparams!(cub, T[0.453639886926089, 0.757300794432156, 0.9833658031003171, 0.2673140133674513, 0.09859868374369533, 0.9313655727780035, 0.9597669540832294, 0.8693869325527526, 0.7389624749052223, 0.5826394788331934, 0.05491312392961762, 0.6373866699608978, 0.9962538574472307, 0.7237541325698926, 0.19371136334839004, 0.543875090517202, 0.31621286243906677, 0.07653701224015233])
      SymCubatures.setweights!(cub, T[0.00041188664492458155, 0.07744252741964325, 0.07978483547877466, 0.022442598130765025, 0.055748397969630784, 0.02612364913831872, 0.047345838243333084, 0.00378591167322386, 0.00494107871623112, 0.005183734332323488, 0.002737379819615485, 0.03279819517102539, 0.041923249456950765, 0.05253038629189004, 0.03478353135937816])
      # #78 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = false,
      #                                 centroid = false,
      #                                 numedge = 4,
      #                                 numS21 = 7,
      #                                 numS111 = 5)
      # SymCubatures.setparams!(cub, T[0.4077750904687781, 0.7603167059782701, 0.09923314655025814, 0.5392503479960836, 0.7336297585443722, 0.8840881622927267, 0.9647641607531148, 0.9597669540832294, 0.8693869325527526, 0.7389624749052223, 0.5826394788331934, 0.5906777979138421, 0.22187382664358296, 0.4240838353085315, 0.9561037136275863, 0.32312514688527244, 0.1922339643728221, 0.3059644597113809, 0.06164364302370151, 0.61117935281228, 0.06810139957582922])
      # SymCubatures.setweights!(cub, T[0.00038560167936041434, 0.051828320319575245, 0.04897897079965172, 0.025976977006482195, 0.013736493536953337, 0.008950065388866491, 0.060619930406823286, 0.04266756023023606, 0.0038903087243377906, 0.003952860653539124, 0.005383423214263016, 0.006285453618608518, 0.0547038473668017, 0.028976866620338368, 0.039091444929817185, 0.026295287780152696, 0.03818188074150067])
      # SymCubatures.setparams!(cub, T[0.926073884304195, 0.6324214560155823, 0.2732070265874975, 0.9800392449549591, 0.9275624142028438, 0.09803698045583892, 0.4678434864713421, 0.9597669540832294, 0.8693869325527526, 0.7389624749052223, 0.5826394788331934, 0.6281478611396281, 0.05988739351469942, 0.5609548008743171, 0.2126550295547661, 0.31652535690496086, 0.08109010595545985, 0.8143850982384193, 0.7648461377280185, 0.932544753465994, 0.30649918776292634])
      # SymCubatures.setweights!(cub, T[0.00043899617260247656, 0.02832593621352893, 0.0284355895772454, 0.05708633451373208, 0.023936844515619497, 0.022118700922270935, 0.02609816009431366, 0.07746633875563608, 0.0036965993872807804, 0.005287955943337329, 0.005332396620054094, 0.00394995932661379, 0.03487470485376105, 0.06094440138731721, 0.03641363868298402, 0.024295609746207345, 0.026584617003303227])
      tol=5e-15
    elseif q <= 16 # 72 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 5,
                                      numedge = 4,
                                      numS111 = 5,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.06522886751875544, 0.22580096183071185, 0.7596981294914541, 0.4709436867165254, 0.9626832123175285, 0.9597669540832294, 0.8693869325527526, 0.7389624749052223, 0.5826394788331934, 0.6837234484596031, 0.07332394938311322, 0.0696195668842623, 0.212933356067202, 0.23253598278247095, 0.45142512187590106, 0.07099634477078105, 0.4271102931595644, 0.7291817098161258, 0.24027678681212955])
      SymCubatures.setweights!(cub, T[0.00046398080499042475, 0.011683737917833769, 0.03748058638601257, 0.08029394321068659, 0.07225832122458112, 0.0350040593727091, 0.002438057356674435, 0.004350159403927255, 0.005775046495650795, 0.006715480702033995, 0.03300813082784586, 0.021283575846091102, 0.05161492134502783, 0.028186978295003252, 0.061368668602672])
      cub_degree = 16
      tol = 1e-14
      #75 nodes
      # cub = SymCubatures.TriSymCub{Float64}(vertices=true, numedge=4,
      #                                       numS21=4, numS111=6)
      # SymCubatures.setparams!
      # (cub, T[0.0768946752469594; 0.6109907336234316; 0.4179369130153705;
      #         0.23221732669622028; 0.5826394788331936; 0.7389624749052223;
      #         0.8693869325527526; 0.9597669540832294; 0.930330150896981;
      #         0.6679157686119799; 0.8075697058065031; 1.1286554029515519;
      #         0.07116357646128006; 0.25084609708056116; 0.20705942876211147;
      #         0.46901560437791967; 0.06341212586405608; 0.504346239131436;
      # 0.7497163430497378; 1.0452430326898021])
      # SymCubatures.setweights!
      # (cub, T[0.0004934174763938973; 0.016280379328615396; 0.03461636692031435;
      #         0.0544553219641883; 0.04016179760262383; 0.005790378245871778;
      #         0.005099447551845056; 0.004407182316786121; 0.002824162676005338;
      #         0.05536224506473959; 0.03347658309313654; 0.025582552500936707;
      #         0.047069550256364515; 0.029816087172401348; 0.050901502809179516])
      # cub_degree = 16
    elseif q<=17 #78 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      centroid = false,
                                      numedge = 4,
                                      numS21 = 6,
                                      numS111 = 5)
      SymCubatures.setparams!(cub, T[0.9609712932934811, 0.041766468223703515, 0.28302191880327227, 0.760185750106773, 0.46292528148197715, 0.173093973905005, 0.9670007152040296, 0.8922417368315723, 0.7826176634981026, 0.6478790677934697, 0.1610470535559693, 0.04528225579552883, 0.7226288131709319, 0.2400281193275206, 0.36438411202861126, 0.07009663480762927, 0.06638724199110738, 0.6541614471267511, 0.20944043731114892, 0.46300947657489105])
      SymCubatures.setweights!(cub, T[0.00028908921124315597, 0.006749278476735112, 0.04010182067917665, 0.006389003512011476, 0.031350458689112376, 0.08284141010739084, 0.07608792396452779, 0.028383320820692025, 0.0008879784296625767, 0.0025339572186700047, 0.005105107386578676, 0.004945455787338109, 0.013263503268280243, 0.06280733563151443, 0.028782576708273472, 0.0346732468195029, 0.04423801935306822])
      #81 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = true,
      #                                 centroid = false,
      #                                 numedge = 4,
      #                                 numS21 = 7,
      #                                 numS111 = 5)
      # SymCubatures.setparams!(cub, T[0.9593342591248476, 0.055020860803980914, 0.28526621062594437, 0.8387490974325792, 0.7575470957605024, 0.4643330654690022, 0.16444412944990114, 0.9670007152040296, 0.8922417368315723, 0.7826176634981026, 0.6478790677934697, 0.18880867265986118, 0.035555254484675405, 0.7184784536560056, 0.23691797352429475, 0.3741094228156754, 0.07603572398114046, 0.062436028204803366, 0.6623692442916642, 0.20752925425214702, 0.4787260323140681])
      # SymCubatures.setweights!(cub, T[0.00039358697853621074, 0.0071635142083468235, 0.04019813807428178, 0.009282438655081281, 0.03869786195664633, 0.011258459094176387, 0.07803814033007511, 0.07672709475642626, 0.02944041155775245, 0.0013436036287747311, 0.0013462041905516386, 0.0055617515345457525, 0.004323946528663438, 0.01238657052502277, 0.05821065906057333, 0.029593585323296128, 0.03325710394842388, 0.04171008578782039])
      cub_degree = 17
      tol = 1e-14
    elseif q <= 18 # 93 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = true,
                                      numS21 = 5,
                                      numedge = 4,
                                      numS111 = 8,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.9707518160195144, 0.20641902821770264, 0.05909234695936303, 0.7110519489326181, 0.39622899824483554, 0.9670007152040296, 0.8922417368315723, 0.7826176634981026, 0.6478790677934697, 0.3906058825830485, 0.19331653245290645, 0.8008189404907055, 0.17109449910705912, 0.24603845950256573, 0.9650206949829071, 0.20142808520129168, 0.5868483166437959, 0.4088076515219679, 0.6463428792626565, 0.4093195068370922, 0.05981030752842776, 0.1973673755323838, 0.06444759127052886, 0.673801874205518, 0.06033239562749046])
      SymCubatures.setweights!(cub, T[0.00031737159521531795, 0.00491592992346147, 0.02867493867862121, 0.029107246299861025, 0.009795965466384334, 0.031357172956706196, 0.05321207918099946, 0.0017927521866819998, 0.0033404931412304374, 0.004059862113462083, 0.004795116727932584, 0.03247349913877226, 0.022030499945312042, 0.0192471893187778, 0.03594121829355437, 0.05970283389294487, 0.02407056230858391, 0.018810685588525493, 0.028378268626930925])
      cub_degree = 18
      tol = 1e-14
    elseif q<=19 #96 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      centroid = false,
                                      numedge = 5,
                                      numS21 = 5,
                                      numS111 = 8)
      SymCubatures.setparams!(cub, T[0.584118160705067, 0.05594736023975348, 0.35206685470664284, 0.8193250037475288, 0.9106302368288487, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 0.36516822326597015, 0.05300184234489419, 0.39755390108675603, 0.17428704397839728, 0.05879006206078798, 0.1815645896189941, 0.17519871987759317, 0.2119731854287182, 0.5725515218226456, 0.35807920419385625, 0.054591434526018896, 0.8428999975969531, 0.17752640292448835, 0.6418878905694481, 0.5915035416317636, 0.05407771456049964])
      SymCubatures.setweights!(cub, T[0.0002316646038971307, 0.058241783304694345, 0.008602602423026792, 0.04069174280772659, 0.05226677627664586, 0.042505709104633006, 0.0014109503198198763, 0.00256960026599889, 0.00309831866176889, 0.003741487391936303, 0.00407918767311633, 0.018331624199230855, 0.034095002356800554, 0.015212277721266073, 0.015179420870583436, 0.04887476434442788, 0.023431653507726604, 0.04029884352862801, 0.021740063231717795])
      
      #99 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = true,
      #                                 midedges = false,
      #                                 centroid = false,
      #                                 numedge = 5,
      #                                 numS21 = 6,
      #                                 numS111 = 8)
      # SymCubatures.setparams!(cub, T[0.3358980598350944, 0.9789706742159475, 0.9289390131831053, 0.7440240615283973, 0.06501518566521776, 0.19142118573089564, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 0.1583978333661437, 0.3969064666813339, 0.042212410750339456, 0.6904634855900522, 0.42850251199567174, 0.04755503506287861, 0.9879382272825736, 0.46840576395362366, 0.9475849438753418, 0.29884780482324186, 0.6542434414374048, 0.14303581650779246, 0.5459090121405945, 0.29780285727201, 0.21229041305026342, 0.05731545867768344])
      # SymCubatures.setweights!(cub, T[0.00021393564288845607, 0.03852805201142253, 0.020863457513064957, 0.03579097887221835, 0.055167566276037854, 0.01162155336267164, 0.02853211993973326, 0.001663406453602736, 0.0025528949149080204, 0.002946876800429615, 0.0029333115372134344, 0.003080267714905865, 0.03252630590418761, 0.019886495371279594, 0.019158063559635667, 0.028048303639912196, 0.036293762745006924, 0.034722432073318175, 0.03654424141681496, 0.01761813939309992])
      cub_degree = 19
      tol = 5e-15
    elseif q <= 20 # 103 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = true,
                                      midedges = false,
                                      numS21 = 5,
                                      numedge = 5,
                                      numS111 = 9,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.5706569630167525, 0.045476786416401474, 0.3386020213960199, 0.8173378833539444, 0.921990277619815, 0.9724496361114411, 0.9096396608220033, 0.8164380765159303, 0.6997654704826745, 0.5682764664274638, 0.33734350166550275, 0.051994927730595886, 0.39122487789147525, 0.17029307283717016, 0.0555956144787866, 0.15691138175383615, 0.6048213057399955, 0.16379237699452232, 0.1617644050989192, 0.21607595265283405, 0.56000817972725, 0.34672951075674346, 0.053470945381353656, 0.8336640057658877, 0.19072411914891063, 0.7807442886366605, 0.5689238612251344, 0.05028584743589604])
      SymCubatures.setweights!(cub, T[0.00020287039679960927, 0.05480860557627294, 0.0059385235222309644, 0.03842886096227663, 0.05547473001542249, 0.0169887933138356, 0.0011773441836911722, 0.0024488282936603527, 0.0030066460137970824, 0.0034778130522919964, 0.004015524716182081, 0.018032136009343752, 0.030456546263333304, 0.013428745082940518, 0.030846728655919257, 0.015872166479982706, 0.0501088680789255, 0.02445659700889463, 0.025248076325159272, 0.02109467351877063, 0.02244868654213111])
      cub_degree = 20
      tol = 1e-14
    else
      error("polynomial degree must be <= 20 (presently)\n")
    end

  else
    # do not include vertices in the cubature
    if q <= 1 #6 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 0,
                                      numedge = 1,
                                      numS111 = 0,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.21132486540518702])
      SymCubatures.setweights!(cub, T[0.3333333333333332])
      cub_degree = 1
      tol = tol
      # # Taken from Chen-Shu, 2017
      # cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=false, midedges=false, 
      #                                 numedge=1, numS21=0, numS111=0)
      # SymCubatures.setparams!(cub,T[0.5 - sqrt(3)/6])
      # SymCubatures.setweights!(cub, T[1/3])
      # cub_degree = 1
    elseif q <= 2 # 7nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 0,
                                      numedge = 1,
                                      numS111 = 0,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.21132486540518702])
      SymCubatures.setweights!(cub, T[0.16666666666666638, 0.9999999999999993])
      cub_degree = 2
      tol = 1e-14
      # cub = SymCubatures.TriSymCub{T}(vertices=false, numedge=1,
      #                                 centroid=true)
      # SymCubatures.setweights!(cub, T[1/6; 1.0])
      # SymCubatures.setparams!(cub, T[0.5*(1 + 1/sqrt(3))])
      # cub_degree = 2
    elseif q <=3 #10 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 0,
                                      numedge = 1,
                                      numS111 = 0,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.1127016653792583])
      SymCubatures.setweights!(cub, T[0.1999999999999999, 0.08333333333333327, 0.8999999999999997])
      cub_degree = 3
      tol = 1e-14
      # # Taken from Chen-Shu, 2017
      # cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=true, midedges=true, 
      #                     numedge=1, numS21=0, numS111=0)
      # SymCubatures.setparams!(cub,T[0.5 - sqrt(15)/10])
      # SymCubatures.setweights!(cub, T[1/5; 1/12; 9/10])
      # cub_degree = 3
    elseif q <= 4 #12 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 1,
                                      numedge = 1,
                                      numS111 = 0,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.40936128314152403, 0.1127016653792583])
      SymCubatures.setweights!(cub, T[0.1168096321294691, 0.44903690872888175, 0.05041006290415779])
      cub_degree = 4
      tol = 1e-14
      # cub = SymCubatures.TriSymCub{T}(vertices=false, numedge=1,
      #                                 midedges=true, numS21=1, numS111=0,
      #                                 centroid=false)
      # SymCubatures.setparams!(cub, T[0.4093612831422611;
      #                                0.5*(1 + sqrt(3/5))])
      # SymCubatures.setweights!(cub, T[0.11680963211922607;
      #                                 0.449036908703613;
      #                                 0.050410062902983145])
      # cub_degree = 4
      # tol = 5e-15
    elseif q <=5 #18 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 0,
                                      numedge = 2,
                                      numS111 = 1,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.06943184420297371, 0.33000947820757187, 0.457538013666131, 0.37414775838255504])
      SymCubatures.setweights!(cub, T[0.03019802974513124, 0.08091308136597997, 0.22222222222222215])
      cub_degree = 5
      tol = 1e-14
      # # Taken from Chen-Shu, 2017
      # cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=false, midedges=false, numedge=2,
      #                   numS21=0, numS111=1)
      # SymCubatures.setparams!(cub,[0.330009478207572; 
      #                             0.0694318442029737; 
      #                             0.3741477583825542; 
      #                             1.168314227951314])
      # SymCubatures.setweights!(cub, T[0.08091308136598; 
      #                                 0.0301980297451312; 
      #                                 0.2222222222222224])
      # cub_degree = 5
      # tol = 5e-15
    elseif q <= 6 #21 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 1,
                                      numedge = 2,
                                      numS111 = 1,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.8478748148896262, 0.06943184420297371, 0.33000947820757187, 0.31541668315764937, 0.21896598857216268])
      SymCubatures.setweights!(cub, T[0.30870082526638737, 0.019432110768375376, 0.053997899611202924, 0.10555291032056137])
      cub_degree = 6
      tol = 1e-14
      # cub = SymCubatures.TriSymCub{T}(vertices=false, numedge=2,
      #                                 midedges=false,
      #                                 numS21=1, numS111=1, centroid=false)
      # SymCubatures.setparams!(cub, T[0.8478748148895112; 
      #                                0.5*(1 + sqrt(3/7 - 2/7*sqrt(6/5)));
      #                                0.5*(1 + sqrt(3/7 + 2/7*sqrt(6/5)));
      #                                0.3154166831592224;
      #                                0.2189659885706149])
      # SymCubatures.setweights!(cub, T[0.30870082526604714;
      #                                 0.0539978996110132;
      #                                 0.01943211076834772;
      #                                 0.10555291032094156])
      # cub_degree = 6
      # tol = 3e-15
    elseif q <= 7 #22 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 2,
                                      numedge = 2,
                                      numS111 = 0,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.8768479048816373, 0.27886746283090746, 0.04691007703066802, 0.23076534494715845])
      SymCubatures.setweights!(cub, T[0.03707416966789951, 0.24947346457954658, 0.21085865924168884, 0.013202630162003234, 0.04106091936085787, 0.18219982239542837])
      cub_degree = 7
      tol = 1e-14
      # # Taken from Chen-Shu, 2017
      # cub = SymCubatures.TriSymCub{T}(vertices=false, centroid=true, midedges=true, numedge=2,
      # numS21=2, numS111=0)
      # SymCubatures.setparams!(cub, T[0.876847904881637;
      #       0.2788674628309072; 
      #       0.230765344947159; 
      #       0.046910077030668])
      # SymCubatures.setweights!(cub, T[0.03707416966789956; 
      #       0.2494734645795472; 
      #       0.2108586592416888; 
      #       0.041060919360858; 
      #       0.01320263016200324;
      #       0.1821998223954268])
      # cub_degree = 7
    elseif q <= 8 #28 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 2,
                                      numedge = 2,
                                      numS111 = 1,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.18847440090927203, 0.4284337169945454, 0.04691007703066802, 0.23076534494715845, 0.16905856872047995, 0.6829175075115606])
      SymCubatures.setweights!(cub, T[0.021995779196974322, 0.106502996693885, 0.11909567136143454, 0.009155740220708433, 0.029322113201317515, 0.13315455373827445, 0.227422215281317])
      cub_degree = 8
      tol = 1e-14
      # cub = SymCubatures.TriSymCub{T}(vertices=false, numedge=2,
      #                                 midedges=true,
      #                                 numS21=0, numS111=2, centroid=true)
      # SymCubatures.setparams!(cub, T[0.5*(1 +(1/3)*sqrt(5 - 2*sqrt(10/7)));
      #                                0.5*(1 +(1/3)*sqrt(5 + 2*sqrt(10/7)));
      #                                0.22099843842186342;
      #                                1.1631912073287645;
      #                                0.24591837943530243;
      #                                0.12484120739275503])
      # SymCubatures.setweights!(cub, T[0.03946097492219484;
      #                                 0.022750835455651795;
      #                                 0.008818681441998117;
      #                                 0.18379350807584505;
      #                                 0.055866234615190226;
      #                                 0.2542415177018296])
      # cub_degree = 8
      # tol=1e-14
    elseif q <= 9 #34 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 1,
                                      numedge = 3,
                                      numS111 = 2,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.4490679047936404, 0.033765242898423975, 0.16939530676686776, 0.38069040695840156, 0.13350211698690437, 0.2551152473178156, 0.6591253866697547, 0.18678645216668188])
      SymCubatures.setweights!(cub, T[0.0952257519287601, 0.006658503841344434, 0.016313739546142424, 0.025165807591009056, 0.06038167775800903, 0.14235832969148426, 0.20905439364578426])
      cub_degree = 9
      tol = 1e-14
      #36 nodes
      # cub = SymCubatures.TriSymCub{T}(vertices = false,
      #                                 midedges = false,
      #                                 numS21 = 2,
      #                                 numedge = 3,
      #                                 numS111 = 2,
      #                                 centroid = false)
      # SymCubatures.setparams!(cub, T[0.5267710694668858, 0.22581307560584146, 0.033765242898423975, 0.16939530676686776, 0.38069040695840156, 0.6556337684088949, 0.1729852501295953, 0.23384053251428352, 0.07010081712018432])
      # SymCubatures.setweights!(cub, T[0.16521286119089668, 0.09190991876959623, 0.006155310695930864, 0.015651771073189685, 0.02308837408072088, 0.13970715125683908, 0.020169336246406353])
      # cub_degree = 9
      # tol = 1e-14
    elseif q <= 10 #39 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 1,
                                      numedge = 3,
                                      numS111 = 3,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.07910043691664528, 0.033765242898423975, 0.16939530676686776, 0.38069040695840156, 0.7326534517622117, 0.16178998830038974, 0.15435874951213008, 0.3673737616984141, 0.47137911926303555, 0.5572798360390397])
      SymCubatures.setweights!(cub, T[0.03388444217225276, 0.0008514906758374474, 0.015608482275332266, 0.021266529714946317, 0.09830017774056278, 0.08553487986793913, 0.09482955197258906])
      cub_degree = 10
      tol = 1e-14
    elseif q <= 11 #42 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 5,
                                      numedge = 3,
                                      numS111 = 1,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.35406264500305373, 0.5349706832314904, 0.14620092105924026, 0.9344133285246975, 0.8608802427946676, 0.025446043828620757, 0.12923440720030277, 0.2970774243113014, 0.11288494321746362, 0.5189980362860592])
      SymCubatures.setweights!(cub, T[0.019568317008650336, 0.09027010441205281, 0.12668034231431224, 0.06200006585228387, 0.06976602848889582, 0.0831816128947378, 0.003678928929350306, 0.012519064458309727, 0.009408240262222357, 0.08199386419798452])
      cub_degree = 11
      tol = 1e-14
    elseif q <= 12 #49 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 3,
                                      numedge = 3,
                                      numS111 = 3,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.944532850510785, 0.061319596747426944, 0.3670399210288128, 0.025446043828620757, 0.12923440720030277, 0.2970774243113014, 0.9643317137065957, 0.6791559176167298, 0.2800879551160052, 0.1108945564420906, 0.11759149845575509, 0.5898408487390895])
      SymCubatures.setweights!(cub, T[0.013033842870623723, 0.06517521728204109, 0.018835815821139735, 0.09104423407752954, 0.0007806066251635613, 0.008577761805886172, 0.012522944584712063, 0.08556960015635551, 0.04947754382659775, 0.0651662877443358, 0.10316420138769361])
      cub_degree = 12
      tol = 1e-14

    elseif q <= 13 #54 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 4,
                                      numedge = 4,
                                      numS111 = 3,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.08465074046995401, 0.9461615506412667, 0.674767791012365, 0.323267179957891, 0.019855071751231856, 0.10166676129318664, 0.2372337950418355, 0.4082826787521751, 0.5887945858454915, 0.10444119220651123, 0.6453371751137852, 0.34019090826916276, 0.09961977082750291, 0.29042446777707714])
      SymCubatures.setweights!(cub, T[0.021591599231498757, 0.06539932723237385, 0.041145418196861044, 0.07713054222104607, 0.0017474343859017836, 0.0059951229414363285, 0.008892253927912585, 0.010775615503604317, 0.05805264312502276, 0.10317694130545241, 0.0420598787031132])
      cub_degree = 13
      tol = 1e-14
    elseif q <= 14 #60 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 6,
                                      numedge = 4,
                                      numS111 = 3,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.3438598694407348, 0.08667709181254661, 0.23438830525749135, 0.5520607432679493, 0.9587852382179086, 0.8565327499763229, 0.019855071751231856, 0.10166676129318664, 0.2372337950418355, 0.4082826787521751, 0.09143726717923181, 0.31045214698586704, 0.24432313634373679, 0.5801444096258322, 0.07777721618947092, 0.6080084479074117])
      SymCubatures.setweights!(cub, T[0.054542527743094635, 0.023228939031591043, 0.02513762015408967, 0.10503584963092831, 0.05101963837324333, 0.07707180717717142, 0.0017455842845511217, 0.006113152456609967, 0.007165616322632858, 0.007968949263384184, 0.03911755531627651, 0.06148029583288167, 0.04172398880193796])
      cub_degree = 14
      tol = 1e-14
    elseif q <= 15 #69 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 4,
                                      numedge = 4,
                                      numS111 = 5,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.07246516124875677, 0.9597333330363812, 0.5578271063863942, 0.25839266543174033, 0.015919880246186957, 0.08198444633668206, 0.19331428364970477, 0.33787328829809554, 0.08945010025247659, 0.7665909769206356, 0.890366171837056, 0.8251504415876004, 0.27617698719951544, 0.5241882285519328, 0.24339513574480046, 0.08051715832948862, 0.4824242225832316, 0.08250669708663373])
      SymCubatures.setweights!(cub, T[0.007975211738078417, 0.01555795199453867, 0.0142287699046545, 0.09767682471952172, 0.04995668055544628, 0.0011984999191526156, 0.004057653525909466, 0.0057771570575940985, 0.007510437420733494, 0.041736293516403324, 0.043023900670943645, 0.07279841199275705, 0.027493765777003604, 0.03703949399671621])
      cub_degree = 15
      tol = 1e-14
    elseif q <= 16 #72 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 5,
                                      numedge = 4,
                                      numS111 = 5,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.04538479805133307, 0.22611045578916183, 0.7597184979885856, 0.4715660213460911, 0.962698662025747, 0.015919880246186957, 0.08198444633668206, 0.19331428364970477, 0.33787328829809554, 0.6642487752933809, 0.07317173680920826, 0.06769959788003582, 0.18359591467251551, 0.23363169305593603, 0.45309708381258046, 0.07154380405671455, 0.3955820133312661, 0.7299488157852585, 0.2402324353785272])
      SymCubatures.setweights!(cub, T[0.0068099595522426995, 0.007902251972551224, 0.03822679813612106, 0.08015956854713617, 0.07195779765864842, 0.03751715583277405, 0.0006453194613776787, 0.0033040868832377697, 0.0051468228469984285, 0.006303392104459153, 0.035057823925930626, 0.019961382833761403, 0.051626832885843146, 0.028924233456013058, 0.06107667308597509])
      cub_degree = 16
      tol = 1e-14

    elseif q <= 17 #81 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 5,
                                      numedge = 5,
                                      numS111 = 6,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.9673783153130373, 0.2285847353426557, 0.05840202475924722, 0.6899436810572654, 0.42617465029864315, 0.013046735741414128, 0.06746831665550773, 0.16029521585048778, 0.2833023029353764, 0.4255628305091844, 0.45997260718861244, 0.21436618225556944, 0.21119339160494885, 0.7433245013095352, 0.41869075536293493, 0.6667211762963206, 0.41159337035640003, 0.06703091561356278, 0.2024719492716725, 0.07034605000725505, 0.673480779947899, 0.06531116254401612])
      SymCubatures.setweights!(cub, T[0.03235009501507383, 0.03894783481393893, 0.010304291352971073, 0.023353742881921363, 0.05363694764402083, 0.000780347117914888, 0.0028284869642989496, 0.004156670349616389, 0.004888093779533072, 0.00532868081300272, 0.047922726227699664, 0.05405338591315981, 0.056272829027598656, 0.026314129806583212, 0.020943399709124497, 0.030548127770838585])
      cub_degree = 17
      tol = 1e-14
    elseif q <= 18 #93 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = false,
                                      numS21 = 5,
                                      numedge = 5,
                                      numS111 = 8,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.9733910151194378, 0.21626534843071651, 0.05150550298873318, 0.736133255549125, 0.37909135549530876, 0.013046735741414128, 0.06746831665550773, 0.16029521585048778, 0.2833023029353764, 0.4255628305091844, 0.3896301312622026, 0.1746624859673528, 0.7730831662450132, 0.08870182260636436, 0.21955949136089495, 0.978215349126951, 0.19128424042417297, 0.5706040919579012, 0.3950618272652804, 0.6260607152152675, 0.39506027304380237, 0.05823520227247069, 0.18823228992001764, 0.06686354946549783, 0.628350432052669, 0.04462834824245546])
      SymCubatures.setweights!(cub, T[0.021617628605668425, 0.030859684307638195, 0.008492717071468023, 0.048302368344678256, 0.05248145434612612, 0.0006756820364337081, 0.0026372964003534435, 0.003926222398471446, 0.003254151514351532, 0.00491785483632988, 0.025721171124103152, 0.020968675195521817, 0.035338100044051406, 0.037915704809379495, 0.057440515681043525, 0.021551248401012964, 0.019465345373257455, 0.0186444391812341])
      cub_degree = 18
      tol = 5e-13
    elseif q <= 19 #96 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 5,
                                      numedge = 5,
                                      numS111 = 8,
                                      centroid = false)
      SymCubatures.setparams!(cub, T[0.5856059294178624, 0.047740969681407466, 0.36491553785928776, 0.818481134242735, 0.9099765582336337, 0.010885670926971514, 0.05646870011595234, 0.13492399721297532, 0.2404519353965941, 0.36522842202382755, 0.344515016089839, 0.05634768668558025, 0.3892623612407014, 0.1816493745951281, 0.0584137152546641, 0.16730180032945222, 0.5835945178396305, 0.3615699266691124, 0.19685560078378486, 0.1865954834733815, 0.05507891514214126, 0.8345201463254764, 0.17940998372119654, 0.635984563672555, 0.5715152111527353, 0.05498532291050753])
      SymCubatures.setweights!(cub, T[0.004144191713042269, 0.05671201687778216, 0.006952954178298502, 0.04182915174250011, 0.049236730600897276, 0.04362659548797906, 0.0005315273623156564, 0.0019591053495906494, 0.002952187933786907, 0.0035504848806637655, 0.003981570473101632, 0.018993243389823684, 0.03526459339555751, 0.014596159925672434, 0.047289681869773156, 0.014035668547812305, 0.025039659732181335, 0.041217429530185055, 0.022671200642619533])
      cub_degree = 19
      tol = 1e-14
    elseif q <= 20 #103 nodes
      cub = SymCubatures.TriSymCub{T}(vertices = false,
                                      midedges = true,
                                      numS21 = 5,
                                      numedge = 5,
                                      numS111 = 9,
                                      centroid = true)
      SymCubatures.setparams!(cub, T[0.5436910855480633, 0.040950582871744044, 0.33887371496698127, 0.8109732389249305, 0.9624177194107271, 0.010885670926971514, 0.05646870011595234, 0.13492399721297532, 0.2404519353965941, 0.36522842202382755, 0.3334872800151244, 0.05388954156501192, 0.3895234197816185, 0.1718277785489534, 0.05577562117440405, 0.15389464625504432, 0.599253467648353, 0.14153780032196547, 0.1643849285186809, 0.21519353620251885, 0.5588940962149358, 0.3266182856250745, 0.04660262406085601, 0.801240166454181, 0.18939648231227538, 0.7890574015599582, 0.5525675668379594, 0.04329675777029787])
      SymCubatures.setweights!(cub, T[0.0043014134434171, 0.051486376804974164, 0.005540917165528219, 0.03694293372504983, 0.05781723615931413, 0.016819298793165075, 0.00044140660013591816, 0.0018162275341354346, 0.002944300406209069, 0.0028771277952727914, 0.003198212554961351, 0.0180917026903521, 0.029958760377674705, 0.013668758583780152, 0.02893920761789538, 0.01556903102368063, 0.04880488827336611, 0.01915230387979439, 0.035450671026311734, 0.017313930955951226, 0.05191629580852894])
      cub_degree = 20
      tol = 5e-14
    else
      error("polynomial degree must be <= 20 (presently)\n")
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
  weights[:] *= 2.0/sum(wts)

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
    # SymCubatures.setweights!(cub, T[0.0015106273303336273,0.060490542374353584,
    #                                 0.004038881996228382, 0.10344930834722398,
    #                                 0.005696088152131421, 0.02424296133613638,
    #                                 0.08113091859465722])
    # SymCubatures.setparams!(cub, T[0.28418700275470193,0.21742832019555544,
    #                                0.25737274681480826,0.45008848310824695])
    SymCubatures.setweights!(cub, T[0.0015106273303336273,0.060490542374353584,
                                    0.004038881996228382, 0.10344930834722398,
                                    0.02424296133613638,0.005696088152131421,
                                    0.08113091859465722])
    SymCubatures.setparams!(cub, T[0.28418700275470193,0.21742832019555544,
                                   0.45008848310824695,0.25737274681480826])
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
* `faceopertype`: the operator type on the facets of the tetrahedron
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right tetrahedron
* `vtx`: vertices for the right tetrahedron

"""
function getTetCubatureDiagE(q::Int, T=Float64; faceopertype::Symbol = :DiagE, tol=10*eps(typeof(real(one(T)))))
  cub_degree = q
  mask = zeros(Int64, (0))

  if faceopertype == :Omega
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
      
      cub_degree = 8
      tol = 1e-14      
    else
      error("polynomial degree must be <= 4 (presently)\n")
    end
  else
    if q <=1 # 6 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = false,
                                      midedges = true,
                                      facecentroid = false,
                                      numedge = 0,
                                      numfaceS21 = 0,
                                      numfaceS111 = 0,
                                      centroid = false,
                                      numS31 = 0,
                                      numS22 = 0,
                                      numS211= 0,
                                      numS1111 =0)
      SymCubatures.setweights!(cub, T[0.22222222222222235])
      cub_degree = 1
      tol = 1e-14
    elseif q <= 2 #7 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = false,
                                      midedges = true,
                                      facecentroid = false,
                                      numedge = 0,
                                      numfaceS21 = 0,
                                      numfaceS111 = 0,
                                      centroid = true,
                                      numS31 = 0,
                                      numS22 = 0,
                                      numS211= 0,
                                      numS1111 =0)
      SymCubatures.setweights!(cub, T[0.1333333333333337, 0.5333333333333319])
      cub_degree = 2
      tol = 1e-14

    # elseif q <= 3
      # #23 nodes based on 7 node facet operator
      # cub = SymCubatures.TetSymCub{T}(vertices = false,
      #                                 midedges = true,
      #                                 facecentroid = true,
      #                                 numedge = 0,
      #                                 numfaceS21 = 1,
      #                                 numfaceS111 = 0,
      #                                 centroid = true,
      #                                 numS31 = 0,
      #                                 numS22 = 0,
      #                                 numS211= 0,
      #                                 numS1111 =0)
      # SymCubatures.setparams!(cub, T[0.22222222222222215])
      # SymCubatures.setweights!(cub, T[0.01225481148524747, 0.02799913080121597, 0.16219856069976812, 0.27502065200818465])
      # cub_degree = 3
      # tol = 1e-14

    elseif q <=4 #23 nodes (for q=3 and q=4)
      cub = SymCubatures.TetSymCub{T}(vertices = true,
                                      midedges = true,
                                      facecentroid = false,
                                      numedge = 0,
                                      numfaceS21 = 1,
                                      numfaceS111 = 0,
                                      centroid = true,
                                      numS31 = 0,
                                      numS22 = 0,
                                      numS211= 0,
                                      numS1111 =0)
      SymCubatures.setparams!(cub, T[0.3771609693928902])
      SymCubatures.setweights!(cub, T[0.0017328144234984297, 0.010017745385041531, 0.07166219974832357, 0.40634920634920585])
      cub_degree = 4
      tol = 1e-14

      # #26 nodes based on 7 node facet operators
      # cub = SymCubatures.TetSymCub{T}(vertices = false,
      #                                 midedges = true,
      #                                 facecentroid = true,
      #                                 numedge = 0,
      #                                 numfaceS21 = 1,
      #                                 numfaceS111 = 0,
      #                                 centroid = false,
      #                                 numS31 = 1,
      #                                 numS22 = 0,
      #                                 numS211= 0,
      #                                 numS1111 =0)
      # SymCubatures.setparams!(cub, T[0.4842953457301756, 0.22222222222222215])
      # SymCubatures.setweights!(cub, T[0.18290740424677793, 0.022508956569022617, 0.01644731543334552, 0.06732054793298498])
      # cub_degree = 4
      # tol = 1e-14
    elseif q <= 5 #44 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = true,
                midedges = false,
                facecentroid = false,
                numedge = 1,
                numfaceS21 = 2,
                numfaceS111 = 0,
                centroid = false,
                numS31 = 1,
                numS22 = 0,
                numS211= 0,
                numS1111 =0)
      SymCubatures.setparams!(cub, T[0.5008941915142769, 0.8506802519794945, 0.23722737279318576, 0.3077459416259917])
      SymCubatures.setweights!(cub, T[0.0015673886232196292, 0.17081759879508043, 0.033441261076507856, 0.01477813407660693, 0.00543005348522964])
      cub_degree = 5
      tol = 1e-14
    elseif q <= 6 #51 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = true,
                midedges = false,
                facecentroid = false,
                numedge = 1,
                numfaceS21 = 2,
                numfaceS111 = 0,
                centroid = true,
                numS31 = 1,
                numS22 = 1,
                numS211= 0,
                numS1111 =0)
      SymCubatures.setparams!(cub, T[0.393148466877178, 0.25786391856319707, 0.8506802519794945, 0.23722737279318576, 0.3077459416259917])
      SymCubatures.setweights!(cub, T[0.0008311992138364736, 0.067074516283967, 0.08700624630786795, 0.025281032406421995, 0.01448700395673796, 0.0046508273850214945, 0.006646628516734983])
      cub_degree = 6
      tol = 1e-14

    elseif q<=7
      #76 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = true,
                                      midedges = true,
                                      facecentroid = true,
                                      numedge = 1,
                                      numfaceS21 = 1,
                                      numfaceS111 = 1,
                                      centroid = false,
                                      numS31 = 2,
                                      numS22 = 1,
                                      numS211= 0,
                                      numS1111 =0)
      # SymCubatures.setparams!(cub, T[0.9465552688408668, 0.326390241974141, 0.7461816901961349, 0.160991838340075, 0.8003678928805428, 0.6058255660767269, 0.215183643569735])
      # SymCubatures.setweights!(cub, T[0.0008255000130940284, 0.01874509335958128, 0.06004170145017678, 0.0019891243235202693, 0.09315724230205255, 0.0030561455241973827, 0.0021889210067103236, 0.013884472504428686, 0.011959453952826883])
      SymCubatures.setparams!(cub, T[0.3649069878051798, 0.9069464284900626, 0.19709340992078445, 0.160991838340075, 0.8003678928805428, 0.6058255660767269, 0.215183643569735])
      SymCubatures.setweights!(cub, T[0.0003131195741790005, 0.06872580049477045, 0.0770182765903434, 0.001297874862616682, 0.05122298827420443, 0.004967708556549984, 0.001966184574009691, 0.013198465506081449, 0.00850236954064094])
      cub_degree = 7
      tol=1e-14
    elseif q<=8
      #89 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = true,
                                      midedges = true,
                                      facecentroid = true,
                                      numedge = 1,
                                      numfaceS21 = 1,
                                      numfaceS111 = 1,
                                      centroid = true,
                                      numS31 = 2,
                                      numS22 = 1,
                                      numS211= 1,
                                      numS1111 =0)
      SymCubatures.setparams!(cub, T[0.26967567337250364, 0.9186565269236643, 0.7878184728172669, 0.16099183834007277, 0.800367892880542, 0.19236824251877926, 1.0726528991802897, 0.6058255660767288, 0.2151836435697346])
      SymCubatures.setweights!(cub, T[0.00017788261247842277, 0.022946869343355818, 0.05842047515941764, 0.0025968334566281045, 0.012039117612706244, 0.00471464641426865, 0.0014485783932065817, 0.041098991084533104, 0.010486947759748596, 0.01215056364118055, 0.051901126953532224])
      cub_degree = 8
      tol=1e-14

    elseif q<=9
      #121 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = true,
                midedges = true,
                facecentroid = true,
                numedge = 1,
                numfaceS21 = 3,
                numfaceS111 = 1,
                centroid = true,
                numS31 = 1,
                numS22 = 1,
                numS211= 2,
                numS1111 =0)
      SymCubatures.setparams!(cub, T[0.2595307420691634, 0.8441512685303726, 0.40849548761177584, 0.1693510205976773, 0.895541492424858, 0.9192858150590572, 0.1362955833435948, 0.5043846222657248, 0.4701004391173425, 0.1604772007159291, 0.5199457984883094, 0.07662224108118755])
      SymCubatures.setweights!(cub, T[0.00015279459961010358, 0.01975401788383916, 0.001314202393676677, 0.021390096706851512, 0.006798995776920046, 0.0039159324894427845, 0.006885934192465157, 0.0005069027474549871, 0.020477522536811944, 0.04175001086337507, 0.0023605549505845433, 0.009102594817477005, 0.060393007434791986])
      cub_degree = 9
      tol=1e-14

    elseif q<=10
      #145 nodes
      cub = SymCubatures.TetSymCub{T}(vertices = true,
                midedges = true,
                facecentroid = true,
                numedge = 1,
                numfaceS21 = 3,
                numfaceS111 = 1,
                centroid = true,
                numS31 = 1,
                numS22 = 1,
                numS211= 4,
                numS1111 =0)
      SymCubatures.setparams!(cub, T[0.5700558294241873, 0.8819914623822424, 0.40849548761177584, 0.1693510205976773, 0.895541492424858, 0.9192858150590572, 0.3674432698076532, 0.12171080589025479, 0.7415577501465243, 0.1346989794996326, 0.10308419269846046, 0.5052883884978763, 0.11223871825337459, 0.203319497342745, 0.5199457984883094, 0.07662224108118755])
      SymCubatures.setweights!(cub, T[8.342292637650742e-5, 0.02995257400472967, 0.001178939330700515, 0.014024455736117564, 0.005665640945854503, 0.0025319848901773656, 0.005490413383261791, 0.0006215155149072649, 0.021063871513121817, 0.02981187131916029, 0.013879210846470981, 0.005258641145889365, 0.0016466522072965662, 0.0074534718610153455, 0.04075764008270213])
      cub_degree = 10
      tol=1e-14
      # #145 nodes
      # cub = SymCubatures.TetSymCub{T}(vertices = true,
      #           midedges = true,
      #           facecentroid = true,
      #           numedge = 1,
      #           numfaceS21 = 3,
      #           numfaceS111 = 1,
      #           centroid = true,
      #           numS31 = 1,
      #           numS22 = 1,
      #           numS211= 4,
      #           numS1111 =0)
      # SymCubatures.setparams!(cub, T[0.5831510926071354, 0.855295628963342, 0.40849548761177584, 0.1693510205976773, 0.895541492424858, 0.9192858150590572, 0.34277718068052015, 0.11217384327295098, 0.7459376222015264, 0.1271800044965199, 0.09162790449408253, 0.5994079435205514, 0.07787143633861693, 0.24259229663883067, 0.5199457984883094, 0.07662224108118755])
      # SymCubatures.setweights!(cub, T[0.00021816263981043935, 0.03921759527783465, 0.0011043100191260003, 0.009765154777214433, 0.005209712411026045, 0.0015442679048136543, 0.00532449820711722, 0.0003501468689026116, 0.022900629369277666, 0.0305175435641846, 0.013142498713694687, 0.006879720547610928, 0.0012142200424160243, 0.006991964941085668, 0.022836161062856575])
      # cub_degree = 10
      # tol=1e-14

      # #139 nodes
      # cub = SymCubatures.TetSymCub{T}(vertices = true,
      #           midedges = true,
      #           facecentroid = true,
      #           numedge = 1,
      #           numfaceS21 = 3,
      #           numfaceS111 = 1,
      #           centroid = true,
      #           numS31 = 1,
      #           numS22 = 0,
      #           numS211= 4,
      #           numS1111 =0)
      # SymCubatures.setparams!(cub, T[0.5803333201819036, 0.40849548761177584, 0.1693510205976773, 0.895541492424858, 0.9192858150590572, 0.330429102391519, 0.10730337252239175, 0.7539038145106844, 0.12433214114387407, 0.09670411438362254, 0.638953493007011, 0.06568839621878922, 0.23886759837771818, 0.5199457984883094, 0.07662224108118755])
      # SymCubatures.setweights!(cub, T[0.0003277864005158543, 0.04215578248651099, 0.000992857927615538, 0.00513620782204975, 0.001120039985381197, 0.004759201791331538, 4.93013126328509e-5, 0.022445330296939237, 0.03376662009972908, 0.015646220106639838, 0.006871141551676767, 0.0013019041677660326, 0.006840093117228722, 0.021307082127670487])
      # cub_degree = 10
      # tol=1e-14

    else
      error("polynomial degree must be <= 10 (presently)\n")
    end
  end

  mask = SymCubatures.getInternalParamMask(cub)
  append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))
  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
  Cubature.solvecubature!(cub, cub_degree, mask, tol=tol, hist=false)
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
function equivalenceconstant(cub::TriSymCub{T}, vtx::Array{T,2}, q::Int) where {T}
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
  A = P'*H*P
  return eigmax(0.5.*(A + A'))
end  

function equivalenceconstant(cub::TetSymCub{T}, vtx::Array{T,2}, q::Int) where {T}
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
  A = P'*H*P
  return eigmax(0.5.*(A + A'))
end

end
