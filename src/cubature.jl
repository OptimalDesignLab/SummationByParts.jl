module Cubature
# routines for constructing cubatures

using ..OrthoPoly
using ..SymCubatures

export tricubature, tetcubature

@doc """
### Cubature.cubatureresidual{TetSymCub{T}}

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
  x, y = SymCubatures.calcnodes(cub, vtx)
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
      P = OrthoPoly.proriolpoly(x, y, i, j)
      F[ptr] += (w.'*P)[1]
      dPdx, dPdy = OrthoPoly.diffproriolpoly(x, y, i, j)
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
  x, y, z = SymCubatures.calcnodes(cub, vtx)
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
        P = OrthoPoly.proriolpoly(x, y, z, i, j, k)
        F[ptr] += (w.'*P)[1]
        dPdx, dPdy, dPdz = OrthoPoly.diffproriolpoly(x, y, z, i, j, k)
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
defined by `cub`, which is a parametric abstract type (see domainsymorbits.jl).

**Inputs**

* `q`: maximum (desired) degree for which the cubature is exact
* `tol`: tolerance with which to solve the accuracy conditions
* `hist`: if true, print the residual-norm convergence history

**In/Outs**

* `cub`: on entry, defines the initial guess for the cubature nodes and weights.
  on exit, defines the nodes and weights that satisfy the desired accuracy.

"""->
function solvecubature!{T}(cub::SymCub{T}, q::Int; tol=eps(T(10)),
                           hist::Bool=false)
  Jac = SymCubatures.calcjacobian(cub)

  # compute accuracy for initial guess 
  F, dF = Cubature.cubatureresidual(cub, q)
  res = norm(F)
  res0 = res
  res_old = res
  hist ? print("calccubature:\n") : nothing
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
  for k = 1:maxiter
    JtJ = Jac.'*dF.'*dF*Jac
    H = JtJ + nu*diagm(diag(JtJ))
    g = -Jac.'*dF.'*F
    dv = H\g

    # update cubature definition and check for convergence
    v += dv
    SymCubatures.setparams!(cub, v[1:cub.numparams])
    SymCubatures.setweights!(cub, v[cub.numparams+1:end])
    F, dF = Cubature.cubatureresidual(cub, q)
    res = norm(F)
    hist ? print("\titer ",k,": res norm = ",res,"\n") : nothing
    if (res < tol)
      return
    end

    # trust-region like update
    if (res > res_old)
      v -= dv
      SymCubatures.setparams!(cub, v[1:cub.numparams])
      SymCubatures.setweights!(cub, v[cub.numparams+1:end])
      F, dF = Cubature.cubatureresidual(cub, q)
      nu *= 4
    else
      nu /= 2
      res_old = res
    end

  end
  error("calccubature failed to find solution in ",maxiter," iterations")
end

@doc """
### Cubature.tricubature{T}

This high-level function computes and returns a symmetric cubature of requested
accuracy on the right triangle.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""->
function tricubature(q::Int, T=Float64; tol=eps(T(10)))
  @assert( q >= 1 && q <= 7 && mod(q,2) == 1 )
  if q == 1
    # P1 (vertices only); 2nd order cubature
    cub = SymCubatures.TriSymCub{T}() 
    SymCubatures.setweights!(cub, T[2/3])
  elseif q == 3
    # P2 + 1 bubble node; 4th order cubature
    cub = SymCubatures.TriSymCub{T}(midedges=true, centroid=true)
    SymCubatures.setweights!(cub, T[1/10, 4/15, 9/10])
  elseif q == 5
    # P3 + 3 bubble nodes; 6th order cubature
    cub = SymCubatures.TriSymCub{T}(numedge=1, numS21=1)
    SymCubatures.setweights!(cub, T[0.02974582604964118,0.09768336246810204,
                                    0.4415541156808217])
    SymCubatures.setparams!(cub, T[0.29346955590904017,0.41469035132718185])
  elseif q == 7
    # P4 + 6 bubble nodes; 8th order cubature
    cub = SymCubatures.TriSymCub{T}(midedges=true, numedge=1, numS21=2)
    SymCubatures.setweights!(cub, T[0.012698412698412695,0.05079365079365077,
                                    0.04285714285714284,0.2023354595827503,
                                    0.3151248578775673])
    SymCubatures.setparams!(cub, T[0.2113248654051872,0.2615831876594899,
                                   0.8495279234516212])
  elseif q == 9
    # P5 + 10 bubble nodes; 10th order cubature
    cub = SymCubatures.TriSymCub{T}(numedge=2, centroid=true, numS21=1,
                                    numS111=1)
    SymCubatures.setweights!(cub, T[0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
    SymCubatures.setparams!(cub, T[0.1 0.25 0.1 0.2 0.6])
  else
    error("polynomial degree must be 1, 3, 5, or 7 (presently)\n")
  end
  vtx = T[-1 -1; 1 -1; -1 1]
  Cubature.solvecubature!(cub, q, tol=tol)
  return cub, vtx
end

@doc """
### Cubature.tetcubature{T}

This high-level function computes and returns a symmetric cubature of requested
accuracy on the right tetrahedron.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `cub`: a symmetric cubature for the right tetrahedron
* `vtx`: vertices for the right tetrahedron

"""->
function tetcubature(q::Int, T=Float64; tol=eps(T(10)))
  @assert( q >= 1 && q <= 7 && mod(q,2) == 1)
  if q == 1
    # P1 (vertices only); 2nd order cubature
    cub = SymCubatures.TetSymCub{T}()
    SymCubatures.setweights!(cub, T[1/3])
  elseif q == 3
    # P2 + 1 bubble node; 4th order cubature
    cub = SymCubatures.TetSymCub{T}(midedges=true, centroid=true)
    SymCubatures.setweights!(cub, T[1/45 4/45 32/45])
  elseif q == 5
    # P3 + 4 bubble nodes; 6th order cubature
    cub = SymCubatures.TetSymCub{T}(facecentroid=true,
                                    numedge=1, numS31=1)
    SymCubatures.setweights!(cub, T[0.004421633248304776 0.0176754534336105 
                                    0.06935370366814568 0.20653163611605146])
    SymCubatures.setparams!(cub, T[0.30480589839889616 0.45720884759834435])
  elseif q == 7
    # P3 + 11 bubble nodes; 8th order cubature
    cub = SymCubatures.TetSymCub{T}(midedges=true, centroid=true, numedge=1,
                                    numfaceS21=1, numS31=1, numS22=1)
    SymCubatures.setweights!(cub, T[0.0015106273303336273,0.004038881996228382,
                                    0.005696088152131421,0.02424296133613638,
                                    0.08113091859465722,0.060490542374353584,
                                    0.10344930834722398])
    SymCubatures.setparams!(cub, T[0.7426272531851917,0.45008848310824695,
                                   0.28418700275470193,0.21742832019555544])
  else
    error("polynomial degree must be 1, 3, 5, or 7 (presently)\n")
  end
  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
  Cubature.solvecubature!(cub, q, tol=tol)
  return cub, vtx
end

end