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
  maxiter = 500
  nu = 100.0
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

This high-level function computes and returns the weights and nodes for a
cubature of requested accuracy on the right-triangle.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `w`: cubature weights
* `x`: cubature node coordinates (x,y)

"""->
function tricubature(q::Int, T=Float64; tol=eps(T(10)))
  @assert( q >= 1 && q <= 7 && mod(q,2) == 1 )
  if q == 1
    # P1 (vertices only); 2nd order cubature
    cub = SymCubatures.TriSymCub{T}() 
    SymCubatures.setweights!(cub, T[0.5])
  elseif q == 3
    # P2 + 1 bubble node; 4th order cubature
    cub = SymCubatures.TriSymCub{T}(midedges=true, centroid=true)
    SymCubatures.setweights!(cub, T[0.5, 0.5, 0.5])
  elseif q == 5
    # P3 + 3 bubble nodes; 6th order cubature
    cub = SymCubatures.TriSymCub{T}(numedge=1, numS21=1)
    SymCubatures.setweights!(cub, T[0.5, 0.5, 0.5])
    SymCubatures.setparams!(cub, T[0.25, 0.25])
  elseif q == 7
    # P4 + 6 bubble nodes; 8th order cubature
    cub = SymCubatures.TriSymCub{T}(midedges=true, numedge=1, numS21=2)
    SymCubatures.setweights!(cub, T[0.5, 0.5, 0.5, 0.5, 0.5])
    SymCubatures.setparams!(cub, T[0.25, 0.25, 0.75])
  else
    error("polynomial degree must be 1, 3, 5, or 7 (presently)\n")
  end
  vtx = T[-1 -1; 1 -1; -1 1]
  Cubature.solvecubature!(cub, q, tol=tol)
  w = SymCubatures.calcweights(cub)
  x = zeros(T, (cub.numnodes, 2) )
  x[:,1], x[:,2] = SymCubatures.calcnodes(cub, vtx)
  return w, x
end

@doc """
### Cubature.tetcubature{T}

This high-level function computes and returns the weights and nodes for a
cubature of requested accuracy on the right-tetrahedron.

**Inputs**

* `q`: maximum degree of polynomial for which the cubature is exact
* `tol`: tolerance with which to solve the cubature

**Outputs**

* `w`: cubature weights
* `x`: cubature node coordinates (x,y,z)

"""->
function tetcubature(q::Int, T=Float64; tol=eps(T(10)))
  @assert( q >= 1 && q <= 5 && mod(q,2) == 1)
  if q == 1
    # P1 (vertices only); 2nd order cubature
    cub = SymCubatures.TetSymCub{T}()
    SymCubatures.setweights!(cub, T[0.5])
  elseif q == 3
    # P2 + 1 bubble node; 4th order cubature
    cub = SymCubatures.TetSymCub{T}(midedges=true, centroid=true)
    SymCubatures.setweights!(cub, T[0.5, 0.5, 0.5])
  elseif q == 5
    # P3 + 4 bubble nodes; 6th order cubature
    cub = SymCubatures.TetSymCub{T}(facecentroid=true,
                                    numedge=2, numS31=1)
    SymCubatures.setweights!(cub, T[0.1 0.1 0.1 0.1 0.1])
    SymCubatures.setparams!(cub, T[1/5 4/5 5/6])
  else
    error("polynomial degree must be 1, 3, or 5 (presently)\n")
  end
  vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
  Cubature.solvecubature!(cub, q, tol=tol)
  w = SymCubatures.calcweights(cub)
  x = zeros(T, (cub.numnodes, 3))
  x[:,1], x[:,2], x[:,3] = SymCubatures.calcnodes(cub, vtx)
  return w, x
end

end