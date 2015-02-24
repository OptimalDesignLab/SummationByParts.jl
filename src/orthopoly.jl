module OrthoPoly
# utilities for working with orthogonal polynomials

@doc """
### OrthoPoly.lglnodes

Computes the Legendre-Gauss-Lobatto (LGL) quadrature nodes and weights and the
interval [-1,1].  The LGL nodes are the zeros of (1-x^2)*P'_N(x), where P_N(x)
denotes the Nth Legendre polynomial.

*Reference*: C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, \"Spectral 
Methods in Fluid Dynamics,\" Section 2.3. Springer-Verlag 1987

**Inputs**

* `N`: highest degree (number of nodes = N+1)
* `T`: number type

**Outputs**

* `x`: the LGL nodes
* `w`: the LGL weights

Julia version adapted from Matlab code written by Greg von Winckel - 04/17/2004
Contact: gregvw@chtm.unm.edu
 
"""->
function lglnodes(N, T=Float64)
N1 = N+1
# Use the Chebyshev-Gauss-Lobatto nodes as an initial guess
x = -cos(pi*[0:N;]/N)
# The Legendre Vandermonde Matrix 
P = zeros(T, (N1,N1))
# Compute P_(N) using the recursion relation; compute its first and second
# derivatives and update x using the Newton-Raphson method.
xold = (T)(2)
while maxabs(real(x-xold)) > eps(real((T)(1)))
  xold = x
  P[:,1] = one(T)
  P[:,2] = x
  for k=2:N
    P[:,k+1]= ((2k-1)*x.*P[:,k]-(k-1)*P[:,k-1])/k
  end
  x = xold - ( x.*P[:,N1]-P[:,N] )./( N1*P[:,N1] )           
end
w = 2./(N*N1*P[:,N1].^2)
return x, w
end

@doc """
### OrthoPoly.jacobipoly{T}

Evaluate a Jacobi polynomial at some points.  Based on JacobiP in Hesthaven and
Warburton's nodal DG book.

**Inputs**

* `x`: points at which to evaluate polynomial
* `alpha`,`beta`: define the type of Jacobi Polynomial (alpha + beta != 1)
*  `N`: polynomial degree

**Outputs**

* `P`: the polynomial evaluated at x

"""->
function jacobipoly{T}(x::Array{T}, alpha::FloatingPoint, beta::FloatingPoint,
                       N::Int)
  @assert( alpha + beta != -1 )
  @assert( alpha > -1 && beta > -1)
  # Initial values P_0(x) and P_1(x)
  gamma0 = ((2^(alpha+beta+1))/(alpha+beta+1))*gamma(alpha+1)*gamma(beta+1)/
  gamma(alpha+beta+1)
  P_0 = ones(x)/sqrt(gamma0)
  if (N == 0)
    size(P_0,1) > size(P_0,2) ? (return P_0) : (return P_0.')
  end
  gamma1 = (alpha+1)*(beta+1)*gamma0/(alpha+beta+3)
  P_1 = 0.5*((alpha+beta+2).*x + (alpha-beta))/sqrt(gamma1)
  if (N == 1)
    size(P_1,1) > size(P_1,2) ? (return P_1) : (return P_1.')
  end
  # Henceforth, P_0 denotes P_{i} and P_1 denotes P_{i+1}
  # repeat value in recurrence
  aold = (2./(2+alpha+beta))*sqrt((alpha+1)*(beta+1)/(alpha+beta+3))
  save = zeros(x)
  for i = 1:N-1
    h1 = 2*i + alpha + beta
    anew = (2./(h1+2))*sqrt((i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/
                            ((h1+1)*(h1+3)))
    bnew = -(alpha^2 - beta^2)/(h1*(h1+2))
    save = P_1
    P_1 = (1/anew).*(-aold.*P_0 + (x-bnew).*P_1)
    P_0 = save
    aold = anew
  end
  size(P_1,1) > size(P_1,2) ? (return P_1) : (return P_1.')
end

@doc """
### OrthoPoly.diffjacobipoly{T}

Evaluate the first derivative of a Jacobi Polynomial at some points.

**Inputs**

* `x`: points at which to evaluate polynomial derivative
* `alpha`,`beta`: define the type of Jacobi Polynomial (alpha + beta != 1)
* `N`: polynomial degree

**Outputs**

* dP - derivative of polynomial evaluated at x

"""->
function diffjacobipoly{T}(x::Array{T}, alpha::FloatingPoint, 
                           beta::FloatingPoint, N::Int)
  @assert( alpha + beta != -1 )
  @assert( alpha > -1 && beta > -1)
  DP_0 = zeros(x)
  if (N == 0)
    size(DP_0,1) > size(DP_0,2) ? (return DP_0) : (return DP_0.')
  end
  gamma0 = ((2^(alpha+beta+1))/(alpha+beta+1))*gamma(alpha+1)*gamma(beta+1)/
  gamma(alpha+beta+1)
  gamma1 = (alpha+1)*(beta+1)*gamma0/(alpha+beta+3)
  DP_1 = ones(x).*0.5*(alpha+beta+2)/sqrt(gamma1)
  if (N == 1)
    size(DP_1,1) > size(DP_1,2) ? (return DP_1) : (return DP_1.')
  end
  # initialize values P_0(x) and P_1(x) for recurrence
  P_0 = ones(x)./sqrt(gamma0)
  P_1 = 0.5*((alpha+beta+2).*x + (alpha-beta))/sqrt(gamma1)
  # repeat value in recurrence
  aold = (2./(2+alpha+beta))*sqrt((alpha+1)*(beta+1)/(alpha+beta+3))
  save = zeros(x)
  for i = 1:N-1
    h1 = 2*i + alpha + beta
    anew = (2./(h1+2))*sqrt((i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/
                            ((h1+1)*(h1+3)))
    bnew = -(alpha^2 - beta^2)/(h1*(h1+2))
    save = DP_1
    DP_1 = (1/anew).*(-aold.*DP_0 + P_1 + (x-bnew).*DP_1)
    DP_0 = save
    save = P_1
    P_1 = (1/anew).*(-aold.*P_0 + (x-bnew).*P_1)
    P_0 = save
    aold = anew
  end
  size(DP_1,1) > size(DP_1,2) ? (return DP_1) : (return DP_1.')
end

@doc """
### OrthoPoly.proriolpoly{T}: method for right triangle

Evaluate Proriol orthogonal polynomial basis function on the right triangle.

**Inputs**

* `x`,`y`: locations at which to evaluate the polynomial
* `i`,`j`: index pair that defines the basis function to evaluate; see Hesthaven
  and Warburton's Nodal DG book, for example, for a reference.

**Outputs**

* `P`: basis function at (`x`,`y`)

"""->
function proriolpoly{T}(x::Array{T}, y::Array{T}, i::Int, j::Int)
  @assert( i >= 0 && j >= 0 )
  xi = zeros(T, size(x))
  for k = 1:length(x)
    y[k] != 1.0 ? xi[k] = 2.0*(1 + x[k])./(1 - y[k]) -1 : xi[k] = -1
  end
  P = sqrt(2).*jacobipoly(xi, 0.0, 0.0, i
                          ).*jacobipoly(y, convert(FloatingPoint, 2*i+1), 0.0, j
                                        ).*((1-y).^(i))
  return P
end

@doc """
### OrthoPoly.proriolpoly{T}: method for a right tetrahedron

Evaluate Proriol orthogonal polynomial basis function on the right tetrahedron.

**Inputs**

* `x`,`y`,`z`: locations at which to evaluate the polynomial
* `i`,`j`,`k`: index triple that defines the basis function to evaluate; see Hesthaven
  and Warburton's Nodal DG book, for example, for a reference.

**Outputs**

* `P`: basis function at (`x`,`y`,`z`)

"""->
function proriolpoly{T}(x::Array{T}, y::Array{T}, z::Array{T}, i::Int, j::Int,
                        k::Int)
  @assert( i >= 0 && j >= 0 && k >= 0 )
  xi = zeros(T, size(x))
  eta = zeros(T, size(y))
  for m = 1:length(x)
    y[m]+z[m] != 0.0 ? xi[m] = -2.0*(1+x[m])./(y[m]+z[m]) - 1 : xi[m] = -1
    z[m] != 1.0 ? eta[m] = 2.0*(1+y[m])./(1-z[m]) - 1 : eta[m] = -1
  end
  P = sqrt(8).*jacobipoly(xi, 0.0, 0.0, i).*
  jacobipoly(eta, convert(FloatingPoint, 2*i+1), 0.0, j).*((1-eta).^(i)).*
  jacobipoly(z, convert(FloatingPoint, 2*i+2*j+2), 0.0, k).*((1-z).^(i+j))
  return P
end

@doc """
### OrthoPoly.diffproriolpoly

Evaluate the derivatives of a Proriol orthogonal polynomial basis function on
the right triangle.

**Inputs**

* `x`,`y`: locations at which to evaluate the derivative
* `i`,`j`: index pair that defines the basis function to differentiate; see
  Hesthaven and Warburton's Nodal DG book, for example, for a reference.

**Outputs**

* `dPdx`,`dPdy`: derivative of basis function at (`x`,`y`)

"""->
function diffproriolpoly{T}(x::Array{T}, y::Array{T}, i::Int, j::Int)
  xi = zeros(T, size(x))
  for k = 1:length(x)
    real(y[k]) != 1.0 ? xi[k] = 2.0*(1 + x[k])./(1 - y[k]) -1 : xi[k] = -1
  end
  dPdx = zeros(T, size(x))
  dPdy = zeros(T, size(y))
  if (i == 0) && (j == 0)
    return dPdx, dPdy
  end
  # compute some terms for reuse
  Jxi = jacobipoly(xi, 0.0, 0.0, i)
  Jy = jacobipoly(y, convert(FloatingPoint, 2*i+1), 0.0, j)
  dJxidxi = diffjacobipoly(xi, 0.0, 0.0, i)
  dJydy = diffjacobipoly(y, convert(FloatingPoint, 2*i+1), 0.0, j)

  if (i > 0)
    # dPdx is only nonzero if i > 0;
    dPdx += sqrt(2).*dJxidxi.*Jy.*2.0.*((1-y).^(i-1))
    dPdy = -Jy.*i.*((1-y).^(i-1))
  end
  # dPdx is now finished, but dPdy needs additional terms
  dPdy += dJydy.*((1-y).^(i))
  dPdy .*= Jxi
  if (i >= 1)
    for k = 1:length(x)
      real(y[k]) != 1.0 ? dPdy[k] += dJxidxi[k]*Jy[k]*2.0*(1+x[k])*((1-y[k])^(i-2)) : 
      dPdy[k] += 0.0
    end
  end
  dPdy *= sqrt(2)
  
  return dPdx, dPdy
end

@doc """
### OrthoPoly.diffproriolpoly{T}

Evaluate the derivatives of a Proriol orthogonal polynomial basis function on
the right tetrahedron.

*Notes*: the derivatives are computed using the complex-step method (since there
 are many outputs and only 3 inputs); therefore, a different method should be
 used for verification of this method.

**Inputs**

* `x`,`y`,`z`: locations at which to evaluate the derivative
* `i`,`j`,`k`: index triple that defines the basis function to differentiate;
  see Hesthaven and Warburton's Nodal DG book, for example, for a reference.

**Outputs**

* `dPdx`,`dPdy`,`dPdz`: derivatives of basis function at (`x`,`y`,`z`)

"""->
function diffproriolpoly{T}(x::Array{T}, y::Array{T}, z::Array{T}, i::Int, 
                            j::Int, k::Int)
  # each node is independent, so use complex step once for each coordinate. Care
  # is needed at the one vertex, where the xi and eta mappings become singular.
  # To avoid problems, directional derivatives are used.
  eps_step = 1e-60
  xc = complex(x,0)
  yc = complex(y,0)
  zc = complex(z,0)
  # compute derivative with respect to z
  zc -= eps_step*im
  Pc = proriolpoly(xc, yc, zc, i, j, k)
  dPdz = -imag(Pc)./eps_step
  # compute dPdy = -(Grad P) dot (0,-1,-1) - dPdz
  yc -= eps_step*im
  Pc = proriolpoly(xc, yc, zc, i, j, k)
  dPdy = -dPdz - imag(Pc)./eps_step
  # compute dPdx = -(Grad P) dot (-1,-1,-1) - dPdz - dPdy
  xc -= eps_step*im
  Pc = proriolpoly(xc, yc, zc, i, j, k)
  dPdx = -dPdz - dPdy - imag(Pc)./eps_step
  return dPdx, dPdy, dPdz
end

end