module OrthoPoly
# utilities for working with orthogonal polynomials

using SpecialFunctions
using LinearAlgebra
# using ..getComplexStep

# getComplexstep(::Type{T}) where {T<:Float32} = 1f-20
# getComplexstep(::Type{T}) where {T<:Float64} = 1e-60
# getComplexstep(::Type{T}) where {T<:Complex{Float64}} = 1e-60
# function getComplexStep(::Type{T}) where {T<:Float32}
#   return 1f-20
# end

# function getComplexStep{T <: Float64}(::Type{T})
#   return 1e-60
# end

# function getComplexStep{T <: Complex}(::Type{T})
#   return 1e-60
# end

"""
### OrthoPoly.lglnodes

Computes the Legendre-Gauss-Lobatto (LGL) quadrature nodes and weights on the
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
 
"""
function lglnodes(N, T=Float64)
  N1 = N+1
  # Use the Chebyshev-Gauss-Lobatto nodes as an initial guess
  x = -cos.(π*[0:N;]/N)
  # The Legendre Vandermonde Matrix 
  P = zeros(T, (N1,N1))
  # Compute P_(N) using the recursion relation; compute its first and second
  # derivatives and update x using the Newton-Raphson method.
  xold = (T)(2)
  iter = 1; maxiter = 100  
  while maximum(abs, real(x.-xold)) > eps(real(Float64((T)(1)))) && iter < maxiter
    iter += 1
    xold = x
    P[:,1] .= one(T)
    P[:,2] = x
    for k=2:N
      P[:,k+1]= ((2k-1)*x.*P[:,k]-(k-1)*P[:,k-1])/k
    end
    x = xold - ( x.*P[:,N1]-P[:,N] )./( N1*P[:,N1] )           
  end
  w = 2.0 ./(N*N1*P[:,N1].^2)
  return x, w
end

"""
### OrthoPoly.lgnodes

Computes the Legendre-Gauss (LG) quadrature nodes and weights on the
interval [-1,1].  The LG nodes are the zeros of P_N(x), where P_N(x)
denotes the Nth Legendre polynomial.

**Inputs**

* `N`: number of nodes
* `T`: number type

**Outputs**

* `x`: the LG nodes
* `w`: the LG weights

Julia version adapted from Matlab code written by Greg von Winckel - 02/25/2004
Contact: gregvw@chtm.unm.edu
 
"""
function lgnodes(N, T=Float64)
  Nm1 = N-1; Np1 = N+1
  N == 1 ? xu = [0.0] : xu = range(-1,stop=1,length=N)
  # initial guess
  x = -cos.((2*[0:Nm1;] .+ 1)*pi/(2*Nm1+2)) .- (0.27/N)*sin.(pi*xu*Nm1/Np1)
  # Legendre-Gauss Vandermonde Matrix and its derivative
  L = zeros(N,Np1)
  Lp = zeros(N)
  # compute the zeros of the Legendre Polynomial using the recursion relation
  # and Newton's method; loop until new points are uniformly within epsilon of
  # old points
  xold = (T)(2)
  iter = 1; maxiter = 30
  while maximum(abs, real(x .- xold)) > 0.1*eps(real(Float64((T)(1)))) && iter < maxiter
    iter += 1
    L[:,1] .= one(T)
    L[:,2] .= x
    for k = 2:N
      L[:,k+1] = ((2*k-1)*x.*L[:,k]-(k-1)*L[:,k-1])/k
    end
    Lp[:] .= (Np1)*(L[:,N] .- x.*L[:,Np1])./(1 .- x.^2)
    xold = x
    x -= L[:,Np1]./Lp
  end
  w = 2.0 ./((1 .- x.^2).*Lp.^2)*(Np1/N)^2
  return x, w
end

"""
### OrthoPoly.jacobipoly{T}

Evaluate a Jacobi polynomial at some points.  Based on JacobiP in Hesthaven and
Warburton's nodal DG book.

**Inputs**

* `x`: points at which to evaluate polynomial
* `alpha`,`beta`: define the type of Jacobi Polynomial (alpha + beta != 1)
* `N`: polynomial degree

**Outputs**

* `P`: the polynomial evaluated at x

"""
function jacobipoly(x::Array{T}, alpha::AbstractFloat, beta::AbstractFloat,
                       N::Int) where {T<:Number}
  @assert( alpha + beta != -1 )
  @assert( alpha > -1 && beta > -1)
  # Initial values P_0(x) and P_1(x)
  gamma0 = ((2^(alpha+beta+1))/(alpha+beta+1))*gamma(alpha+1)*gamma(beta+1)/
  gamma(alpha+beta+1)
  P_0 = ones(size(x))/sqrt(gamma0)
  if (N == 0)
    size(P_0,1) > size(P_0,2) ? (return P_0) : (return P_0')
  end
  gamma1 = (alpha+1)*(beta+1)*gamma0/(alpha+beta+3)
  P_1 = 0.5*((alpha+beta+2).*x .+ (alpha-beta))/sqrt(gamma1)
  if (N == 1)
    size(P_1,1) > size(P_1,2) ? (return P_1) : (return P_1')
  end
  # Henceforth, P_0 denotes P_{i} and P_1 denotes P_{i+1}
  # repeat value in recurrence
  aold = (2.0/(2+alpha+beta))*sqrt((alpha+1)*(beta+1)/(alpha+beta+3))
  save = zeros(size(x))
  for i = 1:N-1
    h1 = 2*i + alpha + beta
    anew = (2.0/(h1+2))*sqrt((i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/
                            ((h1+1)*(h1+3)))
    bnew = -(alpha^2 - beta^2)/(h1*(h1+2))
    save = P_1
    P_1 = (1/anew).*(-aold.*P_0 + (x .- bnew).*P_1)
    P_0 = save
    aold = anew
  end
  size(P_1,1) > size(P_1,2) ? (return P_1) : (return P_1')
end

"""
### OrthoPoly.diffjacobipoly{T}

Evaluate the first derivative of a Jacobi Polynomial at some points.

**Inputs**

* `x`: points at which to evaluate polynomial derivative
* `alpha`,`beta`: define the type of Jacobi Polynomial (alpha + beta != 1)
* `N`: polynomial degree

**Outputs**

* dP - derivative of polynomial evaluated at x

"""
function diffjacobipoly(x::Array{T}, alpha::AbstractFloat, 
                           beta::AbstractFloat, N::Int) where {T<:Number}
  @assert( alpha + beta != -1 )
  @assert( alpha > -1 && beta > -1)
  DP_0 = zeros(size(x))
  if (N == 0)
    size(DP_0,1) > size(DP_0,2) ? (return DP_0) : (return DP_0')
  end
  gamma0 = ((2^(alpha+beta+1))/(alpha+beta+1))*gamma(alpha+1)*gamma(beta+1)/
  gamma(alpha+beta+1)
  gamma1 = (alpha+1)*(beta+1)*gamma0/(alpha+beta+3)
  DP_1 = ones(size(x)).*0.5*(alpha+beta+2)/sqrt(gamma1)
  if (N == 1)
    size(DP_1,1) > size(DP_1,2) ? (return DP_1) : (return DP_1')
  end
  # initialize values P_0(x) and P_1(x) for recurrence
  P_0 = ones(size(x))./sqrt(gamma0)
  P_1 = 0.5*((alpha+beta+2).*x .+ (alpha-beta))/sqrt(gamma1)
  # repeat value in recurrence
  aold = (2.0/(2+alpha+beta))*sqrt((alpha+1)*(beta+1)/(alpha+beta+3))
  save = zeros(size(x))
  for i = 1:N-1
    h1 = 2*i + alpha + beta
    anew = (2.0/(h1+2))*sqrt((i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/
                            ((h1+1)*(h1+3)))
    bnew = -(alpha^2 - beta^2)/(h1*(h1+2))
    save = DP_1
    DP_1 = (1/anew).*(-aold.*DP_0 .+ P_1 .+ (x .- bnew).*DP_1)
    DP_0 = save
    save = P_1
    P_1 = (1/anew).*(-aold.*P_0 + (x .- bnew).*P_1)
    P_0 = save
    aold = anew
  end
  size(DP_1,1) > size(DP_1,2) ? (return DP_1) : (return DP_1')
end

"""
### OrthoPoly.proriolpoly{T}: method for right triangle

Evaluate Proriol orthogonal polynomial basis function on the right triangle.

**Inputs**

* `x`,`y`: locations at which to evaluate the polynomial
* `i`,`j`: index pair that defines the basis function to evaluate; see Hesthaven
  and Warburton's Nodal DG book, for example, for a reference.

**Outputs**

* `P`: basis function at (`x`,`y`)

"""
function proriolpoly(x::Array{T}, y::Array{T}, i::Int, j::Int) where {T<:Number}
  @assert( i >= 0 && j >= 0 ) 
  xi = zeros(T, size(x))
  for k = 1:length(x)
    y[k] != 1.0 ? xi[k] = 2.0*(1 + x[k])./(1 - y[k]) -1 : xi[k] = -1
  end
  P = sqrt(2).*jacobipoly(xi, 0.0, 0.0, i).*jacobipoly(y, convert(AbstractFloat, 2*i+1), 0.0, j).*((1 .- y).^(i))
  return P
end

"""
### OrthoPoly.proriolpoly{T}: method for a right tetrahedron

Evaluate Proriol orthogonal polynomial basis function on the right tetrahedron.

**Inputs**

* `x`,`y`,`z`: locations at which to evaluate the polynomial
* `i`,`j`,`k`: index triple that defines the basis function to evaluate; see Hesthaven
  and Warburton's Nodal DG book, for example, for a reference.

**Outputs**

* `P`: basis function at (`x`,`y`,`z`)

"""
function proriolpoly(x::Array{T}, y::Array{T}, z::Array{T}, i::Int, j::Int,
                        k::Int) where {T<:Number}
  @assert( i >= 0 && j >= 0 && k >= 0 )
  xi = zeros(T, size(x))
  eta = zeros(T, size(y))
  for m = 1:length(x)
    y[m]+z[m] != 0.0 ? xi[m] = -2.0*(1+x[m])./(y[m]+z[m]) - 1 : xi[m] = -1
    z[m] != 1.0 ? eta[m] = 2.0*(1+y[m])./(1-z[m]) - 1 : eta[m] = -1
  end
  P = sqrt(8).*jacobipoly(xi, 0.0, 0.0, i).*
    jacobipoly(eta, convert(AbstractFloat, 2*i+1), 0.0, j).*((1 .- eta).^(i)).*
    jacobipoly(z, convert(AbstractFloat, 2*i+2*j+2), 0.0, k).*((1 .- z).^(i+j))
  return P
end

"""
### OrthoPoly.diffproriolpoly

Evaluate the derivatives of a Proriol orthogonal polynomial basis function on
the right triangle.

**Inputs**

* `x`,`y`: locations at which to evaluate the derivative
* `i`,`j`: index pair that defines the basis function to differentiate; see
  Hesthaven and Warburton's Nodal DG book, for example, for a reference.

**Outputs**

* `dPdx`,`dPdy`: derivative of basis function at (`x`,`y`)

"""
function diffproriolpoly(x::Array{T}, y::Array{T}, i::Int, j::Int) where {T<:Number}
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
  Jy = jacobipoly(y, convert(AbstractFloat, 2*i+1), 0.0, j)
  dJxidxi = diffjacobipoly(xi, 0.0, 0.0, i)
  dJydy = diffjacobipoly(y, convert(AbstractFloat, 2*i+1), 0.0, j)

  if (i > 0)
    # dPdx is only nonzero if i > 0;
    dPdx += sqrt(2).*dJxidxi.*Jy.*2.0.*((1 .- y).^(i-1))
    dPdy = -Jy.*i.*((1 .- y).^(i-1))
  end
  # dPdx is now finished, but dPdy needs additional terms
  dPdy += dJydy.*((1 .- y).^(i))
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

"""
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

"""
function diffproriolpoly(x::Array{T}, y::Array{T}, z::Array{T}, i::Int, j::Int, k::Int) where {T<:Number}
  # each node is independent, so use complex step once for each coordinate. Care
  # is needed at the one vertex, where the xi and eta mappings become singular.
  # To avoid problems, directional derivatives are used.
  eps_step = 1e-60+0im #getComplexStep(T)
  xc = complex.(x,0)
  yc = complex.(y,0)
  zc = complex.(z,0)
  # compute derivative with respect to z
  zc .-= eps_step*im
  Pc = proriolpoly(xc, yc, zc, i, j, k)
  dPdz = -imag(Pc)./eps_step
  # compute dPdy = -(Grad P) dot (0,-1,-1) - dPdz
  yc .-= eps_step*im
  Pc = proriolpoly(xc, yc, zc, i, j, k)
  dPdy = -dPdz - imag(Pc)./eps_step
  # compute dPdx = -(Grad P) dot (-1,-1,-1) - dPdz - dPdy
  xc .-= eps_step*im
  Pc = proriolpoly(xc, yc, zc, i, j, k)
  dPdx = -dPdz - dPdy - imag(Pc)./eps_step
  return dPdx, dPdy, dPdz
end

"""
### OrthoPoly.vandermonde{T}

Evaluate the Vandermonde matrix using the Proriol polynomials on the right triangle.

**Inputs**

* `p`: the maximum total degree of the polynomial
* `x`,`y`: locations at which to evaluate the derivative
* `compute_grad`: indicates whether to compute the gradient of the Vandermonde matrix

**Outputs**

* `V`: the Vandermonde matrix 
* `Vdx`,`Vdy`: derivatives of the Vandermonde matrix at (`x`,`y`)

"""
function vandermonde(p::Int, x::Array{T}, y::Array{T}; compute_grad::Bool=true) where {T}
  nnodes = length(x)
  num_eq = convert(Int, (p+1)*(p+2)/2)

  # loop over orthogonal polynomials of degree r <= p 
  V = zeros(T, (nnodes, num_eq))
  Vdx = zeros(T, (nnodes, num_eq))
  Vdy = zeros(T, (nnodes, num_eq))
  ptr = 1
  for r = 0:p
    for j = 0:r
      i = r-j
      P = OrthoPoly.proriolpoly(x, y, i, j)
      V[:,ptr] = P
      if compute_grad
        dPdx, dPdy = OrthoPoly.diffproriolpoly(x, y, i, j)
        Vdx[:,ptr] = dPdx
        Vdy[:,ptr] = dPdy
      end
      ptr += 1
    end
  end
  return V, Vdx, Vdy
end

"""
### OrthoPoly.vandermonde{T}

Evaluate the Vandermonde matrix using the Proriol polynomials on the right tetrahedron.

**Inputs**

* `p`: the maximum total degree of the polynomial
* `x`,`y`,`z`: locations at which to evaluate the derivative
* `compute_grad`: indicates whether to compute the gradient of the Vandermonde matrix

**Outputs**

* `V`: the Vandermonde matrix 
* `Vdx`,`Vdy`,`Vdz`: derivatives of the Vandermonde matrix at (`x`,`y`,`z`)

"""
function vandermonde(p::Int, x::Array{T}, y::Array{T}, z::Array{T}; compute_grad::Bool=true) where {T}
  nnodes = length(x)
  num_eq = convert(Int, (p+1)*(p+2)*(p+3)/6)

  # loop over orthogonal polynomials of degree r <= p 
  V = zeros(T, (nnodes, num_eq))
  Vdx = zeros(T, (nnodes, num_eq))
  Vdy = zeros(T, (nnodes, num_eq))
  Vdz = zeros(T, (nnodes, num_eq))
  ptr = 1
  for r = 0:p
    for k = 0:r
      for j = 0:r-k
        i = r-j-k
        P = OrthoPoly.proriolpoly(x, y, z, i, j, k)
        V[:,ptr] = P
        if compute_grad
          dPdx, dPdy, dPdz = OrthoPoly.diffproriolpoly(x, y, z, i, j, k)
          Vdx[:,ptr] = dPdx
          Vdy[:,ptr] = dPdy
          Vdz[:,ptr] = dPdz
        end
        ptr += 1
      end
    end
  end
  return V, Vdx, Vdy, Vdz
end

"""
### OrthoPoly.vandermonde_monomial{T}

Evaluate the Vandermonde matrix using monomials on the line.

**Inputs**

* `p`: the maximum total degree of the polynomial
* `x`: locations at which to evaluate the derivative
* `compute_grad`: indicates whether to compute the gradient of the Vandermonde matrix
* `compute_integ`: indicates whether to compute the integral of the monomial basis functions

**Outputs**

* `V`: the Vandermonde matrix 
* `Vdx`: derivatives of the Vandermonde matrix at (`x`)
* `Vinteg`: the integral of each basis function in the Vandermonde matrix

"""
function vandermonde_monomial(p::Int, x::Array{T}; compute_grad::Bool=true, compute_integ::Bool=true) where {T}
  nnodes = length(x)
  num_eq = p+1
  V = zeros(T, (nnodes, num_eq))
  Vdx = zeros(T, (nnodes, num_eq))
  Vinteg = zeros(Float64, num_eq)

  for indx=1:num_eq
    a = indx-1;
    for n=1:nnodes 
      V[n, indx] = x[n]^a
      if compute_grad
        if a == 0
          Vdx[n, indx] = 0
        else 
          Vdx[n, indx] = a * x[n]^(a-1) 
        end
      end
    end

    if compute_integ # assumes simplex on the interval [0,1]
      Vinteg[indx] = factorial(big(convert(Int,a)))/factorial(big(convert(Int,1.0+a)))
    end
  end

  return V, Vdx, Vinteg
end

"""
### OrthoPoly.vandermonde_monomial{T}

Evaluate the Vandermonde matrix using monomials on the right triangle.

**Inputs**

* `p`: the maximum total degree of the polynomial
* `x`,`y`: locations at which to evaluate the derivative
* `compute_grad`: indicates whether to compute the gradient of the Vandermonde matrix
* `compute_integ`: indicates whether to compute the integral of the monomial basis functions

**Outputs**

* `V`: the Vandermonde matrix 
* `Vdx`,`Vdy`: derivatives of the Vandermonde matrix at (`x`,`y`)
* `Vinteg`: the integral of each basis function in the Vandermonde matrix

"""
function vandermonde_monomial(p::Int, x::Array{T}, y::Array{T}; compute_grad::Bool=true, 
  compute_integ::Bool=true) where {T}
  nnodes = length(x)
  num_eq = convert(Int, (p+1)*(p+2)/2)
  V = zeros(T, (nnodes, num_eq))
  Vdx = zeros(T, (nnodes, num_eq))
  Vdy = zeros(T, (nnodes, num_eq))
  Vinteg = zeros(Float64, num_eq)

  for i=0:p
    for j=0:i
      indx = convert(Int,(i*(i+1))/2) + j + 1
      a  = i-j
      b  = j

      for n=1:nnodes 
        V[n, indx] = x[n]^a * y[n]^b
        if compute_grad
          if a == 0
            Vdx[n, indx] = 0
          else 
            Vdx[n, indx] = a * x[n]^(a-1) * y[n]^b
          end

          if b == 0
            Vdy[n, indx]= 0
          else
            Vdy[n, indx] = b * x[n]^a * y[n]^(b-1) 
          end
        end
      end

      if compute_integ # assumes simplex on the interval [0,1]
        Vinteg[indx] = factorial(big(convert(Int,a)))*factorial(big(convert(Int,b)))/factorial(big(convert(Int,2.0+a+b)))
      end

    end
  end

  return V, Vdx, Vdy, Vinteg
end

"""
### OrthoPoly.vandermonde_monomial{T}

Evaluate the Vandermonde matrix using monomials on the right tetrahedron.

**Inputs**

* `p`: the maximum total degree of the polynomial
* `x`,`y`,`z`: locations at which to evaluate the derivative
* `compute_grad`: indicates whether to compute the gradient of the Vandermonde matrix
* `compute_integ`: indicates whether to compute the integral of the monomial basis functions

**Outputs**

* `V`: the Vandermonde matrix 
* `Vdx`,`Vdy`, `Vdz`: derivatives of the Vandermonde matrix at (`x`,`y`,`z`)
* `Vinteg`: the integral of each basis function in the Vandermonde matrix

"""
function vandermonde_monomial(p::Int, x::Array{T}, y::Array{T}, z::Array{T}; 
  compute_grad::Bool=true, compute_integ::Bool=true) where {T}
  nnodes = length(x)
  num_eq = convert(Int, (p+1)*(p+2)*(p+3)/6)
  V = zeros(T, (nnodes, num_eq))
  Vdx = zeros(T, (nnodes, num_eq))
  Vdy = zeros(T, (nnodes, num_eq))
  Vdz = zeros(T, (nnodes, num_eq))
  Vinteg = zeros(Float64, num_eq)

  for i=0:p
    for j=0:i
      for k=0:j
        indx = convert(Int,(j*(j+1))/2 + (i*(i+1)*(i+2))/6)+ k +1;
        a = i-j;
        b = j-k;
        c = k;
        for n=1:nnodes
          V[n,indx] = x[n]^a * y[n]^b * z[n]^c;

          if compute_grad
            if a == 0
                Vdx[n, indx] = 0.0;
            else 
                Vdx[n, indx] = a * x[n]^(a-1) * y[n]^b * z[n]^c;
            end

            if b == 0
                Vdy[n, indx]= 0.0;
            else
                Vdy[n, indx] = b * x[n]^a * y[n]^(b-1) * z[n]^c; 
            end

            if c == 0
                Vdz[n, indx]= 0.0;
            else
                Vdz[n, indx] = c * x[n]^a * y[n]^(b) * z[n]^(c-1); 
            end
          end
        end
        if compute_integ # assumes simplex on the interval [0,1]
          Vinteg[indx] = factorial(big(convert(Int,a)))*factorial(big(convert(Int,b)))*factorial(big(convert(Int,c)))/factorial(big(convert(Int,3.0+a+b+c)))
        end
      end
    end
  end

  return V, Vdx, Vdy, Vdz, Vinteg
end

"""
### OrthoPoly.vandermonde_arnoldi{T}

Evaluate the Vandermonde with Arnoldi on the line.

**Inputs**

* `p`: the maximum total degree of the polynomial
* `x`: locations at which to evaluate the derivative
* `compute_grad`: indicates whether to compute the gradient of the Vandermonde matrix

**Outputs**

* `V`: the Vandermonde matrix 
* `Vdx`: derivatives of the Vandermonde matrix at (`x`)
* `Hes`: the Hessenberg matrix 

"""
function vandermonde_arnoldi(p::Int, x::Array{T}; compute_grad::Bool=true) where {T}
  nnodes = length(x)
  num_eq = p+1
  V = zeros(T, (nnodes, num_eq))
  V[:,1] .= 1.0
  Vdx = zeros(T, (nnodes, num_eq))
  Hes = zeros(num_eq,num_eq-1)

  for k = 1:num_eq-1
    v = x .* V[:,k]
    for j = 1:k
      Hes[j,k] = (Matrix(V[:,j]') * v)[1]/nnodes
      v = v - Hes[j,k]*V[:,j]
    end
    Hes[k+1,k] = norm(v)/sqrt(nnodes)
    V[:,k+1] = v/Hes[k+1,k]
    if compute_grad
      Vdx[:,k+1] = (x.*Vdx[:,k] - Vdx[:,1:k]*Hes[1:k,k] + V[:,k])/Hes[k+1,k]
    end
  end

  return V,Vdx,Hes
end

"""
### OrthoPoly.vandermonde_arnoldi{T}

Evaluate the Vandermonde with Arnoldi on the right triangle.

**Inputs**

* `p`: the maximum total degree of the polynomial
* `x`,`y`: locations at which to evaluate the derivative
* `compute_grad`: indicates whether to compute the gradient of the Vandermonde matrix

**Outputs**

* `V`: the Vandermonde matrix 
* `Vdx`, `Vdy`: derivatives of the Vandermonde matrix at (`x`,`y`)
* `Hes`: the Hessenberg matrix 

"""
function vandermonde_arnoldi(p::Int, x::Array{T}, y::Array{T}; compute_grad::Bool=true) where {T}
  nnodes = length(x)
  num_eq = convert(Int, (p+1)*(p+2)/2)
  V = zeros(T, (nnodes, num_eq))
  V[:,1] .= 1.0
  Vdx = zeros(T, (nnodes, num_eq))
  Vdy = zeros(T, (nnodes, num_eq))
  Hes = zeros(num_eq,num_eq)
  Hes[1,1] = 1.0

  # get list of indices for multivariate polynomials
  ab = zeros(2,num_eq)
  for i=0:p
    for j=0:i
        indx = convert(Int,i*(i+1)/2) + j + 1
        a = i-j
        b = j
        ab[:,indx] = [a;b]
    end 
  end

  # find the smallest k such that there indx(k)+ e_k = indx(k+1)
  ks = []
  coords = []
  max_indx = minimum([num_eq,nnodes])-1
  for indx=1:max_indx #size(ab,2)-1
    for k = 1:size(ab,2)
        if ab[:,indx+1]-ab[:,k]==[1;0]
          push!(coords, 1)
          push!(ks, k) 
          break;
        elseif ab[:,indx+1]-ab[:,k]==[0;1]
          push!(coords, 2)
          push!(ks, k) 
          break;
        end
    end
  end

  xy = [x y]
  dexi = zeros(nnodes)
  deta = zeros(nnodes)
  for indx=1:max_indx #size(ab,2)-1
    exi = xy[:,coords[indx]]
    eta = copy(exi)

    if coords[indx]==1
        dexi = ones(nnodes)
        deta = zeros(nnodes)
    elseif coords[indx]==2
        dexi = zeros(nnodes)
        deta = ones(nnodes)
    end
    
    v = diagm(exi)*V[:,ks[indx]]
    for t = 1:2
        s = 1/nnodes .* V[:,1:indx]'*v
        v = v - V[:,1:indx]*s
        Hes[1:indx,indx+1] =  Hes[1:indx,indx+1] + s
    end

    Hes[indx+1,indx+1] = 1/sqrt(nnodes).*norm(v,2)
    V[:,indx+1] = v ./Hes[indx+1,indx+1];

    if compute_grad
      s1 = 1/nnodes * Matrix(V[:,1:indx]')* diagm(exi)*V[:,ks[indx]]
      s2 =  s1 - 1/nnodes *Matrix(V[:,1:indx]')*V[:,1:indx]*s1
      Vdx[:,indx+1] = (diagm(dexi)*V[:,ks[indx]]
                      + diagm(exi)*Vdx[:,ks[indx]]- Vdx[:,1:indx]*s1
                      - Vdx[:,1:indx]*s2) ./Hes[indx+1,indx+1]
      Vdy[:,indx+1] = (diagm(deta)*V[:,ks[indx]] 
                      + diagm(eta)*Vdy[:,ks[indx]]- Vdy[:,1:indx]*s1
                      - Vdy[:,1:indx]*s2) ./Hes[indx+1,indx+1]
    end
    
    # #another way to compute Q(:,indx+1) (more important form to obtain derivatives)
    # V[:,indx+1] = (diagm(eta)*V[:,ks[indx]] - V[:,1:indx]*s1 - V[:,1:indx]*s2)./Hes[indx+1,indx+1]
  end
  for indx=max_indx+1:num_eq
    a = ab[1,indx]
    b = ab[2,indx]

    # P = OrthoPoly.proriolpoly(x, y, convert(Int,a), convert(Int,b))
    # V[:,indx] = P
    # if compute_grad
    #   dPdx, dPdy = OrthoPoly.diffproriolpoly(x, y, convert(Int,a), convert(Int,b))
    #   Vdx[:,indx] = dPdx
    #   Vdy[:,indx] = dPdy
    # end
    # indx += 1

    for n=1:nnodes 
      V[n, indx] = x[n]^a * y[n]^b
      if compute_grad
        if a == 0
          Vdx[n, indx] = 0
        else 
          Vdx[n, indx] = a * x[n]^(a-1) * y[n]^b
        end

        if b == 0
          Vdy[n, indx]= 0
        else
          Vdy[n, indx] = b * x[n]^a * y[n]^(b-1) 
        end
      end
    end
  end

  return V,Vdx,Vdy,Hes
end

"""
### OrthoPoly.vandermonde_arnoldi{T}

Evaluate the Vandermonde with Arnoldi on the right tetrahedron.

**Inputs**

* `p`: the maximum total degree of the polynomial
* `x`,`y`,`z`: locations at which to evaluate the derivative
* `compute_grad`: indicates whether to compute the gradient of the Vandermonde matrix

**Outputs**

* `V`: the Vandermonde matrix 
* `Vdx`,`Vdy`,`Vdz`: derivatives of the Vandermonde matrix at (`x`,`y`,`z`)
* `Hes`: the Hessenberg matrix 

"""
function vandermonde_arnoldi(p::Int, x::Array{T}, y::Array{T}, z::Array{T}; compute_grad::Bool=true) where {T}
  nnodes = length(x)
  num_eq = convert(Int, (p+1)*(p+2)*(p+3)/6)
  V = zeros(T, (nnodes, num_eq))
  V[:,1] .= 1.0
  Vdx = zeros(T, (nnodes, num_eq))
  Vdy = zeros(T, (nnodes, num_eq))
  Vdz = zeros(T, (nnodes, num_eq))
  Hes = zeros(num_eq,num_eq)
  Hes[1,1] = 1.0

  # get list of indices for multivariate polynomials
  abc = zeros(3,num_eq)
  for i=0:p
    for j=0:i
      for k=0:j
        indx = convert(Int,(j*(j+1))/2 + (i*(i+1)*(i+2))/6) + k + 1
        a = i-j
        b = j-k
        c = k
        abc[:,indx] = [a;b;c]
      end
    end 
  end

  # find the smallest k such that there indx(k)+ e_k = indx(k+1)
  ks = []
  coords = []
  for indx=1:size(abc,2)-1
    for k = 1:size(abc,2)
        if abc[:,indx+1]-abc[:,k]==[1;0;0]
          push!(coords, 1)
          push!(ks, k) 
          break;
        elseif abc[:,indx+1]-abc[:,k]==[0;1;0]
          push!(coords, 2)
          push!(ks, k) 
          break;
        elseif abc[:,indx+1]-abc[:,k]==[0;0;1]
          push!(coords, 3)
          push!(ks, k) 
          break;
        end
    end
  end

  xyz = [x y z]
  dexi = zeros(nnodes)
  deta = zeros(nnodes)
  dzeta = zeros(nnodes)
  for indx=1:size(abc,2)-1
    exi = xyz[:,coords[indx]]
    eta = copy(exi)
    zeta = copy(exi)

    if coords[indx]==1
      dexi = ones(nnodes)
      deta = zeros(nnodes)
      dzeta = zeros(nnodes)
    elseif coords[indx]==2
      dexi = zeros(nnodes)
      deta = ones(nnodes)
      dzeta = zeros(nnodes)
    elseif coords[indx]==3
      dexi = zeros(nnodes)
      deta = zeros(nnodes)
      dzeta = ones(nnodes)
    end
    
    v = diagm(exi)*V[:,ks[indx]]
    for t = 1:2
        s = 1/nnodes .* V[:,1:indx]'*v
        v = v - V[:,1:indx]*s
        Hes[1:indx,indx+1] =  Hes[1:indx,indx+1] + s
    end

    Hes[indx+1,indx+1] = 1/sqrt(nnodes).*norm(v,2)
    V[:,indx+1] = v ./Hes[indx+1,indx+1];

    if compute_grad
      s1 = 1/nnodes * Matrix(V[:,1:indx]')* diagm(exi)*V[:,ks[indx]]
      s2 =  s1 - 1/nnodes *Matrix(V[:,1:indx]')*V[:,1:indx]*s1
      Vdx[:,indx+1] = (diagm(dexi)*V[:,ks[indx]]
                      + diagm(exi)*Vdx[:,ks[indx]]- Vdx[:,1:indx]*s1
                      - Vdx[:,1:indx]*s2) ./Hes[indx+1,indx+1]
      Vdy[:,indx+1] = (diagm(deta)*V[:,ks[indx]] 
                      + diagm(eta)*Vdy[:,ks[indx]]- Vdy[:,1:indx]*s1
                      - Vdy[:,1:indx]*s2) ./Hes[indx+1,indx+1]
      Vdz[:,indx+1] = (diagm(dzeta)*V[:,ks[indx]] 
                      + diagm(zeta)*Vdz[:,ks[indx]]- Vdz[:,1:indx]*s1
                      - Vdz[:,1:indx]*s2) ./Hes[indx+1,indx+1]
    end

    # #another way to compute Q(:,indx+1) (more important form to obtain derivatives)
    # V[:,indx+1] = (diagm(eta)*V[:,ks[indx]] - V[:,1:indx]*s1 - V[:,1:indx]*s2)./Hes[indx+1,indx+1]
  end
  return V,Vdx,Vdy,Vdz,Hes
end

"""
### OrthoPoly.vandermonde_full{T}

Evaluate the Vandermonde and complements it with the nullspace to create a full rank matrix.

**Inputs**

* `p`: the maximum total degree of the polynomial
* `x`: locations at which to evaluate the derivative

**Outputs**

* `Vfull`: the full rank Vandermonde matrix

"""
function vandermonde_full(p::Int, x::Array{T}) where {T}
  dim = size(x,1)
  num_nodes = size(x,2)
  num_eq = binomial(dim+p, dim)
  V = zeros(num_nodes, num_eq)
  if dim == 2
    V, _, _ = vandermonde(p,x[1,:],x[2,:])
  elseif dim == 3
    V, _, _ = vandermonde(p,x[1,:],x[2,:],x[3,:])
  else
    Error("Unsupported dimenstion")
  end

  Vnull = nullspace(V')
  Vfull = [V Vnull]

  return Vfull
end

end #end of OrthoPoly class