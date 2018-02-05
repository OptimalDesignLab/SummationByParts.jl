facts("Testing OrthoPoly Module...") do

  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing OrthoPoly.lglnodes for DataType "string($T)) do
        x, w = OrthoPoly.lglnodes(2, ($T))
        @fact x --> roughly(($T)[-1, 0, 1], atol=10*eps(typeof(real(one($T)))) )
        @fact w --> roughly(($T)[1/3, 4/3, 1/3],
                            atol=10*eps(typeof(real(one($T)))) )
        x, w = OrthoPoly.lglnodes(3, ($T))
        @fact x --> roughly(($T)[-1, -1/sqrt(5), 1/sqrt(5), 1],
                            atol=10*eps(typeof(real(one($T)))) )
        @fact w --> roughly(($T)[1/6, 5/6, 5/6, 1/6],
                            atol=10*eps(typeof(real(one($T)))) )
        x, w = OrthoPoly.lglnodes(4, ($T))
        @fact x --> roughly(($T)[-1, -sqrt(3/7), 0, sqrt(3/7), 1], 
                            atol=10*eps(typeof(real(one($T)))) )
        @fact w --> roughly(($T)[1/10, 49/90, 32/45, 49/90, 1/10],
                            atol=10*eps(typeof(real(one($T)))) )
      end
    end
  end

  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing OrthoPoly.lgnodes for DataType "string($T)) do
        x, w = OrthoPoly.lgnodes(2, ($T))
        @fact x --> roughly(($T)[-1/sqrt(3), 1/sqrt(3)],
                            atol=10*eps(typeof(real(one($T)))) )
        @fact w --> roughly(($T)[1, 1],
                            atol=10*eps(typeof(real(one($T)))) )
        x, w = OrthoPoly.lgnodes(3, ($T))
        @fact x --> roughly(($T)[-sqrt(3/5), 0, sqrt(3/5)],
                            atol=10*eps(typeof(real(one($T)))) )
        @fact w --> roughly(($T)[5/9, 8/9, 5/9],
                            atol=10*eps(typeof(real(one($T)))) )
        x, w = OrthoPoly.lgnodes(4, ($T))
        @fact x --> roughly(($T)[-sqrt(3/7 + (2/7)*sqrt(6/5)),
                                 -sqrt(3/7 - (2/7)*sqrt(6/5)),
                                 sqrt(3/7 - (2/7)*sqrt(6/5)),
                                 sqrt(3/7 + (2/7)*sqrt(6/5))], 
                            atol=10*eps(typeof(real(one($T)))) )
        @fact w --> roughly(($T)[(18-sqrt(30))/36, (18+sqrt(30))/36,
                                 (18+sqrt(30))/36, (18-sqrt(30))/36],
                            atol=10*eps(typeof(real(one($T)))) )
      end
    end
  end

  for T = (Float32, Float64, Complex64, Complex128)
    @eval begin
      context("Testing OrthoPoly.jacobipoly for DataType "string($T)) do
        # compare against 5th degree Legendre polynomial
        x = ($T)[-1:0.1:1;]
        L = OrthoPoly.jacobipoly(x, 0.0, 0.0, 5)
        L ./= L[21] # necessary for standardization
        @fact L --> roughly( (63.*x.^5 -70.*x.^3 + 15.*x)/8,
                            atol=10*eps(typeof(real(one($T)))) )
      end
    end
  end

  context("Testing OrthoPoly.diffjacobipoly") do
    alpha = 0.3
    beta = 0.65
    x = Float64[-1:0.1:1;]
    eps_step = 1e-60
    for N = 1:10
      dP = OrthoPoly.diffjacobipoly(x, alpha, beta, N)
      xc = complex.(x, eps_step)
      Pc = OrthoPoly.jacobipoly(xc, alpha, beta, N)
      dP_cmplx = imag(Pc)/eps_step
      @fact dP --> roughly(dP_cmplx, atol=1e-15)
    end
  end

  context("Testing OrthoPoly.proriolpoly (Triangle method)") do
    # Integrate the product of the polynomials, and verify that they are
    # orthonormal

    # compute Gauss-Legendre quadrature points for degenerate square
    x_gl = Float64[-sqrt(5 + 2*sqrt(10/7))/3, -sqrt(5 - 2*sqrt(10/7))/3, 0,
            sqrt(5 - 2*sqrt(10/7))/3, sqrt(5 + 2*sqrt(10/7))/3]
    w_gl = Float64[(322-13*sqrt(70))/900, (322+13*sqrt(70))/900, 128/225,
                   (322+13*sqrt(70))/900, (322-13*sqrt(70))/900]
    ptr = 1
    x_q = zeros((5*5))
    y_q = zeros(x_q)
    w_q = zeros(x_q)
    for j = 1:5
      for i = 1:5
        y_q[ptr] = x_gl[j]
        x_q[ptr] = 0.5*(1-y_q[ptr])*x_gl[i] - 0.5*(y_q[ptr]+1)
        w_q[ptr] = w_gl[j]*0.5*(1-y_q[ptr])*w_gl[i]
        ptr += 1
      end
    end
    # evaluate Proriol orthogonal polynomials up to degree q
    q = 4
    N = convert(Int, (q+1)*(q+2)/2)
    P = zeros((5*5, N))
    ptr = 1
    for r = 0:q
      for j = 0:r
        i = r-j
        P[:,ptr] = OrthoPoly.proriolpoly(x_q, y_q, i, j)
        ptr += 1
      end
    end
    # verify that integral inner products produce identity
    innerprod = P.'*diagm(w_q)*P
    @fact innerprod --> roughly(eye(N), atol=1e-14)
  end

  context("Testing OrthoPoly.proriolpoly (Tetrahedron method)") do
    # Integrate the product of the polynomials, and verify that they are
    # orthonormal

    # compute Gauss-Legendre quadrature points for a degenerate cube
    x_gl = Float64[-sqrt(5 + 2*sqrt(10/7))/3, -sqrt(5 - 2*sqrt(10/7))/3, 0,
            sqrt(5 - 2*sqrt(10/7))/3, sqrt(5 + 2*sqrt(10/7))/3]
    w_gl = Float64[(322-13*sqrt(70))/900, (322+13*sqrt(70))/900, 128/225,
                   (322+13*sqrt(70))/900, (322-13*sqrt(70))/900]
    ptr = 1
    x_q = zeros((5*5*5))
    y_q = zeros(x_q)
    z_q = zeros(x_q)
    w_q = zeros(x_q)
    for k = 1:5
      for j = 1:5
        for i = 1:5
          z_q[ptr] = x_gl[k]
          y_q[ptr] = 0.5*(1 + x_gl[j])*(1 - z_q[ptr]) - 1
          x_q[ptr] = -1 - 0.5*(1 + x_gl[i])*(y_q[ptr] + z_q[ptr])
          w_q[ptr] = -0.25*(1 - z_q[ptr])*(y_q[ptr] + z_q[ptr])*
          w_gl[i]*w_gl[j]*w_gl[k]
          ptr += 1
        end
      end
    end
    # evaluate Proriol orthogonal polynomials up to degree q
    q = 3
    N = convert(Int, (q+1)*(q+2)*(q+3)/6)
    P = zeros((5*5*5, N))
    ptr = 1
    for r = 0:q
      for k = 0:r
        for j = 0:r-k
          i = r-j-k
          P[:,ptr] = OrthoPoly.proriolpoly(x_q, y_q, z_q, i, j, k)
          ptr += 1
        end
      end
    end
    # verify that integral inner products produce identity
    innerprod = P.'*diagm(w_q)*P
    @fact innerprod --> roughly(eye(N), atol=1e-14)
  end

  context("Testing OrthoPoly.diffproriolpoly") do
    x_lgl = Float64[-1 -sqrt(3/7) 0 sqrt(3/7) 1]
    ptr = 1
    x = zeros((5*5))
    y = zeros(x)
    for j = 1:5
      for i = 1:5
        y[ptr] = x_lgl[j]
        x[ptr] = 0.5*(1-y[ptr])*x_lgl[i] - 0.5*(y[ptr]+1)
        ptr += 1
      end
    end
    
    eps_step = 1e-60
    q = 10
    N = convert(Int, (q+1)*(q+2)/2)
    for r = 0:q
      for j = 0:r
        i = r-j
        dPdx, dPdy = OrthoPoly.diffproriolpoly(x, y, i, j)
        xc = complex.(x, 0)
        yc = complex.(y, -eps_step)
        Pc = OrthoPoly.proriolpoly(xc, yc, i, j)
        @fact dPdy --> roughly(-imag(Pc)/eps_step, atol=1e-14)
        xc -= eps_step*im
        Pc = OrthoPoly.proriolpoly(xc, yc, i, j)
        @fact dPdx --> roughly(-imag(Pc)/eps_step - dPdy, atol=1e-12)
      end
    end
  end

end
