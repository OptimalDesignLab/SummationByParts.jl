facts("Testing SummationByParts Module (utils.jl file)...") do

  context("Testing SummationByParts.buildinterpolation (TriSymCub method)") do
    # this checks that polynomials of total degree d are reconstructed accurately
    numpoints = 3
    for d = 1:4
      sbp = TriSBP{Float64}(degree=d)
      x = 2.*rand(2,numpoints) - 1.0
      R = SummationByParts.buildinterpolation(sbp, x)
      xsbp = SymCubatures.calcnodes(sbp.cub, sbp.vtx)
      # loop over all monomials
      for r = 0:d
        for j = 0:r
          i = r-j
          u = vec((x[1,:].^i).*(x[2,:].^j))
          usbp = vec((xsbp[1,:].^i).*(xsbp[2,:].^j))
          uinterp = R*usbp
          @fact uinterp --> roughly(u, atol=1e-14)
        end
      end
    end
  end

  context("Testing SummationByParts.buildinterpolation (TriSymCub method, internal=true)") do
    # this checks that polynomials of total degree d are reconstructed accurately
    numpoints = 3
    for d = 1:4
      sbp = TriSBP{Float64}(degree=d, reorder=false, internal=true)
      x = 2.*rand(2,numpoints) - 1.0
      R = SummationByParts.buildinterpolation(sbp, x)
      xsbp = SymCubatures.calcnodes(sbp.cub, sbp.vtx)
      # loop over all monomials
      for r = 0:d
        for j = 0:r
          i = r-j
          u = vec((x[1,:].^i).*(x[2,:].^j))
          usbp = vec((xsbp[1,:].^i).*(xsbp[2,:].^j))
          uinterp = R*usbp
          @fact uinterp --> roughly(u, atol=1e-14)
        end
      end
    end
  end

end