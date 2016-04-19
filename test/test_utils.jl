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

  context("Testing SummationByParts.permuteface!") do
    # test 2D version
    permvec = [2, 1, 3, 4]
    vals = rand(5, 4)
    vals2 = copy(vals)  # permute a copy, leaving the original unchanged
    workarr = zeros(vals)

    SummationByParts.permuteface!(permvec, workarr, vals2)

    for i=1:4
      for j=1:5
        @fact vals[j, i] --> roughly(vals2[j, permvec[i] ], atol=1e-13)
      end
    end

    # test 1D version
    vals = rand(4)
    vals2 = copy(vals)
    workarr = zeros(vals)
    SummationByParts.permuteface!(permvec, workarr, vals2)

    for i=1:4
      @fact vals[i] --> roughly(vals2[permvec[i]], atol=1e-13)
    end

  end  # end context(testing permuteface!)

  context("Testing SummationByParts.permuteinterface!") do

    # test 3D version
    nfaces = 2
    ndofpernode = 5

    # create an SBP with non-trivial number of face nodes
    sbp = TriSBP{Float64}(degree = 3, reorder=false, internal=true)
    ref_verts = [-1. 1 -1; -1 -1 1]
    sbpface = TriFace{Float64}(sbp.degree, sbp.cub, ref_verts.')

    # create some data
    q = rand(ndofpernode, sbpface.numnodes, nfaces)
    q2 = copy(q)

    ifaces = Array(Interface, 2)
    for i=1:nfaces
      ifaces[i] = Interface(1,1,1,1,1)
    end

    SummationByParts.permuteinterface!(sbpface, ifaces, q2)

    for iface=1:nfaces
      for j=1:sbpface.numnodes
        for k=1:ndofpernode
          @fact q[k, j, iface] --> roughly(q2[k, sbpface.nbrperm[j, 1], iface], atol=1e-13)
        end
      end
    end

    # test 2D version
    q = rand(sbpface.numnodes, nfaces)
    q2 = copy(q)
    SummationByParts.permuteinterface!(sbpface, ifaces, q2)
    for iface=1:nfaces
      for j=1:sbpface.numnodes
          @fact q[j, iface] --> roughly(q2[sbpface.nbrperm[j, 1], iface], atol=1e-13)
      end
    end

  end  # end context(Testing permuteinterface!)
        
end
