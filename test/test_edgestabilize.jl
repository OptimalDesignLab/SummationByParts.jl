facts("Testing SummationByParts Module (edge-stabilization methods)...") do

  context("Testing SummationByParts.edgestabilize! (TriFace, scalar field method)") do
    # build a two element grid, and verify that polynomials of degree p or less
    # vanish when edgestabilize is applied
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, [-1. -1.; 1. -1.; -1. 1.])
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x = zeros(Float64, (2,sbp.numnodes,2))
      x[:,:,1] = SummationByParts.SymCubatures.calcnodes(sbp.cub, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SummationByParts.SymCubatures.calcnodes(sbp.cub, vtx)
      dξdx = zeros(Float64, (2,2,sbpface.numnodes,2,1))
      jac = zeros(Float64, (sbpface.numnodes,2,1))
      ifaces = Array(Interface, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      mappingjacobian!(sbpface, ifaces, x, dξdx, jac)
      # build the normal vector to the face; (<dξdx,dηdx> , <dηdx,dηdx>)
      dirvec = zeros(Float64, (2,sbpface.numnodes,2,1))
      for LR = 1:2
        for i = 1:sbpface.numnodes
          dirvec[1,i,LR,1] = dξdx[1,1,i,LR,1]*dξdx[2,1,i,LR,1] + 
          dξdx[1,2,i,LR,1]*dξdx[2,2,i,LR,1]
          dirvec[2,i,LR,1] = dξdx[2,1,i,LR,1]*dξdx[2,1,i,LR,1] +
          dξdx[2,2,i,LR,1]*dξdx[2,2,i,LR,1]
        end
      end
      # scaling factor of 1
      tau = ones(Float64, (sbpface.numnodes,1))
      u = zeros(Float64, (sbp.numnodes,2))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          res = zeros(u)
          edgestabilize!(sbpface, ifaces, dirvec, tau, u, res)
          @fact res --> roughly(zeros(res), atol=1e-10) # !!! note the size of atol
        end
      end
    end
  end

end
