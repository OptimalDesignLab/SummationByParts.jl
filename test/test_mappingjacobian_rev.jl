facts("Testing SummationByParts Module (reverse-diff of mapping Jacobian)...") do
#=
  context("Testing calcMappingJacobian_rev! (TriSBP method)") do
    # build a curvilinear Lagrangian element, and verify against randomly
    # perturbed Lagrangian nodes using complex step
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p, internal=false)

      function mapping(ξ)
        return [0.5*(ξ[1]+1) + (0.5*(ξ[1]+1))^p;
                0.5*(ξ[2]+1) + (0.5*(ξ[2]+1))^p + (0.5*(ξ[1]+1))^(p-1)]
      end
      function diffmapping(ξ)
        return [(0.5 + 0.5*p*(0.5*(ξ[1]+1))^(p-1)) 0.0;
                (0.5*(p-1)*(0.5*(ξ[1]+1))^max(p-2,0))
                (0.5 + 0.5*p*(0.5*(ξ[2]+1))^(p-1))]
      end      

      numdof = div((p+2)*(p+3),2)
      # set the coordinates of the reference and mapped nodes of the Lagrange
      # element
      xref = zeros(2,numdof)
      xlag = zeros(2,numdof,1)
      ptr = 1
      for r = 0:p+1
        for j = 0:r
          i = r-j
          xref[1,ptr] = 2*i/(p+1) - 1
          xref[2,ptr] = 2*j/(p+1) - 1
          xlag[:,ptr,1] = mapping(xref[:,ptr])
          ptr += 1
        end
      end
      # compute the SBP nodes and the mapping Jacobian
      xsbp = zeros(2,sbp.numnodes,1)
      dξdx = zeros(2,2,sbp.numnodes,1)
      jac = zeros(sbp.numnodes,1)
      calcMappingJacobian!(sbp, p+1, xref, xlag, xsbp, dξdx, jac)

      # These arrays are meant to mimic the derivative of the residual with
      # respect to xsbp, dξdx, and jac
      dRdxsbp = randn(2,sbp.numnodes,1)
      dRddξdx = randn(2,2,sbp.numnodes,1)
      dRdjac = randn(sbp.numnodes,1)
      
      # evaluate complex-step gradient to get dRdxlag
      xsbp_c = zeros(Complex{Float64}, (2,sbp.numnodes,1))
      dξdx_c = zeros(Complex{Float64}, (2,2,sbp.numnodes,1))
      jac_c = zeros(Complex{Float64}, (sbp.numnodes,1))
      xref_c = complex(xref, 0.0)
      dRdxlag_cmplx = zeros(2,numdof,1)
      eps_cmplx = 1e-60
      for di = 1:2
        for nd = 1:numdof
          xlag_c = complex(xlag, 0.0)
          xlag_c[di,nd,1] += complex(0.0, eps_cmplx)
          # compute the complex perturbed SBP nodes and the mapping Jacobian
          calcMappingJacobian!(sbp, p+1, xref_c, xlag_c, xsbp_c, dξdx_c, jac_c)
          dRdxlag_cmplx[di,nd,1] +=
            (dot(vec(dRdxsbp[:,:,1]),vec(imag(xsbp_c))) +
             dot(vec(dRddξdx[:,:,:,1]),vec(imag(dξdx_c))) +
             dot(vec(dRdjac[:,1]),vec(imag(jac_c))))
        end
      end
      scale!(dRdxlag_cmplx,1/eps_cmplx)

      # use reverse mode to get dRdxlag
      dRdxlag = zeros(2,numdof,1)
      Eone_bar = zeros(sbp.numnodes,2,1)
      calcMappingJacobian_rev!(sbp, p+1, xref, xlag, dRdxlag, dRdxsbp, dRddξdx,
                               dRdjac, Eone_bar)

      @fact dRdxlag --> roughly(dRdxlag_cmplx, atol=1e-14)
    end
  end

  context("Testing calcMappingJacobian_rev! (TetSBP method)") do
    # build a curvilinear Lagrangian element, and verify against randomly
    # perturbed Lagrangian nodes using complex step
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p, internal=false)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      function mapping(ξ)
        x = 0.5*(ξ[1]+1) + (0.5*(ξ[1]+1))^(p+1)
        y = 0.5*(ξ[2]+1) + (0.5*(ξ[2]+1))^(p+1)
        z = 0.5*(ξ[3]+1) + (0.5*(ξ[3]+1))^(p+1)
        fac = 0.5
        return [fac*x - fac*y; fac*x + fac*y; z]
      end
      # set the coordinates of the Lagrangian face nodes in reference and
      # physical space
      numdof = div((p+2)*(p+3),2)
      xref = zeros(2,numdof)
      ptr = 1
      for r = 0:p+1
        for j = 0:r
          i = r-j
          xref[1,ptr] = 2*i/(p+1) - 1.0
          xref[2,ptr] = 2*j/(p+1) - 1.0
          ptr += 1
        end
      end
      xlag = zeros(3,numdof,4)
      for i = 1:numdof
        xlag[:,i,1] = mapping([xref[1,i]; xref[2,i]; -1.0])
        xlag[:,i,2] = mapping([xref[2,i]; -1.0; xref[1,i]])
        xlag[:,i,3] = mapping([xref[2,i]; xref[1,i];
                               -1.0 - xref[1,i] - xref[2,i]])
        xlag[:,i,4] = mapping([-1.0; xref[1,i]; xref[2,i]]) 
      end
      # get the SBP face nodes and the normal vector
      xsbp = zeros(3,sbpface.numnodes,4)
      nrm = zeros(3,sbpface.numnodes,4)
      E = zeros(sbp.numnodes,sbp.numnodes,3)
      for f = 1:4
        facenormal!(sbpface, p+1, xref, sview(xlag,:,:,f),
                    sview(xsbp,:,:,f), sview(nrm,:,:,f))
        for di = 1:3
          # for the given Lagrangian nodes, the face-normal is inward pointing,
          # so subtract to reverse sign
          E[sbpface.perm[:,f],sbpface.perm[:,f],di] -= 
          sbpface.interp*diagm(sbpface.wface.*vec(nrm[di,:,f]))*sbpface.interp.'
        end
      end
      Eone = zeros(sbp.numnodes,3,1)
      Eone = reshape(sum(E, 2), (sbp.numnodes,3,1))
      # now set the coordinates of the Lagrangian element nodes in reference and
      # physical space
      numdof = binomial(p+1+3,3)
      xref = zeros(3,numdof)
      ptr = 1
      for r = 0:(p+1)
        for k = 0:r
          for j = 0:r-k
            i = r-j-k
            xref[1,ptr] = 2*i/(p+1) - 1.0
            xref[2,ptr] = 2*j/(p+1) - 1.0
            xref[3,ptr] = 2*k/(p+1) - 1.0
            ptr += 1
          end
        end
      end
      xlag = zeros(3,numdof,1)
      for i = 1:numdof
        xlag[:,i,1] = mapping(xref[:,i])
      end
      # compute the SBP nodes and the mapping Jacobian
      xsbp = zeros(3,sbp.numnodes,1)
      dξdx = zeros(3,3,sbp.numnodes,1)
      jac = zeros(sbp.numnodes,1)
      calcMappingJacobian!(sbp, p+1, xref, xlag, xsbp, dξdx, jac, Eone)

      # These arrays are meant to mimic the derivative of the residual with
      # respect to xsbp, dξdx, and jac
      dRdxsbp = randn(3,sbp.numnodes,1)
      dRddξdx = randn(3,3,sbp.numnodes,1)
      dRdjac = randn(sbp.numnodes,1)
      
      # evaluate complex-step gradient to get dRdxlag
      xsbp_c = zeros(Complex{Float64}, (3,sbp.numnodes,1))
      dξdx_c = zeros(Complex{Float64}, (3,3,sbp.numnodes,1))
      jac_c = zeros(Complex{Float64}, (sbp.numnodes,1))
      xref_c = complex(xref, 0.0)                
      Eone_c = complex(Eone, 0.0)
      xlag_c = complex(xlag, 0.0)      
      dRdxlag_cmplx = zeros(3,numdof,1)
      eps_cmplx = 1e-60
      for di = 1:3
        for nd = 1:numdof
          xlag_c[di,nd,1] += complex(0.0, eps_cmplx)
          # compute the complex perturbed SBP nodes and the mapping Jacobian
          calcMappingJacobian!(sbp, p+1, xref_c, xlag_c, xsbp_c, dξdx_c, jac_c,
                               Eone_c)
          xlag_c[di,nd,1] -= complex(0.0, eps_cmplx)
          dRdxlag_cmplx[di,nd,1] +=
            (dot(vec(dRdxsbp[:,:,1]),vec(imag(xsbp_c))) +
             dot(vec(dRddξdx[:,:,:,1]),vec(imag(dξdx_c))) +
             dot(vec(dRdjac[:,1]),vec(imag(jac_c))))
        end
      end
      scale!(dRdxlag_cmplx,1/eps_cmplx)      

      # evaluate complex-step gradient to get dRdEone
      dRdEone_cmplx = zeros(sbp.numnodes,3,1)
      for di = 1:3
        for i = 1:sbp.numnodes
          Eone_c[i,di,1] += complex(0.0, eps_cmplx)
          # compute the complex perturbed SBP nodes and the mapping Jacobian
          calcMappingJacobian!(sbp, p+1, xref_c, xlag_c, xsbp_c, dξdx_c, jac_c,
                               Eone_c)
          Eone_c[i,di,1] -= complex(0.0, eps_cmplx)
          dRdEone_cmplx[i,di,1] +=
            (dot(vec(dRdxsbp[:,:,1]),vec(imag(xsbp_c))) +
             dot(vec(dRddξdx[:,:,:,1]),vec(imag(dξdx_c))) +
             dot(vec(dRdjac[:,1]),vec(imag(jac_c))))
        end
      end
      scale!(dRdEone_cmplx,1/eps_cmplx)      
      
      # use reverse mode to get dRdxlag and dRdEone
      dRdxlag = zeros(3,numdof,1)
      dRdEone = zeros(sbp.numnodes,3,1)
      calcMappingJacobian_rev!(sbp, p+1, xref, xlag, dRdxlag, dRdxsbp, dRddξdx,
                               dRdjac, dRdEone)
      @fact dRdxlag --> roughly(dRdxlag_cmplx, atol=1e-14)
      @fact dRdEone --> roughly(dRdEone_cmplx, atol=1e-14)
    end
  end
=#
  context("Testing mappingjacobian_rev! (TriSBP method)") do
    # build one element grid, and verify against randomly
    # perturbed SBP nodes using complex step
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      vtx = 2.0.*sbp.vtx
      x = zeros(Float64, (2,sbp.numnodes,1))
      x[:,:,1] = calcnodes(sbp, vtx)
      x_bar = zeros(x)
      dξdx_bar = randn(2,2,sbp.numnodes,1)
      jac_bar = randn(sbp.numnodes,1)

      x_bar_cs = zeros(x)
      x_cmplx = complex(x, zeros(x))
      dξdx_cmplx = zeros(Complex128, (2,2,sbp.numnodes,1))
      jac_cmplx = zeros(Complex128, (sbp.numnodes,1))
      ceps = 1e-60
      for i = 1:sbp.numnodes
        for di = 1:2
          x_cmplx[di,i,1] += complex(0.0, ceps)
          mappingjacobian!(sbp, x_cmplx, dξdx_cmplx, jac_cmplx)
          x_cmplx[di,i,1] -= complex(0.0, ceps)
          x_bar_cs[di,i,1] += (dot(vec(dξdx_bar),vec(imag(dξdx_cmplx))) +
                               dot(vec(jac_bar),vec(imag(jac_cmplx))))
        end
      end
      scale!(x_bar_cs, 1/ceps)

      # this needs to come at the end, because dξdx_bar is modified
      mappingjacobian_rev!(sbp, x, x_bar, dξdx_bar, jac_bar)
      
      @fact x_bar --> roughly(x_bar_cs, atol=1e-14)      
    end
  end

  context("Testing mappingjacobian_rev! (TetSBP method)") do
    # build one element grid, and verify against randomly
    # perturbed SBP nodes using complex step
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      vtx = [0. 0. 0.; 2. 0. 0.; 0. 2. 0.; 0. 0. 2.]
      x = zeros(Float64, (3,sbp.numnodes,1))
      x[:,:,1] = calcnodes(sbp, vtx)
      x_bar = zeros(x)
      dξdx_bar = randn(3,3,sbp.numnodes,1)
      jac_bar = randn(sbp.numnodes,1)

      x_bar_cs = zeros(x)
      x_cmplx = complex(x, zeros(x))
      dξdx_cmplx = zeros(Complex128, (3,3,sbp.numnodes,1))
      jac_cmplx = zeros(Complex128, (sbp.numnodes,1))
      ceps = 1e-60
      for i = 1:sbp.numnodes
        for di = 1:3
          x_cmplx[di,i,1] += complex(0.0, ceps)
          mappingjacobian!(sbp, x_cmplx, dξdx_cmplx, jac_cmplx)
          x_cmplx[di,i,1] -= complex(0.0, ceps)
          x_bar_cs[di,i,1] += (dot(vec(dξdx_bar),vec(imag(dξdx_cmplx))) +
                               dot(vec(jac_bar),vec(imag(jac_cmplx))))
        end
      end
      scale!(x_bar_cs, 1/ceps)

      # this needs to come at the end, because dξdx_bar may be modified
      mappingjacobian_rev!(sbp, x, x_bar, dξdx_bar, jac_bar)
      
      @fact x_bar --> roughly(x_bar_cs, atol=1e-14)      
    end
  end

end
