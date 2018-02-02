facts("Testing SummationByParts Module (mapping Jacobian methods)...") do

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing calcMappingJacobian! ("string($TSBP)" method)") do
        # build a curvilinear Lagrangian element, and verify components of the
        # Jacobian and its determinant
        for p = 1:4
          sbp = ($TSBP)(degree=p)

          function mapping(ξ)
            return [0.5*(ξ[1]+1) + (0.5*(ξ[1]+1))^p;
                    0.5*(ξ[2]+1) + (0.5*(ξ[2]+1))^p + (0.5*(ξ[1]+1))^(p-1)]
          end
          function diffmapping(ξ)
            return [(0.5 + 0.5*p*(0.5*(ξ[1]+1))^(p-1)) 0.0;
                    (0.5*(p-1)*(0.5*(ξ[1]+1))^max(p-2,0)) (0.5 + 0.5*p*(0.5*(ξ[2]+1))^(p-1))]
          end      
          
          # function mapping(ξ)
          #   return [(0.5*(ξ[1]+1))^p; (0.5*(ξ[2]+1))^p + (0.5*(ξ[1]+1))^(p-1)]
          # end
          # function diffmapping(ξ)
          #   return [(0.5*p*(0.5*(ξ[1]+1))^(p-1)) 0.0;
          #           (0.5*(p-1)*(0.5*(ξ[1]+1))^max(p-2,0)) (0.5*p*(0.5*(ξ[2]+1))^(p-1))]
          # end

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

          x = calcnodes(sbp)
          for i = 1:sbp.numnodes
            dxdxi = diffmapping(x[:,i])
            @fact dxdxi[2,2] --> roughly(dξdx[1,1,i,1], atol=5e-14)
            @fact dxdxi[1,2] --> roughly(-dξdx[1,2,i,1], atol=5e-14)
            @fact dxdxi[2,1] --> roughly(-dξdx[2,1,i,1], atol=5e-14)
            @fact dxdxi[1,1] --> roughly(dξdx[2,2,i,1], atol=5e-14)
            @fact dxdxi[1,1,1]*dxdxi[2,2,1] - dxdxi[1,2,1]*dxdxi[2,1,1] -->
            roughly(1./jac[i], atol=5e-14)
          end
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      context("Testing calcMappingJacobian! ("string($TSBP)" method)") do
        # build a curvilinear Lagrangian element, and verify metric invariants are
        # satisfied
        for p = 1:4
          sbp = ($TSBP)(degree=p)
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
            facenormal!(sbpface, p+1, xref, view(xlag,:,:,f),
                        view(xsbp,:,:,f), view(nrm,:,:,f))
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
          # verify the metric invariants
          Qt = [sbp.Q[:,:,1].' sbp.Q[:,:,2].' sbp.Q[:,:,3].']
          metrics = zeros(3*sbp.numnodes)
          for di = 1:3
            for di2 = 1:3
              for i = 1:sbp.numnodes      
                metrics[i + (di2-1)*sbp.numnodes] = dξdx[di2,di,i]
              end
            end
            res = Qt*metrics - Eone[:,di]
            @fact res --> roughly(zeros(sbp.numnodes), atol=1e-14)
          end
        end
      end
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing calcMappingJacobianElement! ("string($TSBP)" method)") do
        # build a curvilinear Lagrangian element, and verify components of the
        # Jacobian and its determinant
        for p = 1:4
          sbp = ($TSBP)(degree=p)

          function mapping(ξ)
            return [0.5*(ξ[1]+1) + (0.5*(ξ[1]+1))^p;
                    0.5*(ξ[2]+1) + (0.5*(ξ[2]+1))^p + (0.5*(ξ[1]+1))^(p-1)]
          end
          function diffmapping(ξ)
            return [(0.5 + 0.5*p*(0.5*(ξ[1]+1))^(p-1)) 0.0;
                    (0.5*(p-1)*(0.5*(ξ[1]+1))^max(p-2,0)) (0.5 + 0.5*p*(0.5*(ξ[2]+1))^(p-1))]
          end      
          
          numdof = div((p+2)*(p+3),2)
          # set the coordinates of the reference and mapped nodes of the Lagrange
          # element
          xref = zeros(2,numdof)
          xlag = zeros(xref)
          ptr = 1
          for r = 0:p+1
            for j = 0:r
              i = r-j
              xref[1,ptr] = 2*i/(p+1) - 1
              xref[2,ptr] = 2*j/(p+1) - 1
              xlag[:,ptr] = mapping(xref[:,ptr])
              ptr += 1
            end
          end
          # compute the SBP nodes and the mapping Jacobian
          xsbp = zeros(2,sbp.numnodes)
          dξdx = zeros(2,2,sbp.numnodes)
          jac = zeros(sbp.numnodes)
          calcMappingJacobianElement!(sbp, p+1, xref, xlag, xsbp, dξdx, jac)
          
          x = calcnodes(sbp)
          for i = 1:sbp.numnodes
            dxdxi = diffmapping(x[:,i])
            @fact dxdxi[2,2] --> roughly(dξdx[1,1,i], atol=5e-14)
            @fact dxdxi[1,2] --> roughly(-dξdx[1,2,i], atol=5e-14)
            @fact dxdxi[2,1] --> roughly(-dξdx[2,1,i], atol=5e-14)
            @fact dxdxi[1,1] --> roughly(dξdx[2,2,i], atol=5e-14)
            @fact dxdxi[1,1]*dxdxi[2,2] - dxdxi[1,2]*dxdxi[2,1] -->
            roughly(1./jac[i], atol=5e-14)
          end
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      context("Testing calcMappingJacobianElement! ("string($TSBP)" method)") do
        # build a curvilinear Lagrangian element, and verify metric invariants are
        # satisfied
        for p = 1:4
          sbp = ($TSBP)(degree=p)
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
            facenormal!(sbpface, p+1, xref, view(xlag,:,:,f),
                        view(xsbp,:,:,f), view(nrm,:,:,f))
            for di = 1:3
              # for the given Lagrangian nodes, the face-normal is inward pointing,
              # so subtract to reverse sign
              E[sbpface.perm[:,f],sbpface.perm[:,f],di] -= 
              sbpface.interp*diagm(sbpface.wface.*vec(nrm[di,:,f]))*sbpface.interp.'
            end
          end
          Eone = zeros(sbp.numnodes,3)
          Eone = reshape(sum(E, 2), (sbp.numnodes,3))
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
          xlag = zeros(xref)
          for i = 1:numdof
            xlag[:,i] = mapping(xref[:,i])
          end
          # compute the SBP nodes and the mapping Jacobian
          xsbp = zeros(3,sbp.numnodes)
          dξdx = zeros(3,3,sbp.numnodes)
          jac = zeros(sbp.numnodes)
          calcMappingJacobianElement!(sbp, p+1, xref, xlag, xsbp, dξdx, jac, Eone)
          # verify the metric invariants
          Qt = [sbp.Q[:,:,1].' sbp.Q[:,:,2].' sbp.Q[:,:,3].']
          metrics = zeros(3*sbp.numnodes)
          for di = 1:3
            for di2 = 1:3
              for i = 1:sbp.numnodes      
                metrics[i + (di2-1)*sbp.numnodes] = dξdx[di2,di,i]
              end
            end
            res = Qt*metrics - Eone[:,di]
            @fact res --> roughly(zeros(sbp.numnodes), atol=1e-14)
          end
        end
      end
    end
  end

  # context("Testing metric invariants (TetSBP type)") do
  #   # build one element grid, and verify components of the Jacobian and its
  #   # determinant
  #   mapdegree = [1;2] 
  #   for p = 1:2 # metric invariants only satisfied for p=1 and p=2
  #     sbp = TetSBP{Float64}(degree=p)
  #     vtx = [0. 0. 0.; 2. 1. 0.; -1. 2. 0.; 0.5 0.5 2.]
  #     xi = zeros(Float64, (3,sbp.numnodes,1))
  #     xi[:,:,1] = calcnodes(sbp, vtx)
  #     for r = 0:p
  #       for k = 0:r
  #         for j = 0:r-k
  #           i = r-j-k
  #           x = deepcopy(xi)
  #           x += (xi[1,:,:].^i).*(xi[2,:,:].^j).*(xi[3,:,:].^k)
  #           dξdx = zeros(Float64, (3,3,sbp.numnodes,1))
  #           jac = zeros(Float64, (sbp.numnodes,1))
  #           mappingjacobian!(sbp, x, dξdx, jac)
  
  #           invariant = zeros(Float64, (3,sbp.numnodes,1))
  #           differentiate!(sbp, 1, view(dξdx,1,:,:,:), invariant)
  #           differentiate!(sbp, 2, view(dξdx,2,:,:,:), invariant)
  #           differentiate!(sbp, 3, view(dξdx,3,:,:,:), invariant)
  #           @fact invariant --> roughly(zeros(Float64, (3,sbp.numnodes,1)), atol=1e-13)
  #         end
  #       end
  #     end
  #   end
  # end

  context("Testing mappingjacobian! (TriFace method)") do
    # build a two element grid, and verify components of the Jacobian and its
    # determinant
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
      # verify on element 1
      @fact vec(dξdx[1,1,:,1,1]) --> roughly(zeros(sbpface.numnodes), atol=5e-13)
      @fact vec(dξdx[1,2,:,1,1]) --> roughly(2.0*ones(sbpface.numnodes), atol=5e-13)
      @fact vec(dξdx[2,1,:,1,1]) --> roughly(-2.0*ones(sbpface.numnodes), atol=5e-13)
      @fact vec(dξdx[2,2,:,1,1]) --> roughly(-2.0*ones(sbpface.numnodes), atol=5e-13)
      @fact vec(jac[:,1,1]) --> roughly(4.0*ones(sbpface.numnodes), atol=5e-13)
      # verify on element 2
      @fact vec(dξdx[1,1,:,2,1]) --> roughly(zeros(sbpface.numnodes), atol=5e-13)
      @fact vec(dξdx[1,2,:,2,1]) --> roughly(-2.0*ones(sbpface.numnodes), atol=5e-13)
      @fact vec(dξdx[2,1,:,2,1]) --> roughly(2.0*ones(sbpface.numnodes), atol=5e-13)
      @fact vec(dξdx[2,2,:,2,1]) --> roughly(2.0*ones(sbpface.numnodes), atol=5e-13)
      @fact vec(jac[:,2,1]) --> roughly(4.0*ones(sbpface.numnodes), atol=5e-13)
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      context("Testing mappingjacobian! ("string($TSBP)" method)") do
        # build a two element grid, and verify components of the Jacobian and its
        # determinant
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          vtx = [0. 0.; 1. 0.; 0. 1.]
          x = zeros(Float64, (2,sbp.numnodes,2))
          x[:,:,1] = calcnodes(sbp, vtx)
          vtx = [1. 0.; 1. 1.; 0. 1.]
          x[:,:,2] = calcnodes(sbp, vtx)
          dξdx = zeros(Float64, (2,2,sbp.numnodes,2))
          jac = zeros(Float64, (sbp.numnodes,2))
          mappingjacobian!(sbp, x, dξdx, jac)
          # verify on element 1
          @fact vec(dξdx[1,1,:,1]) --> roughly(0.5*ones(sbp.numnodes), atol=5e-13)
          @fact vec(dξdx[1,2,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-13)
          @fact vec(dξdx[2,2,:,1]) --> roughly(0.5*ones(sbp.numnodes), atol=5e-13)
          @fact vec(dξdx[2,1,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-13)
          @fact vec(jac[:,1]) --> roughly(4.0*ones(sbp.numnodes), atol=5e-13)
          # verify on element 2
          @fact vec(dξdx[1,1,:,2]) --> roughly(0.5*ones(sbp.numnodes), atol=5e-13)
          @fact vec(dξdx[1,2,:,2]) --> roughly(0.5*ones(sbp.numnodes), atol=5e-13)
          @fact vec(dξdx[2,2,:,2]) --> roughly(zeros(sbp.numnodes), atol=5e-13)
          @fact vec(dξdx[2,1,:,2]) --> roughly(-0.5*ones(sbp.numnodes), atol=5e-13)
          @fact vec(jac[:,2]) --> roughly(4.0*ones(sbp.numnodes), atol=5e-13)
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      context("Testing mappingjacobian! ("string($TSBP)" method)") do
        # build one element grid, and verify components of the Jacobian and its
        # determinant
        for p = 1:4
          sbp = ($TSBP)(degree=p)
          vtx = [0. 0. 0.; 2. 0. 0.; 0. 2. 0.; 0. 0. 2.]
          x = zeros(Float64, (3,sbp.numnodes,1))
          x[:,:,1] = calcnodes(sbp, vtx)
          dξdx = zeros(Float64, (3,3,sbp.numnodes,1))
          jac = zeros(Float64, (sbp.numnodes,1))
          mappingjacobian!(sbp, x, dξdx, jac)
          # dxi/dx = (1,0,0)
          @fact vec(dξdx[1,1,:,1]) --> roughly(ones(sbp.numnodes), atol=5e-12)
          @fact vec(dξdx[1,2,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-12)
          @fact vec(dξdx[1,3,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-12)
          # deta/dx = (0,1,0)
          @fact vec(dξdx[2,1,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-12)
          @fact vec(dξdx[2,2,:,1]) --> roughly(ones(sbp.numnodes), atol=5e-12)
          @fact vec(dξdx[2,3,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-12)
          # dzeta/dx = (0,1,0)
          @fact vec(dξdx[3,1,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-12)
          @fact vec(dξdx[3,2,:,1]) --> roughly(zeros(sbp.numnodes), atol=5e-12)
          @fact vec(dξdx[3,3,:,1]) --> roughly(ones(sbp.numnodes), atol=5e-12)
          # jac = 1
          @fact vec(jac[:,1]) --> roughly(ones(sbp.numnodes), atol=5e-12)
        end
      end
    end
  end

end
