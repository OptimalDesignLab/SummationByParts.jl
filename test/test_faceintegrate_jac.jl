facts("Testing SummationByParts Module (Jacobian of face integration methods)...") do
  
  for TSBP = ((getLineSegSBPLobbato,getLineSegFace,1),
              (getLineSegSBPLegendre,getLineSegFace,1),
              (getTriSBPGamma,TriFace{Float64},2),
              (getTriSBPOmega,TriFace{Float64},2),
              (getTriSBPDiagE,getTriFaceForDiagE,2),
              (getTetSBPGamma,TetFace{Float64},3),
              (getTetSBPOmega,TetFace{Float64},3),
              (getTetSBPDiagE,getTetFaceForDiagE,3))
    @eval begin
      context("Testing interiorFaceIntegrate_jac! ("string($TSBP[1])" scalar field method)") do
        for p = 1:4
          # create a two element mesh with random orientation
          sbp = ($TSBP[1])(degree=p)
          sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          ifaces = Array(Interface, 1)
          if size(sbp.Q,3) == 1
            ifaces[1] = Interface(1,2,rand(1:2),rand(1:2),1)
          elseif size(sbp.Q,3) == 2
            ifaces[1] = Interface(1,2,rand(1:3),rand(1:3),1)
          else
            ifaces[1] = Interface(1,2,rand(1:4),rand(1:4),rand(1:3))
          end

          u = randn(sbp.numnodes,2)
          
          # get the Jacobians directly
          uface = zeros(sbpface.numnodes,2)
          interiorFaceInterpolate!(sbpface, ifaces[1],
                                   view(u,:,ifaces[1].elementL),
                                   view(u,:,ifaces[1].elementR),
                                   view(uface,:,ifaces[1].elementL),
                                   view(uface,:,ifaces[1].elementR))
          flux_jac1 = zeros(Float64, (sbpface.numnodes))
          flux_jac2 = zeros(Float64, (sbpface.numnodes))
          for i = 1:sbpface.numnodes
            # flux = cos(uL)*sin(uR)
            flux_jac1[i] = -sin(uface[i,1])*sin(uface[i,2])
            flux_jac2[i] = cos(uface[i,1])*cos(uface[i,2])
          end
          dR1du1 = zeros(Float64, (sbp.numnodes, sbp.numnodes))
          dR1du2 = zeros(dR1du1)
          dR2du1 = zeros(dR1du1)
          dR2du2 = zeros(dR1du1)
          interiorFaceIntegrate_jac!(sbpface, ifaces[1],
                                     flux_jac1, flux_jac2,
                                     dR1du1, dR1du2, dR2du1, dR2du2)

          # get the Jacobians using complex step
          u_c = complex(u)
          uface_c = complex(uface)
          flux_c = zeros(Complex128, (sbpface.numnodes))
          res_c = zeros(u_c)
          dR1du1_cmplx = zeros(dR1du1)
          dR1du2_cmplx = zeros(dR1du1)
          dR2du1_cmplx = zeros(dR1du1)
          dR2du2_cmplx = zeros(dR1du1)
          ceps = 1e-60
          for e = 1:2
            for i = 1:sbp.numnodes
              u_c[i,e] += complex(0.0, ceps)
              # interpolate to face
              interiorFaceInterpolate!(sbpface, ifaces[1],
                                       view(u_c,:,ifaces[1].elementL),
                                       view(u_c,:,ifaces[1].elementR),
                                       view(uface_c,:,ifaces[1].elementL),
                                       view(uface_c,:,ifaces[1].elementR))
              # compute the (complexified) flux
              for j = 1:sbpface.numnodes
                flux_c[j] = cos(uface_c[j,1])*sin(uface_c[j,2])
              end
              # get residuals
              fill!(res_c, zero(Complex128))
              interiorFaceIntegrate!(sbpface, ifaces[1], flux_c,
                                     view(res_c,:,ifaces[1].elementL),
                                     view(res_c,:,ifaces[1].elementR))
              if e == 1
                dR1du1_cmplx[:,i] = imag(res_c[:,1])/ceps
                dR2du1_cmplx[:,i] = imag(res_c[:,2])/ceps
              else
                dR1du2_cmplx[:,i] = imag(res_c[:,1])/ceps
                dR2du2_cmplx[:,i] = imag(res_c[:,2])/ceps
              end
              u_c[i,e] -= complex(0.0, ceps)
            end
          end

          # check for equality
          @fact dR1du1 --> roughly(dR1du1_cmplx, atol=1e-15)
          @fact dR1du2 --> roughly(dR1du2_cmplx, atol=1e-15)
          @fact dR2du1 --> roughly(dR2du1_cmplx, atol=1e-15)
          @fact dR2du2 --> roughly(dR2du2_cmplx, atol=1e-15)          
        end
      end
    end
  end

  for TSBP = ((getLineSegSBPLobbato,getLineSegFace,1),
              (getLineSegSBPLegendre,getLineSegFace,1),
              (getTriSBPGamma,TriFace{Float64},2),
              (getTriSBPOmega,TriFace{Float64},2),
              (getTriSBPDiagE,getTriFaceForDiagE,2),
              (getTetSBPGamma,TetFace{Float64},3),
              (getTetSBPOmega,TetFace{Float64},3),
              (getTetSBPDiagE,getTetFaceForDiagE,3))
    @eval begin
      context("Testing interiorFaceIntegrate_jac! ("string($TSBP[1])" vector field method)") do
        for p = 1:4
          # create a two element mesh with random orientation
          sbp = ($TSBP[1])(degree=p)
          sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          ifaces = Array(Interface, 1)
          if size(sbp.Q,3) == 1
            ifaces[1] = Interface(1,2,rand(1:2),rand(1:2),1)
          elseif size(sbp.Q,3) == 2
            ifaces[1] = Interface(1,2,rand(1:3),rand(1:3),1)
          else
            ifaces[1] = Interface(1,2,rand(1:4),rand(1:4),rand(1:3))
          end

          u = randn(2,sbp.numnodes,2)
          
          # get the Jacobians directly
          uface = zeros(2,sbpface.numnodes,2)
          interiorFaceInterpolate!(sbpface, ifaces[1],
                                   view(u,:,:,ifaces[1].elementL),
                                   view(u,:,:,ifaces[1].elementR),
                                   view(uface,:,:,ifaces[1].elementL),
                                   view(uface,:,:,ifaces[1].elementR))
          flux_jac1 = zeros(Float64, (2,2,sbpface.numnodes))
          flux_jac2 = zeros(Float64, (2,2,sbpface.numnodes))
          for i = 1:sbpface.numnodes
            # u1 = 0.5*(uL[1] + uR[1])
            # u2 = (2*uL[2]*uR[2]/(uL[2] + uR[2])
            # flux = [u1*cos(u2); u1*u1*sin(u2)]
            u1 = 0.5*(uface[1,i,1] + uface[1,i,2])
            u2 = (2.0*uface[2,i,1]*uface[2,i,2])/(uface[2,i,1] + uface[2,i,2])
            flux_jac1[1,1,i] = 0.5*cos(u2)
            flux_jac1[1,2,i] = -u1*sin(u2)*(u2/uface[2,i,1]
                                            - u2/(uface[2,i,1] + uface[2,i,2]))
            flux_jac1[2,1,i] = 2.0*u1*sin(u2)*0.5
            flux_jac1[2,2,i] = u1*u1*cos(u2)*(u2/uface[2,i,1]
                                              - u2/(uface[2,i,1] + uface[2,i,2]))
            
            flux_jac2[1,1,i] = 0.5*cos(u2)
            flux_jac2[1,2,i] = -u1*sin(u2)*(u2/uface[2,i,2]
                                            - u2/(uface[2,i,1] + uface[2,i,2]))
            flux_jac2[2,1,i] = 2.0*u1*sin(u2)*0.5
            flux_jac2[2,2,i] = u1*u1*cos(u2)*(u2/uface[2,i,2]
                                              - u2/(uface[2,i,1] + uface[2,i,2]))
          end
          dR1du1 = zeros(Float64, (2, 2, sbp.numnodes, sbp.numnodes))
          dR1du2 = zeros(dR1du1)
          dR2du1 = zeros(dR1du1)
          dR2du2 = zeros(dR1du1)
          interiorFaceIntegrate_jac!(sbpface, ifaces[1],
                                     flux_jac1, flux_jac2,
                                     dR1du1, dR1du2, dR2du1, dR2du2)

          # get the Jacobians using complex step
          u_c = complex(u)
          uface_c = complex(uface)
          flux_c = zeros(Complex128, (2,sbpface.numnodes))
          res_c = zeros(u_c)
          dR1du1_cmplx = zeros(dR1du1)
          dR1du2_cmplx = zeros(dR1du1)
          dR2du1_cmplx = zeros(dR1du1)
          dR2du2_cmplx = zeros(dR1du1)
          ceps = 1e-60
          for e = 1:2
            for i = 1:sbp.numnodes
              for p = 1:2
                u_c[p,i,e] += complex(0.0, ceps)
                # interpolate to face
                interiorFaceInterpolate!(sbpface, ifaces[1],
                                         view(u_c,:,:,ifaces[1].elementL),
                                         view(u_c,:,:,ifaces[1].elementR),
                                         view(uface_c,:,:,ifaces[1].elementL),
                                         view(uface_c,:,:,ifaces[1].elementR))
                # compute the (complexified) flux
                for j = 1:sbpface.numnodes
                  # u1 = 0.5*(uL[1] + uR[1])
                  # u2 = (2*uL[2]*uR[2]/(uL[2] + uR[2])
                  # flux = [u1*cos(u2); u1*u1*sin(u2)]
                  u1 = 0.5*(uface_c[1,j,1] + uface_c[1,j,2])
                  u2 = ((2.0*uface_c[2,j,1]*uface_c[2,j,2])/
                        (uface_c[2,j,1] + uface_c[2,j,2]))
                  flux_c[1,j] = u1*cos(u2)
                  flux_c[2,j] = u1*u1*sin(u2)
                end
                # get residuals
                fill!(res_c, zero(Complex128))
                interiorFaceIntegrate!(sbpface, ifaces[1], flux_c,
                                       view(res_c,:,:,ifaces[1].elementL),
                                       view(res_c,:,:,ifaces[1].elementR))
                if e == 1
                  dR1du1_cmplx[:,p,:,i] = imag(res_c[:,:,1])/ceps
                  dR2du1_cmplx[:,p,:,i] = imag(res_c[:,:,2])/ceps
                else
                  dR1du2_cmplx[:,p,:,i] = imag(res_c[:,:,1])/ceps
                  dR2du2_cmplx[:,p,:,i] = imag(res_c[:,:,2])/ceps
                end
                u_c[p,i,e] -= complex(0.0, ceps)
              end
            end
          end

          # check for equality
          @fact dR1du1 --> roughly(dR1du1_cmplx, atol=1e-15)
          @fact dR1du2 --> roughly(dR1du2_cmplx, atol=1e-15)
          @fact dR2du1 --> roughly(dR2du1_cmplx, atol=1e-15)
          @fact dR2du2 --> roughly(dR2du2_cmplx, atol=1e-15)          
        end
      end
    end
  end
  
end
     
