@testset "Testing SummationByParts Module (Jacobian of face integration methods)..." begin

  for TSBP = ((getLineSegSBPLobbato,getLineSegFace,1),
              (getLineSegSBPLegendre,getLineSegFace,1),
              (getTriSBPGamma,TriFace{Float64},2),
              (getTriSBPOmega,TriFace{Float64},2),
              (getTriSBPDiagE,getTriFaceForDiagE,2),
              (getTetSBPGamma,TetFace{Float64},3),
              (getTetSBPOmega,TetFace{Float64},3),
              (getTetSBPDiagE,getTetFaceForDiagE,3))
    @eval begin
      @testset "Testing boundaryFaceIntegrate_jac! ($(string($TSBP[1])) scalar field method)" begin
        for p = 1:4
          # evaluate residual based on randomly selected face and random state
          # sbp = ($TSBP[1])(degree=p)
          # sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          if $TSBP == (getTetSBPDiagE,getTetFaceForDiagE,3)
            sbp = ($TSBP[1])(degree=p, faceopertype=:Omega)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
          else
            sbp = ($TSBP[1])(degree=p)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          end
          face = 0
          if size(sbp.Q,3) == 1
            face = rand(1:2)
          elseif size(sbp.Q,3) == 2
            face = rand(1:3)
          else
            face = rand(1:4)
          end
          u = randn(sbp.numnodes)
          
          # get the Jacobian directly
          uface = zeros(sbpface.numnodes)
          boundaryFaceInterpolate!(sbpface, face, u, uface)
          flux_jac = zeros(Float64, (sbpface.numnodes))
          for i = 1:sbpface.numnodes
            # flux = 2*u*cos(u*u)
            flux_jac[i] = -2.0*uface[i]*sin(uface[i]*uface[i])
          end
          dRdu = zeros(Float64, (sbp.numnodes, sbp.numnodes))
          boundaryFaceIntegrate_jac!(sbpface, face, flux_jac, dRdu)

          # get the Jacobian using complex step
          u_c = complex(u)
          uface_c = complex(uface)
          flux_c = zeros(ComplexF64, (sbpface.numnodes))
          res_c = zeros(ComplexF64, size(u_c))
          dRdu_cmplx = zeros(size(dRdu))
          ceps = 1e-60
          for i = 1:sbp.numnodes
            u_c[i] += complex(0.0, ceps)
            # interpolate to face
            boundaryFaceInterpolate!(sbpface, face, u_c, uface_c)
            # compute the (complexified) flux
            for j = 1:sbpface.numnodes
              flux_c[j] = cos(uface_c[j]*uface_c[j])
            end
            # get residuals
            fill!(res_c, zero(ComplexF64))
            boundaryFaceIntegrate!(sbpface, face, flux_c, res_c)
            dRdu_cmplx[:,i] = imag(res_c[:,1])/ceps
            u_c[i] -= complex(0.0, ceps)
          end

          # check for equality
          @test ≈(dRdu, dRdu_cmplx, atol=1e-13)
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
      @testset "Testing boundaryFaceIntegrate_jac! ($(string($TSBP[1])) vector field method)" begin
        for p = 1:4
          # evaluate residual based on randomly selected face and random state
          # sbp = ($TSBP[1])(degree=p)
          # sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          if $TSBP == (getTetSBPDiagE,getTetFaceForDiagE,3)
            sbp = ($TSBP[1])(degree=p, faceopertype=:Omega)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
          else
            sbp = ($TSBP[1])(degree=p)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          end
          face = 0
          if size(sbp.Q,3) == 1
            face = rand(1:2)
          elseif size(sbp.Q,3) == 2
            face = rand(1:3)
          else
            face = rand(1:4)
          end
          u = randn(2,sbp.numnodes)
          
          # get the Jacobian directly
          uface = zeros(2, sbpface.numnodes)
          boundaryFaceInterpolate!(sbpface, face, u, uface)
          flux_jac = zeros(Float64, (2, 2, sbpface.numnodes))
          for i = 1:sbpface.numnodes
            # flux = [u1*cos(u2); u1*u1*sin(u2)]
            flux_jac[1,1,i] = cos(uface[2,i])
            flux_jac[1,2,i] = -uface[1,i]*sin(uface[2,i])
            flux_jac[2,1,i] = 2.0*uface[1,i]*sin(uface[2,i])
            flux_jac[2,2,i] = uface[1,i]*uface[1,i]*cos(uface[2,i])
          end
          dRdu = zeros(Float64, (2, 2, sbp.numnodes, sbp.numnodes))
          boundaryFaceIntegrate_jac!(sbpface, face, flux_jac, dRdu)

          # get the Jacobian using complex step
          u_c = complex(u)
          uface_c = complex(uface)
          flux_c = zeros(ComplexF64, (2, sbpface.numnodes))
          res_c = zeros(ComplexF64, size(u_c))
          dRdu_cmplx = zeros(size(dRdu))
          ceps = 1e-60
          for i = 1:sbp.numnodes
            for p = 1:2
              u_c[p,i] += complex(0.0, ceps)
              # interpolate to face
              boundaryFaceInterpolate!(sbpface, face, u_c, uface_c)
              # compute the (complexified) flux
              for j = 1:sbpface.numnodes
                # flux = [u1*cos(u2); u1*u1*sin(u2)]
                flux_c[1,j] = uface_c[1,j]*cos(uface_c[2,j])
                flux_c[2,j] = uface_c[1,j]*uface_c[1,j]*sin(uface_c[2,j])
              end
              # get residuals
              fill!(res_c, zero(ComplexF64))
              boundaryFaceIntegrate!(sbpface, face, flux_c, res_c)
              dRdu_cmplx[:,p,:,i] = imag(res_c[:,:])/ceps
              u_c[p,i] -= complex(0.0, ceps)
            end
          end

          # check for equality
          @test ≈(dRdu, dRdu_cmplx, atol=1e-12)
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
      @testset "Testing interiorFaceIntegrate_jac! ($(string($TSBP[1])) scalar field method)" begin
        for p = 1:4
          # create a two element mesh with random orientation
          # sbp = ($TSBP[1])(degree=p)
          # sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          if $TSBP == (getTetSBPDiagE,getTetFaceForDiagE,3)
            sbp = ($TSBP[1])(degree=p, faceopertype=:Omega)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
          else
            sbp = ($TSBP[1])(degree=p)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          end
          ifaces = Array{Interface}(undef, 1)
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
          dR1du2 = zeros(size(dR1du1))
          dR2du1 = zeros(size(dR1du1))
          dR2du2 = zeros(size(dR1du1))
          interiorFaceIntegrate_jac!(sbpface, ifaces[1],
                                     flux_jac1, flux_jac2,
                                     dR1du1, dR1du2, dR2du1, dR2du2)

          # get the Jacobians using complex step
          u_c = complex(u)
          uface_c = complex(uface)
          flux_c = zeros(ComplexF64, (sbpface.numnodes))
          res_c = zeros(ComplexF64, size(u_c))
          dR1du1_cmplx = zeros(size(dR1du1))
          dR1du2_cmplx = zeros(size(dR1du1))
          dR2du1_cmplx = zeros(size(dR1du1))
          dR2du2_cmplx = zeros(size(dR1du1))
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
              fill!(res_c, zero(ComplexF64))
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
          @test ≈(dR1du1, dR1du1_cmplx, atol=1e-13)
          @test ≈(dR1du2, dR1du2_cmplx, atol=1e-13)
          @test ≈(dR2du1, dR2du1_cmplx, atol=1e-13)
          @test ≈(dR2du2, dR2du2_cmplx, atol=1e-13)          
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
      @testset "Testing interiorFaceIntegrate_jac! ($(string($TSBP[1])) vector field method)" begin
        for p = 1:4
          # create a two element mesh with random orientation
          # sbp = ($TSBP[1])(degree=p)
          # sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          if $TSBP == (getTetSBPDiagE,getTetFaceForDiagE,3)
            sbp = ($TSBP[1])(degree=p, faceopertype=:Omega)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx, faceopertype=:Omega)
          else
            sbp = ($TSBP[1])(degree=p)
            sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
          end
          ifaces = Array{Interface}(undef, 1)
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
          dR1du2 = zeros(size(dR1du1))
          dR2du1 = zeros(size(dR1du1))
          dR2du2 = zeros(size(dR1du1))
          interiorFaceIntegrate_jac!(sbpface, ifaces[1],
                                     flux_jac1, flux_jac2,
                                     dR1du1, dR1du2, dR2du1, dR2du2)

          # get the Jacobians using complex step
          u_c = complex(u)
          uface_c = complex(uface)
          flux_c = zeros(ComplexF64, (2,sbpface.numnodes))
          res_c = zeros(ComplexF64, size(u_c))
          dR1du1_cmplx = zeros(size(dR1du1))
          dR1du2_cmplx = zeros(size(dR1du1))
          dR2du1_cmplx = zeros(size(dR1du1))
          dR2du2_cmplx = zeros(size(dR1du1))
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
                fill!(res_c, zero(ComplexF64))
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
          @test ≈(dR1du1, dR1du1_cmplx, atol=1e-8)
          @test ≈(dR1du2, dR1du2_cmplx, atol=1e-8)
          @test ≈(dR2du1, dR2du1_cmplx, atol=1e-8)
          @test ≈(dR2du2, dR2du2_cmplx, atol=1e-8)          
        end
      end
    end
  end
  
end
     
