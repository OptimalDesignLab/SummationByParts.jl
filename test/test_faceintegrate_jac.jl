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
      context("Testing boundaryFaceIntegrate_jac! ("string($TSBP[1])" vector field method)") do
        for p = 1:4
          # evaluate residual based on randomly selected face and random state
          sbp = ($TSBP[1])(degree=p)
          sbpface = ($TSBP[2])(p, sbp.cub, sbp.vtx)
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
          dFdu = zeros(Float64, (2, 2, sbp.numnodes, sbp.numnodes))
          dFdu_face = zeros(Float64, (2, 2, sbpface.numnodes, sbp.numnodes))
          flux_jac = zeros(Float64, (2, 2, sbp.numnodes))
          for i = 1:sbp.numnodes
            # nominal flux is ( u[2]*sin(u[1]), u[1]*cos(u[2]) )
            flux_jac[1,1,i] = cos(u[1,i])*u[2,i]
            flux_jac[1,2,i] = sin(u[1,i])
            flux_jac[2,1,i] = cos(u[2,i])
            flux_jac[2,2,i] = -u[1,i]*sin(u[2,i])
          end
          # just differentiate flux with respect to di=1
          di = 1
          trans = false
          weakDifferentiateElement_jac!(sbp, di, flux_jac, dFdu, Add(),
                                        trans)
          boundaryFaceInterpolate_jac!(sbpface, face, dFdu, dFdu_face)
          fill!(dFdu, 0.0)
          boundaryFaceIntegrate_jac!(sbpface, face, dFdu_face, dFdu)

          # get the Jacobian using complex step
          u_c = complex(u)
          flux_c = zeros(Complex128, (2, sbp.numnodes))
          dfdx_c = zeros(Complex128, (2, sbp.numnodes))
          dfdx_face_c = zeros(Complex128, (2, sbpface.numnodes))
          dFdu_cmplx = zeros(dFdu)
          ceps = 1e-60
          for i = 1:sbp.numnodes
            for p = 1:2
              u_c[p,i] += complex(0.0, ceps)
              # compute complex flux
              for j = 1:sbp.numnodes
                # nominal flux is ( u[2]*sin(u[1]), u[1]*cos(u[2]) )
                flux_c[1,j] = u_c[2,j]*sin(u_c[1,j])
                flux_c[2,j] = u_c[1,j]*cos(u_c[2,j])
              end
              fill!(dfdx_c, zero(Complex128))
              weakDifferentiateElement!(sbp, di, flux_c, dfdx_c, Add(), trans)
              # interpolate to face
              fill!(dfdx_face_c, zero(Complex128))
              boundaryFaceInterpolate!(sbpface, face, dfdx_c, dfdx_face_c)
              # integrate and apply R^T
              fill!(dfdx_c, zero(Complex128))
              boundaryFaceIntegrate!(sbpface, face, dfdx_face_c, dfdx_c)
              # get residuals
              dFdu_cmplx[:,p,:,i] = imag(dfdx_c[:,:])/ceps
              u_c[p,i] -= complex(0.0, ceps)
            end
          end
          # check for equality
          @fact dFdu --> roughly(dFdu_cmplx, atol=1e-15)

          # test 5D method consistency with 4D (integrate)
          dim = size(sbp.Q, 3)
          facejac_4d = rand(2, 2, sbpface.numnodes, sbp.numnodes, dim)
          facejac_5d = zeros(2, 2, dim, sbpface.numnodes, sbp.numnodes)
          resjac_4d = zeros(2, 2, sbp.numnodes, sbp.numnodes, dim)
          resjac_5d = zeros(2, 2, dim, sbp.numnodes, sbp.numnodes
                            )
          for q=1:sbp.numnodes
            for p=1:sbpface.numnodes
              for d=1:dim
                for j=1:2
                  for i=1:2
                    facejac_5d[i, j, d, p, q] = facejac_4d[i, j, p, q, d]
                  end
                end
              end
            end
          end

          op = SummationByParts.Subtract()
          for include_quad in [true, false]
            for d=1:dim
              facejac_d = sview(facejac_4d, :, :, :, :, d)
              resjac_d = sview(resjac_4d, :, :, :, :, d)
              boundaryFaceIntegrate_jac!(sbpface, face, facejac_d,
                               resjac_d, op, include_quadrature=include_quad)
            end

            boundaryFaceIntegrate_jac!(sbpface, face, facejac_5d,
                            resjac_5d, op, include_quadrature=include_quad)


            for d=1:dim
              @fact maximum(abs.(resjac_5d[:, :, d, :, :] - resjac_4d[:, :, :, :, d])) --> roughly(0.0, atol=1e-13)
            end
          end  # end include_quad

          # test 5D method consistency with 4D (interpolate)
          dim = size(sbp.Q, 3)
          facejac_4d = zeros(2, 2, sbpface.numnodes, sbp.numnodes, dim)
          facejac_5d = zeros(2, 2, dim, sbpface.numnodes, sbp.numnodes)
          resjac_4d = rand(2, 2, sbp.numnodes, sbp.numnodes, dim)
          resjac_5d = zeros(2, 2, dim, sbp.numnodes, sbp.numnodes)
 
          for q=1:sbp.numnodes
            for p=1:sbp.numnodes
              for d=1:dim
                for j=1:2
                  for i=1:2
                    resjac_5d[i, j, d, p, q] = resjac_4d[i, j, p, q, d]
                  end
                end
              end
            end
          end

          for d=1:dim
            facejac_d = sview(facejac_4d, :, :, :, :, d)
            resjac_d = sview(resjac_4d, :, :, :, :, d)
            boundaryFaceInterpolate_jac!(sbpface, face, resjac_d, facejac_d)
          end

          boundaryFaceInterpolate_jac!(sbpface, face, resjac_5d, facejac_5d)

          for d=1:dim
            @fact maximum(abs.(facejac_5d[:, :, d, :, :] - facejac_4d[:, :, :, :, d])) --> roughly(0.0, atol=1e-13)
          end
 

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
          ifaces = Array{Interface}(1)
          if size(sbp.Q,3) == 1
            ifaces[1] = Interface(1,2,rand(1:2),rand(1:2),1)
          elseif size(sbp.Q,3) == 2
            ifaces[1] = Interface(1,2,rand(1:3),rand(1:3),1)
          else
            ifaces[1] = Interface(1,2,rand(1:4),rand(1:4),rand(1:3))
          end
          u = randn(2,sbp.numnodes,2)

          # get Jacobians directly
          dFdu = zeros(Float64, (2, 2, sbp.numnodes, sbp.numnodes, 2))
          dFdu_face = zeros(Float64, (2, 2, sbpface.numnodes, sbp.numnodes, 2))
          flux_jac = zeros(Float64, (2, 2, sbp.numnodes, 2))
          for e = 1:2
            for i = 1:sbp.numnodes
              # nominal flux is ( u[2]*sin(u[1]), u[1]*cos(u[2]) )
              flux_jac[1,1,i,e] = cos(u[1,i,e])*u[2,i,e]
              flux_jac[1,2,i,e] = sin(u[1,i,e])
              flux_jac[2,1,i,e] = cos(u[2,i,e])
              flux_jac[2,2,i,e] = -u[1,i,e]*sin(u[2,i,e])
            end
          end
          # just differentiate flux with respect to di=1
          di = 1
          trans = false
          for e = 1:2
            weakDifferentiateElement_jac!(sbp, di, view(flux_jac,:,:,:,e),
                                          view(dFdu,:,:,:,:,e), Add(), trans)
          end
          interiorFaceInterpolate_jac!(sbpface, ifaces[1],
                                       view(dFdu,:,:,:,:,1),
                                       view(dFdu,:,:,:,:,2),
                                       view(dFdu_face,:,:,:,:,1),
                                       view(dFdu_face,:,:,:,:,2))
          fill!(dFdu, 0.0)
          interiorFaceIntegrate_jac!(sbpface, ifaces[1],
                                     view(dFdu_face,:,:,:,:,1),
                                     view(dFdu_face,:,:,:,:,1), # use left side only
                                     view(dFdu,:,:,:,:,1),
                                     view(dFdu,:,:,:,:,2),
                                     Add(), include_quadrature=true)
          
          # get the Jacobian using complex step
          u_c = complex(u)
          flux_c = zeros(Complex128, (2, sbp.numnodes, 2))
          dfdx_c = zeros(Complex128, (2, sbp.numnodes, 2))
          dfdx_face_c = zeros(Complex128, (2, sbpface.numnodes, 2))
          dFdu_cmplx = zeros(dFdu)
          ceps = 1e-60
          for i = 1:sbp.numnodes
            for p = 1:2
              fill!(dfdx_c, zero(Complex128))
              for e = 1:2
                u_c[p,i,e] += complex(0.0, ceps)
                # compute complex flux
                for j = 1:sbp.numnodes
                  # nominal flux is ( u[2]*sin(u[1]), u[1]*cos(u[2]) )
                  flux_c[1,j,e] = u_c[2,j,e]*sin(u_c[1,j,e])
                  flux_c[2,j,e] = u_c[1,j,e]*cos(u_c[2,j,e])
                end
                weakDifferentiateElement!(sbp, di, view(flux_c,:,:,e),
                                          view(dfdx_c,:,:,e), Add(), trans)
              end
              # interpolate to face
              fill!(dfdx_face_c, zero(Complex128))
              interiorFaceInterpolate!(sbpface, ifaces[1],
                                       view(dfdx_c,:,:,1),
                                       view(dfdx_c,:,:,2),
                                       view(dfdx_face_c,:,:,1),
                                       view(dfdx_face_c,:,:,2))
              # integrate and move back to elements
              fill!(dfdx_c, zero(Complex128))
              interiorFaceIntegrate!(sbpface, ifaces[1],
                                     view(dfdx_face_c,:,:,:1),
                                     view(dfdx_c,:,:,1),
                                     view(dfdx_c,:,:,2),
                                     Add())
              # get Jacobian entries and reset perturbed variable
              for e = 1:2
                dFdu_cmplx[:,p,:,i,e] = imag(dfdx_c[:,:,e])/ceps
                u_c[p,i,e] -= complex(0.0, ceps)
              end
            end
          end
          # check for equality
          @fact dFdu --> roughly(dFdu_cmplx, atol=1e-15)

          # test 5D method consistency with 4D method (faceIntegrate)
          dim = size(sbp.Q, 3)
          facejacL_4d = rand(2, 2, sbpface.numnodes, sbp.numnodes, dim)
          facejacR_4d = rand(2, 2, sbpface.numnodes, sbp.numnodes, dim)
          facejacL_5d = zeros(2, 2, dim, sbpface.numnodes, sbp.numnodes)
          facejacR_5d = zeros(2, 2, dim, sbpface.numnodes, sbp.numnodes)
          resjacL_4d = zeros(2, 2, sbp.numnodes, sbp.numnodes, dim)
          resjacR_4d = zeros(2, 2, sbp.numnodes, sbp.numnodes, dim)
          resjacL_5d = zeros(2, 2, dim, sbp.numnodes, sbp.numnodes)
          resjacR_5d = zeros(2, 2, dim, sbp.numnodes, sbp.numnodes)

          for q=1:sbp.numnodes
            for p=1:sbpface.numnodes
              for d=1:dim
                for j=1:2
                  for i=1:2
                    facejacL_5d[i, j, d, p, q] = facejacL_4d[i, j, p, q, d]
                    facejacR_5d[i, j, d, p, q] = facejacR_4d[i, j, p, q, d]
                  end
                end
              end
            end
          end

          op = SummationByParts.Subtract()
          for include_quad in [true, false]
            for d=1:dim
              facejacL_d = sview(facejacL_4d, :, :, :, :, d)
              facejacR_d = sview(facejacR_4d, :, :, :, :, d)
              resjacL_d = sview(resjacL_4d, :, :, :, :, d)
              resjacR_d = sview(resjacR_4d, :, :, :, :, d)
              interiorFaceIntegrate_jac!(sbpface, ifaces[1], facejacL_d,
                    facejacR_d, resjacL_d, resjacR_d, op, include_quadrature=include_quad)
            end

            interiorFaceIntegrate_jac!(sbpface, ifaces[1], facejacL_5d,
                    facejacR_5d, resjacL_5d, resjacR_5d, op, include_quadrature=include_quad)


            for d=1:dim
              @fact maximum(abs.(resjacL_5d[:, :, d, :, :] - resjacL_4d[:, :, :, :, d])) --> roughly(0.0, atol=1e-13)
              @fact maximum(abs.(resjacR_5d[:, :, d, :, :] - resjacR_4d[:, :, :, :, d])) --> roughly(0.0, atol=1e-13)
            end
          end  # end include_quad

          # test 5D method consistency with 4D method (faceIntegrate)
          facejacL_4d = zeros(2, 2, sbpface.numnodes, sbp.numnodes, dim)
          facejacR_4d = zeros(2, 2, sbpface.numnodes, sbp.numnodes, dim)
          facejacL_5d = zeros(2, 2, dim, sbpface.numnodes, sbp.numnodes)
          facejacR_5d = zeros(2, 2, dim, sbpface.numnodes, sbp.numnodes)
          resjacL_4d = rand(2, 2, sbp.numnodes, sbp.numnodes, dim)
          resjacR_4d = rand(2, 2, sbp.numnodes, sbp.numnodes, dim)
          resjacL_5d = zeros(2, 2, dim, sbp.numnodes, sbp.numnodes)
          resjacR_5d = zeros(2, 2, dim, sbp.numnodes, sbp.numnodes)

          for q=1:sbp.numnodes
            for p=1:sbp.numnodes
              for d=1:dim
                for j=1:2
                  for i=1:2
                    resjacL_5d[i, j, d, p, q] = resjacL_4d[i, j, p, q, d]
                    resjacR_5d[i, j, d, p, q] = resjacR_4d[i, j, p, q, d]
                  end
                end
              end
            end
          end

          for d=1:dim
            facejacL_d = sview(facejacL_4d, :, :, :, :, d)
            facejacR_d = sview(facejacR_4d, :, :, :, :, d)
            resjacL_d = sview(resjacL_4d, :, :, :, :, d)
            resjacR_d = sview(resjacR_4d, :, :, :, :, d)
            interiorFaceInterpolate_jac!(sbpface, ifaces[1], resjacL_d,
                  resjacR_d, facejacL_d, facejacR_d)
          end

            interiorFaceInterpolate_jac!(sbpface, ifaces[1], resjacL_5d,
                    resjacR_5d, facejacL_5d, facejacR_5d)

          for d=1:dim
            @fact maximum(abs.(facejacL_5d[:, :, d, :, :] - facejacL_4d[:, :, :, :, d])) --> roughly(0.0, atol=1e-13)
            @fact maximum(abs.(facejacR_5d[:, :, d, :, :] - facejacR_4d[:, :, :, :, d])) --> roughly(0.0, atol=1e-13)
          end
             
        end
      end
    end
  end
  
end
