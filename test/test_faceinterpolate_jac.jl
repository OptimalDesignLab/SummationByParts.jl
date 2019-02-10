facts("Testing SummationByParts Module (Jacobian of face interpolation methods)...") do

  for TSBP = ((getLineSegSBPLobbato,getLineSegFace,1),
              (getLineSegSBPLegendre,getLineSegFace,1),
              (getTriSBPGamma,TriFace{Float64},2),
              (getTriSBPOmega,TriFace{Float64},2),
              (getTriSBPDiagE,getTriFaceForDiagE,2),
              (getTetSBPGamma,TetFace{Float64},3),
              (getTetSBPOmega,TetFace{Float64},3),
              (getTetSBPDiagE,getTetFaceForDiagE,3))
    @eval begin
      context("Testing boundaryFaceInterpolate_jac! ("string($TSBP[1])" vector field method)") do
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

          # get the Jacobian using complex step
          u_c = complex(u)
          flux_c = zeros(Complex128, (2, sbp.numnodes))
          dfdx_c = zeros(Complex128, (2, sbp.numnodes))
          dfdx_face_c = zeros(Complex128, (2, sbpface.numnodes))
          dFdu_face_cmplx = zeros(dFdu_face)
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
              # get residuals
              dFdu_face_cmplx[:,p,:,i] = imag(dfdx_face_c[:,:])/ceps
              u_c[p,i] -= complex(0.0, ceps)
            end
          end
          # check for equality
          @fact dFdu_face --> roughly(dFdu_face_cmplx, atol=1e-15)
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
      context("Testing interiorFaceInterpolate_jac! ("string($TSBP[1])" vector field method)") do
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

          # get the Jacobian using complex step
          u_c = complex(u)
          flux_c = zeros(Complex128, (2, sbp.numnodes, 2))
          dfdx_c = zeros(Complex128, (2, sbp.numnodes, 2))
          dfdx_face_c = zeros(Complex128, (2, sbpface.numnodes, 2))
          dFdu_face_cmplx = zeros(dFdu_face)
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
              # get Jacobian entries and reset perturbed variable
              for e = 1:2
                dFdu_face_cmplx[:,p,:,i,e] = imag(dfdx_face_c[:,:,e])/ceps
                u_c[p,i,e] -= complex(0.0, ceps)
              end
            end
          end
          # check for equality
          @fact dFdu_face --> roughly(dFdu_face_cmplx, atol=1e-15)          
          
        end
      end
    end
  end
  
end
