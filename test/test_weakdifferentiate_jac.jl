@testset "Testing SummationByParts Module (Jacobian of weak differentiate methods)..." begin

 for TSBP = ((getLineSegSBPLobbato,1), (getLineSegSBPLegendre,1),
             (getTriSBPGamma,2), (getTriSBPOmega,2), (getTriSBPDiagE,2),
             (getTetSBPGamma,3), (getTetSBPOmega,3), (getTetSBPDiagE,3))
    @eval begin
      @testset "Testing weakDifferentiateElement_jac! ($(string($TSBP[1])) scalar field method)" begin
        for p = 1:4
          for trans in [true, false]
            sbp = ($TSBP[1])(degree=p)
            u = rand(Float64, (sbp.numnodes))

            # get the Jacobian directly
            dRdu = zeros(Float64, (sbp.numnodes, sbp.numnodes))
            flux_jac = zeros(Float64, (sbp.numnodes))
            for i = 1:sbp.numnodes
              flux_jac[i] = cos(u[i])
            end
            for di = 1:size(sbp.Q,3)
              weakDifferentiateElement_jac!(sbp, di, flux_jac, dRdu, Add(), trans)
            end
            
            # get the Jacobian using complex step
            dRdu_cmplx = zeros(ComplexF64, size(dRdu))
            u_c = complex(u)
            flux_c = zeros(ComplexF64, size(u_c))
            res_c = zeros(ComplexF64, size(u_c))
            ceps = 1e-60
            for i = 1:sbp.numnodes
              u_c[i] += complex(0.0, ceps)
              flux_c = sin.(u_c)
              fill!(res_c, zero(ComplexF64))
              for di = 1:size(sbp.Q,3)
                weakDifferentiateElement!(sbp, di, flux_c, res_c, Add(), trans)
              end
              dRdu_cmplx[:,i] = imag(res_c)./ceps
              u_c[i] -= complex(0.0, ceps)
            end

            # check for equality
            @test ≈(dRdu, dRdu_cmplx, atol=1e-15)
          end          
        end
      end
    end
  end

  for TSBP = ((getLineSegSBPLobbato,1), (getLineSegSBPLegendre,1),
              (getTriSBPGamma,2), (getTriSBPOmega,2), (getTriSBPDiagE,2),
              (getTetSBPGamma,3), (getTetSBPOmega,3), (getTetSBPDiagE,3))
    @eval begin
      @testset "Testing weakDifferentiateElement_jac! ($(string($TSBP[1])) vector field method)" begin
        for p = 1:4
          for trans in [true, false]
            sbp = ($TSBP[1])(degree=p)
            u = rand(Float64, (2,sbp.numnodes))

            # get the Jacobian directly
            dRdu = zeros(Float64, (2, 2, sbp.numnodes, sbp.numnodes))
            flux_jac = zeros(Float64, (2, 2, sbp.numnodes))
            for i = 1:sbp.numnodes
              # nominal flux is ( u[2]*sin(u[1]), u[1]*cos(u[2]) )
              flux_jac[1,1,i] = cos(u[1,i])*u[2,i]
              flux_jac[1,2,i] = sin(u[1,i])
              flux_jac[2,1,i] = cos(u[2,i])
              flux_jac[2,2,i] = -u[1,i]*sin(u[2,i])
            end
            for di = 1:size(sbp.Q,3)
              weakDifferentiateElement_jac!(sbp, di, flux_jac, dRdu, Add(), trans)
            end
            
            # get the Jacobian using complex step
            dRdu_cmplx = zeros(ComplexF64, size(dRdu))
            u_c = complex(u)
            flux_c = zeros(ComplexF64, size(u_c))
            res_c = zeros(ComplexF64, size(u_c))
            ceps = 1e-60
            for i = 1:sbp.numnodes
              for j = 1:2
                u_c[j,i] += complex(0.0, ceps)
                for k = 1:sbp.numnodes
                  flux_c[1,k] = u_c[2,k]*sin(u_c[1,k])
                  flux_c[2,k] = u_c[1,k]*cos(u_c[2,k])
                end
                fill!(res_c, zero(ComplexF64))
                for di = 1:size(sbp.Q,3)
                  weakDifferentiateElement!(sbp, di, flux_c, res_c, Add(), trans)
                end
                dRdu_cmplx[:,j,:,i] = imag(res_c)/ceps
                u_c[j,i] -= complex(0.0, ceps)
              end
            end

            # check for equality
            @test ≈(dRdu, dRdu_cmplx, atol=1e-15)
          end          
        end
      end
    end
  end
  
end
