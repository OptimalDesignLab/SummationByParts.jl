@testset "Testing SummationByParts Module (directional differentiate methods)..." begin

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      @testset "Testing directionalDifferentiateElement! ($(string($TSBP)) scalar field method)" begin
        # build a single element grid, define u = x+y, and verify that Ddir = 2
        for p = 1:4 
          sbp = ($TSBP)(degree=p)
          x = zeros(Float64, (2,sbp.numnodes))
          vtx = [-1. -1.; -1. 1.; 1. -1.]
          x[:,:] = calcnodes(sbp, vtx)
          u = ones(Float64, (sbp.numnodes))
          u = vec(x[1,:] + x[2,:])
          dir = [1.;1.]
          for i = 1:sbp.numnodes
            Ddir = directionalDifferentiateElement!(sbp, dir, u, i)
            @test ≈(Ddir, 2.0, atol=1e-12)
          end
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      @testset "Testing directionalDifferentiateElement! ($(string($TSBP)) scalar field method)" begin
        # build a single element grid, define u = x+y+z, and verify that Ddir = 3.0
        for p = 1:4 
          sbp = ($TSBP)(degree=p)
          vtx = [-1. -1. -1.; 1. -1. -1.; -1. 1. -1.; -1. -1. 1.]
          x = zeros(Float64, (3,sbp.numnodes))
          x[:,:] = calcnodes(sbp, vtx)
          u = ones(Float64, (sbp.numnodes))
          u = vec(x[1,:] + x[2,:] + x[3,:])
          dir = [1.;1.;1.]
          for i = 1:sbp.numnodes
            Ddir = directionalDifferentiateElement!(sbp, dir, u, i)
            @test ≈(Ddir, 3.0, atol=1e-10)
          end
        end
      end
    end
  end

  for TSBP = (getTriSBPGamma, getTriSBPOmega, getTriSBPDiagE)
    @eval begin
      @testset "Testing directionalDifferentiateElement! ($(string($TSBP)) vector field method)" begin
        # build a single element grid, define u = x+y, and verify that Ddir = 2
        for p = 1:4 
          sbp = ($TSBP)(degree=p)
          x = zeros(Float64, (2,sbp.numnodes))
          vtx = [-1. -1.; -1. 1.; 1. -1.]
          x[:,:] = calcnodes(sbp, vtx)
          u = ones(Float64, (2, sbp.numnodes))
          u[1,:] = vec(x[1,:] + x[2,:])
          u[2,:] = -u[1,:]
          dir = [1.;1.]
          Ddir = zeros(Float64, size(u,1))
          for i = 1:sbp.numnodes
            fill!(Ddir, 0.0)
            directionalDifferentiateElement!(sbp, dir, u, i, Ddir)
            @test ≈(Ddir[1], 2.0, atol=1e-12)
            @test ≈(Ddir[2], -2.0, atol=1e-12)
          end
        end
      end
    end
  end

  for TSBP = (getTetSBPGamma, getTetSBPOmega, getTetSBPDiagE)
    @eval begin
      @testset "Testing directionalDifferentiateElement! ($(string($TSBP)) vector field method)" begin
        # build a single element grid, define u = x+y+z, and verify that Ddir = 3.0
        for p = 1:4 
          sbp = ($TSBP)(degree=p)
          vtx = [-1. -1. -1.; 1. -1. -1.; -1. 1. -1.; -1. -1. 1.]
          x = zeros(Float64, (3,sbp.numnodes))
          x[:,:] = calcnodes(sbp, vtx)
          u = ones(Float64, (2, sbp.numnodes))
          u[1,:] = vec(x[1,:] + x[2,:] + x[3,:])
          u[2,:] = -u[1,:]
          dir = [1.;1.;1.]
          Ddir = zeros(Float64, size(u,1))
          for i = 1:sbp.numnodes
            fill!(Ddir, 0.0)
            directionalDifferentiateElement!(sbp, dir, u, i, Ddir)
            @test ≈(Ddir[1], 3.0, atol=1e-10)
            @test ≈(Ddir[2], -3.0, atol=1e-10)
          end
        end
      end
    end
  end
      
end
