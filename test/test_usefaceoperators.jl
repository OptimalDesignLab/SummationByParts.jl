facts("Testing SummationByParts Module (usefaceoperators.jl file)...") do

  context("Testing SummationByParts.boundaryintegrate! (TriSBP, scalar field method)") do
    # build a two element grid, and verify the accuracy of the boundary integration
    function bndryflux(u, dξdx, nrm)
      return u*sum(nrm.'*dξdx)
    end
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      x = zeros(Float64, (2,sbp.numnodes,2))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = calcnodes(sbp, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = calcnodes(sbp, vtx)
      dξdx = zeros(Float64, (2,2,sbp.numnodes,2))
      jac = zeros(Float64, (sbp.numnodes,2))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)

      u = zeros(Float64, (sbp.numnodes,2))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          res = zeros(u)
          flux = zeros(sbp.numfacenodes, size(bndryfaces,1))
          for k = 1:size(bndryfaces,1)
            kB = bndryfaces[k].element
            for l = 1:sbp.numfacenodes
              lB = sbp.facenodes[l, bndryfaces[k].face]
              flux[l,k] = bndryflux(u[lB,kB], dξdx[:,:,lB,kB],
                                    sbp.facenormal[:,bndryfaces[k].face])
            end
          end
          boundaryintegrate!(sbp, bndryfaces, flux, res)
          exact = 0.0
          if i == 0 && j != 0
            exact = 1.
          elseif j == 0 && i != 0
            exact = 1.
          elseif i != 0 && j != 0
            exact = 1/(j+1) + 1/(i+1)
          end
          #println("i,j = ",i,",",j,": i+j = ",i+j)
          @fact sum(res) --> roughly(exact, atol=1e-13)
        end
      end
    end
  end

  context("Testing SummationByParts.boundaryintegrate! (TetSBP, scalar field method)") do
    # build a two element grid, and verify the accuracy of the boundary integration
    function bndryflux(u, dξdx, nrm)
      return u*sum(nrm.'*dξdx)
    end
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      x = zeros(Float64, (3,sbp.numnodes,4))
      vtx = Float64[0 0 0; 1 0 0; 0 1 0; 0 0 1]
      x[:,:,1] = calcnodes(sbp, vtx)
      vtx = Float64[1 0 1; 0 0 1; 1 1 1; 1 0 0]
      x[:,:,2] = calcnodes(sbp, vtx)
      vtx = Float64[1 1 0; 0 1 0; 1 0 0; 1 1 1]
      x[:,:,3] = calcnodes(sbp, vtx)
      vtx = Float64[0 1 1; 1 1 1; 0 0 1; 0 1 0]
      x[:,:,4] = calcnodes(sbp, vtx)

      dξdx = zeros(Float64, (3,3,sbp.numnodes,4))
      jac = zeros(Float64, (sbp.numnodes,4))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, 12)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(1,4)
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,3)
      bndryfaces[6] = Boundary(2,4)
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,3)
      bndryfaces[9] = Boundary(3,4)
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,3)
      bndryfaces[12] = Boundary(4,4)

      u = zeros(Float64, (sbp.numnodes,4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            res = zeros(u)
            flux = zeros(sbp.numfacenodes, size(bndryfaces,1))
            for bindex = 1:size(bndryfaces,1)
              kB = bndryfaces[bindex].element
              for l = 1:sbp.numfacenodes
                lB = sbp.facenodes[l, bndryfaces[bindex].face]
                flux[l,bindex] = bndryflux(u[lB,kB], dξdx[:,:,lB,kB],
                                           sbp.facenormal[:,bndryfaces[bindex].face])
              end
            end
            boundaryintegrate!(sbp, bndryfaces, flux, res)
            exact = 0.0
            if !(i == j == k == 0)
              i != 0 ? exact += 1/((j+1)*(k+1)) : nothing
              j != 0 ? exact += 1/((i+1)*(k+1)) : nothing
              k != 0 ? exact += 1/((i+1)*(j+1)) : nothing
            end
            #println("i,j,k = ",i,",",j,",",k," : i+j+k = ",i+j+k)
            @fact sum(res) --> roughly(exact, atol=1e-13)
          end
        end
      end
    end
  end

  context("Testing SummationByParts.boundaryintegrate! (TriSBP, vector field method)") do
    # build a two element grid, and verify the accuracy of the boundary integration
    function bndryflux{T}(u::AbstractArray{T,1}, dξdx::AbstractArray{T,2}, 
                          nrm::AbstractArray{T,1}, flux::AbstractArray{T,1})
      tmp = sum(nrm.'*dξdx)
      for field = 1:size(u,1)
        flux[field] = u[field]*tmp
      end
    end

    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      x = zeros(Float64, (2,sbp.numnodes,2))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = calcnodes(sbp, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = calcnodes(sbp, vtx)
      dξdx = zeros(Float64, (2,2,sbp.numnodes,2))
      jac = zeros(Float64, (sbp.numnodes,2))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)

      u = zeros(Float64, (1,sbp.numnodes,2))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          res = zeros(u)
          flux = zeros(1, sbp.numfacenodes, size(bndryfaces,1))
          for k = 1:size(bndryfaces,1)
            kB = bndryfaces[k].element
            for l = 1:sbp.numfacenodes
              lB = sbp.facenodes[l, bndryfaces[k].face]
              bndryflux(u[:,lB,kB], dξdx[:,:,lB,kB],
                        sbp.facenormal[:,bndryfaces[k].face],
                        view(flux,:,l,k))
            end
          end
          boundaryintegrate!(sbp, bndryfaces, flux, res)
          exact = 0.0
          if i == 0 && j != 0
            exact = 1.
          elseif j == 0 && i != 0
            exact = 1.
          elseif i != 0 && j != 0
            exact = 1/(j+1) + 1/(i+1)
          end
          #println("i,j = ",i,",",j,": i+j = ",i+j)
          @fact sum(res) --> roughly(exact, atol=1e-13)
        end
      end
    end
  end

  context("Testing SummationByParts.boundaryintegrate! (TetSBP, vector field method)") do
    # build a two element grid, and verify the accuracy of the boundary integration
    function bndryflux{T}(u::AbstractArray{T,1}, dξdx::AbstractArray{T,2}, 
                          nrm::AbstractArray{T,1}, flux::AbstractArray{T,1})
      tmp = sum(nrm.'*dξdx)
      for field = 1:size(u,1)
        flux[field] = u[field]*tmp
      end
    end
    for p = 1:4
      sbp = TetSBP{Float64}(degree=p)
      x = zeros(Float64, (3,sbp.numnodes,4))
      vtx = Float64[0 0 0; 1 0 0; 0 1 0; 0 0 1]
      x[:,:,1] = calcnodes(sbp, vtx)
      vtx = Float64[1 0 1; 0 0 1; 1 1 1; 1 0 0]
      x[:,:,2] = calcnodes(sbp, vtx)
      vtx = Float64[1 1 0; 0 1 0; 1 0 0; 1 1 1]
      x[:,:,3] = calcnodes(sbp, vtx)
      vtx = Float64[0 1 1; 1 1 1; 0 0 1; 0 1 0]
      x[:,:,4] = calcnodes(sbp, vtx)

      dξdx = zeros(Float64, (3,3,sbp.numnodes,4))
      jac = zeros(Float64, (sbp.numnodes,4))
      mappingjacobian!(sbp, x, dξdx, jac)
      bndryfaces = Array(Boundary, 12)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(1,4)
      bndryfaces[4] = Boundary(2,1)
      bndryfaces[5] = Boundary(2,3)
      bndryfaces[6] = Boundary(2,4)
      bndryfaces[7] = Boundary(3,1)
      bndryfaces[8] = Boundary(3,3)
      bndryfaces[9] = Boundary(3,4)
      bndryfaces[10] = Boundary(4,1)
      bndryfaces[11] = Boundary(4,3)
      bndryfaces[12] = Boundary(4,4)

      u = zeros(Float64, (1,sbp.numnodes,4))
      for d = 0:p
        for k = 0:d
          for j = 0:d-k
            i = d-j-k
            u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j).*(x[3,:,:].^k)
            res = zeros(u)
            flux = zeros(1,sbp.numfacenodes, size(bndryfaces,1))
            for bindex = 1:size(bndryfaces,1)
              kB = bndryfaces[bindex].element
              for l = 1:sbp.numfacenodes
                lB = sbp.facenodes[l, bndryfaces[bindex].face]
                bndryflux(u[:,lB,kB], dξdx[:,:,lB,kB],
                          sbp.facenormal[:,bndryfaces[bindex].face],
                          view(flux,:,l,bindex))
              end
            end
            boundaryintegrate!(sbp, bndryfaces, flux, res)
            exact = 0.0
            if !(i == j == k == 0)
              i != 0 ? exact += 1/((j+1)*(k+1)) : nothing
              j != 0 ? exact += 1/((i+1)*(k+1)) : nothing
              k != 0 ? exact += 1/((i+1)*(j+1)) : nothing
            end
            #println("i,j,k = ",i,",",j,",",k," : i+j+k = ",i+j+k)
            @fact sum(res) --> roughly(exact, atol=1e-13)
          end
        end
      end
    end
  end

  context("Testing SummationByParts.interiorfaceinterpolate! (TriSBP, scalar field method)") do
    # build a two element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p, reorder=false)
      sbpface = TriFace{Float64}(degree=p, faceonly=true)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,1))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      ifaces = Array(Interface, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      u = zeros(Float64, (sbp.numnodes,2))
      uface = zeros(Float64, (sbpface.numnodes, 2, 1))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          interiorfaceinterpolate!(sbpface, ifaces, u, uface)
          # check that interpolation from left and right elements is exact
          @fact uface[:,1,1] --> roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @fact uface[sbpface.nbrperm[:,1],2,1] --> 
          roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
        end
      end
    end
  end

  context("Testing SummationByParts.interiorfaceinterpolate! (TriSBP, vector field method)") do
    # build a two element grid and verify that interiorfaceinterpolate
    # interpolates all polynomials of degree p exactly
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p, reorder=false)
      sbpface = TriFace{Float64}(degree=p, faceonly=true)
      x = zeros(Float64, (2,sbp.numnodes,2))
      xf = zeros(Float64, (2,sbpface.numnodes,1))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = SymCubatures.calcnodes(sbp.cub, vtx)
      xf[:,:,1] = SymCubatures.calcnodes(sbpface.cub, vtx[[2;3],:])
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = SymCubatures.calcnodes(sbp.cub, vtx)
      ifaces = Array(Interface, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      u = zeros(Float64, (2,sbp.numnodes,2))
      uface = zeros(Float64, (2,sbpface.numnodes, 2, 1))
      for d = 0:p
        for j = 0:d
          i = d-j
          u[1,:,:] = (x[1,:,:].^i).*(x[2,:,:].^j)
          u[2,:,:] = 2.0.*(x[1,:,:].^i).*(x[2,:,:].^j)
          interiorfaceinterpolate!(sbpface, ifaces, u, uface)
          # check that interpolation from left and right elements is exact
          @fact vec(uface[1,:,1,1]) --> roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @fact vec(uface[2,:,1,1]) --> roughly(2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @fact vec(uface[1,sbpface.nbrperm[:,1],2,1]) --> 
          roughly(vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
          @fact vec(uface[2,sbpface.nbrperm[:,1],2,1]) --> 
          roughly(2.0.*vec((xf[1,:,:].^i).*(xf[2,:,:].^j)), atol=1e-13)
        end
      end
    end
  end

  context("Testing SummationByParts.interiorfaceintegrate! (TriSBP, scalar field method)") do
    function fluxfunc(uL, uR, dξdxL, dξdxR, jacL, jacR, αL, αR, nrmL, nrmR)
      flux = sum(nrmL.'*dξdxL)
      flux >= zero(flux) ? flux *= uL : flux *= uR
      return flux
    end
    function bndryflux(u, dξdx, nrm)
      return u*sum(nrm.'*dξdx)
    end
    # build a two element grid and verify that interiorfaceintegrate does
    # nothing when given a continuous linear field
    for p = 1:4
      sbp = TriSBP{Float64}(degree=p)
      x = zeros(Float64, (2,sbp.numnodes,2))
      vtx = [0. 0.; 1. 0.; 0. 1.]
      x[:,:,1] = calcnodes(sbp, vtx)
      vtx = [1. 0.; 1. 1.; 0. 1.]
      x[:,:,2] = calcnodes(sbp, vtx)
      dξdx = zeros(Float64, (2,2,sbp.numnodes,2))
      jac = zeros(Float64, (sbp.numnodes,2))
      mappingjacobian!(sbp, x, dξdx, jac)
      α = zeros(dξdx)
      for k = 1:2
        for i = 1:sbp.numnodes
          for di1 = 1:2
            for di2 = 1:2
              α[di1,di2,i,k] = (dξdx[di1,1,i,k].*dξdx[di2,1,i,k] + 
                                dξdx[di1,2,i,k].*dξdx[di2,2,i,k])*jac[i,k]
            end
          end
        end
      end
      ifaces = Array(Interface, 1)
      ifaces[1] = Interface(1,2,2,3,1)
      bndryfaces = Array(Boundary, 4)
      bndryfaces[1] = Boundary(1,1)
      bndryfaces[2] = Boundary(1,3)
      bndryfaces[3] = Boundary(2,1)
      bndryfaces[4] = Boundary(2,2)
      u = zeros(Float64, (sbp.numnodes,2))
      u[:,:] = x[1,:,:] + x[2,:,:]
      res = zeros(u)
      Fξ = zeros(u)
      Fη = zeros(u)
      for k = 1:2
        for i = 1:sbp.numnodes
          Fξ[i,k] = u[i,k]*(dξdx[1,1,i,k] + dξdx[1,2,i,k])
          Fη[i,k] = u[i,k]*(dξdx[2,1,i,k] + dξdx[2,2,i,k])
        end
      end
      weakdifferentiate!(sbp, 1, Fξ, res, trans=true)
      weakdifferentiate!(sbp, 2, Fη, res, trans=true)
      res *= -1.0

      # compute boundary fluxes and integrate over boundary
      bflux = zeros(sbp.numfacenodes,size(bndryfaces,1))
      for findex = 1:size(bndryfaces,1)
        for i = 1:sbp.numfacenodes
          iB = sbp.facenodes[i, bndryfaces[findex].face]
          bflux[i,findex] = 
          bndryflux(u[iB,bndryfaces[findex].element], 
                    view(dξdx,:,:,iB,bndryfaces[findex].element),
                    view(sbp.facenormal,:,bndryfaces[findex].face))
        end
      end
      boundaryintegrate!(sbp, bndryfaces, bflux, res)

      # compute interior fluxes and integrate over interior faces
      flux = zeros(sbp.numfacenodes,size(ifaces,1))
      for findex = 1:size(ifaces,1)
        for i = 1:sbp.numfacenodes
          iL = sbp.facenodes[i, ifaces[findex].faceL]
          iR = sbp.facenodes[sbp.numfacenodes-i+1, ifaces[findex].faceR]
          flux[i,findex] = 
          fluxfunc(u[iL,ifaces[findex].elementL], u[iR,ifaces[findex].elementR],
                   view(dξdx,:,:,iL,ifaces[findex].elementL),
                   view(dξdx,:,:,iR,ifaces[findex].elementR),
                   jac[iL,ifaces[findex].elementL], jac[iR,ifaces[findex].elementR],
                   view(α,:,:,iL,ifaces[findex].elementL),
                   view(α,:,:,iR,ifaces[findex].elementR),
                   view(sbp.facenormal,:,ifaces[findex].faceL),
                   view(sbp.facenormal,:,ifaces[findex].faceR))
        end
      end
      interiorfaceintegrate!(sbp, ifaces, flux, res)
      for k = 1:2
        for i = 1:sbp.numnodes
          res[i,k] /= (sbp.w[i]/jac[i,k])
        end
      end
      @fact res --> roughly(2.0*ones(res), atol=1e-11)
    end
  end 

end