facts("Testing SummationByParts Module (reverse diff of face-normal methods)...") do

  context("Testing SummationByParts.calcFaceNormals_rev! (TriFace method)") do
    # build a curvilinear element, differentiate xsbp and nrm with using complex
    # step, and then compare with reverse mode
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
      function mapping(ξ)
        x = 1 - (1 - 0.5*(ξ[1]+1))^(p+1)
        y = 1 - (1 - 0.5*(ξ[2]+1))^(p+1)
        fac = 1/sqrt(2)
        return [fac*x - fac*y; fac*x + fac*y]
      end
      # set the coordinates of the Lagrangian nodes in reference and physical
      # space
      xref = zeros(1,p+2)
      for i = 0:(p+1)
        xref[1,i+1] = 2*i/(p+1) - 1.0
      end
      xlag = zeros(2,p+2,3)
      for i = 0:(p+1)
        xlag[:,i+1,1] = mapping([xref[1,i+1]; -1.0])
        xlag[:,i+1,2] = mapping([xref[1,p-i+2]; xref[1,i+1]])
        xlag[:,i+1,3] = mapping([-1.0; xref[1,p-i+2]])
      end
      # get the SBP nodes and the normal vector
      xsbp = zeros(2,sbpface.numnodes,3)
      nrm = zeros(2,sbpface.numnodes,3)
      calcFaceNormals!(sbpface, p+1, xref, xlag, xsbp, nrm)

      # set vector that multiplies from the left
      xsbp_bar = rand(size(xsbp))
      nrm_bar = rand(size(nrm))
      xlag_bar_cmplx = zeros(size(xlag))
      
      # differentiate with respect to the Lagrangian nodes using complex step
      xlag_cmplx = complex(xlag, 0.0)
      xsbp_cmplx = complex(xsbp, 0.0)
      nrm_cmplx = complex(nrm, 0.0)
      ceps = 1e-60
      for f = 1:3
        for i = 0:(p+1)
          for di = 1:2
            xlag_cmplx[di,i+1,f] += complex(0.0, ceps)
            calcFaceNormals!(sbpface, p+1, xref, xlag_cmplx, xsbp_cmplx,
                             nrm_cmplx)
            xlag_bar_cmplx[di,i+1,f] = (sum(xsbp_bar.*imag(xsbp_cmplx)) +
                                        sum(nrm_bar.*imag(nrm_cmplx)))./ceps
            xlag_cmplx[di,i+1,f] -= complex(0.0, ceps)
          end
        end
      end
      xlag_bar = zeros(size(xlag))
      calcFaceNormals_rev!(sbpface, p+1, xref, xlag, xlag_bar, xsbp_bar, nrm_bar)
      @fact xlag_bar --> roughly(xlag_bar_cmplx, atol=1e-15)
    end
  end

  context("Testing SummationByParts.calcFaceNormals_rev! (TriSparseFace method)") do
    # build a curvilinear element, differentiate xsbp and nrm with using complex
    # step, and then compare with reverse mode
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      function mapping(ξ)
        x = 1 - (1 - 0.5*(ξ[1]+1))^(p+1)
        y = 1 - (1 - 0.5*(ξ[2]+1))^(p+1)
        fac = 1/sqrt(2)
        return [fac*x - fac*y; fac*x + fac*y]
      end
      # set the coordinates of the Lagrangian nodes in reference and physical
      # space
      xref = zeros(1,p+2)
      for i = 0:(p+1)
        xref[1,i+1] = 2*i/(p+1) - 1.0
      end
      xlag = zeros(2,p+2,3)
      for i = 0:(p+1)
        xlag[:,i+1,1] = mapping([xref[1,i+1]; -1.0])
        xlag[:,i+1,2] = mapping([xref[1,p-i+2]; xref[1,i+1]])
        xlag[:,i+1,3] = mapping([-1.0; xref[1,p-i+2]])
      end
      # get the SBP nodes and the normal vector
      xsbp = zeros(2,sbpface.numnodes,3)
      nrm = zeros(2,sbpface.numnodes,3)
      calcFaceNormals!(sbpface, p+1, xref, xlag, xsbp, nrm)

      # set vector that multiplies from the left
      xsbp_bar = rand(size(xsbp))
      nrm_bar = rand(size(nrm))
      xlag_bar_cmplx = zeros(size(xlag))
      
      # differentiate with respect to the Lagrangian nodes using complex step
      xlag_cmplx = complex(xlag, 0.0)
      xsbp_cmplx = complex(xsbp, 0.0)
      nrm_cmplx = complex(nrm, 0.0)
      ceps = 1e-60
      for f = 1:3
        for i = 0:(p+1)
          for di = 1:2
            xlag_cmplx[di,i+1,f] += complex(0.0, ceps)
            calcFaceNormals!(sbpface, p+1, xref, xlag_cmplx, xsbp_cmplx,
                             nrm_cmplx)
            xlag_bar_cmplx[di,i+1,f] = (sum(xsbp_bar.*imag(xsbp_cmplx)) +
                                        sum(nrm_bar.*imag(nrm_cmplx)))./ceps
            xlag_cmplx[di,i+1,f] -= complex(0.0, ceps)
          end
        end
      end
      xlag_bar = zeros(size(xlag))
      calcFaceNormals_rev!(sbpface, p+1, xref, xlag, xlag_bar, xsbp_bar, nrm_bar)
      @fact xlag_bar --> roughly(xlag_bar_cmplx, atol=1e-15)
    end
  end

  context("Testing SummationByParts.calcFaceNormals_rev! (TetFace method)") do
    # build a curvilinear element, differentiate xsbp and nrm with using complex
    # step, and then compare with reverse mode
    for p = 1:4
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      function mapping(ξ)
        x = 1 - (1 - 0.5*(ξ[1]+1))^(p+1)
        y = 1 - (1 - 0.5*(ξ[2]+1))^(p+1)
        z = 1 - (1 - 0.5*(ξ[3]+1))^(p+1)
        fac = 1/sqrt(2)
        return [fac*x - fac*y; fac*x + fac*y; z]
      end
      # set the coordinates of the Lagrangian nodes in reference and physical
      # space
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
      # get the SBP nodes and the normal vector
      xsbp = zeros(3,sbpface.numnodes,4)
      nrm = zeros(3,sbpface.numnodes,4)
      calcFaceNormals!(sbpface, p+1, xref, xlag, xsbp, nrm)

      # set vector that multiplies from the left
      xsbp_bar = rand(size(xsbp))
      nrm_bar = rand(size(nrm))
      xlag_bar_cmplx = zeros(size(xlag))
      
      # differentiate with respect to the Lagrangian nodes using complex step
      xlag_cmplx = complex(xlag, 0.0)
      xsbp_cmplx = complex(xsbp, 0.0)
      nrm_cmplx = complex(nrm, 0.0)
      ceps = 1e-60
      for f = 1:4
        for i = 1:numdof
          for di = 1:3
            xlag_cmplx[di,i,f] += complex(0.0, ceps)
            calcFaceNormals!(sbpface, p+1, xref, xlag_cmplx, xsbp_cmplx,
                             nrm_cmplx)
            xlag_bar_cmplx[di,i,f] = (sum(xsbp_bar.*imag(xsbp_cmplx)) +
                                        sum(nrm_bar.*imag(nrm_cmplx)))./ceps
            xlag_cmplx[di,i,f] -= complex(0.0, ceps)
          end
        end
      end
      xlag_bar = zeros(size(xlag))
      calcFaceNormals_rev!(sbpface, p+1, xref, xlag, xlag_bar, xsbp_bar, nrm_bar)
      @fact xlag_bar --> roughly(xlag_bar_cmplx, atol=1e-15)
    end
  end

  context("Testing SummationByParts.calcFaceNormals_rev! (TetSparseFace method)") do
    # build a curvilinear element, differentiate xsbp and nrm with using complex
    # step, and then compare with reverse mode
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx)
      function mapping(ξ)
        x = 1 - (1 - 0.5*(ξ[1]+1))^(p+1)
        y = 1 - (1 - 0.5*(ξ[2]+1))^(p+1)
        z = 1 - (1 - 0.5*(ξ[3]+1))^(p+1)
        fac = 1/sqrt(2)
        return [fac*x - fac*y; fac*x + fac*y; z]
      end
      # set the coordinates of the Lagrangian nodes in reference and physical
      # space
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
      # get the SBP nodes and the normal vector
      xsbp = zeros(3,sbpface.numnodes,4)
      nrm = zeros(3,sbpface.numnodes,4)
      calcFaceNormals!(sbpface, p+1, xref, xlag, xsbp, nrm)

      # set vector that multiplies from the left
      xsbp_bar = rand(size(xsbp))
      nrm_bar = rand(size(nrm))
      xlag_bar_cmplx = zeros(size(xlag))
      
      # differentiate with respect to the Lagrangian nodes using complex step
      xlag_cmplx = complex(xlag, 0.0)
      xsbp_cmplx = complex(xsbp, 0.0)
      nrm_cmplx = complex(nrm, 0.0)
      ceps = 1e-60
      for f = 1:4
        for i = 1:numdof
          for di = 1:3
            xlag_cmplx[di,i,f] += complex(0.0, ceps)
            calcFaceNormals!(sbpface, p+1, xref, xlag_cmplx, xsbp_cmplx,
                             nrm_cmplx)
            xlag_bar_cmplx[di,i,f] = (sum(xsbp_bar.*imag(xsbp_cmplx)) +
                                        sum(nrm_bar.*imag(nrm_cmplx)))./ceps
            xlag_cmplx[di,i,f] -= complex(0.0, ceps)
          end
        end
      end
      xlag_bar = zeros(size(xlag))
      calcFaceNormals_rev!(sbpface, p+1, xref, xlag, xlag_bar, xsbp_bar, nrm_bar)
      @fact xlag_bar --> roughly(xlag_bar_cmplx, atol=1e-15)
    end
  end

  context("Testing SummationByParts.facenormal_rev! (TriFace method)") do
    # build a curvilinear element, differentiate xsbp and nrm with using complex
    # step, and then compare with reverse mode
    for p = 1:4
      sbp = getTriSBPGamma(degree=p)
      sbpface = TriFace{Float64}(p, sbp.cub, sbp.vtx)
      function mapping(ξ)
        x = 1 - (1 - 0.5*(ξ[1]+1))^(p+1)
        y = 1 - (1 - 0.5*(ξ[2]+1))^(p+1)
        fac = 1/sqrt(2)
        return [fac*x - fac*y; fac*x + fac*y]
      end
      # set the coordinates of the Lagrangian nodes in reference and physical
      # space
      xref = zeros(1,p+2)
      for i = 0:(p+1)
        xref[1,i+1] = 2*i/(p+1) - 1.0
      end
      xlag = zeros(2,p+2,3)
      for i = 0:(p+1)
        xlag[:,i+1,1] = mapping([xref[1,i+1]; -1.0])
        xlag[:,i+1,2] = mapping([xref[1,p-i+2]; xref[1,i+1]])
        xlag[:,i+1,3] = mapping([-1.0; xref[1,p-i+2]])
      end
      # get the SBP nodes and the normal vector
      xsbp = zeros(2,sbpface.numnodes,3)
      nrm = zeros(2,sbpface.numnodes,3)
      for f = 1:3
        facenormal!(sbpface, p+1, xref, view(xlag,:,:,f),
                    view(xsbp,:,:,f), view(nrm,:,:,f))
      end
      
      # set vector that multiplies from the left
      xsbp_bar = rand(size(xsbp))
      nrm_bar = rand(size(nrm))
      xlag_bar_cmplx = zeros(size(xlag))
      
      # differentiate with respect to the Lagrangian nodes using complex step
      xlag_cmplx = complex(xlag, 0.0)
      xsbp_cmplx = complex(xsbp, 0.0)
      nrm_cmplx = complex(nrm, 0.0)
      ceps = 1e-60
      for f = 1:3
        for i = 0:(p+1)
          for di = 1:2
            xlag_cmplx[di,i+1,f] += complex(0.0, ceps)
            facenormal!(sbpface, p+1, xref, view(xlag_cmplx,:,:,f),
                        view(xsbp_cmplx,:,:,f), view(nrm_cmplx,:,:,f))
            xlag_bar_cmplx[di,i+1,f] = (sum(xsbp_bar[:,:,f].*
                                            imag(xsbp_cmplx[:,:,f])) +
                                        sum(nrm_bar[:,:,f].*
                                            imag(nrm_cmplx[:,:,f])))./ceps
            xlag_cmplx[di,i+1,f] -= complex(0.0, ceps)
          end
        end
      end
      xlag_bar = zeros(size(xlag))
      for f = 1:3
        facenormal_rev!(sbpface, p+1, xref, view(xlag,:,:,f),
                        view(xlag_bar,:,:,f), view(xsbp_bar,:,:,f),
                        view(nrm_bar,:,:,f))
      end
      @fact xlag_bar --> roughly(xlag_bar_cmplx, atol=1e-15)
    end
  end

  context("Testing SummationByParts.facenormal_rev! (TriSparseFace method)") do
    # build a curvilinear element, differentiate xsbp and nrm with using complex
    # step, and then compare with reverse mode
    for p = 1:4
      sbp = getTriSBPDiagE(degree=p)
      sbpface = getTriFaceForDiagE(p, sbp.cub, sbp.vtx)
      function mapping(ξ)
        x = 1 - (1 - 0.5*(ξ[1]+1))^(p+1)
        y = 1 - (1 - 0.5*(ξ[2]+1))^(p+1)
        fac = 1/sqrt(2)
        return [fac*x - fac*y; fac*x + fac*y]
      end
      # set the coordinates of the Lagrangian nodes in reference and physical
      # space
      xref = zeros(1,p+2)
      for i = 0:(p+1)
        xref[1,i+1] = 2*i/(p+1) - 1.0
      end
      xlag = zeros(2,p+2,3)
      for i = 0:(p+1)
        xlag[:,i+1,1] = mapping([xref[1,i+1]; -1.0])
        xlag[:,i+1,2] = mapping([xref[1,p-i+2]; xref[1,i+1]])
        xlag[:,i+1,3] = mapping([-1.0; xref[1,p-i+2]])
      end
      # get the SBP nodes and the normal vector
      xsbp = zeros(2,sbpface.numnodes,3)
      nrm = zeros(2,sbpface.numnodes,3)
      for f = 1:3
        facenormal!(sbpface, p+1, xref, view(xlag,:,:,f),
                    view(xsbp,:,:,f), view(nrm,:,:,f))
      end
      
      # set vector that multiplies from the left
      xsbp_bar = rand(size(xsbp))
      nrm_bar = rand(size(nrm))
      xlag_bar_cmplx = zeros(size(xlag))
      
      # differentiate with respect to the Lagrangian nodes using complex step
      xlag_cmplx = complex(xlag, 0.0)
      xsbp_cmplx = complex(xsbp, 0.0)
      nrm_cmplx = complex(nrm, 0.0)
      ceps = 1e-60
      for f = 1:3
        for i = 0:(p+1)
          for di = 1:2
            xlag_cmplx[di,i+1,f] += complex(0.0, ceps)
            facenormal!(sbpface, p+1, xref, view(xlag_cmplx,:,:,f),
                        view(xsbp_cmplx,:,:,f), view(nrm_cmplx,:,:,f))
            xlag_bar_cmplx[di,i+1,f] = (sum(xsbp_bar[:,:,f].*
                                            imag(xsbp_cmplx[:,:,f])) +
                                        sum(nrm_bar[:,:,f].*
                                            imag(nrm_cmplx[:,:,f])))./ceps
            xlag_cmplx[di,i+1,f] -= complex(0.0, ceps)
          end
        end
      end
      xlag_bar = zeros(size(xlag))
      for f = 1:3
        facenormal_rev!(sbpface, p+1, xref, view(xlag,:,:,f),
                        view(xlag_bar,:,:,f), view(xsbp_bar,:,:,f),
                        view(nrm_bar,:,:,f))
      end
      @fact xlag_bar --> roughly(xlag_bar_cmplx, atol=1e-15)
    end
  end

  context("Testing SummationByParts.facenormal_rev! (TetFace method)") do
    # build a curvilinear element, differentiate xsbp and nrm with using complex
    # step, and then compare with reverse mode
    for p = 1:4
      sbp = getTetSBPGamma(degree=p)
      sbpface = TetFace{Float64}(p, sbp.cub, sbp.vtx)
      function mapping(ξ)
        x = 1 - (1 - 0.5*(ξ[1]+1))^(p+1)
        y = 1 - (1 - 0.5*(ξ[2]+1))^(p+1)
        z = 1 - (1 - 0.5*(ξ[3]+1))^(p+1)
        fac = 1/sqrt(2)
        return [fac*x - fac*y; fac*x + fac*y; z]
      end
      # set the coordinates of the Lagrangian nodes in reference and physical
      # space
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
      # get the SBP nodes and the normal vector
      xsbp = zeros(3,sbpface.numnodes,4)
      nrm = zeros(3,sbpface.numnodes,4)
      for f = 1:4
        facenormal!(sbpface, p+1, xref, view(xlag,:,:,f), view(xsbp,:,:,f),
                    view(nrm,:,:,f))
      end
      # set vector that multiplies from the left
      xsbp_bar = rand(size(xsbp))
      nrm_bar = rand(size(nrm))
      xlag_bar_cmplx = zeros(size(xlag))
      
      # differentiate with respect to the Lagrangian nodes using complex step
      xlag_cmplx = complex(xlag, 0.0)
      xsbp_cmplx = complex(xsbp, 0.0)
      nrm_cmplx = complex(nrm, 0.0)
      ceps = 1e-60
      for f = 1:4
        for i = 1:numdof
          for di = 1:3
            xlag_cmplx[di,i,f] += complex(0.0, ceps)
            facenormal!(sbpface, p+1, xref, view(xlag_cmplx,:,:,f),
                        view(xsbp_cmplx,:,:,f), view(nrm_cmplx,:,:,f))
            xlag_bar_cmplx[di,i,f] = (sum(xsbp_bar[:,:,f].*
                                          imag(xsbp_cmplx[:,:,f])) +
                                      sum(nrm_bar[:,:,f].*
                                          imag(nrm_cmplx[:,:,f])))./ceps
            xlag_cmplx[di,i,f] -= complex(0.0, ceps)
          end
        end
      end
      xlag_bar = zeros(size(xlag))
      for f = 1:4
        facenormal_rev!(sbpface, p+1, xref, view(xlag,:,:,f),
                        view(xlag_bar,:,:,f), view(xsbp_bar,:,:,f),
                        view(nrm_bar,:,:,f))
      end
      @fact xlag_bar --> roughly(xlag_bar_cmplx, atol=1e-15)
    end
  end

  context("Testing SummationByParts.facenormal_rev! (TetSparseFace method)") do
    # build a curvilinear element, differentiate xsbp and nrm with using complex
    # step, and then compare with reverse mode
    for p = 1:4
      sbp = getTetSBPDiagE(degree=p)
      sbpface = getTetFaceForDiagE(p, sbp.cub, sbp.vtx)
      function mapping(ξ)
        x = 1 - (1 - 0.5*(ξ[1]+1))^(p+1)
        y = 1 - (1 - 0.5*(ξ[2]+1))^(p+1)
        z = 1 - (1 - 0.5*(ξ[3]+1))^(p+1)
        fac = 1/sqrt(2)
        return [fac*x - fac*y; fac*x + fac*y; z]
      end
      # set the coordinates of the Lagrangian nodes in reference and physical
      # space
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
      # get the SBP nodes and the normal vector
      xsbp = zeros(3,sbpface.numnodes,4)
      nrm = zeros(3,sbpface.numnodes,4)
      for f = 1:4
        facenormal!(sbpface, p+1, xref, view(xlag,:,:,f), view(xsbp,:,:,f),
                    view(nrm,:,:,f))
      end
      # set vector that multiplies from the left
      xsbp_bar = rand(size(xsbp))
      nrm_bar = rand(size(nrm))
      xlag_bar_cmplx = zeros(size(xlag))
      
      # differentiate with respect to the Lagrangian nodes using complex step
      xlag_cmplx = complex(xlag, 0.0)
      xsbp_cmplx = complex(xsbp, 0.0)
      nrm_cmplx = complex(nrm, 0.0)
      ceps = 1e-60
      for f = 1:4
        for i = 1:numdof
          for di = 1:3
            xlag_cmplx[di,i,f] += complex(0.0, ceps)
            facenormal!(sbpface, p+1, xref, view(xlag_cmplx,:,:,f),
                        view(xsbp_cmplx,:,:,f), view(nrm_cmplx,:,:,f))
            xlag_bar_cmplx[di,i,f] = (sum(xsbp_bar[:,:,f].*
                                          imag(xsbp_cmplx[:,:,f])) +
                                      sum(nrm_bar[:,:,f].*
                                          imag(nrm_cmplx[:,:,f])))./ceps
            xlag_cmplx[di,i,f] -= complex(0.0, ceps)
          end
        end
      end
      xlag_bar = zeros(size(xlag))
      for f = 1:4
        facenormal_rev!(sbpface, p+1, xref, view(xlag,:,:,f),
                        view(xlag_bar,:,:,f), view(xsbp_bar,:,:,f),
                        view(nrm_bar,:,:,f))
      end
      @fact xlag_bar --> roughly(xlag_bar_cmplx, atol=1e-15)
    end
  end
  
end
