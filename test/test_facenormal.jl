facts("Testing SummationByParts Module (face-normal methods)...") do

  context("Testing SummationByParts.calcFaceNormals! (TriFace method)") do
    # build a curvilinear element, and verify that the geometric conservation
    # law holds
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
      divfree = zeros(2)
      for f = 1:3
        divfree[1] += dot(vec(nrm[1,:,f]),sbpface.wface)
        divfree[2] += dot(vec(nrm[2,:,f]),sbpface.wface)
      end
      @fact divfree[1] --> roughly(0.0, atol=1e-15)
      @fact divfree[2] --> roughly(0.0, atol=1e-15)
    end
  end

  context("Testing SummationByParts.calcFaceNormals! (TriSparseFace method)") do
    # build a curvilinear element, and verify that the geometric conservation
    # law holds
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
      divfree = zeros(2)
      for f = 1:3
        divfree[1] += dot(vec(nrm[1,:,f]),sbpface.wface)
        divfree[2] += dot(vec(nrm[2,:,f]),sbpface.wface)
      end
      @fact divfree[1] --> roughly(0.0, atol=1e-15)
      @fact divfree[2] --> roughly(0.0, atol=1e-15)
    end
  end

  context("Testing SummationByParts.calcFaceNormals! (TetFace method)") do
    # build a curvilinear element, and verify that the geometric conservation
    # law holds
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
      divfree = zeros(3)
      for f = 1:4
        divfree[1] += dot(vec(nrm[1,:,f]),sbpface.wface)
        divfree[2] += dot(vec(nrm[2,:,f]),sbpface.wface)
        divfree[3] += dot(vec(nrm[3,:,f]),sbpface.wface)
      end
      @fact divfree[1] --> roughly(0.0, atol=5e-15)
      @fact divfree[2] --> roughly(0.0, atol=5e-15)
      @fact divfree[3] --> roughly(0.0, atol=5e-15)
    end
  end

  context("Testing SummationByParts.calcFaceNormals! (TetSparseFace method)") do
    # build a curvilinear element, and verify that the geometric conservation
    # law holds
    for p = 1:2
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
      divfree = zeros(3)
      for f = 1:4
        divfree[1] += dot(vec(nrm[1,:,f]),sbpface.wface)
        divfree[2] += dot(vec(nrm[2,:,f]),sbpface.wface)
        divfree[3] += dot(vec(nrm[3,:,f]),sbpface.wface)
      end
      @fact divfree[1] --> roughly(0.0, atol=5e-15)
      @fact divfree[2] --> roughly(0.0, atol=5e-15)
      @fact divfree[3] --> roughly(0.0, atol=5e-15)
    end
  end

  context("Testing SummationByParts.facenormal! (TriFace method)") do
    # build a curvilinear element, and verify that the geometric conservation
    # law holds
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
      divfree = zeros(2)
      for f = 1:3
        facenormal!(sbpface, p+1, xref, sview(xlag,:,:,f),
                    sview(xsbp,:,:,f), sview(nrm,:,:,f))
        divfree[1] += dot(vec(nrm[1,:,f]),sbpface.wface)
        divfree[2] += dot(vec(nrm[2,:,f]),sbpface.wface)
      end
      @fact divfree[1] --> roughly(0.0, atol=1e-15)
      @fact divfree[2] --> roughly(0.0, atol=1e-15)
    end
  end

  context("Testing SummationByParts.facenormal! (TriSparseFace method)") do
    # build a curvilinear element, and verify that the geometric conservation
    # law holds
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
      divfree = zeros(2)
      for f = 1:3
        facenormal!(sbpface, p+1, xref, sview(xlag,:,:,f),
                    sview(xsbp,:,:,f), sview(nrm,:,:,f))
        divfree[1] += dot(vec(nrm[1,:,f]),sbpface.wface)
        divfree[2] += dot(vec(nrm[2,:,f]),sbpface.wface)
      end
      @fact divfree[1] --> roughly(0.0, atol=1e-15)
      @fact divfree[2] --> roughly(0.0, atol=1e-15)
    end
  end

  context("Testing SummationByParts.facenormal! (TetFace method)") do
    # build a curvilinear element, and verify that the geometric conservation
    # law holds
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
      divfree = zeros(3)
      for f = 1:4
        facenormal!(sbpface, p+1, xref, sview(xlag,:,:,f), 
                    sview(xsbp,:,:,f), sview(nrm,:,:,f))
        divfree[1] += dot(vec(nrm[1,:,f]),sbpface.wface)
        divfree[2] += dot(vec(nrm[2,:,f]),sbpface.wface)
        divfree[3] += dot(vec(nrm[3,:,f]),sbpface.wface)
      end
      @fact divfree[1] --> roughly(0.0, atol=1e-15)
      @fact divfree[2] --> roughly(0.0, atol=1e-15)
      @fact divfree[3] --> roughly(0.0, atol=1e-15)
    end
  end

  context("Testing SummationByParts.facenormal! (TetSparseFace method)") do
    # build a curvilinear element, and verify that the geometric conservation
    # law holds
    for p = 1:2
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
      divfree = zeros(3)
      for f = 1:4
        facenormal!(sbpface, p+1, xref, sview(xlag,:,:,f), 
                    sview(xsbp,:,:,f), sview(nrm,:,:,f))
        divfree[1] += dot(vec(nrm[1,:,f]),sbpface.wface)
        divfree[2] += dot(vec(nrm[2,:,f]),sbpface.wface)
        divfree[3] += dot(vec(nrm[3,:,f]),sbpface.wface)
      end
      @fact divfree[1] --> roughly(0.0, atol=1e-15)
      @fact divfree[2] --> roughly(0.0, atol=1e-15)
      @fact divfree[3] --> roughly(0.0, atol=1e-15)
    end
  end

end
