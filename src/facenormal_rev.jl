# This file contains the reverse-mode version of the methods in
# facenormal.jl

for TFACE = (TriFace, TriSparseFace)
  @eval begin
    """
    ### SummationByParts.calcFaceNormals_rev!

    This is the reverse differentiated version of calcFaceNormals!.  See
    facenormal.jl for further details of the primal method.  This function is
    differentiated with respect to the primal version's `xsbp` and `nrm`
    variables.

    **Note**: `xlag` must be provided, but is only needed for the `TetFace` method.

    **Inputs**

    * `sbpface`: an SBP face operator type
    * `mapdegree`: the polynomial degree of the mapping
    * `xref`: Lagrangian nodes in reference space; [coord, Lag node]
    * `xlag`: Lagrangian nodes in physical space; [coord, Lag node, face]
    * `xsbp_bar`: multiplies d(xsbp)/d(xlag) from left; [coord, sbp node, face]
    * `nrm_bar`: multiplies d(nrm)/d(xlag) from the left; [component, sbp node, face]
      
    **In/Outs**
      
    * `xlag_bar`: result of vector Jacobian product; [coord, Lag node, face]
      
    """
    function calcFaceNormals_rev!{Tsbp,Tmsh}(sbpface::($TFACE){Tsbp},
                                             mapdegree::Int,
                                             xref::AbstractMatrix,
                                             xlag::AbstractArray{Tmsh,3},
                                             xlag_bar::AbstractArray{Tmsh,3},
                                             xsbp_bar::AbstractArray{Tmsh,3},
                                             nrm_bar::AbstractArray{Tmsh,3})
      @assert( size(xlag,1) == size(xlag_bar,1) == size(xsbp_bar,1) ==
               size(nrm_bar,1) == 2 )
      @assert( size(xref,1) == 1 )
      @assert( size(xsbp_bar,2) == size(nrm_bar,2) )
      numdof = (mapdegree+1)
      @assert( size(xlag,2) == size(xlag_bar,2) == size(xref,2) == numdof )
      @assert( size(xlag,3) == size(xlag_bar,3) == size(xsbp_bar,3) ==
               size(nrm_bar,3) )
      # find the inverse of the Vandermonde matrix
      V = zeros(Tmsh, (numdof,numdof) )
      for i = 0:mapdegree
        V[:,i+1] = OrthoPoly.jacobipoly(vec(xref[1,:]), 0.0, 0.0, i)
      end
      Vinv = inv(V)
      # get the SBP nodes in reference space in order to find the orthogonal
      # polynomials and their derivatives at these nodes
      x = SymCubatures.calcnodes(sbpface.cub, sbpface.vtx)
      P = zeros(Tmsh, (sbpface.numnodes, numdof))
      dPdξ = zeros(Tmsh, (sbpface.numnodes, numdof))
      for i = 0:mapdegree
        P[:,i+1] = OrthoPoly.jacobipoly(vec(x[1,:]), 0.0, 0.0, i)
        dPdξ[:,i+1] = OrthoPoly.diffjacobipoly(vec(x[1,:]), 0.0, 0.0, i)
      end
      coeff_bar = zeros(Tmsh, (numdof,2))
      # reverse sweep over each face
      for f = 1:size(xlag_bar,3)
        # compute coeff_bar = xsbp_bar*d(xsbp)/d(coeff) + nrm_bar*d(nrm)/d(coeff)
        fill!(coeff_bar, zero(Tmsh))
        for i = 1:numdof      
          for nd = 1:sbpface.numnodes
            # xsbp[1,nd,f] += coeff[i,1]*P[nd,i]
            coeff_bar[i,1] += xsbp_bar[1,nd,f]*P[nd,i]
            # xsbp[2,nd,f] += coeff[i,2]*P[nd,i]
            coeff_bar[i,2] += xsbp_bar[2,nd,f]*P[nd,i]
            # nrm[1,nd,f] += coeff[i,2]*dPdξ[nd,i]
            coeff_bar[i,2] += nrm_bar[1,nd,f]*dPdξ[nd,i]
            # nrm[2,nd,f] -= coeff[i,1]*dPdξ[nd,i]
            coeff_bar[i,1] -= nrm_bar[2,nd,f]*dPdξ[nd,i]
          end
        end
        # compute xlag_bar = coeff_bar*d(xlag)/d(coeff)
        for i = 1:numdof
          #coeff[i,1] = zero(Tmsh)
          #coeff[i,2] = zero(Tmsh)
          for j = 1:numdof
            #coeff[i,1] += Vinv[i,j]*xlag[1,j,f]
            xlag_bar[1,j,f] += Vinv[i,j]*coeff_bar[i,1]
            #coeff[i,2] += Vinv[i,j]*xlag[2,j,f]
            xlag_bar[2,j,f] += Vinv[i,j]*coeff_bar[i,2]
          end
        end
      end
    end
  end
end

for TFACE = (TetFace, TetSparseFace)
  @eval begin
    function calcFaceNormals_rev!{Tsbp,Tmsh}(sbpface::($TFACE){Tsbp},
                                             mapdegree::Int,
                                             xref::AbstractMatrix,
                                             xlag::AbstractArray{Tmsh,3},
                                             xlag_bar::AbstractArray{Tmsh,3},
                                             xsbp_bar::AbstractArray{Tmsh,3},
                                             nrm_bar::AbstractArray{Tmsh,3})
      @assert( size(xlag,1) == size(xlag_bar,1) == size(xsbp_bar,1) ==
               size(nrm_bar,1) == 3 )
      @assert( size(xref,1) == 2 )
      @assert( size(xsbp_bar,2) == size(nrm_bar,2) )
      numdof = binomial(mapdegree+2,2)
      @assert( size(xlag,2) == size(xlag_bar,2) == size(xref,2) == numdof )
      @assert( size(xlag,3) == size(xlag_bar,3) ==size(xsbp_bar,3) ==
               size(nrm_bar,3) )
      # find the inverse of the Vandermonde matrix
      V = zeros(Tmsh, (numdof,numdof) )
      ptr = 1
      for r = 0:mapdegree
        for j = 0:r
          i = r-j
          V[:,ptr] = OrthoPoly.proriolpoly(vec(xref[1,:]), vec(xref[2,:]), i, j)
          ptr += 1
        end
      end
      Vinv = inv(V)
      # get the SBP nodes in reference space in order to find the orthogonal
      # polynomials and their derivatives at these nodes
      x = SymCubatures.calcnodes(sbpface.cub, sbpface.vtx)
      P = zeros(Tmsh, (sbpface.numnodes, numdof))
      dPdξ = zeros(Tmsh, (sbpface.numnodes, numdof))
      dPdη = zeros(Tmsh, (sbpface.numnodes, numdof))
      ptr = 1
      for r = 0:mapdegree
        for j = 0:r
          i = r-j
          P[:,ptr] = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
          dPdξ[:,ptr], dPdη[:,ptr] =
            OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
          ptr += 1
        end
      end
      coeff = zeros(Tmsh, (numdof,3))
      dxdξ = zeros(Tmsh, (3,2,sbpface.numnodes))
      coeff_bar = zeros(Tmsh, (numdof,3))
      dxdξ_bar = zeros(Tmsh, (3,2,sbpface.numnodes))
      # loop over each face...
      for f = 1:size(xlag_bar,3)
        # unlike faces in 2D, the normals of the faces in 3D depend nonlinearly on
        # the tangent vector dxdξ, so this vector must be recomputed first
        
        # find the coefficents of the polynomial mapping using xlag and Vinv
        for di = 1:3
          for i = 1:numdof
            coeff[i,di] = zero(Tmsh)
            for j = 1:numdof
              coeff[i,di] += Vinv[i,j]*xlag[di,j,f]
            end
          end
        end
        # compute the tangent vectors at the SBP nodes
        fill!(dxdξ, zero(Tmsh))
        for i = 1:numdof
          for di = 1:3
            for nd = 1:sbpface.numnodes
              dxdξ[di,1,nd] += coeff[i,di]*dPdξ[nd,i]
              dxdξ[di,2,nd] += coeff[i,di]*dPdη[nd,i]
            end
          end
        end
        # start reverse sweep
        fill!(dxdξ_bar, zero(Tmsh))
        # compute dxdxi_bar = nrm_bar*d(nrm)/d(dxdxi)
        for di = 1:3
          it1 = mod(di,3)+1
          it2 = mod(di+1,3)+1
          for i = 1:sbpface.numnodes
            # nrm[di,i,f] = dxdξ[it1,1,i]*dxdξ[it2,2,i] - dxdξ[it2,1,i]*dxdξ[it1,2,i]
            dxdξ_bar[it1,1,i] += nrm_bar[di,i,f]*dxdξ[it2,2,i]
            dxdξ_bar[it2,2,i] += nrm_bar[di,i,f]*dxdξ[it1,1,i]
            dxdξ_bar[it2,1,i] -= nrm_bar[di,i,f]*dxdξ[it1,2,i]
            dxdξ_bar[it1,2,i] -= nrm_bar[di,i,f]*dxdξ[it2,1,i]
          end
        end
        # compute coeff_bar = xsbp_bar*d(xsbp)/d(coeff) + dxdξ_bar*d(dxdξ)/d(coeff)
        fill!(coeff_bar, zero(Tmsh))
        for i = 1:numdof
          for di = 1:3
            for nd = 1:sbpface.numnodes
              # xsbp[di,nd,f] += coeff[i,di]*P[nd,i]
              coeff_bar[i,di] += xsbp_bar[di,nd,f]*P[nd,i]
              # dxdξ[di,1,nd] += coeff[i,di]*dPdξ[nd,i]
              coeff_bar[i,di] += dxdξ_bar[di,1,nd]*dPdξ[nd,i]
              # dxdξ[di,2,nd] += coeff[i,di]*dPdη[nd,i]
              coeff_bar[i,di] += dxdξ_bar[di,2,nd]*dPdη[nd,i]
            end
          end
        end
        # compute xlag_bar = coeff_bar*d(coeff)/d(xlag)
        for di = 1:3
          for i = 1:numdof
            for j = 1:numdof
              #coeff[i,di] += Vinv[i,j]*xlag[di,j,f]
              xlag_bar[di,j,f] += Vinv[i,j]*coeff_bar[i,di]
            end
          end
        end    
      end
    end
  end
end

for TFACE = (TriFace, TriSparseFace)
  @eval begin
    """
    ### SummationByParts.facenormal_rev!

    This is the reverse differentiated version of facenormal!.  See
    facenormal.jl for further details of the primal method.  This function is
    differentiated with respect to the primal version's `xsbp` and `nrm` variables.

    **Note**: `xlag` must be provided, but is only needed for the `TetFace` method.

    **Inputs**

    * `sbpface`: an SBP face operator type
    * `mapdegree`: the polynomial degree of the mapping
    * `xref`: Lagrangian nodes in reference space; [coord, Lagrangian node]
    * `xlag`: Lagrangian nodes in physical space; [coord, Lagragnian node]
    * `xsbp_bar`: multiplies d(xsbp)/d(xlag) from left; [coord, sbp node]
    * `nrm_bar`: multiplies d(nrm)/d(xlag) from the left; [component, sbp node]

    **In/Outs**

    * `xlag_bar`: result of vector Jacobian product; [coord, Lagrangian node]

    """
    function facenormal_rev!{Tsbp,Tmsh}(sbpface::($TFACE){Tsbp},
                                        mapdegree::Int,
                                        xref::AbstractArray{Tmsh,2},
                                        xlag::AbstractArray{Tmsh,2},
                                        xlag_bar::AbstractArray{Tmsh,2},
                                        xsbp_bar::AbstractArray{Tmsh,2},
                                        nrm_bar::AbstractArray{Tmsh,2})
      @assert( size(xlag,1) == size(xlag_bar,1) == size(xsbp_bar,1) ==
               size(nrm_bar,1) == 2 )
      @assert( size(xref,1) == 1 )
      @assert( size(xsbp_bar,2) == size(nrm_bar,2) )
      numdof = (mapdegree+1)
      @assert( size(xlag,2) == size(xlag_bar,2) == size(xref,2) == numdof )
      # find the Vandermonde matrix
      V = zeros(Tmsh, (numdof,numdof) )
      for i = 0:mapdegree
        V[:,i+1] = OrthoPoly.jacobipoly(vec(xref[1,:]), 0.0, 0.0, i)
      end
      # compute the mapped SBP nodes and the analytical normal at sbp nodes
      x = SymCubatures.calcnodes(sbpface.cub, sbpface.vtx) # <-- SBP nodes, ref. spc
      # compute coeff_bar = xsbp_bar*d(xsbp)/d(coeff) + nrm_bar*d(nrm)/d(coeff)
      coeff_bar = zeros(Tmsh, (numdof,2))
      for i = 0:mapdegree
        P = OrthoPoly.jacobipoly(vec(x[1,:]), 0.0, 0.0, i)
        dPdξ = OrthoPoly.diffjacobipoly(vec(x[1,:]), 0.0, 0.0, i)
        for nd = 1:sbpface.numnodes
          # xsbp[1,nd] += coeff[i+1,1]*P[nd]
          coeff_bar[i+1,1] += xsbp_bar[1,nd]*P[nd]
          # xsbp[2,nd] += coeff[i+1,2]*P[nd]
          coeff_bar[i+1,2] += xsbp_bar[2,nd]*P[nd]
          # nrm[1,nd] += coeff[i+1,2]*dPdξ[nd]
          coeff_bar[i+1,2] += nrm_bar[1,nd]*dPdξ[nd]
          # nrm[2,nd] -= coeff[i+1,1]*dPdξ[nd]
          coeff_bar[i+1,1] -= nrm_bar[2,nd]*dPdξ[nd]
        end
      end
      xlag_bar[:,:] += ((V.')\coeff_bar).'
    end
  end
end

for TFACE = (TetFace, TetSparseFace)
  @eval begin
    function facenormal_rev!{Tsbp,Tmsh}(sbpface::($TFACE){Tsbp},
                                        mapdegree::Int,
                                        xref::AbstractArray{Tmsh,2},
                                        xlag::AbstractArray{Tmsh,2},
                                        xlag_bar::AbstractArray{Tmsh,2},
                                        xsbp_bar::AbstractArray{Tmsh,2},
                                        nrm_bar::AbstractArray{Tmsh,2})
      @assert( size(xlag,1) == size(xlag_bar,1) == size(xsbp_bar,1) ==
               size(nrm_bar,1) == 3 )
      @assert( size(xref,1) == 2 )
      @assert( size(xsbp_bar,2) == size(nrm_bar,2) )
      numdof = binomial(mapdegree+2,2)
      @assert( size(xlag,2) == size(xlag_bar,2) == size(xref,2) == numdof )
      # find the polynomial mapping using xlag
      V = zeros(Tmsh, (numdof,numdof) )
      ptr = 1
      for r = 0:mapdegree
        for j = 0:r
          i = r-j
          V[:,ptr] = OrthoPoly.proriolpoly(vec(xref[1,:]), vec(xref[2,:]), i, j)
          ptr += 1
        end
      end
      coeff = zeros(Tmsh, (numdof,3))
      coeff = V\(xlag.')
      # compute the mapped SBP nodes and the tangent vectors at sbp nodes
      x = SymCubatures.calcnodes(sbpface.cub, sbpface.vtx)
      # unlike faces in 2D, the normals of the faces in 3D depend nonlinearly on
      # the tangent vector dxdξ, so this vector must be recomputed first
      dxdξ = zeros(Tmsh, (3,2,sbpface.numnodes))
      ptr = 1
      for r = 0:mapdegree
        for j = 0:r
          i = r-j
          P = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
          dPdξ, dPdη = OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
          for di = 1:3
            for nd = 1:sbpface.numnodes
              dxdξ[di,1,nd] += coeff[ptr,di]*dPdξ[nd]
              dxdξ[di,2,nd] += coeff[ptr,di]*dPdη[nd]
            end
          end
          ptr += 1
        end
      end

      # Start the reverse sweep

      # compute dxdξ_bar = nrm_bar*d(nrm)/d(dxdξ)
      dxdξ_bar = zeros(Tmsh, (3,2,sbpface.numnodes))
      for di = 1:3
        it1 = mod(di,3)+1
        it2 = mod(di+1,3)+1
        for i = 1:sbpface.numnodes
          # nrm[di,i,f] = dxdξ[it1,1,i]*dxdξ[it2,2,i] - dxdξ[it2,1,i]*dxdξ[it1,2,i]
          dxdξ_bar[it1,1,i] += nrm_bar[di,i]*dxdξ[it2,2,i]
          dxdξ_bar[it2,2,i] += nrm_bar[di,i]*dxdξ[it1,1,i]
          dxdξ_bar[it2,1,i] -= nrm_bar[di,i]*dxdξ[it1,2,i]
          dxdξ_bar[it1,2,i] -= nrm_bar[di,i]*dxdξ[it2,1,i]
        end
      end
      # compute coeff_bar = xsbp_bar*d(xsbp)/d(coeff) + dxdξ_bar*d(dxdξ)/d(coeff)
      coeff_bar = zeros(Tmsh, (numdof,3))
      ptr = 1
      for r = 0:mapdegree
        for j = 0:r
          i = r-j
          P = OrthoPoly.proriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
          dPdξ, dPdη = OrthoPoly.diffproriolpoly(vec(x[1,:]), vec(x[2,:]), i, j)
          for di = 1:3
            for nd = 1:sbpface.numnodes
              # xsbp[di,nd] += coeff[i,di]*P[nd]
              coeff_bar[ptr,di] += xsbp_bar[di,nd]*P[nd]
              # dxdξ[di,1,nd] += coeff[ptr,di]*dPdξ[nd,i]
              coeff_bar[ptr,di] += dxdξ_bar[di,1,nd]*dPdξ[nd]
              # dxdξ[di,2,nd] += coeff[ptr,di]*dPdη[nd,i]
              coeff_bar[ptr,di] += dxdξ_bar[di,2,nd]*dPdη[nd]
            end
          end
          ptr += 1
        end
      end
      # compute xlag_bar = coeff_bar*d(coeff)/d(xlag)
      xlag_bar[:,:] += ((V.')\coeff_bar).'
    end
  end
end
