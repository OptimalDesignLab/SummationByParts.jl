# This file gathers together methods related to computing the Jacobian of
# face-based integration terms.

"""
### SummationByParts.boundaryFaceIntegrate_jac!

This function applies the transposed interpolation operator to the given
face-based Jacobian, `dfdu_face`, to get the element-based Jacobian, `dfdu`.

Note that `dfdu` is **added to**, so it may need to be initialized to zero for
expected behavior.

**Currently only available for vector fields.**  The index range for the arrays
is `dfdu[1:nvar,1:nvar,1:n,1:n]`, and `dfdu_face[1:nvar,1:nvar,1:nf,1:n]`, where
`nvar` is the number of state variables at each node, `n` is the number of 
element nodes, and `nf` is the number of face nodes.

**Inputs**

* `sbpface`: an SBP AbstractFace type
* `face`: the relevant face (index) of the element
* `dfdu_face`: Jacobian at face nodes of some function w.r.t. the element state
* `±`: PlusFunctor to add to fun, MinusFunctor to subract
* `include_quadrature`: apply face quadrature weight if true, do not if false

**In/Outs**

* `dfdu`: Jacobian afer applying transposed interpolation operator

"""
function boundaryFaceIntegrate_jac!(sbpface::DenseFace{Tsbp},
                                    face::Integer,
                                    dfdu_face::AbstractArray{Tjac,4},
                                    dfdu::AbstractArray{Tjac,4},
                                    (±)::UnaryFunctor=Add();
                                    include_quadrature::Bool=true
                                    ) where {Tsbp,Tjac}
  @asserts_enabled begin
    @assert( size(dfdu_face,3) == size(sbpface.interp,2) )
    @assert( size(sbpface.interp,1) <= size(dfdu,3) == size(dfdu,4) == 
             size(dfdu_face,4) )
    @assert( size(dfdu,1) == size(dfdu,2) == size(dfdu_face,1) ==
             size(dfdu_face,2) )
  end
  if include_quadrature
    # loop over the volume variables that we are differentiating w.r.t.
    for k = 1:size(dfdu,4)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        for j = 1:sbpface.stencilsize
          for q = 1:size(dfdu,2)
            for p = 1:size(dfdu,1)
              dfdu[p,q,sbpface.perm[j,face],k] += ±(sbpface.interp[j,i]*
                                                    sbpface.wface[i]*
                                                    dfdu_face[p,q,i,k])
            end
          end
        end
      end
    end
  else
    # do not include quadrature, just apply R^T
    # loop over the volume variables that we are differentiating w.r.t.
    for k = 1:size(dfdu,4)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        for j = 1:sbpface.stencilsize
          for q = 1:size(dfdu,2)
            for p = 1:size(dfdu,1)
              dfdu[p,q,sbpface.perm[j,face],k] += ±(sbpface.interp[j,i]*
                                                    dfdu_face[p,q,i,k])
            end
          end
        end
      end
    end
  end
end

function boundaryFaceIntegrate_jac!(sbpface::SparseFace{Tsbp},
                                    face::Integer,
                                    dfdu_face::AbstractArray{Tjac,4},
                                    dfdu::AbstractArray{Tjac,4},
                                    (±)::UnaryFunctor=Add();
                                    include_quadrature::Bool=true
                                    ) where {Tsbp,Tjac}
  @asserts_enabled begin
    @assert( size(dfdu_face,3) == sbpface.numnodes )
    @assert( size(dfdu,3) == size(dfdu,4) == size(dfdu_face,4) )
    @assert( size(dfdu,1) == size(dfdu,2) == size(dfdu_face,1) ==
             size(dfdu_face,2) )
  end
  if include_quadrature
    # loop over the volume variables that we are differentiating w.r.t.
    for k = 1:size(dfdu,4)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        for q = 1:size(dfdu,2)
          for p = 1:size(dfdu,1)
            dfdu[p,q,sbpface.perm[i,face],k] += ±(sbpface.wface[i]*
                                                  dfdu_face[p,q,i,k])
          end
        end
      end
    end
  else
    # do not include face quadrature weight, just apply R^T
    for k = 1:size(dfdu,4)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        for q = 1:size(dfdu,2)
          for p = 1:size(dfdu,1)
            dfdu[p,q,sbpface.perm[i,face],k] += ±(dfdu_face[p,q,i,k])
          end
        end
      end
    end
  end
end

function boundaryFaceIntegrate_jac!(sbpface::DenseFace{Tsbp},
                                    face::Integer,
                                    dfdu_face::AbstractArray{Tjac,5},
                                    dfdu::AbstractArray{Tjac,5},
                                    (±)::UnaryFunctor=Add();
                                    include_quadrature::Bool=true
                                    ) where {Tsbp,Tjac}
  @asserts_enabled begin
    @assert( size(dfdu_face,4) == size(sbpface.interp,2) )
    @assert( size(sbpface.interp,1) <= size(dfdu,4) == size(dfdu,5) == 
             size(dfdu_face,5) )
    @assert( size(dfdu,3) == size(dfdu_face,3) )
    @assert( size(dfdu,1) == size(dfdu,2) == size(dfdu_face,1) ==
             size(dfdu_face,2) )
  end
  if include_quadrature
    # loop over the volume variables that we are differentiating w.r.t.
    for k = 1:size(dfdu,5)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        for j = 1:sbpface.stencilsize
          for d = 1:size(dfdu,3)
            for q = 1:size(dfdu,2)
              for p = 1:size(dfdu,1)
                dfdu[p,q,d,sbpface.perm[j,face],k] += ±(sbpface.interp[j,i]*
                                                        sbpface.wface[i]*
                                                        dfdu_face[p,q,d,i,k])
              end
            end
          end
        end
      end
    end
  else
    # do not include quadrature, just apply R^T
    # loop over the volume variables that we are differentiating w.r.t.
    for k = 1:size(dfdu,5)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        for j = 1:sbpface.stencilsize
          for d = 1:size(dfdu,3)
            for q = 1:size(dfdu,2)
              for p = 1:size(dfdu,1)
                dfdu[p,q,d,sbpface.perm[j,face],k] += ±(sbpface.interp[j,i]*
                                                        dfdu_face[p,q,d,i,k])
              end
            end
          end
        end
      end
    end
  end
end

function boundaryFaceIntegrate_jac!(sbpface::SparseFace{Tsbp},
                                    face::Integer,
                                    dfdu_face::AbstractArray{Tjac,5},
                                    dfdu::AbstractArray{Tjac,5},
                                    (±)::UnaryFunctor=Add();
                                    include_quadrature::Bool=true
                                    ) where {Tsbp,Tjac}
  @asserts_enabled begin
    @assert( size(dfdu_face,4) == sbpface.numnodes )
    @assert( size(dfdu,4) == size(dfdu,5) == size(dfdu_face,5) )
    @assert( size(dfdu,3) == size(dfdu_face,3) )
    @assert( size(dfdu,1) == size(dfdu,2) == size(dfdu_face,1) ==
             size(dfdu_face,2) )
  end
  if include_quadrature
    # loop over the volume variables that we are differentiating w.r.t.
    for k = 1:size(dfdu,5)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        for d = 1:size(dfdu,3)
          for q = 1:size(dfdu,2)
            for p = 1:size(dfdu,1)
              dfdu[p,q,d,sbpface.perm[i,face],k] += ±(sbpface.wface[i]*
                                                      dfdu_face[p,q,d,i,k])
            end
          end
        end
      end
    end
  else
    # do not include face quadrature weight, just apply R^T
    for k = 1:size(dfdu,5)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        for d = 1:size(dfdu,3)
          for q = 1:size(dfdu,2)
            for p = 1:size(dfdu,1)
              dfdu[p,q,d,sbpface.perm[i,face],k] += ±(dfdu_face[p,q,d,i,k])
            end
          end
        end
      end
    end
  end
end

"""
### SummationByParts.interiorFaceIntegrate_jac!

This function applies the transposed interpolation operator to the given
face-based Jacobians, `dfduL_face` and `dfduR_face`, to get the element-based
Jacobians, `dfduL` and `dfduR`.

Note that `dfduL` and `dfduR` are **added to**, so they may need to be
initialized to zero for expected behavior.

**Currently only available for vector fields.** The index range for the arrays
is `dfdu*[1:nvar,1:nvar,1:n,1:n]`, and `dfdu*_face[1:nvar,1:nvar,1:nf,1:n]`,
where `nvar` is the number of state variables at each node, `n` is the number of
element nodes, and `nf` is the number of face nodes.

**Inputs**

* `sbpface`: an SBP AbstractFace type
* `face`: the relevant face (index) of the element
* `dfduL_face`: Jacobian of function w.r.t. the left state
* `dfduR_face`: Jacobian of function w.r.t. the right state
* `±`: PlusFunctor to add to fun, MinusFunctor to subract
* `include_quadrature`: apply face quadrature weight if true, do not if false

**In/Outs**

* `dfduL`: Jacobian w.r.t. left state after applying transposed interpolation 
* `dfduR`: Jacobian w.r.t. right state after applying transposed interpolation

"""
function interiorFaceIntegrate_jac!(sbpface::DenseFace{Tsbp},
                                    iface::Interface,
                                    dfduL_face::AbstractArray{Tjac,4},
                                    dfduR_face::AbstractArray{Tjac,4},
                                    dfduL::AbstractArray{Tjac,4},
                                    dfduR::AbstractArray{Tjac,4},
                                    (±)::UnaryFunctor=Add();
                                    include_quadrature::Bool=true
                                    ) where {Tsbp,Tjac}
#=
  println("size(interp, 1) = ", size(sbpface.interp, 1))
  println("size(dfduL, 3) = ", size(dfduL, 3))
  println("size(dfduL, 4) = ", size(dfduL, 4))
  println("size(dfduR, 3) = ", size(dfduR, 4))
  println("size(dfduL_face, 4) = ", size(dfdu_face, 4))
  =#
  @asserts_enabled begin
    @assert( size(dfduL_face,3) == size(dfduR_face,3) == size(sbpface.interp,2) )
    @assert( size(sbpface.interp,1) <= size(dfduL,3) == size(dfduL,4) ==
             size(dfduR,3) == size(dfduR,4) == size(dfduL_face,4) ==
             size(dfduR_face,4) )
    @assert( size(dfduL,1) == size(dfduL,2) == size(dfduR,1) == size(dfduR,2) ==
             size(dfduL_face,1) == size(dfduL_face,2) == size(dfduR_face,1) ==
             size(dfduR_face,2) )
  end
  if include_quadrature 
    # loop over the volume variables that we are differentiating w.r.t.
    for k = 1:size(dfduL,4)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,iface.orient]
        # loop over the interpolation stencil
        for j = 1:sbpface.stencilsize
          for q = 1:size(dfduL,2)
            for p = 1:size(dfduL,1)
              dfduL[p,q,sbpface.perm[j,iface.faceL],k] +=
                ±(sbpface.interp[j,i]*sbpface.wface[i]*dfduL_face[p,q,i,k])
              dfduR[p,q,sbpface.perm[j,iface.faceR],k] -=
                ±(sbpface.interp[j,iR]*sbpface.wface[iR]*dfduR_face[p,q,i,k])
            end
          end
        end
      end
    end
  else
    # do not apply the quadrature, just apply R^T
    # loop over the volume variables that we are differentiating w.r.t.
    for k = 1:size(dfduL,4)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,iface.orient]
        # loop over the interpolation stencil
        for j = 1:sbpface.stencilsize
          for q = 1:size(dfduL,2)
            for p = 1:size(dfduL,1)
              dfduL[p,q,sbpface.perm[j,iface.faceL],k] +=
                ±(sbpface.interp[j,i]*dfduL_face[p,q,i,k])
              dfduR[p,q,sbpface.perm[j,iface.faceR],k] -=
                ±(sbpface.interp[j,iR]*dfduR_face[p,q,i,k])
            end
          end
        end
      end
    end
  end
end

function interiorFaceIntegrate_jac!(sbpface::SparseFace{Tsbp},
                                    iface::Interface,
                                    dfduL_face::AbstractArray{Tjac,4},
                                    dfduR_face::AbstractArray{Tjac,4},
                                    dfduL::AbstractArray{Tjac,4},
                                    dfduR::AbstractArray{Tjac,4},
                                    (±)::UnaryFunctor=Add();
                                    include_quadrature::Bool=true
                                    ) where {Tsbp,Tjac}
  @asserts_enabled begin
    @assert( size(dfduL_face,3) == size(dfduR_face,3) == sbpface.numnodes )
    @assert( size(dfduL,3) == size(dfduL,4) == size(dfduR,3) == size(dfduR,4) ==
             size(dfduL_face,4) == size(dfduR_face,4) )
    @assert( size(dfduL,1) == size(dfduL,2) == size(dfduR,1) == size(dfduR,2) ==
             size(dfduL_face,1) == size(dfduL_face,2) == size(dfduR_face,1) ==
             size(dfduR_face,2) )
  end
  if include_quadrature
    # loop over the volume variables that we are differentiating w.r.t.
    for k = 1:size(dfduL,4)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,iface.orient]
        for q = 1:size(dfduL,2)
          for p = 1:size(dfduL,1)
            dfduL[p,q,sbpface.perm[i,iface.faceL],k] +=
              ±(sbpface.wface[i]*dfduL_face[p,q,i,k])
            dfduR[p,q,sbpface.perm[iR,iface.faceR],k] -=
              ±(sbpface.wface[i]*dfduR_face[p,q,i,k])
          end
        end
      end
    end
  else
    # do not include the quadrature weights, just apply R^T
    # loop over the volume variables that we are differentiating w.r.t.
    for k = 1:size(dfduL,4)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,iface.orient]
        for q = 1:size(dfduL,2)
          for p = 1:size(dfduL,1)
            dfduL[p,q,sbpface.perm[i,iface.faceL],k] += ±(dfduL_face[p,q,i,k])
            dfduR[p,q,sbpface.perm[iR,iface.faceR],k] -= ±(dfduR_face[p,q,i,k])
          end
        end
      end
    end
  end
end

function interiorFaceIntegrate_jac!(sbpface::DenseFace{Tsbp},
                                    iface::Interface,
                                    dfduL_face::AbstractArray{Tjac,5},
                                    dfduR_face::AbstractArray{Tjac,5},
                                    dfduL::AbstractArray{Tjac,5},
                                    dfduR::AbstractArray{Tjac,5},
                                    (±)::UnaryFunctor=Add();
                                    include_quadrature::Bool=true
                                    ) where {Tsbp,Tjac}
  @asserts_enabled begin
    @assert( size(dfduL_face,4) == size(dfduR_face,4) == size(sbpface.interp,2) )
    @assert( size(sbpface.interp,1) <= size(dfduL,4) == size(dfduL,5) ==
             size(dfduR,4) == size(dfduR,5) == size(dfduL_face,5) ==
             size(dfduR_face,5) )
    @assert( size(dfduL,3) == size(dfduR,3) == size(dfduL_face,3) ==
             size(dfduR_face,3) )
    @assert( size(dfduL,1) == size(dfduL,2) == size(dfduR,1) == size(dfduR,2) ==
             size(dfduL_face,1) == size(dfduL_face,2) == size(dfduR_face,1) ==
             size(dfduR_face,2) )
  end
  if include_quadrature 
    # loop over the volume variables that we are differentiating w.r.t.
    for k = 1:size(dfduL,5)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,iface.orient]
        # loop over the interpolation stencil
        for j = 1:sbpface.stencilsize
          for d = 1:size(dfduL,3)
            for q = 1:size(dfduL,2)
              for p = 1:size(dfduL,1)
                dfduL[p,q,d,sbpface.perm[j,iface.faceL],k] +=
                  ±(sbpface.interp[j,i]*sbpface.wface[i]*dfduL_face[p,q,d,i,k])
                dfduR[p,q,d,sbpface.perm[j,iface.faceR],k] -=
                  ±(sbpface.interp[j,iR]*sbpface.wface[iR]*dfduR_face[p,q,d,i,k])
              end
            end
          end
        end
      end
    end
  else
    # do not apply the quadrature, just apply R^T
    # loop over the volume variables that we are differentiating w.r.t.
    for k = 1:size(dfduL,5)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,iface.orient]
        # loop over the interpolation stencil
        for j = 1:sbpface.stencilsize
          for d = 1:size(dfduL,3)
            for q = 1:size(dfduL,2)
              for p = 1:size(dfduL,1)
                dfduL[p,q,d,sbpface.perm[j,iface.faceL],k] +=
                  ±(sbpface.interp[j,i]*dfduL_face[p,q,d,i,k])
                dfduR[p,q,d,sbpface.perm[j,iface.faceR],k] -=
                  ±(sbpface.interp[j,iR]*dfduR_face[p,q,d,i,k])
              end
            end
          end
        end
      end
    end
  end
end

function interiorFaceIntegrate_jac!(sbpface::SparseFace{Tsbp},
                                    iface::Interface,
                                    dfduL_face::AbstractArray{Tjac,5},
                                    dfduR_face::AbstractArray{Tjac,5},
                                    dfduL::AbstractArray{Tjac,5},
                                    dfduR::AbstractArray{Tjac,5},
                                    (±)::UnaryFunctor=Add();
                                    include_quadrature::Bool=true
                                    ) where {Tsbp,Tjac}
  @asserts_enabled begin
    @assert( size(dfduL_face,4) == size(dfduR_face,4) == sbpface.numnodes )
    @assert( size(dfduL,4) == size(dfduL,5) == size(dfduR,4) == size(dfduR,5) ==
             size(dfduL_face,5) == size(dfduR_face,5) )
    @assert( size(dfduL,3) == size(dfduR,3) == size(dfduL_face,3) ==
             size(dfduR_face,3) )    
    @assert( size(dfduL,1) == size(dfduL,2) == size(dfduR,1) == size(dfduR,2) ==
             size(dfduL_face,1) == size(dfduL_face,2) == size(dfduR_face,1) ==
             size(dfduR_face,2) )
  end
  if include_quadrature
    # loop over the volume variables that we are differentiating w.r.t.
    for k = 1:size(dfduL,5)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,iface.orient]
        for d = 1:size(dfduL,3)
          for q = 1:size(dfduL,2)
            for p = 1:size(dfduL,1)
              dfduL[p,q,d,sbpface.perm[i,iface.faceL],k] +=
                ±(sbpface.wface[i]*dfduL_face[p,q,d,i,k])
              dfduR[p,q,d,sbpface.perm[iR,iface.faceR],k] -=
                ±(sbpface.wface[i]*dfduR_face[p,q,d,i,k])
            end
          end
        end
      end
    end
  else
    # do not include the quadrature weights, just apply R^T
    # loop over the volume variables that we are differentiating w.r.t.
    for k = 1:size(dfduL,5)
      # loop over the face nodes
      for i = 1:sbpface.numnodes
        iR = sbpface.nbrperm[i,iface.orient]
        for d = 1:size(dfduL,3)
          for q = 1:size(dfduL,2)
            for p = 1:size(dfduL,1)
              dfduL[p,q,d,sbpface.perm[i,iface.faceL],k] += ±(dfduL_face[p,q,d,i,k])
              dfduR[p,q,d,sbpface.perm[iR,iface.faceR],k] -= ±(dfduR_face[p,q,d,i,k])
            end
          end
        end
      end
    end
  end
end


  
