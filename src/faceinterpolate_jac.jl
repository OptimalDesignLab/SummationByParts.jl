# This file gathers together methods related to computing the Jacobian of
# face-based interpolation terms.

# Boundary methods
#------------------------------------------------------------------------------
# 4D methods

"""
### SummationByParts.boundaryFaceInterpolate_jac!

This function applies the interpolation operator to the given element-based
Jacobian, `dfdu`, to get the face-based Jacobian, `dfdu_face`.

Note that `dfdu_face` is **added to**, so it may need to be initialized to zero
for expected behavior. 

**Currently only available for vector fields.**  The index range for the arrays
is `dfdu[1:nvar,1:nvar,1:n,1:n]`, and `dfdu_face[1:nvar,1:nvar,1:nf,1:n]`, where
`nvar` is the number of state variables at each node, `n` is the number of 
element nodes, and `nf` is the number of face nodes. Methods are also
available for `dfdu[1:nbar,1:nvar,1:n]` (`dfdu_face` remains the same size)

**Inputs**

* `sbpface`: an SBP AbstractFace type
* `face`: the relevant face (index) of the element
* `dfdu`: Jacobian of some function w.r.t. the element state

**In/Outs**

* `dfdu_face`: Jacobian of interpolated function w.r.t. the state

"""
function boundaryFaceInterpolate_jac!(sbpface::DenseFace{Tsbp},
                                      face::Integer,
                                      dfdu::AbstractArray{Tsol,4},
                                      dfdu_face::AbstractArray{Tsol,4}
                                      ) where {Tsbp,Tsol}
  @asserts_enabled begin
    @assert( size(dfdu_face,3) == size(sbpface.interp,2) )
    @assert( size(sbpface.interp,1) <= size(dfdu,3) == size(dfdu,4) == 
             size(dfdu_face,4) )
    @assert( size(dfdu,1) == size(dfdu,2) == size(dfdu_face,1) ==
             size(dfdu_face,2) )
  end
  # loop over the volume variables that we are differentiating w.r.t.
  for k = 1:size(dfdu,4)
    # loop over the face nodes
    for i = 1:sbpface.numnodes
      for j = 1:sbpface.stencilsize
        for q = 1:size(dfdu,2)
          for p = 1:size(dfdu,1)
            dfdu_face[p,q,i,k] += (sbpface.interp[j,i]*
                                   dfdu[p,q,sbpface.perm[j,face],k])
          end
        end
      end
    end
  end
end

function boundaryFaceInterpolate_jac!(sbpface::SparseFace{Tsbp},
                                      face::Integer,
                                      dfdu::AbstractArray{Tsol,4},
                                      dfdu_face::AbstractArray{Tsol,4}
                                      ) where {Tsbp,Tsol}
  @asserts_enabled begin
    @assert( size(dfdu_face,3) == sbpface.numnodes )
    @assert( size(dfdu,3) == size(dfdu,4) == size(dfdu_face,4) )
    @assert( size(dfdu,1) == size(dfdu,2) == size(dfdu_face,1) ==
             size(dfdu_face,2) )
  end
  # loop over the volume variables that we are differentiating w.r.t.
  for k = 1:size(dfdu,4)
    # loop over the face nodes
    for i = 1:sbpface.numnodes
      for q = 1:size(dfdu,2)
        for p = 1:size(dfdu,1)
          dfdu_face[p,q,i,k] = dfdu[p,q,sbpface.perm[i,face],k]
        end
      end
    end
  end
end

#------------------------------------------------------------------------------
# 3D methods

function boundaryFaceInterpolate_jac!(sbpface::DenseFace{Tsbp},
                                      face::Integer,
                                      dfdu::AbstractArray{Tsol,3},
                                      dfdu_face::AbstractArray{Tsol,4}
                                      ) where {Tsbp,Tsol}
  @asserts_enabled begin
    @assert( size(dfdu_face,3) == size(sbpface.interp,2) )
    @assert( size(sbpface.interp,1) <= size(dfdu,3)  == size(dfdu_face,4) )
    @assert( size(dfdu,1) == size(dfdu,2) == size(dfdu_face,1) ==
             size(dfdu_face,2) )
  end
  # loop over the volume variables that we are differentiating w.r.t.
    # loop over the face nodes
    for i = 1:sbpface.numnodes
      for j = 1:sbpface.stencilsize
        perm = sbpface.perm[j, face]
        for q = 1:size(dfdu,2)
          for p = 1:size(dfdu,1)
            dfdu_face[p,q,i,perm] += (sbpface.interp[j,i]*
                                   dfdu[p,q,perm])
          end
        end
      end
    end
end


function boundaryFaceInterpolate_jac!(sbpface::SparseFace{Tsbp},
                                      face::Integer,
                                      dfdu::AbstractArray{Tsol,3},
                                      dfdu_face::AbstractArray{Tsol,4}
                                      ) where {Tsbp,Tsol}
  @asserts_enabled begin
    @assert( size(dfdu_face,3) == sbpface.numnodes )
    @assert( size(dfdu,3) == size(dfdu_face,4) )
    @assert( size(dfdu,1) == size(dfdu,2) == size(dfdu_face,1) ==
             size(dfdu_face,2) )
  end
  # loop over the face nodes
  for i = 1:sbpface.numnodes
    perm = sbpface.perm[i, face]
    for q = 1:size(dfdu,2)
      for p = 1:size(dfdu,1)
        dfdu_face[p,q,i,perm] = dfdu[p,q,perm]
      end
    end
  end
end

#------------------------------------------------------------------------------
# 5D methods

function boundaryFaceInterpolate_jac!(sbpface::DenseFace{Tsbp},
                                      face::Integer,
                                      dfdu::AbstractArray{Tsol,5},
                                      dfdu_face::AbstractArray{Tsol,5}
                                      ) where {Tsbp,Tsol}
  @asserts_enabled begin
    @assert( size(dfdu_face,4) == size(sbpface.interp,2) )
    @assert( size(sbpface.interp,1) <= size(dfdu,4) == size(dfdu,5) == 
             size(dfdu_face,5) )
    @assert( size(dfdu,3) == size(dfdu_face,3) )
    @assert( size(dfdu,1) == size(dfdu,2) == size(dfdu_face,1) ==
             size(dfdu_face,2) )
  end
  # loop over the volume variables that we are differentiating w.r.t.
  for k = 1:size(dfdu,5)
    # loop over the face nodes
    for i = 1:sbpface.numnodes
      for j = 1:sbpface.stencilsize
        for d = 1:size(dfdu,3)
          for q = 1:size(dfdu,2)
            for p = 1:size(dfdu,1)
              dfdu_face[p,q,d,i,k] += (sbpface.interp[j,i]*
                                       dfdu[p,q,d,sbpface.perm[j,face],k])
            end
          end
        end
      end
    end
  end
end

function boundaryFaceInterpolate_jac!(sbpface::SparseFace{Tsbp},
                                      face::Integer,
                                      dfdu::AbstractArray{Tsol,5},
                                      dfdu_face::AbstractArray{Tsol,5}
                                      ) where {Tsbp,Tsol}
  @asserts_enabled begin
    @assert( size(dfdu_face,4) == sbpface.numnodes )
    @assert( size(dfdu,4) == size(dfdu,5) == size(dfdu_face,5) )
    @assert( size(dfdu,3) == size(dfdu_face,3) )
    @assert( size(dfdu,1) == size(dfdu,2) == size(dfdu_face,1) ==
             size(dfdu_face,2) )
  end
  # loop over the volume variables that we are differentiating w.r.t.
  for k = 1:size(dfdu,5)
    # loop over the face nodes
    for i = 1:sbpface.numnodes
      for d = 1:size(dfdu,3)
        for q = 1:size(dfdu,2)
          for p = 1:size(dfdu,1)
            dfdu_face[p,q,d,i,k] = dfdu[p,q,d,sbpface.perm[i,face],k]
          end
        end
      end
    end
  end
end

# Interface methods
#------------------------------------------------------------------------------
# 4D methods

"""
### SummationByParts.interiorFaceInterpolate_jac!

This function applies the interpolation operator to the given element-based
Jacobian, `dfdu`, to get the face-based Jacobian, `dfdu_face`.

Note that `dfduL_face` and `dfduR_face` are **added to**, so it may need to be
initialized to zero for expected behavior.

**Currently only available for vector fields.** The index range for the arrays
is `dfdu*[1:nvar,1:nvar,1:n,1:n]`, and `dfdu*_face[1:nvar,1:nvar,1:nf,1:n]`,
where `nvar` is the number of state variables at each node, `n` is the number of
left/right element nodes, and `nf` is the number of face nodes. Methods are also
available for `dfdu[1:nbar,1:nvar,1:n]` (`dfdu_face` remains the same size)


**Inputs**

* `sbpface`: an SBP AbstractFace type
* `face`: the relevant face (index) of the element
* `dfduL`: Jacobian of some function w.r.t. the left element state
* `dfduR`: Jacobian of some function w.r.t. the right element state

**In/Outs**

* `dfduL_face`: Jacobian of interpolated function w.r.t. the left state
* `dfduR_face`: Jacobian of interpolated function w.r.t. the right state

"""
function interiorFaceInterpolate_jac!(sbpface::DenseFace{Tsbp},
                                      iface::Interface,
                                      dfduL::AbstractArray{Tsol,4},
                                      dfduR::AbstractArray{Tsol,4},
                                      dfduL_face::AbstractArray{Tsol,4},
                                      dfduR_face::AbstractArray{Tsol,4}
                                      ) where {Tsbp,Tsol}
  @asserts_enabled begin
    @assert( size(dfduL_face,3) == size(dfduR_face,3) == size(sbpface.interp,2) )
    @assert( size(sbpface.interp,1) <= size(dfduL,3) == size(dfduL,4) ==
             size(dfduR,3) == size(dfduR,4) == size(dfduL_face,4) ==
             size(dfduR_face,4) )
    @assert( size(dfduL,1) == size(dfduL,2) == size(dfduR,1) == size(dfduR,2) ==
             size(dfduL_face,1) == size(dfduL_face,2) == size(dfduR_face,1) ==
             size(dfduR_face,2) )
  end
  # loop over the volume variables that we are differentiating w.r.t.
  for k = 1:size(dfduL,4)
    # loop over the face nodes
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,iface.orient]
      # loop over the interpolation stencil
      for j = 1:sbpface.stencilsize
        for q = 1:size(dfduL,2)
          for p = 1:size(dfduL,1)
            dfduL_face[p,q,i,k] += (sbpface.interp[j,i]*
                                    dfduL[p,q,sbpface.perm[j,iface.faceL],k])
            dfduR_face[p,q,i,k] += (sbpface.interp[j,iR]*
                                    dfduR[p,q,sbpface.perm[j,iface.faceR],k])
          end
        end
      end
    end
  end
end




function interiorFaceInterpolate_jac!(sbpface::SparseFace{Tsbp},
                                      iface::Interface,
                                      dfduL::AbstractArray{Tsol,4},
                                      dfduR::AbstractArray{Tsol,4},
                                      dfduL_face::AbstractArray{Tsol,4},
                                      dfduR_face::AbstractArray{Tsol,4}
                                      ) where {Tsbp,Tsol}
  @asserts_enabled begin
    @assert( size(dfduL_face,3) == size(dfduR_face,3) == sbpface.numnodes )
    @assert( size(dfduL,3) == size(dfduL,4) == size(dfduR,3) == size(dfduR,4) ==
             size(dfduL_face,4) == size(dfduR_face,4) )
    @assert( size(dfduL,1) == size(dfduL,2) == size(dfduR,1) == size(dfduR,2) ==
             size(dfduL_face,1) == size(dfduL_face,2) == size(dfduR_face,1) ==
             size(dfduR_face,2) )
  end
  # loop over the volume variables that we are differentiating w.r.t.
  for k = 1:size(dfduL,4)
    # loop over the face nodes
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,iface.orient]
      for q = 1:size(dfduL,2)
        for p = 1:size(dfduL,1)
          dfduL_face[p,q,i,k] += dfduL[p,q,sbpface.perm[i,iface.faceL],k]
          dfduR_face[p,q,i,k] += dfduR[p,q,sbpface.perm[iR,iface.faceR],k]
        end
      end
    end
  end
end

#------------------------------------------------------------------------------
# 3D methods

function interiorFaceInterpolate_jac!(sbpface::DenseFace{Tsbp},
                                      iface::Interface,
                                      dfduL::AbstractArray{Tsol,3},
                                      dfduR::AbstractArray{Tsol,3},
                                      dfduL_face::AbstractArray{Tsol,4},
                                      dfduR_face::AbstractArray{Tsol,4}
                                      ) where {Tsbp,Tsol}
  @asserts_enabled begin
    @assert( size(dfduL_face,3) == size(dfduR_face,3) == size(sbpface.interp,2) )
    @assert( size(sbpface.interp,1) <= size(dfduL,3)  ==
             size(dfduR,3) == size(dfduL_face,4) == size(dfduR_face,4) )
    @assert( size(dfduL,1) == size(dfduL,2) == size(dfduR,1) == size(dfduR,2) ==
             size(dfduL_face,1) == size(dfduL_face,2) == size(dfduR_face,1) ==
             size(dfduR_face,2) )
  end
  # loop over the volume variables that we are differentiating w.r.t.
  # loop over the face nodes
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    # loop over the interpolation stencil
    for j = 1:sbpface.stencilsize
      permL = sbpface.perm[j, iface.faceL]
      permR = sbpface.perm[j, iface.faceR]
      for q = 1:size(dfduL,2)
        for p = 1:size(dfduL,1)
          dfduL_face[p,q,i,permL] += (sbpface.interp[j,i]*
                                  dfduL[p,q,permL])
          dfduR_face[p,q,i,permR] += (sbpface.interp[j,iR]*
                                  dfduR[p,q,permR])
        end
      end
    end
  end
end




function interiorFaceInterpolate_jac!(sbpface::SparseFace{Tsbp},
                                      iface::Interface,
                                      dfduL::AbstractArray{Tsol,3},
                                      dfduR::AbstractArray{Tsol,3},
                                      dfduL_face::AbstractArray{Tsol,4},
                                      dfduR_face::AbstractArray{Tsol,4}
                                      ) where {Tsbp,Tsol}
  @asserts_enabled begin
    @assert( size(dfduL_face,3) == size(dfduR_face,3) == sbpface.numnodes )
    @assert( size(dfduL,3) == size(dfduR,3) == size(dfduL_face,4) == 
             size(dfduR_face,4) )
    @assert( size(dfduL,1) == size(dfduL,2) == size(dfduR,1) == size(dfduR,2) ==
             size(dfduL_face,1) == size(dfduL_face,2) == size(dfduR_face,1) ==
             size(dfduR_face,2) )
  end
  # loop over the face nodes
  for i = 1:sbpface.numnodes
    iR = sbpface.nbrperm[i,iface.orient]
    permL = sbpface.perm[i, iface.faceL]
    permR = sbpface.perm[iR, iface.faceR]
    for q = 1:size(dfduL,2)
      for p = 1:size(dfduL,1)
        dfduL_face[p,q,i,permL] += dfduL[p,q,permL]
        dfduR_face[p,q,i,permR] += dfduR[p,q,permR]
      end
    end
  end
end

#------------------------------------------------------------------------------
# 5D methods

function interiorFaceInterpolate_jac!(sbpface::DenseFace{Tsbp},
                                      iface::Interface,
                                      dfduL::AbstractArray{Tsol,5},
                                      dfduR::AbstractArray{Tsol,5},
                                      dfduL_face::AbstractArray{Tsol,5},
                                      dfduR_face::AbstractArray{Tsol,5}
                                      ) where {Tsbp,Tsol}

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
              dfduL_face[p,q,d,i,k] += (sbpface.interp[j,i]*
                                        dfduL[p,q,d,sbpface.perm[j,iface.faceL],k])
              dfduR_face[p,q,d,i,k] += (sbpface.interp[j,iR]*
                                        dfduR[p,q,d,sbpface.perm[j,iface.faceR],k])
            end
          end
        end
      end
    end
  end
end

function interiorFaceInterpolate_jac!(sbpface::SparseFace{Tsbp},
                                      iface::Interface,
                                      dfduL::AbstractArray{Tsol,5},
                                      dfduR::AbstractArray{Tsol,5},
                                      dfduL_face::AbstractArray{Tsol,5},
                                      dfduR_face::AbstractArray{Tsol,5}
                                      ) where {Tsbp,Tsol}
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
  # loop over the volume variables that we are differentiating w.r.t.
  for k = 1:size(dfduL,5)
    # loop over the face nodes
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,iface.orient]
      for d = 1:size(dfduL,3)
        for q = 1:size(dfduL,2)
          for p = 1:size(dfduL,1)
            dfduL_face[p,q,d,i,k] += dfduL[p,q,d,sbpface.perm[i,iface.faceL],k]
            dfduR_face[p,q,d,i,k] += dfduR[p,q,d,sbpface.perm[iR,iface.faceR],k]
          end
        end
      end
    end
  end
end
