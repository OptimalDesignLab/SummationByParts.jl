# This file gathers together functions related to edge stabilization

"""
### SummationByParts.edgestabilize!

Applies edge stabilization to a given field, differentiating in the direction
specified by `dirvec`, and scaling by the `tau` field.

**Inputs**

* `sbpface`: an SBP face operator type
* `ifaces`: list of element interfaces stored as an array of `Interface`s
* `dirvec`: direction to differentiate in [xi coord, face node, L/R, face] format
* `tau`: scaling term in [face node, face] format
* `u`: field being stablized in [vol node, element] format
* `±`: PlusFunctor to add to res, MinusFunctor to subract

**In/Outs**

* `res`: where the result is stored in [vol node, element] format

"""
function edgestabilize!(sbpface::AbstractFace{Tsbp},
                                        ifaces::Array{Interface},
                                        dirvec::Array{Tmsh,4},
                                        tau::Array{Tmsh,2},
                                        u::Array{Tsol,2},
                                        res::Array{Tsol,2},
                                        (±)::UnaryFunctor=Add()) where {Tsbp,Tmsh,Tsol}
  @assert( size(u) == size(res) )
  @assert( size(ifaces,1) == size(dirvec,4) == size(tau,2) )
  @assert( size(dirvec,2) == size(tau,1) == sbpface.numnodes )
  @assert( size(dirvec,3) == 2 )
  dudξ = zeros(Tsol, (2))
  left = 1
  right = 1
  for (findex, face) in enumerate(ifaces)
    for i = 1:sbpface.numnodes
      iR = sbpface.nbrperm[i,face.orient]
      Du = zero(Tsol)
      for di = 1:2
        fill!(dudξ, zero(Tsol))
        # compute the derivatives in the ξ[di] direction
        for j = 1:sbpface.dstencilsize
          dudξ[left] += sbpface.deriv[j,i,di]*u[sbpface.dperm[j,face.faceL],
                                                face.elementL]
          dudξ[right] += sbpface.deriv[j,iR,di]*u[sbpface.dperm[j,face.faceR],
                                                  face.elementR]
        end
        # contract with direction vector
        Du += dudξ[left]*dirvec[di,i,left,findex] + 
        dudξ[right]*dirvec[di,iR,right,findex]
      end
      # scale by tau and face cubature
      Du *= sbpface.wface[i]*tau[i,findex]
      # now apply transposed jump derivative
      for di = 1:2
        dudξ[left] = Du*dirvec[di,i,left,findex]
        dudξ[right] = Du*dirvec[di,iR,right,findex]
        for j = 1:sbpface.dstencilsize
          res[sbpface.dperm[j,face.faceL],face.elementL] +=
            ±(sbpface.deriv[j,i,di]*dudξ[left])
          res[sbpface.dperm[j,face.faceR],face.elementR] += 
          ±(sbpface.deriv[j,iR,di]*dudξ[right])
        end
      end
    end
  end
end
