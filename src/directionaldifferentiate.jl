# This file gathers together functions related to directional differentiation
# using the SBP operators

"""
### SummationByParts.directionalDifferentiateElement!

Performs a directional derivative (in reference space) at a given node.  The
input field `u` is for **a single element**, not a collection of elements.

**WARNING**: In the case of a vector field u, the directional derivative is
  added to the output Ddir; the user must zero this before.

**Inputs**

* `sbp`: an SBP operator type
* `dir`: a direction vector for the directional derivative
* `u`: the field that is being differentiated (either a scalar or vector)
* `i`: index of the node at which the derivative is desired

**Returns or In/Outs**

* `Ddir`: derivative of `u` in direction `dir`

"""
function directionalDifferentiateElement!(sbp::AbstractSBP{Tsbp},
                                                          dir::Array{Tmsh,1},
                                                          u::AbstractArray{Tsol,1},
                                                          i::Int) where {Tsbp,Tmsh,Tsol}
  @asserts_enabled begin
    @assert( size(sbp.Q, 3) == size(dir,1) )
  end

  Ddir = zero(Tmsh)*zero(Tsol)
  for di = 1:size(sbp.Q, 3)
    tmp = zero(Tmsh)*zero(Tsol)
    for j = 1:sbp.numnodes
      tmp += sbp.Q[i,j,di]*u[j]
    end
    Ddir += dir[di]*tmp/sbp.w[i]
  end
  return Ddir
end

function directionalDifferentiateElement!(sbp::AbstractSBP{Tsbp},
                                                               dir::Array{Tmsh,1}, 
                                                               u::AbstractArray{Tsol,2},
                                                               i::Int,
                                                               Ddir::Array{Tres,1}) where {Tsbp,Tmsh,Tsol,Tres}
  @asserts_enabled begin
    @assert( size(sbp.Q, 3) == size(dir,1) )
  end

  for di = 1:size(sbp.Q, 3)
    tmp = zeros(size(Ddir))
    for j = 1:sbp.numnodes
      for field = 1:size(u,1)
        tmp[field] += sbp.Q[i,j,di]*u[field,j]
      end
    end
    for field = 1:size(u,1)
      Ddir[field] += dir[di]*tmp[field]/sbp.w[i]
    end
  end
end
