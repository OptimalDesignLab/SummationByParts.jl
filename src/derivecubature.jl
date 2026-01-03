# This file contains a collection of functions to derive quadrature rules.
using SummationByParts
using SummationByParts.OrthoPoly
using SummationByParts.Cubature
using SummationByParts.SymCubatures
using SummationByParts.Optimizer
using DataStructures

"""
### SummationByParts.deriveTriCubatureOmega

This function derives quadrature rules for SBP-Omega operators
on the triangle.

**Inputs**

* `q`: the degree of the operator
* `vertices` : if true, vertices are present in the set of nodes
* `midedges` : if true, edge midpoints are present in set of nodes
* `centroid` : if true, centroid is present in set of nodes
* `numedge` : number of unique edge parameters
* `numS21` : number of S21 orbits (vertex to opposite face)
* `numS111` : number of S111 orbits
* `delta1`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin <= 0.1
* `delta2`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin > 0.1
* `verbose`: print out iteration results
* `xinit`: initial guess of the prameter and/or the parameter and weights
* `xinit_sym_group`: list of the symmetry group ordering provided in xinit

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function deriveTriCubatureOmega(;q::Int=1,
                                vertices::Bool=false, 
                                midedges::Bool=false, 
                                numS21::Int=0, 
                                numedge::Int=0, 
                                numS111::Int=0, 
                                centroid::Bool=false,
                                delta1=1e-2,
                                delta2=1e-2,
                                verbose::Bool=false,
                                xinit_sym_group=[],
                                xinit=[],
                                T=Float64)

    delta1 = T(delta1)
    delta2 = T(delta2)

    @assert(vertices==false && midedges==false && numedge==0, "Unsupported type of symmetry group provided.")

    if xinit_sym_group==[]
        if vertices==true
            push!(xinit_sym_group,"vertices")
        end
        if midedges==true
            push!(xinit_sym_group,"midedges")
        end
        if numS21 != 0
            push!(xinit_sym_group,"numS21")
        end
        if numedge != 0
            push!(xinit_sym_group,"numedge")
        end
        if numS111 != 0
            push!(xinit_sym_group,"numS111")
        end
        if centroid==true
            push!(xinit_sym_group,"centroid")
        end
    end

    tol = Cubature.default_tol(T)
    vtx = T[-1 -1; 1 -1; -1 1]
    cub = SymCubatures.TriSymCub{T}(vertices=vertices, 
                                    midedges=midedges, 
                                    numS21=numS21, 
                                    numedge=numedge,
                                    numS111=numS111, 
                                    centroid=centroid)
    
    # find the indices of the parameters we want to solve for 
    mask = 1:(cub.numparams+cub.numweights)

    # compute the number of parameters and weights
    numparams = cub.numparams
    numweights = cub.numweights

    # If initial guess is provided, sort them to match the ordering of symmetry groups considered in the main code
    if length(xinit)==numparams
        x_w = T(0.1) .* ones(T, numweights, 1)
        xinit = collect(Iterators.flatten(collect(Iterators.flatten(vcat(xinit, x_w)))))
    end

    sym_group = ["vertices", "midedges", "numS21", "numedge", "numS111", "centroid"]
    param_dict = OrderedDict{String, Vector{T}}()
    weight_dict = OrderedDict{String, Vector{T}}()    
    if xinit != []
        p_loc = [0]
        w_loc = [numparams]
        for s in xinit_sym_group
            if s == sym_group[1]
                push!(w_loc, w_loc[end]+1)
            elseif s ==sym_group[2]
                push!(w_loc, w_loc[end]+1)
            elseif s ==sym_group[3]
                push!(p_loc, p_loc[end]+numS21)
                push!(w_loc, w_loc[end]+numS21)        
            elseif s ==sym_group[4]
                push!(p_loc, p_loc[end]+numedge)
                push!(w_loc, w_loc[end]+numedge)
            elseif s==sym_group[5]
                push!(p_loc, p_loc[end]+2*numS111)
                push!(w_loc, w_loc[end]+numS111)
            elseif s==sym_group[6]
                push!(w_loc, w_loc[end]+1)
            end
        end

        p_cnt = 1
        for i=eachindex(xinit_sym_group)
            if (xinit_sym_group[i] ∉ ["vertices", "midedges", "centroid"])
                param_dict[xinit_sym_group[i]] = xinit[p_loc[p_cnt]+1:p_loc[p_cnt+1]]
                p_cnt += 1
            end
        end

        for i=2:length(xinit_sym_group)+1
            weight_dict[xinit_sym_group[i-1]] = xinit[w_loc[i-1]+1:w_loc[i]]
        end

        param_sorted = OrderedDict(key => param_dict[key] for key in sym_group if haskey(param_dict, key))
        weight_sorted = OrderedDict(key => weight_dict[key] for key in sym_group if haskey(weight_dict, key))
        xinit_param = values(param_sorted)
        xinit_weight = values(weight_sorted)

        xinit = collect(Iterators.flatten(collect(Iterators.flatten(hcat(xinit_param, xinit_weight)))))
    end 
    # solve for quadrature rule
    Cubature.solvecubature!(cub, q, mask, tol=tol, hist=true, verbose=verbose, xinit=xinit, delta1=delta1, delta2=delta2)

    return cub, vtx
end

"""
### SummationByParts.deriveTriCubatureGamma

This function derives quadrature rules for SBP-Gamma operators
on the triangle.

**Inputs**

* `q`: the degree of the operator
* `vertices` : if true, vertices are present in the set of nodes
* `midedges` : if true, edge midpoints are present in set of nodes
* `centroid` : if true, centroid is present in set of nodes
* `numedge` : number of unique edge parameters
* `numS21` : number of S21 orbits (vertex to opposite face)
* `numS111` : number of S111 orbits
* `delta1`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin <= 0.1
* `delta2`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin > 0.1
* `verbose`: print out iteration results
* `xinit`: initial guess of the prameter and/or the parameter and weights
* `xinit_sym_group`: list of the symmetry group ordering provided in xinit

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function deriveTriCubatureGamma(;q::Int=1,
                                vertices::Bool=false, 
                                midedges::Bool=false, 
                                numS21::Int=0, 
                                numedge::Int=0, 
                                numS111::Int=0, 
                                centroid::Bool=false,
                                delta1=1e-2,
                                delta2=1e-2,
                                verbose::Bool=false,
                                xinit_sym_group=[], 
                                xinit=[],
                                T=Float64)

    delta1 = T(delta1)
    delta2 = T(delta2)

    if xinit_sym_group==[]
        if vertices==true
            push!(xinit_sym_group,"vertices")
        end
        if midedges==true
            push!(xinit_sym_group,"midedges")
        end
        if numS21 != 0
            push!(xinit_sym_group,"numS21")
        end
        if numedge != 0
            push!(xinit_sym_group,"numedge")
        end
        if numS111 != 0
            push!(xinit_sym_group,"numS111")
        end
        if centroid==true
            push!(xinit_sym_group,"centroid")
        end
    end

    tol = Cubature.default_tol(T)
    vtx = T[-1 -1; 1 -1; -1 1]
    cub = SymCubatures.TriSymCub{T}(vertices=vertices, 
                                    midedges=midedges, 
                                    numS21=numS21, 
                                    numedge=numedge,
                                    numS111=numS111, 
                                    centroid=centroid)

    # find the indices of the parameters we want to solve for 
    mask = 1:(cub.numparams+cub.numweights)

    # compute the number of parameters and weights
    numparams = cub.numparams
    numweights = cub.numweights

    # If initial guess is provided, sort them to match the ordering of symmetry groups considered in the main code
    if length(xinit)==numparams
        x_w = T(0.1) .* ones(T, numweights, 1)
        xinit = collect(Iterators.flatten(collect(Iterators.flatten(vcat(xinit, x_w)))))
    end
    sym_group = ["vertices", "midedges", "numS21", "numedge", "numS111", "centroid"]
    if xinit != []
        p_loc = [0]
        w_loc = [numparams]
        for s in xinit_sym_group
            if s == sym_group[1]
                push!(w_loc, w_loc[end]+1)
            elseif s ==sym_group[2]
                push!(w_loc, w_loc[end]+1)
            elseif s ==sym_group[3]
                push!(p_loc, p_loc[end]+numS21)
                push!(w_loc, w_loc[end]+numS21)        
            elseif s ==sym_group[4]
                push!(p_loc, p_loc[end]+numedge)
                push!(w_loc, w_loc[end]+numedge)
            elseif s==sym_group[5]
                push!(p_loc, p_loc[end]+2*numS111)
                push!(w_loc, w_loc[end]+numS111)
            elseif s==sym_group[6]
                push!(w_loc, w_loc[end]+1)
            end
        end

        param_dict = OrderedDict{String, Vector{T}}()
        weight_dict = OrderedDict{String, Vector{T}}()
        p_cnt = 1
        for i=eachindex(xinit_sym_group)
            if (xinit_sym_group[i] ∉ ["vertices", "midedges", "centroid"])
                param_dict[xinit_sym_group[i]] = xinit[p_loc[p_cnt]+1:p_loc[p_cnt+1]]
                p_cnt += 1
            end
        end

        for i=2:length(xinit_sym_group)+1
            weight_dict[xinit_sym_group[i-1]] = xinit[w_loc[i-1]+1:w_loc[i]]
        end

        param_sorted = OrderedDict(key => param_dict[key] for key in sym_group if haskey(param_dict, key))
        weight_sorted = OrderedDict(key => weight_dict[key] for key in sym_group if haskey(weight_dict, key))
        xinit_param = values(param_sorted)
        xinit_weight = values(weight_sorted)

        xinit = collect(Iterators.flatten(collect(Iterators.flatten(hcat(xinit_param, xinit_weight)))))
    else
        xinit = 0.1 .* ones(numparams+numweights, 1)
    end

    # solve for quadrature rule
    Cubature.solvecubature!(cub, q, mask, tol=tol, hist=true, verbose=verbose, xinit=xinit, delta1=delta1, delta2=delta2)

    return cub, vtx
end

"""
### SummationByParts.deriveTriCubatureDiagE

This function derives quadrature rules for SBP-Omega operators
on the triangle.

**Inputs**

* `q`: the degree of the operator
* `vertices` : if true, vertices are present in the set of nodes
* `midedges` : if true, edge midpoints are present in set of nodes
* `centroid` : if true, centroid is present in set of nodes
* `numedge` : number of unique edge parameters
* `numS21` : number of S21 orbits (vertex to opposite face)
* `numS111` : number of S111 orbits
* `delta1`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin <= 0.1
* `delta2`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin > 0.1
* `verbose`: print out iteration results
* `xinit`: initial guess of the prameter and/or the parameter and weights
* `xinit_sym_group`: list of the symmetry group ordering provided in xinit
* `xedge`: parameters for the edge symmetry groups
* `xedge_sym_group`: list of the symmetry group ordering provided in xedge

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function deriveTriCubatureDiagE(;q::Int=1,
                                vertices::Bool=false, 
                                midedges::Bool=false, 
                                numS21::Int=0, 
                                numedge::Int=0, 
                                numS111::Int=0, 
                                centroid::Bool=false,
                                delta1=1e-2,
                                delta2=1e-2,
                                verbose::Bool=false,
                                xinit_sym_group=[], 
                                xinit=[],
                                xedge_sym_group=[], 
                                xedge=[],
                                T=Float64)

    delta1 = T(delta1)
    delta2 = T(delta2)

    if xinit_sym_group==[]
        if vertices==true
            push!(xinit_sym_group,"vertices")
        end
        if midedges==true
            push!(xinit_sym_group,"midedges")
        end
        if numS21 != 0
            push!(xinit_sym_group,"numS21")
        end
        if numedge != 0
            push!(xinit_sym_group,"numedge")
        end
        if numS111 != 0
            push!(xinit_sym_group,"numS111")
        end
        if centroid==true
            push!(xinit_sym_group,"centroid")
        end
    end

    if xedge_sym_group==[]
        if vertices==true
            push!(xedge_sym_group,"vertices")
        end
        if midedges==true
            push!(xedge_sym_group,"midedges")
        end
        if numedge != 0
            push!(xedge_sym_group,"numedge")
        end
    end

    tol = Cubature.default_tol(T)
    vtx = T[-1 -1; 1 -1; -1 1]
    cub = SymCubatures.TriSymCub{T}(vertices=vertices, 
                                    midedges=midedges, 
                                    numS21=numS21, 
                                    numedge=numedge,
                                    numS111=numS111, 
                                    centroid=centroid)

    # find the indices of the parameters we want to solve for 
    mask = SymCubatures.getInternalParamMask(cub)
    append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights))  

    # compute the number of parameters and weights
    numparams = cub.numparams
    numweights = cub.numweights

    # If initial guess is provided, sort them to match the ordering of symmetry groups considered in the main code
    if length(xinit)==numparams
        x_w = T(0.1) .* ones(T, numweights, 1)
        xinit = collect(Iterators.flatten(collect(Iterators.flatten(vcat(xinit, x_w)))))
    end
    if xinit == []
        xinit = 0.1 .* ones(numparams+numweights, 1)
    end
    sym_group = ["vertices", "midedges", "numS21", "numedge", "numS111", "centroid"]
    param_dict = OrderedDict{String, Vector{T}}()
    weight_dict = OrderedDict{String, Vector{T}}()
    if xinit != []
        p_loc = [0]
        w_loc = [numparams]
        for s in xinit_sym_group
            if s == sym_group[1]
                push!(w_loc, w_loc[end]+1)
            elseif s ==sym_group[2]
                push!(w_loc, w_loc[end]+1)
            elseif s ==sym_group[3]
                push!(p_loc, p_loc[end]+numS21)
                push!(w_loc, w_loc[end]+numS21)        
            elseif s ==sym_group[4]
                push!(p_loc, p_loc[end]+numedge)
                push!(w_loc, w_loc[end]+numedge)
            elseif s==sym_group[5]
                push!(p_loc, p_loc[end]+2*numS111)
                push!(w_loc, w_loc[end]+numS111)
            elseif s==sym_group[6]
                push!(w_loc, w_loc[end]+1)
            end
        end

        p_cnt = 1
        for i=eachindex(xinit_sym_group)
            if (xinit_sym_group[i] ∉ ["vertices", "midedges", "centroid"])
                param_dict[xinit_sym_group[i]] = xinit[p_loc[p_cnt]+1:p_loc[p_cnt+1]]
                p_cnt += 1
            end
        end

        for i=2:length(xinit_sym_group)+1
            weight_dict[xinit_sym_group[i-1]] = xinit[w_loc[i-1]+1:w_loc[i]]
        end

        param_sorted = OrderedDict(key => param_dict[key] for key in sym_group if haskey(param_dict, key))
        weight_sorted = OrderedDict(key => weight_dict[key] for key in sym_group if haskey(weight_dict, key))
        xinit_param = values(param_sorted)
        xinit_weight = values(weight_sorted)

        xinit = collect(Iterators.flatten(collect(Iterators.flatten(hcat(xinit_param, xinit_weight)))))
    end

    # If edge parameters are provided, sort them to match the ordering of symmetry groups considered in the main code
    if xedge != []
        p_loc = [0]
        for s in xedge_sym_group
            if s ==sym_group[4]
                push!(p_loc, p_loc[end]+numedge)
            end
        end

        param_dict_edge = OrderedDict{String, Vector{T}}()
        p_cnt = 1
        for i=eachindex(xedge_sym_group)
            if (xedge_sym_group[i] ∉ ["vertices", "midedges", "centroid"])
                param_dict_edge[xedge_sym_group[i]] = xedge[p_loc[p_cnt]+1:p_loc[p_cnt+1]]
                p_cnt += 1
            end
        end

        for i=eachindex(xedge_sym_group)
            if (xedge_sym_group[i] ∉ ["vertices", "midedges", "centroid"])
                param_dict[xedge_sym_group[i]] = param_dict_edge[xedge_sym_group[i]]
            end
        end

        param_sorted = OrderedDict(key => param_dict[key] for key in sym_group if haskey(param_dict, key))
        weight_sorted = OrderedDict(key => weight_dict[key] for key in sym_group if haskey(weight_dict, key))
        xinit_param = values(param_sorted)
        xinit_weight = values(weight_sorted)

        xinit = collect(Iterators.flatten(collect(Iterators.flatten(hcat(xinit_param, xinit_weight)))))
    end

    # solve for quadrature rule
    Cubature.solvecubature!(cub, q, mask, tol=tol, hist=true, verbose=verbose, xinit=xinit, delta1=delta1, delta2=delta2)

    return cub, vtx
end

"""
### SummationByParts.deriveTetCubatureOmega

This function derives quadrature rules for SBP-Omega operators
on the tetrahedron.

**Inputs**

* `q`: the degree of the operator
* `vertices` : if true, vertices are present in the set of nodes
* `midedges` : if true, edge midpoints are present in set of nodes
* `centroid` : if true, centroid is present in set of nodes
* `facecentroid` : if true, face centroids are present in the set of nodes
* `numedge` : number of unique edge parameters
* `numfaceS21` : number of S21 face orbits (same tri orbit on face)
* `numfaceS111` : number of S111 face orbits (same tri orbit on face)
* `numS31` : number of S31 orbits (vertex to opposite face)
* `numS22` : number of S22 orbits
* `numS211`: number of S211 orbits
* `numS1111`: number of S1111 orbits
* `delta1`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin <= 0.1
* `delta2`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin > 0.1
* `verbose`: print out iteration results
* `xinit`: initial guess of the prameter and/or the parameter and weights
* `xinit_sym_group`: list of the symmetry group ordering provided in xinit

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function deriveTetCubatureOmega(;q::Int=1,
                                vertices::Bool=false, 
                                numS31::Int=0,
                                midedges::Bool=false, 
                                numS22::Int=0,
                                numfaceS21::Int=0, 
                                numedge::Int=0, 
                                numS211::Int=0,
                                numfaceS111::Int=0, 
                                facecentroid::Bool=false,
                                numS1111::Int=0,
                                centroid::Bool=false,
                                delta1=1e-2,
                                delta2=1e-2,
                                verbose::Bool=false,
                                xinit_sym_group=[],
                                xinit=[],
                                T=Float64)

    delta1 = T(delta1)
    delta2 = T(delta2)

    @assert(vertices==false && midedges==false && numedge==0, 
            numfaceS21==0 && numfaceS111==0 && facecentroid==false,
            "Unsupported type of symmetry group provided.")

    if xinit_sym_group==[]
        if vertices==true
            push!(xinit_sym_group,"vertices")
        end
        if numS31!=0
            push!(xinit_sym_group,"numS31")
        end
        if midedges==true
            push!(xinit_sym_group,"midedges")
        end
        if numS22 != 0
            push!(xinit_sym_group,"numS22")
        end
        if numfaceS21 != 0
            push!(xinit_sym_group,"numfaceS21")
        end
        if numedge != 0
            push!(xinit_sym_group,"numedge")
        end
        if numS211 != 0
            push!(xinit_sym_group,"numS211")
        end
        if numfaceS111 != 0
            push!(xinit_sym_group,"numfaceS111")
        end
        if facecentroid != 0
            push!(xinit_sym_group,"facecentroid")
        end
        if numS1111 != 0
            push!(xinit_sym_group,"numS1111")
        end
        if centroid==true
            push!(xinit_sym_group,"centroid")
        end
    end
    tol = Cubature.default_tol(T)
    vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
    cub = SymCubatures.TetSymCub{T}(vertices=vertices, 
                                    numS31=numS31,
                                    midedges=midedges,
                                    numS22=numS22,
                                    numfaceS21=numfaceS21,
                                    numedge=numedge, 
                                    numS211=numS211,
                                    numfaceS111=numfaceS111,
                                    facecentroid=facecentroid,
                                    numS1111=numS1111,
                                    centroid=centroid)
    
    # find the indices of the parameters we want to solve for 
    mask = 1:(cub.numparams+cub.numweights)

    # compute the number of parameters and weights
    numparams = cub.numparams
    numweights = cub.numweights

    # If initial guess is provided, sort them to match the ordering of symmetry groups considered in the main code
    if length(xinit)==numparams
        x_w = T(0.1) .* ones(T, numweights, 1)
        xinit = collect(Iterators.flatten(collect(Iterators.flatten(vcat(xinit, x_w)))))
    end
    sym_group = ["vertices", "numS31", "midedges", "numS22", "numfaceS21", "numedge", 
                 "numS211", "numfaceS111", "facecentroid", "numS1111", "centroid"]
    param_dict = OrderedDict{String, Vector{T}}()
    weight_dict = OrderedDict{String, Vector{T}}()
    if xinit != []
        p_loc = [0]
        w_loc = [numparams]
        for s in xinit_sym_group
            if s == sym_group[1]
                push!(w_loc, w_loc[end]+1)
            elseif s ==sym_group[2]
                push!(p_loc, p_loc[end]+numS31)
                push!(w_loc, w_loc[end]+numS31)
            elseif s ==sym_group[3]
                push!(w_loc, w_loc[end]+1)        
            elseif s ==sym_group[4]
                push!(p_loc, p_loc[end]+numS22)
                push!(w_loc, w_loc[end]+numS22)
            elseif s==sym_group[5]
                push!(p_loc, p_loc[end]+numfaceS21)
                push!(w_loc, w_loc[end]+numfaceS21)
            elseif s==sym_group[6]
                push!(p_loc, p_loc[end]+numedge)
                push!(w_loc, w_loc[end]+numedge)
            elseif s==sym_group[7]
                push!(p_loc, p_loc[end]+2*numS211)
                push!(w_loc, w_loc[end]+numS211)
            elseif s==sym_group[8]
                push!(p_loc, p_loc[end]+2*numfaceS111)
                push!(w_loc, w_loc[end]+numfaceS111)
            elseif s==sym_group[9]
                push!(w_loc, w_loc[end]+1)
            elseif s==sym_group[10]
                push!(p_loc, p_loc[end]+3*numS1111)
                push!(w_loc, w_loc[end]+numS1111)
            elseif s==sym_group[11]
                push!(w_loc, w_loc[end]+1)
            end
        end

        p_cnt = 1
        for i=eachindex(xinit_sym_group)
            if (xinit_sym_group[i] ∉ ["vertices", "midedges", "facecentroid", "centroid"])
                param_dict[xinit_sym_group[i]] = xinit[p_loc[p_cnt]+1:p_loc[p_cnt+1]]
                p_cnt += 1
            end
        end

        for i=2:length(xinit_sym_group)+1
            weight_dict[xinit_sym_group[i-1]] = xinit[w_loc[i-1]+1:w_loc[i]]
        end

        param_sorted = OrderedDict(key => param_dict[key] for key in sym_group if haskey(param_dict, key))
        weight_sorted = OrderedDict(key => weight_dict[key] for key in sym_group if haskey(weight_dict, key))
        xinit_param = values(param_sorted)
        xinit_weight = values(weight_sorted)

        xinit = collect(Iterators.flatten(collect(Iterators.flatten(hcat(xinit_param, xinit_weight)))))
    end 
    # solve for quadrature rule
    Cubature.solvecubature!(cub, q, mask, tol=tol, hist=true, verbose=verbose, xinit=xinit, delta1=delta1, delta2=delta2)

    return cub, vtx
end

"""
### SummationByParts.deriveTetCubatureGamma

This function derives quadrature rules for SBP-Gamma operators
on the tetrahedron.

**Inputs**

* `q`: the degree of the operator
* `vertices` : if true, vertices are present in the set of nodes
* `midedges` : if true, edge midpoints are present in set of nodes
* `centroid` : if true, centroid is present in set of nodes
* `facecentroid` : if true, face centroids are present in the set of nodes
* `numedge` : number of unique edge parameters
* `numfaceS21` : number of S21 face orbits (same tri orbit on face)
* `numfaceS111` : number of S111 face orbits (same tri orbit on face)
* `numS31` : number of S31 orbits (vertex to opposite face)
* `numS22` : number of S22 orbits
* `numS211`: number of S211 orbits
* `numS1111`: number of S1111 orbits
* `delta1`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin <= 0.1
* `delta2`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin > 0.1
* `verbose`: print out iteration results
* `xinit`: initial guess of the prameter and/or the parameter and weights
* `xinit_sym_group`: list of the symmetry group ordering provided in xinit

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function deriveTetCubatureGamma(;q::Int=1,
                                vertices::Bool=false, 
                                numS31::Int=0,
                                midedges::Bool=false, 
                                numS22::Int=0,
                                numfaceS21::Int=0, 
                                numedge::Int=0, 
                                numS211::Int=0,
                                numfaceS111::Int=0, 
                                facecentroid::Bool=false,
                                numS1111::Int=0,
                                centroid::Bool=false,
                                delta1=1e-2,
                                delta2=1e-2,
                                verbose::Bool=false,
                                xinit_sym_group=[],
                                xinit=[],
                                T=Float64)

    delta1 = T(delta1)
    delta2 = T(delta2)

    if xinit_sym_group==[]
        if vertices==true
            push!(xinit_sym_group,"vertices")
        end
        if numS31!=0
            push!(xinit_sym_group,"numS31")
        end
        if midedges==true
            push!(xinit_sym_group,"midedges")
        end
        if numS22 != 0
            push!(xinit_sym_group,"numS22")
        end
        if numfaceS21 != 0
            push!(xinit_sym_group,"numfaceS21")
        end
        if numedge != 0
            push!(xinit_sym_group,"numedge")
        end
        if numS211 != 0
            push!(xinit_sym_group,"numS211")
        end
        if numfaceS111 != 0
            push!(xinit_sym_group,"numfaceS111")
        end
        if facecentroid != 0
            push!(xinit_sym_group,"facecentroid")
        end
        if numS1111 != 0
            push!(xinit_sym_group,"numS1111")
        end
        if centroid==true
            push!(xinit_sym_group,"centroid")
        end
    end
    tol = Cubature.default_tol(T)
    vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
    cub = SymCubatures.TetSymCub{T}(vertices=vertices, 
                                    numS31=numS31,
                                    midedges=midedges,
                                    numS22=numS22,
                                    numfaceS21=numfaceS21,
                                    numedge=numedge, 
                                    numS211=numS211,
                                    numfaceS111=numfaceS111,
                                    facecentroid=facecentroid,
                                    numS1111=numS1111,
                                    centroid=centroid)
    
    # find the indices of the parameters we want to solve for 
    mask = 1:(cub.numparams+cub.numweights)

    # compute the number of parameters and weights
    numparams = cub.numparams
    numweights = cub.numweights

    # If initial guess is provided, sort them to match the ordering of symmetry groups considered in the main code
    if length(xinit)==numparams
        x_w = T(0.1) .* ones(T, numweights, 1)
        xinit = collect(Iterators.flatten(collect(Iterators.flatten(vcat(xinit, x_w)))))
    end
    sym_group = ["vertices", "numS31", "midedges", "numS22", "numfaceS21", "numedge", 
                 "numS211", "numfaceS111", "facecentroid", "numS1111", "centroid"]
    param_dict = OrderedDict{String, Vector{T}}()
    weight_dict = OrderedDict{String, Vector{T}}()
    if xinit != []
        p_loc = [0]
        w_loc = [numparams]
        for s in xinit_sym_group
            if s == sym_group[1]
                push!(w_loc, w_loc[end]+1)
            elseif s ==sym_group[2]
                push!(p_loc, p_loc[end]+numS31)
                push!(w_loc, w_loc[end]+numS31)
            elseif s ==sym_group[3]
                push!(w_loc, w_loc[end]+1)        
            elseif s ==sym_group[4]
                push!(p_loc, p_loc[end]+numS22)
                push!(w_loc, w_loc[end]+numS22)
            elseif s==sym_group[5]
                push!(p_loc, p_loc[end]+numfaceS21)
                push!(w_loc, w_loc[end]+numfaceS21)
            elseif s==sym_group[6]
                push!(p_loc, p_loc[end]+numedge)
                push!(w_loc, w_loc[end]+numedge)
            elseif s==sym_group[7]
                push!(p_loc, p_loc[end]+2*numS211)
                push!(w_loc, w_loc[end]+numS211)
            elseif s==sym_group[8]
                push!(p_loc, p_loc[end]+2*numfaceS111)
                push!(w_loc, w_loc[end]+numfaceS111)
            elseif s==sym_group[9]
                push!(w_loc, w_loc[end]+1)
            elseif s==sym_group[10]
                push!(p_loc, p_loc[end]+3*numS1111)
                push!(w_loc, w_loc[end]+numS1111)
            elseif s==sym_group[11]
                push!(w_loc, w_loc[end]+1)
            end
        end

        p_cnt = 1
        for i=eachindex(xinit_sym_group)
            if (xinit_sym_group[i] ∉ ["vertices", "midedges", "facecentroid", "centroid"])
                param_dict[xinit_sym_group[i]] = xinit[p_loc[p_cnt]+1:p_loc[p_cnt+1]]
                p_cnt += 1
            end
        end

        for i=2:length(xinit_sym_group)+1
            weight_dict[xinit_sym_group[i-1]] = xinit[w_loc[i-1]+1:w_loc[i]]
        end

        param_sorted = OrderedDict(key => param_dict[key] for key in sym_group if haskey(param_dict, key))
        weight_sorted = OrderedDict(key => weight_dict[key] for key in sym_group if haskey(weight_dict, key))
        xinit_param = values(param_sorted)
        xinit_weight = values(weight_sorted)

        xinit = collect(Iterators.flatten(collect(Iterators.flatten(hcat(xinit_param, xinit_weight)))))
    end 
    # solve for quadrature rule
    Cubature.solvecubature!(cub, q, mask, tol=tol, hist=true, verbose=verbose, xinit=xinit, delta1=delta1, delta2=delta2)

    return cub, vtx
end

"""
### SummationByParts.deriveTetCubatureDiagE

This function derives quadrature rules for SBP diagonal-E operators
on the tetrahedron.

**Inputs**

* `q`: the degree of the operator
* `vertices` : if true, vertices are present in the set of nodes
* `midedges` : if true, edge midpoints are present in set of nodes
* `centroid` : if true, centroid is present in set of nodes
* `facecentroid` : if true, face centroids are present in the set of nodes
* `numedge` : number of unique edge parameters
* `numfaceS21` : number of S21 face orbits (same tri orbit on face)
* `numfaceS111` : number of S111 face orbits (same tri orbit on face)
* `numS31` : number of S31 orbits (vertex to opposite face)
* `numS22` : number of S22 orbits
* `numS211`: number of S211 orbits
* `numS1111`: number of S1111 orbits
* `delta1`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin <= 0.1
* `delta2`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin > 0.1
* `verbose`: print out iteration results
* `xinit`: initial guess of the prameter and/or the parameter and weights
* `xinit_sym_group`: list of the symmetry group ordering provided in xinit
* `xedge`: parameters for the edge symmetry groups
* `xedge_sym_group`: list of the symmetry group ordering provided in xedge

**Outputs**

* `cub`: a symmetric cubature for the right triangle
* `vtx`: vertices for the right triangle

"""
function deriveTetCubatureDiagE(;q::Int=1,
                                vertices::Bool=false, 
                                numS31::Int=0,
                                midedges::Bool=false, 
                                numS22::Int=0,
                                numfaceS21::Int=0, 
                                numedge::Int=0, 
                                numS211::Int=0,
                                numfaceS111::Int=0, 
                                facecentroid::Bool=false,
                                numS1111::Int=0,
                                centroid::Bool=false,
                                delta1=1e-2,
                                delta2=1e-2,
                                verbose::Bool=false,
                                xinit_sym_group=[],
                                xinit=[],
                                xedge_sym_group=[],
                                xedge=[],
                                T=Float64)

    delta1 = T(delta1)
    delta2 = T(delta2)

    if xinit_sym_group==[]
        if vertices==true
            push!(xinit_sym_group,"vertices")
        end
        if numS31!=0
            push!(xinit_sym_group,"numS31")
        end
        if midedges==true
            push!(xinit_sym_group,"midedges")
        end
        if numS22 != 0
            push!(xinit_sym_group,"numS22")
        end
        if numfaceS21 != 0
            push!(xinit_sym_group,"numfaceS21")
        end
        if numedge != 0
            push!(xinit_sym_group,"numedge")
        end
        if numS211 != 0
            push!(xinit_sym_group,"numS211")
        end
        if numfaceS111 != 0
            push!(xinit_sym_group,"numfaceS111")
        end
        if facecentroid != 0
            push!(xinit_sym_group,"facecentroid")
        end
        if numS1111 != 0
            push!(xinit_sym_group,"numS1111")
        end
        if centroid==true
            push!(xinit_sym_group,"centroid")
        end
    end

    if xedge_sym_group==[]
        if vertices==true
            push!(xedge_sym_group,"vertices")
        end
        if midedges==true
            push!(xedge_sym_group,"midedges")
        end
        if numfaceS21 != 0
            push!(xedge_sym_group,"numfaceS21")
        end
        if numedge != 0
            push!(xedge_sym_group,"numedge")
        end
        if numfaceS111 != 0
            push!(xedge_sym_group,"numfaceS111")
        end
        if facecentroid != 0
            push!(xedge_sym_group,"facecentroid")
        end
    end

    tol = Cubature.default_tol(T)
    vtx = T[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
    cub = SymCubatures.TetSymCub{T}(vertices=vertices, 
                                    numS31=numS31,
                                    midedges=midedges,
                                    numS22=numS22,
                                    numfaceS21=numfaceS21,
                                    numedge=numedge, 
                                    numS211=numS211,
                                    numfaceS111=numfaceS111,
                                    facecentroid=facecentroid,
                                    numS1111=numS1111,
                                    centroid=centroid)
    
    # find the indices of the parameters we want to solve for 
    mask = SymCubatures.getInternalParamMask(cub)
    append!(mask, (cub.numparams+1):(cub.numparams+cub.numweights)) 

    # compute the number of parameters and weights
    numparams = cub.numparams
    numweights = cub.numweights

    # If initial guess is provided, sort them to match the ordering of symmetry groups considered in the main code
    if length(xinit)==numparams
        x_w = T(0.1) .* ones(T, numweights, 1)
        xinit = collect(Iterators.flatten(collect(Iterators.flatten(vcat(xinit, x_w)))))
    end
    if xinit == []
        xinit = 0.1 .* ones(numparams+numweights, 1)
    end
    sym_group = ["vertices", "numS31", "midedges", "numS22", "numfaceS21", "numedge", 
                 "numS211", "numfaceS111", "facecentroid", "numS1111", "centroid"]
    param_dict = OrderedDict{String, Vector{T}}()
    weight_dict = OrderedDict{String, Vector{T}}()
    if xinit != []
        p_loc = [0]
        w_loc = [numparams]
        for s in xinit_sym_group
            if s == sym_group[1]
                push!(w_loc, w_loc[end]+1)
            elseif s ==sym_group[2]
                push!(p_loc, p_loc[end]+numS31)
                push!(w_loc, w_loc[end]+numS31)
            elseif s ==sym_group[3]
                push!(w_loc, w_loc[end]+1)        
            elseif s ==sym_group[4]
                push!(p_loc, p_loc[end]+numS22)
                push!(w_loc, w_loc[end]+numS22)
            elseif s==sym_group[5]
                push!(p_loc, p_loc[end]+numfaceS21)
                push!(w_loc, w_loc[end]+numfaceS21)
            elseif s==sym_group[6]
                push!(p_loc, p_loc[end]+numedge)
                push!(w_loc, w_loc[end]+numedge)
            elseif s==sym_group[7]
                push!(p_loc, p_loc[end]+2*numS211)
                push!(w_loc, w_loc[end]+numS211)
            elseif s==sym_group[8]
                push!(p_loc, p_loc[end]+2*numfaceS111)
                push!(w_loc, w_loc[end]+numfaceS111)
            elseif s==sym_group[9]
                push!(w_loc, w_loc[end]+1)
            elseif s==sym_group[10]
                push!(p_loc, p_loc[end]+3*numS1111)
                push!(w_loc, w_loc[end]+numS1111)
            elseif s==sym_group[11]
                push!(w_loc, w_loc[end]+1)
            end
        end

        p_cnt = 1
        for i=eachindex(xinit_sym_group)
            if (xinit_sym_group[i] ∉ ["vertices", "midedges", "facecentroid", "centroid"])
                param_dict[xinit_sym_group[i]] = xinit[p_loc[p_cnt]+1:p_loc[p_cnt+1]]
                p_cnt += 1
            end
        end

        for i=2:length(xinit_sym_group)+1
            weight_dict[xinit_sym_group[i-1]] = xinit[w_loc[i-1]+1:w_loc[i]]
        end

        param_sorted = OrderedDict(key => param_dict[key] for key in sym_group if haskey(param_dict, key))
        weight_sorted = OrderedDict(key => weight_dict[key] for key in sym_group if haskey(weight_dict, key))
        xinit_param = values(param_sorted)
        xinit_weight = values(weight_sorted)

        xinit = collect(Iterators.flatten(collect(Iterators.flatten(hcat(xinit_param, xinit_weight)))))
    end 

    # If edge parameters are provided, sort them to match the ordering of symmetry groups considered in the main code
    if xedge != []
        p_loc = [0]
        for s in xedge_sym_group
            if s ==sym_group[5]
                push!(p_loc, p_loc[end]+numfaceS21)        
            elseif s ==sym_group[6]
                push!(p_loc, p_loc[end]+numedge)
            elseif s==sym_group[8]
                push!(p_loc, p_loc[end]+2*numfaceS111)
            end
        end

        param_dict_edge = OrderedDict{String, Vector{T}}()
        p_cnt = 1
        for i=eachindex(xedge_sym_group)
            if (xedge_sym_group[i] ∉ ["vertices", "midedges", "facecentroid","centroid"])
                param_dict_edge[xedge_sym_group[i]] = xedge[p_loc[p_cnt]+1:p_loc[p_cnt+1]]
                p_cnt += 1
            end
        end

        for i=eachindex(xedge_sym_group)
            if (xedge_sym_group[i] ∉ ["vertices", "midedges","facecentroid","centroid"])
                param_dict[xedge_sym_group[i]] = param_dict_edge[xedge_sym_group[i]]
            end
        end

        param_sorted = OrderedDict(key => param_dict[key] for key in sym_group if haskey(param_dict, key))
        weight_sorted = OrderedDict(key => weight_dict[key] for key in sym_group if haskey(weight_dict, key))
        xinit_param = values(param_sorted)
        xinit_weight = values(weight_sorted)

        xinit = collect(Iterators.flatten(collect(Iterators.flatten(hcat(xinit_param, xinit_weight)))))
    end

    # solve for quadrature rule
    Cubature.solvecubature!(cub, q, mask, tol=tol, hist=true, verbose=verbose, xinit=xinit, delta1=delta1, delta2=delta2)

    return cub, vtx
end
