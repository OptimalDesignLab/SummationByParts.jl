using SummationByParts
using SummationByParts.Cubature
using SummationByParts.SymCubatures
using Printf

"""
### SummationByParts.write_tri_data

This function writes to a file the SBP diagonal-E quadruature data on the right triangle.

**Inputs**

* `quad_degree`: the degree of the operator
* `vertices`: indicates whether to include vertices (LG or LGL type quadrature)
* `print_only`: indicates whether to print to screen or to file

"""
function write_tri_data(;quad_degree=nothing, vertices::Bool=true, print_only::Bool=false)

    quad_start=1
    quad_end=20
    if !isnothing(quad_degree)
        @assert(isa(quad_degree, Int) && quad_degree>=1 && quad_degree<=20, "Quadrature degree should be an integer in [1,20].")
        quad_start=quad_degree
        quad_end=quad_degree
    end

    for q = quad_start:quad_end
        cub, vtx = Cubature.getTriCubatureDiagE(q, Float64, vertices=vertices)
        vert_node = cub.vertices
        midedges_node = cub.midedges
        numS21_node = cub.numS21
        numedge_node = cub.numedge
        numS111_node = cub.numS111
        centroid_node = cub.centroid
        params = cub.params
        weights = cub.weights
        numnodes = cub.numnodes

        if vertices
            face_oper = "lgl"
        else
            face_oper = "lg"
        end
        dir = pwd() 
        path_compact = joinpath(dir, "quadrature_data/tri/compact")
        path_expanded = joinpath(dir, "quadrature_data/tri/expanded")
        file_cmp = string("tri_$face_oper","_q$q","_n$numnodes","_cmp.dat")
        file_ext = string("tri_$face_oper","_q$q","_n$numnodes","_ext.dat")
        file_compact = joinpath(path_compact,file_cmp)
        file_expanded = joinpath(path_expanded,file_ext)

        if !print_only
            pad_len = 10
            pad_len2 = 50
            open(file_compact,"w") do file_compact 
                write(file_compact, string(rpad("vertices",pad_len), rpad("midedges",pad_len), rpad("numS21",pad_len), rpad("numedge",pad_len), rpad("numS111",pad_len),rpad("centroid",pad_len),"\n"))
                write(file_compact, string(rpad("$vert_node",pad_len), rpad("$midedges_node",pad_len), rpad("$numS21_node",pad_len), rpad("$numedge_node",pad_len), 
                                rpad("$numS111_node",pad_len),rpad("$centroid_node",pad_len),"\n\n"))
                
                kp = 0
                kw = 0
                if cub.vertices
                    write(file_compact, "\nVertices \n")
                end
                if vert_node
                    w_vert = weights[1]
                    write(file_compact, string(rpad("[0.0]",pad_len2), rpad("$w_vert",pad_len2),"\n"))
                    kw+=1
                end

                if cub.midedges
                    write(file_compact, "\nMidedges \n")
                end
                if midedges_node
                    w_mid = weights[kw+1]
                    write(file_compact, string(rpad("[0.5]",pad_len2), rpad("$w_mid",pad_len2),"\n"))
                    kw+=1
                end

                if cub.numS21 != 0
                    write(file_compact, "\nS21 \n")
                end
                cnt=0
                for i=kw+1:kw+numS21_node
                    w_S21 = weights[i]
                    par1 = params[kp+1+cnt]./2
                    p_S21 = string("[", "$par1","]")
                    write(file_compact, string(rpad(p_S21,pad_len2), rpad("$w_S21",pad_len2),"\n"))
                    if (2*par1)>=1.0
                        error("Invalid quadrature, problem with S21 orbit.")
                    end
                    cnt+=1
                end
                kw+=numS21_node
                kp+=numS21_node 
                
                if cub.numedge != 0
                    write(file_compact, "\nEdge \n")
                end
                cnt=0
                for i=kw+1:kw+numedge_node
                    w_edge = weights[i]
                    par1 = params[kp+1+cnt]
                    p_edge = string("[","$par1","]")
                    write(file_compact, string(rpad(p_edge,pad_len2), rpad("$w_edge",pad_len2),"\n"))
                    if (par1)>=1.0
                        error("Invalid quadrature, problem with edge orbit.")
                    end
                    cnt+=1
                end
                kw+=numedge_node
                kp+=numedge_node
                
                if cub.numS111 != 0
                    write(file_compact, "\nS111 \n")
                end
                cnt=0
                for i=kw+1:kw+numS111_node
                    w_S111 = weights[i]
                    par1 = params[kp+1+cnt]./2
                    par2 = params[kp+2+cnt]./2
                    p_S111 = string("[","$par1",", ","$par2","]")
                    write(file_compact, string(rpad(p_S111,pad_len2), rpad("$w_S111",pad_len2),"\n"))
                    cnt+=2
                    if (par1+par2)>=1.0
                        error("Invalid quadrature, problem with S111 orbit.")
                    end
                end
                kw+=numS111_node
                kp+=2*numS111_node

                if cub.centroid
                    write(file_compact, "\nCentroid \n")
                end
                if centroid_node
                    w_cent = weights[end]
                    par1 = 1.0/3.0
                    p_cent = string("[", "$par1","]")
                    write(file_compact, string(rpad(p_cent,pad_len2), rpad("$w_cent",pad_len2),"\n"))
                    kw+=1
                end
            end
            open(file_expanded,"w") do file_expanded
                pad_len = 12
                write(file_expanded, string(rpad("centroid",pad_len),rpad("vertices",pad_len), rpad("midedges",pad_len), rpad("numS21",pad_len), rpad("numedge",pad_len), 
                                    rpad("numS111",pad_len), "\n"))
                write(file_expanded, string(rpad("$centroid_node (1)",pad_len),rpad("$vert_node (3)",pad_len), rpad("$midedges_node (3)",pad_len), rpad("$numS21_node (3)",pad_len), rpad("$numedge_node (6)",pad_len), 
                                    rpad("$numS111_node (6)",pad_len), "\n\n"))
                
                xy = SymCubatures.calcnodes(cub, vtx)
                ws = SymCubatures.calcweights(cub)

                pad_len = 25
                for i=1:size(xy,2)
                    x = round(xy[1,i], digits=16)
                    y = round(xy[2,i], digits=16)
                    w = round(ws[i], digits=16)
                    write(file_expanded, string(rpad(@sprintf("%.16f", x),pad_len), rpad(@sprintf("%.16f", y),pad_len), rpad(@sprintf("%.16f", w),pad_len),"\n"))
                end
            end
        else
            xy = SymCubatures.calcnodes(cub, vtx)
            ws = SymCubatures.calcweights(cub)
            
            println(" ")
            println(" ")
            println("-----------------------------------------------------------")
            println("      Triangle quadrature (q = $q, nnodes = $numnodes)  ")
            println("-----------------------------------------------------------")
            for i = 1:size(xy)[1]
                println(join(xy[i,:], ','))
                println(" ")
            end
            println(join(ws',','))
            println("------------------------------------------------------------")
        end
    end

end

"""
### SummationByParts.write_tet_data

This function writes to a file the SBP diagonal-E quadruature data on the right tetrahedron.

**Inputs**

* `quad_degree`: the degree of the operator
* `print_only`: indicates whether to print to screen or to file

"""
function write_tet_data(;quad_degree=nothing, print_only::Bool=false)

    quad_start=1
    quad_end=10
    if !isnothing(quad_degree)
        @assert(isa(quad_degree, Int) && quad_degree>=1 && quad_degree<=20, "Quadrature degree should be an integer in [1,20].")
        quad_start=quad_degree
        quad_end=quad_degree
    end

    for q = quad_start:quad_end
        cub, vtx = Cubature.getTetCubatureDiagE(q)
        vert_node = cub.vertices
        numS31_node =  cub.numS31
        midedges_node = cub.midedges
        numS22_node = cub.numS22
        numfaceS21_node = cub.numfaceS21
        numedge_node = cub.numedge
        numS211_node = cub.numS211
        numS1111_node = cub.numS1111
        facetcentroid_node = cub.facecentroid
        numfaceS111_node = cub.numfaceS111
        centroid_node = cub.centroid

        params = cub.params
        weights = cub.weights
        numnodes = cub.numnodes

        dir = pwd() 
        path_compact = joinpath(dir, "quadrature_data/tet/compact")
        path_expanded = joinpath(dir, "quadrature_data/tet/expanded")
        file_cmp = string("tet","_q$q","_n$numnodes","_cmp.dat")
        file_ext = string("tet","_q$q","_n$numnodes","_ext.dat")
        file_compact = joinpath(path_compact,file_cmp)
        file_expanded = joinpath(path_expanded,file_ext)

        if !print_only
            pad_len = 15
            pad_len2 = 50
            open(file_compact,"w") do file_compact 
                write(file_compact, string(rpad("vertices",pad_len), rpad("numS31",pad_len), rpad("midedges",pad_len), 
                        rpad("numS22",pad_len), rpad("numfaceS21",pad_len), rpad("numedge",pad_len), "\n"))
                write(file_compact, string(rpad("$vert_node (4)",pad_len), rpad("$numS31_node (4)",pad_len), rpad("$midedges_node (6)",pad_len), 
                        rpad("$numS22_node (6)",pad_len), rpad("$numfaceS21_node (12)",pad_len), rpad("$numedge_node (12)",pad_len),"\n"))
                write(file_compact, string(rpad("numS211",pad_len), rpad("numS1111",pad_len), rpad("facecentroid",pad_len), rpad("numfaceS111",pad_len),rpad("centroid",pad_len),"\n"))
                write(file_compact, string(rpad("$numS211_node (12)",pad_len),rpad("$numS1111_node (24)",pad_len), rpad("$facetcentroid_node (4)",pad_len), rpad("$numfaceS111_node (24)",pad_len), rpad("$centroid_node (1)",pad_len),"\n\n"))
                
                kp = 0
                kw = 0
                if cub.vertices != 0
                    write(file_compact, "\nVertices \n")
                end
                if vert_node
                    w = weights[1]
                    write(file_compact, string(rpad("[0.0]",pad_len2), rpad("$w",pad_len2),"\n"))
                    kw+=1
                end

                if cub.numS31 != 0
                    write(file_compact, "\nS31 \n")
                end
                cnt=0
                for i=kw+1:kw+numS31_node
                    w = weights[i]
                    par1 = params[kp+1+cnt]./3
                    par = string("[", "$par1","]")
                    write(file_compact, string(rpad(par,pad_len2), rpad("$w",pad_len2),"\n"))
                    if (3*par1)>=1.0
                        error("Invalid quadrature, problem with S31 orbit.")
                    end
                    cnt+=1
                end
                kw+=numS31_node
                kp+=numS31_node 

                if cub.midedges != 0
                    write(file_compact, "\nMidedges \n")
                end
                if midedges_node
                    w = weights[kw+1]
                    write(file_compact, string(rpad("[0.5]",pad_len2), rpad("$w",pad_len2),"\n"))
                    kw+=1
                end

                if cub.numS22 != 0
                    write(file_compact, "\nS22 \n")
                end
                cnt=0
                for i=kw+1:kw+numS22_node
                    w = weights[i]
                    par1 = params[kp+1+cnt]./2
                    par = string("[", "$par1","]")
                    write(file_compact, string(rpad(par,pad_len2), rpad("$w",pad_len2),"\n"))
                    if (2*par1)>=1.0
                        error("Invalid quadrature, problem with S22 orbit.")
                    end
                    cnt+=1
                end
                kw+=numS22_node
                kp+=numS22_node

                if cub.numfaceS21 != 0
                    write(file_compact, "\nFaceS21 \n")
                end
                cnt=0
                for i=kw+1:kw+numfaceS21_node
                    w = weights[i]
                    par1 = params[kp+1+cnt]./2
                    par = string("[", "$par1","]")
                    write(file_compact, string(rpad(par,pad_len2), rpad("$w",pad_len2),"\n"))
                    if (2*par1)>=1.0
                        error("Invalid quadrature, problem with S21 orbit.")
                    end
                    cnt+=1
                end
                kw+=numfaceS21_node
                kp+=numfaceS21_node 
                
                if cub.numedge != 0
                    write(file_compact, "\nEdge \n")
                end
                cnt=0
                for i=kw+1:kw+numedge_node
                    w = weights[i]
                    par1 = params[kp+1+cnt]
                    par = string("[","$par1","]")
                    write(file_compact, string(rpad(par,pad_len2), rpad("$w",pad_len2),"\n"))
                    if (par1)>=1.0
                        error("Invalid quadrature, problem with edge orbit.")
                    end
                    cnt+=1
                end
                kw+=numedge_node
                kp+=numedge_node
                
                if cub.numS211 != 0
                    write(file_compact, "\nS211 \n")
                end
                cnt=0
                for i=kw+1:kw+numS211_node
                    w = weights[i]
                    par1 = params[kp+1+cnt]./2
                    par2 = params[kp+2+cnt]./2
                    par = string("[","$par1",", ","$par2","]")
                    write(file_compact, string(rpad(par,pad_len2), rpad("$w",pad_len2),"\n"))
                    cnt+=2
                    if (2*par1+par2)>=1.0
                        error("Invalid quadrature, problem with S211 orbit.")
                    end
                end
                kw+=numS211_node
                kp+=2*numS211_node
                
                if cub.numfaceS111 != 0
                    write(file_compact, "\nFaceS111 \n")
                end
                cnt=0
                for i=kw+1:kw+numfaceS111_node
                    w_faceS111 = weights[i]
                    par1 = params[kp+1+cnt]./2
                    par2 = params[kp+2+cnt]./2
                    p_faceS111 = string("[","$par1",", ","$par2","]")
                    write(file_compact, string(rpad(p_faceS111,pad_len2), rpad("$w_faceS111",pad_len2),"\n"))
                    cnt+=2
                    if (par1+par2)>=1.0
                        error("Invalid quadrature, problem with S111 orbit.")
                    end
                end
                kw+=numfaceS111_node
                kp+=2*numfaceS111_node
                
                if cub.facecentroid
                    write(file_compact, "\nFacecentroid \n")
                end
                if facetcentroid_node
                    w = weights[kw+1]
                    par1 = 1.0/3.0
                    par = string("[", "$par1","]")
                    write(file_compact, string(rpad(par,pad_len2), rpad("$w",pad_len2),"\n"))
                    kw+=1
                end

                if cub.numS1111 != 0
                    write(file_compact, "\nS1111 \n")
                end
                cnt=0
                for i=kw+1:kw+numS1111_node
                    w = weights[i]
                    par1 = params[kp+1+cnt]./2
                    par2 = params[kp+2+cnt]./2
                    par3 = params[kp+3+cnt]./2
                    par = string("[","$par1",", ","$par2",", ","$par3","]")
                    write(file_compact, string(rpad(par,pad_len2), rpad("$w",pad_len2),"\n"))
                    cnt+=3
                    if (par1+par2+par3)>=1.0
                        error("Invalid quadrature, problem with S1111 orbit.")
                    end
                end
                kw+=numS1111_node
                kp+=3*numS1111_node

                if cub.centroid
                    write(file_compact, "\nCentroid \n")
                end
                if centroid_node
                    w = weights[end]
                    par1 = 1.0/4.0
                    par = string("[", "$par1","]")
                    write(file_compact, string(rpad(par,pad_len2), rpad("$w",pad_len2),"\n"))
                    kw+=1
                end
            end
            open(file_expanded,"w") do file_expanded
                pad_len = 15
                write(file_expanded, string(rpad("vertices",pad_len), rpad("numS31",pad_len), rpad("midedges",pad_len), 
                        rpad("numS22",pad_len), rpad("numfaceS21",pad_len), rpad("numedge",pad_len), "\n"))
                write(file_expanded, string(rpad("$vert_node (4)",pad_len), rpad("$numS31_node (4)",pad_len), rpad("$midedges_node (6)",pad_len), 
                        rpad("$numS22_node (6)",pad_len), rpad("$numfaceS21_node (12)",pad_len), rpad("$numedge_node (12)",pad_len),"\n"))
                write(file_expanded, string(rpad("numS211",pad_len), rpad("numS1111",pad_len), rpad("facecentroid",pad_len), rpad("numfaceS111",pad_len),rpad("centroid",pad_len),"\n"))
                write(file_expanded, string(rpad("$numS211_node (12)",pad_len),rpad("$numS1111_node (24)",pad_len), rpad("$facetcentroid_node (4)",pad_len), rpad("$numfaceS111_node (24)",pad_len), rpad("$centroid_node (1)",pad_len),"\n\n"))
                
                
                xyz = SymCubatures.calcnodes(cub, vtx)
                ws = SymCubatures.calcweights(cub)

                pad_len = 25
                for i=1:size(xyz,2)
                    x = round(xyz[1,i], digits=16)
                    y = round(xyz[2,i], digits=16)
                    z = round(xyz[3,i], digits=16)
                    w = round(ws[i], digits=16)
                    write(file_expanded, string(rpad(@sprintf("%.16f", x),pad_len), rpad(@sprintf("%.16f", y),pad_len), 
                    rpad(@sprintf("%.16f", z),pad_len), rpad(@sprintf("%.16f", w),pad_len),"\n"))
                end
            end
        else
            xyz = SymCubatures.calcnodes(cub, vtx)
            ws = SymCubatures.calcweights(cub)

            println(" ")
            println(" ")
            println("-------------------------------------------------------------------")
            println("       Tetrahedron quadrature  (q = $q, nnodes = $numnodes) ")
            println("-------------------------------------------------------------------")
            for i = 1:size(xyz)[1]
                println(join(xyz[i,:], ','))
                println(" ")
            end
            println(join(ws',','))
            println("-------------------------------------------------------------------")
        end
        
        # add facet data for the triangles of the tet
        cub, vtx = Cubature.getTriCubatureForTetFaceDiagE(q)
        vert_node = cub.vertices
        midedges_node = cub.midedges
        numS21_node = cub.numS21
        numedge_node = cub.numedge
        numS111_node = cub.numS111
        centroid_node = cub.centroid
        params = cub.params
        weights = cub.weights
        numnodes = cub.numnodes
        qf = convert(Int,ceil(q/2)*2)

        if !print_only
            pad_len = 10
            pad_len2 = 50
            open(file_compact,"a") do file_compact 
                write(file_compact, string("\n\n"))
                write(file_compact, string("===========================================================================================\n"))
                write(file_compact, string("   Facet quadrature data (degree $qf quadrature on a triangle with $numnodes nodes)\n"))
                write(file_compact, string("===========================================================================================\n"))
                write(file_compact, string(rpad("vertices",pad_len), rpad("midedges",pad_len), rpad("numS21",pad_len), rpad("numedge",pad_len), rpad("numS111",pad_len),
                                rpad("centroid",pad_len),"\n"))
                write(file_compact, string(rpad("$vert_node",pad_len), rpad("$midedges_node",pad_len), rpad("$numS21_node",pad_len), rpad("$numedge_node",pad_len), 
                                rpad("$numS111_node",pad_len), rpad("$centroid_node",pad_len),"\n\n"))
                
                kp = 0
                kw = 0
                if cub.vertices
                    write(file_compact, "\nVertices \n")
                end
                if vert_node
                    w_vert = weights[kw+1]
                    write(file_compact, string(rpad("[0.0]",pad_len2), rpad("$w_vert",pad_len2),"\n"))
                    kw+=1
                end

                if cub.midedges
                    write(file_compact, "\nMidedges \n")
                end
                if midedges_node
                    w_mid = weights[kw+1]
                    write(file_compact, string(rpad("[0.5]",pad_len2), rpad("$w_mid",pad_len2),"\n"))
                    kw+=1
                end

                if cub.numS21 != 0
                    write(file_compact, "\nS21 \n")
                end
                cnt=0
                for i=kw+1:kw+numS21_node
                    w_S21 = weights[i]
                    par1 = params[kp+1+cnt]./2
                    p_S21 = string("[", "$par1","]")
                    write(file_compact, string(rpad(p_S21,pad_len2), rpad("$w_S21",pad_len2),"\n"))
                    if (2*par1)>=1.0
                        error("Invalid quadrature, problem with S21 orbit.")
                    end
                    cnt+=1
                end
                kw+=numS21_node
                kp+=numS21_node 
                
                if cub.numedge != 0
                    write(file_compact, "\nEdge \n")
                end
                cnt=0
                for i=kw+1:kw+numedge_node
                    w_edge = weights[i]
                    par1 = params[kp+1+cnt]
                    p_edge = string("[","$par1","]")
                    write(file_compact, string(rpad(p_edge,pad_len2), rpad("$w_edge",pad_len2),"\n"))
                    if (par1)>=1.0
                        error("Invalid quadrature, problem with edge orbit.")
                    end
                    cnt+=1
                end
                kw+=numedge_node
                kp+=numedge_node
                
                if cub.numS111 != 0
                    write(file_compact, "\nS111 \n")
                end
                cnt=0
                for i=kw+1:kw+numS111_node
                    w_S111 = weights[i]
                    par1 = params[kp+1+cnt]./2
                    par2 = params[kp+2+cnt]./2
                    p_S111 = string("[","$par1",", ","$par2","]")
                    write(file_compact, string(rpad(p_S111,pad_len2), rpad("$w_S111",pad_len2),"\n"))
                    cnt+=2
                    if (par1+par2)>=1.0
                        error("Invalid quadrature, problem with S111 orbit.")
                    end
                end
                kw+=numS111_node
                kp+=2*numS111_node

                if cub.centroid
                    write(file_compact, "\nCentroid \n")
                end
                if centroid_node
                    w_cent = weights[end]
                    par1 = 1.0/3.0
                    p_cent = string("[", "$par1","]")
                    write(file_compact, string(rpad(p_cent,pad_len2), rpad("$w_cent",pad_len2),"\n"))
                    kw+=1
                end
            end
            open(file_expanded,"a") do file_expanded
                write(file_expanded, string("\n\n"))
                write(file_expanded, string("===================================================================================\n"))
                write(file_expanded, string("   Facet quadrature data (degree $qf quadrature on a triangle with $numnodes nodes)\n"))
                write(file_expanded, string("===================================================================================\n"))
                pad_len = 12
                write(file_expanded, string(rpad("vertices",pad_len), rpad("midedges",pad_len), rpad("numS21",pad_len), rpad("numedge",pad_len), rpad("numS111",pad_len),
                                rpad("centroid",pad_len),"\n"))
                write(file_expanded, string(rpad("$vert_node (3)",pad_len), rpad("$midedges_node (3)",pad_len), rpad("$numS21_node (3)",pad_len), rpad("$numedge_node (6)",pad_len), 
                                rpad("$numS111_node (6)",pad_len), rpad("$centroid_node (1)",pad_len),"\n\n"))
                
                xy = SymCubatures.calcnodes(cub, vtx)
                ws = SymCubatures.calcweights(cub)

                pad_len = 25
                for i=1:size(xy,2)
                    x = round(xy[1,i], digits=16)
                    y = round(xy[2,i], digits=16)
                    w = round(ws[i], digits=16)
                    write(file_expanded, string(rpad(@sprintf("%.16f", x),pad_len), rpad(@sprintf("%.16f", y),pad_len), rpad(@sprintf("%.16f", w),pad_len),"\n"))
                end
            end
        else
            xy = SymCubatures.calcnodes(cub, vtx)
            ws = SymCubatures.calcweights(cub)
            println(" ")
            println("------------------------------------------------------------------")
            println("      Triangle quadrature for tet (q = $qf, nnodes = $numnodes) ")
            println("------------------------------------------------------------------")
            for i = 1:size(xy)[1]
                println(join(xy[i,:], ','))
                println(" ")
            end
            println(join(ws',','))
            println("-------------------------------------------------------------------")
        end
    end

end

# write_tri_data(quad_degree=nothing, vertices=true, print_only=true)
# write_tet_data(quad_degree=nothing,print_only=true)