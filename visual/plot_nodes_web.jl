using Plots; gr()
using DelimitedFiles

"""
### SummationByParts.plotly_tri_nodes

This function plots quadruature data on the right triangle using Plot.jl.

**Inputs**

* `q`: the degree of the operator
* `n`: the number of quadrature nodes
* `facet_type`: the facet quadrature rule type (LG or LGL)
* `x`: the quadrature points in Cartesian coordinates

"""
function plot_tri_nodes(;q::Int=1,n::Int=-1,facet_type::String="lg",x=[])
    dir=""
    path_file=""

    current_file_path = @__FILE__
    current_dir = dirname(current_file_path)
    dir = joinpath(dirname(current_dir),"quadrature_data/tri/expanded/")

    xvol = []
    yvol = []
    if x==[]
        if n==-1     
            path_file = string("tri_$facet_type","_q$q")

            # Get a list of files in the directory
            files = readdir(dir)
            # Filter files based on the provided half name
            matching_files = filter(x -> occursin(path_file, x), files)
            
            ns=[]
            for i = 1:length(matching_files)
                file_name = matching_files[i]
                split_name = split(file_name,"_n")
                split_name2 = split(split_name[2],"_")
                push!(ns,parse(Int,split_name2[1]))
            end
            n = minimum(ns)
            path_file = string("tri_$facet_type","_q$q","_n$n","_ext.dat")
        else
            path_file = string("tri_$facet_type","_q$q","_n$n","_ext.dat")
        end
        path = joinpath(dir,path_file)

        lines = readdlm(path)
        xvol = convert(Array{Float64},lines[3:end,1])
        yvol = convert(Array{Float64},lines[3:end,2])
    else
        xvol = x[1,:]
        yvol = x[2,:]
    end
    xvol = convert(Array{Float64},xvol)
    yvol = convert(Array{Float64},yvol)

    xvert = Array([-1,1,-1,-1])
    yvert = Array([-1,-1,1,-1])
    xfacet = []
    yfacet = []
    for i = 1:length(xvol)
        if (abs(xvol[i] + yvol[i]-0.0) < 1e-13 || abs(xvol[i]+1.0)<1e-13 || abs(yvol[i]+1.0)<1e-13)
            push!(xfacet,xvol[i])
            push!(yfacet,yvol[i])
        end
    end

    Plots.plot(xvert,yvert, seriestype=:path, linestyle=:solid, lc="black", lw=1.5, label="")
    Plots.scatter!(xvol,yvol,m = (3.5, :black),label="volume nodes") 
    label_facet=""
    if xfacet!=[]
        label_facet="volume nodes on facet"
    end
    p=Plots.scatter!(xfacet, yfacet; #aspect_ratio=:equal, grid=false,  xaxis=false, yaxis=false,
    ms=4.5, markercolor=:transparent, markerstrokecolor=:red, markerstrokewidth = 2,  label=label_facet, 
    legend=true, legendposition=:top, legendfontsize=8)

    display(p)

    return 
end

"""
### SummationByParts.plotly_tet_nodes

This function plots quadruature data on the right tetrahedron using Plot.jl.

**Inputs**

* `q`: the degree of the operator
* `n`: the number of quadrature nodes

"""
function plot_tet_nodes(;q::Int=1,n::Int=-1,x=[])
    dir=""
    path_file=""

    current_file_path = @__FILE__
    current_dir = dirname(current_file_path)
    dir = joinpath(dirname(current_dir),"quadrature_data/tet/expanded/")

    xvol = []
    yvol = []
    zvol = []
    if x==[]
        if n==-1
            path_file = string("tet","_q$q")
            
            # Get a list of files in the directory
            files = readdir(dir)
            # Filter files based on the provided half name
            matching_files = filter(x -> occursin(path_file, x), files)
            
            ns=[]
            for i = 1:length(matching_files)
                file_name = matching_files[i]
                split_name = split(file_name,"_n")
                split_name2 = split(split_name[2],"_")
                push!(ns,parse(Int,split_name2[1]))
            end
            n = minimum(ns)
            path_file = string("tet","_q$q","_n$n","_ext.dat")
        else
            path_file = string("tet","_q$q","_n$n","_ext.dat")
        end
        path = joinpath(dir,path_file)

        lines = readdlm(path)
        xvol = convert(Array{Float64},lines[5:n+4,1])
        yvol = convert(Array{Float64},lines[5:n+4,2])
        zvol = convert(Array{Float64},lines[5:n+4,3])

    else
        xvol = x[1,:]
        yvol = x[2,:]
        zvol = x[3,:]
    end

    xvert = Array([-1,1,-1,-1,-1,-1,1,-1])
    yvert = Array([-1,-1,1,-1,-1,1,-1,-1])
    zvert = Array([-1,-1,-1,-1,1,-1,-1,1])
    xfacet = []
    yfacet = []
    zfacet = []

    for i = 1:length(xvol)
        if (abs(xvol[i] + yvol[i] + zvol[i]+1.0) < 1e-13 || abs(xvol[i]+1.0)<1e-13 || abs(yvol[i]+1.0)<1e-13 || abs(zvol[i]+1.0)<1e-13)
            push!(xfacet,xvol[i])
            push!(yfacet,yvol[i])
            push!(zfacet,zvol[i])
        end
    end
    
    Plots.plot(xvert, yvert, zvert, linestyle=:solid, lc="black", lw=1.5, label="")
    Plots.scatter!(xvol, yvol, zvol, m = (3.5, :black),label="volume nodes")
    label_facet=""
    if xfacet!=[]
        label_facet="volume nodes on facet"
    end
    p=Plots.scatter!(xfacet, yfacet, zfacet; #aspect_ratio=:equal, grid=false,  xaxis=false, yaxis=false, zaxis=false,
    ms=3.5, markercolor=:transparent, markerstrokecolor=:red, markerstrokewidth = 2,  label=label_facet, 
    legend=true, legendposition=:top, legendfontsize=8)

    display(p)

    return 
end

# plot_tet_nodes(q=10,n=145)
plot_tet_nodes(q=3)
# plot_tri_nodes(q=20,facet_type="lg")