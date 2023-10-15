using PlotlyJS
using DelimitedFiles

"""
### SummationByParts.plotly_tri_nodes

This function plots quadruature data on the right triangle using PlotlyJS.jl.

**Inputs**

* `q`: the degree of the operator
* `n`: the number of quadrature nodes
* `facet_type`: the facet quadrature rule type (LG or LGL)
* `x`: the quadrature points in Cartesian coordinates
* `save_fig`: indicates whether to save the plots or not

"""
function plotly_tri_nodes(;q::Int=1,n::Int=-1,facet_type::String="lg",x=[],save_fig::Bool=false)
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
            path_file = string("tri_$facet_type","_q$q","_n$n","_ext.txt")
        else
            path_file = string("tri_$facet_type","_q$q","_n$n","_ext.txt")
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

    p_vert = scatter(x=xvert, y=yvert, mode="lines", line_color="gray", line_width=2)
    p_vol_nodes = scatter(x=xvol, y=yvol, mode="markers", marker_size=8, marker_color="black")
    p_facet_nodes = scatter(x=xfacet, y=yfacet, mode="markers", marker_symbol="circle-open", marker_color="red", marker_size=9)

    layout = Layout(axis = attr(title = "",showgrid = false,showticklabels = false,showbackground = false),
                    yaxis = attr(title = "",showgrid = false,showticklabels = false,showbackground = false),
                    template="simple_white",xaxis_visible=false,yaxis_visible=false,
                    margin=attr(l=0, r=0, b=0, t=0),showlegend=false, aspectratio = Dict(:x => 1, :y => 1))

    p = plot([p_vert,p_vol_nodes,p_facet_nodes], layout)
    display(p)

    if save_fig == true
        path = joinpath(dirname(current_dir), "visual/tri/$facet_type/")
        file_name = string("tri_$facet_type","_q$q","_n$n",".png")
        file= joinpath(path,file_name)
        savefig(p,file,width=400,height=400,scale=1)
    end

    return 
end

"""
### SummationByParts.plotly_tet_nodes

This function plots quadruature data on the right tetrahedron using PlotlyJS.jl.

**Inputs**

* `q`: the degree of the operator
* `n`: the number of quadrature nodes
* `save_fig`: indicates whether to save the plots or not

"""
function plotly_tet_nodes(;q::Int=1,n::Int=-1,save_fig::Bool=false)
    dir=""
    path_file=""

    current_file_path = @__FILE__
    current_dir = dirname(current_file_path)
    dir = joinpath(dirname(current_dir),"quadrature_data/tet/expanded/")

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
        path_file = string("tet","_q$q","_n$n","_ext.txt")
    else
        path_file = string("tet","_q$q","_n$n","_ext.txt")
    end
    path = joinpath(dir,path_file)

    lines = readdlm(path)
    xvert = Array([-1,1,-1,-1,-1,-1,1,-1])
    yvert = Array([-1,-1,1,-1,-1,1,-1,-1])
    zvert = Array([-1,-1,-1,-1,1,-1,-1,1])
    xvol = convert(Array{Float64},lines[5:n+4,1])
    yvol = convert(Array{Float64},lines[5:n+4,2])
    zvol = convert(Array{Float64},lines[5:n+4,3])

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
    
    vertices = [[-1.0, -1.0, -1.0], 
                [1.0, -1.0, -1.0],  
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0]]

    # Define the faces of the tetrahedron (as indices of vertices)
    faces = [[1, 2, 3],
             [1, 2, 4],
             [2, 3, 4],
             [1, 3, 4]]

    # Define colors for each face
    colors = ["yellow", "rgba(0,0,0,0)", "blue", "green"] #"rgba(0,0,0,0)"
    # colors = ["gray", "rgba(0,0,0,0)", "blue", "gray"] #"rgba(0,0,0,0)"
    traces = mesh3d(x=[v[1] for v in vertices],
                    y=[v[2] for v in vertices],
                    z=[v[3] for v in vertices],
                    i=[f[1] - 1 for f in faces],  # Adjust indices to start from 0
                    j=[f[2] - 1 for f in faces],
                    k=[f[3] - 1 for f in faces],
                    opacity=0.3,  # Adjust opacity for visibility
                    facecolor=colors)

    #-------- Option to color each facet separately ---------
    # trace1 = PlotlyJS.mesh3d(x=[v[1] for v in vertices],
    #         y=[v[2] for v in vertices],
    #         z=[v[3] for v in vertices],
    #         i=[[f[1] - 1 for f in faces][1]],  
    #         j=[[f[2] - 1 for f in faces][1]],
    #         k=[[f[3] - 1 for f in faces][1]],
    #         opacity=0.9,  # Adjust opacity for visibility
    #         color=colors[1])
    # trace2 = PlotlyJS.mesh3d(x=[v[1] for v in vertices],
    #         y=[v[2] for v in vertices],
    #         z=[v[3] for v in vertices],
    #         i=[[f[1] - 1 for f in faces][2]],  
    #         j=[[f[2] - 1 for f in faces][2]],
    #         k=[[f[3] - 1 for f in faces][2]],
    #         opacity=0.3,  # Adjust opacity for visibility
    #         color=colors[2])
    # trace3 = PlotlyJS.mesh3d(x=[v[1] for v in vertices],
    #         y=[v[2] for v in vertices],
    #         z=[v[3] for v in vertices],
    #         i=[[f[1] - 1 for f in faces][3]],  
    #         j=[[f[2] - 1 for f in faces][3]],
    #         k=[[f[3] - 1 for f in faces][3]],
    #         opacity=0.3,  # Adjust opacity for visibility
    #         color=colors[3])
    # trace4 = PlotlyJS.mesh3d(x=[v[1] for v in vertices],
    #         y=[v[2] for v in vertices],
    #         z=[v[3] for v in vertices],
    #         i=[[f[1] - 1 for f in faces][4]],  
    #         j=[[f[2] - 1 for f in faces][4]],
    #         k=[[f[3] - 1 for f in faces][4]],
    #         opacity=0.9,  # Adjust opacity for visibility
    #         color=colors[4])

    # traces = [trace1, trace2, trace3, trace4]
    #----------------

    p_vert = scatter(x=xvert, y=yvert, z=zvert, type="scatter3d", mode="lines", line_color="black", line_width=2)
    p_vol_nodes = scatter(x=xvol, y=yvol, z=zvol, type="scatter3d", mode="markers", 
                          marker=attr(color="black",size=4,symbol="circle"))
    p_facet_nodes = scatter(x=xfacet, y=yfacet, z=zfacet, type="scatter3d", mode="markers", 
                            marker=attr(color="red",size=5,symbol="circle-open"))
    layout = Layout(scene = attr(template="simple_white",xaxis_visible=false,yaxis_visible=false,zaxis_visible=false,       
                            showlegend=false), showlegend=false)

    # p = plot([p_vert,p_vol_nodes,p_facet_nodes, trace1, trace2, trace3, trace4],layout)
    p = plot([p_vert,p_vol_nodes,p_facet_nodes,traces],layout)
    camera = attr(eye=attr(x=1.5, y=-1.5, z=0.2))
    relayout!(p, scene_camera=camera)
    
    display(p)

    if save_fig == true
        path = joinpath(dirname(current_dir), "visual/tet/")
        file_name = string("tet","_q$q","_n$n",".png")
        file= joinpath(path,file_name)
        # savefig(p,file)
        savefig(p,file,width=3000,height=3000,scale=4)
    end

    return 
end

# Save all figures
plot_all = false
save_fig = false
if plot_all==true
    for i=1:20  
        plotly_tri_nodes(q=i,facet_type="lg",save_fig=save_fig)
        plotly_tri_nodes(q=i,facet_type="lgl",save_fig=save_fig)
    end
    for i=1:10
        plotly_tet_nodes(q=i,save_fig=save_fig)
    end
end
# plotly_tet_nodes(q=10,n=145,save_fig=save_fig)
# plotly_tri_nodes(q=10,facet_type="lg",save_fig=save_fig)
