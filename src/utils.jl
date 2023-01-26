export get_edge_plot, get_arrow_plot, get_node_plot, get_nlabel_plot, get_elabel_plot

"Get the `EdgePlot` subplot from a `GraphPlot`."
function get_edge_plot(gp::GraphPlot)
    p = gp.plots[1]
    @assert p isa EdgePlot
    return p
end

"Get the scatter plot of the arrow heads from a `GraphPlot`."
function get_arrow_plot(gp::GraphPlot)
    p = gp.plots[2]
    @assert p isa Scatter
    @assert p.marker[] == Arrow#'➤'
    return p
end

"Get the scatter plot of the nodes from a `GraphPlot`."
function get_node_plot(gp::GraphPlot)
    p = gp.plots[3]
    @assert p isa Scatter
    @assert p[1][] == gp[:node_pos][]
    return p
end

"Get the text plot of the node labels from a `GraphPlot`."
get_nlabel_plot(gp::GraphPlot) = _get_label_plot(gp, gp.nlabels[])
"Get the text plot of the edge labels from a `GraphPlot`."
get_elabel_plot(gp::GraphPlot) = _get_label_plot(gp, gp.elabels[])

function _get_label_plot(gp::GraphPlot, labels)
    ps = filter(p -> p isa Makie.Text && p[:text][] == labels, gp.plots)
    if isempty(ps)
        return nothing
    elseif length(ps) == 1
        return ps[1]
    else
        error("Could not determine plot $ps")
    end
end

"""
    getattr(o::Observable, idx, default=nothing)

If observable wraps an AbstractVector or AbstractDict return
the value at idx. If dict has no key idx rerturn default.
Else return the one and only element.
"""
function getattr(o::Observable, idx, default=nothing)
    if o[] isa AbstractVector && !isa(o[], Point)
        return o[][idx]
    elseif o[] isa AbstractDict
        return get(o[], idx, default)
    else
        return o[] === nothing ? default : o[]
    end
end

"""
    Pointf(p::Point{N, T})

Convert Point{N, T} or NTuple{N, T} to Point{N, Float32}.
"""
Pointf(p::Union{Point{N,T}, NTuple{N,T}}) where {N,T} = Point{N, Float32}(p)
Pointf(p::StaticVector{N, T}) where {N,T} = Point{N, Float32}(p)
Pointf(p::Vararg{T,N}) where {N,T} = Point{N, Float32}(p)
Pointf(p::Vector{T}) where {T} = Point{length(p), Float32}(p)

"""
    align_to_dir(align::Tuple{Symbol, Symbol})

Given a tuple of alignment (i.e. `(:left, :bottom)`) return a normalized
2d vector which points in the direction of the offset.
"""
function align_to_dir(align::Tuple{Symbol, Symbol})
    halign, valign = align

    x = 0.0
    if halign === :left
        x = 1.0
    elseif halign === :right
        x = -1.0
    end

    y = 0.0
    if valign === :top
        y = -1.0
    elseif valign === :bottom
        y = 1.0
    end
    norm = x==y==0.0 ? 1 : sqrt(x^2 + y^2)
    return Point2f(x/norm, y/norm)
end

"""
    plot_controlpoints!(ax::Axis, gp::GraphPlot)
    plot_controlpoints!(ax::Axis, path::BezierPath)

Add all the bezier controlpoints of graph plot or a single
path to the axis `ax`.
"""
function plot_controlpoints!(ax::Axis, gp::GraphPlot)
    ep = get_edge_plot(gp)
    ep.plots[1] isa BezierSegments || return
    ep = ep.plots[1]
    paths = ep[:paths][]

    for (i, p) in enumerate(paths)
        p isa Line && continue
        color = getattr(gp.edge_color, i)
        plot_controlpoints!(ax, p; color)
    end
end

function plot_controlpoints!(ax::Axis, p::BezierPath; color=:black)
    for (j, c) in enumerate(p.commands)
        if c isa CurveTo
            segs = [p.commands[j-1].p, c.c1, c.p, c.c2]
            linesegments!(ax, segs; color, linestyle=:dot)
            scatter!(ax, [c.c1, c.c2]; color)
        end
    end
end

"""
    move_arrows_to_nodes!(ax::Axis, gp::GraphPlot; t=0.9)

Moves arrowheads to the surface of the each destination node.
Only supported for markers of type `Circle`.

Call this function only after all changes have been made to the plot.
"""
function move_arrows_to_nodes!(ax::Axis, gp::GraphPlot; t=1)
    if gp.node_marker[] !== Circle
        error("`move_arrows_to_nodes! is only supported for plots that use `Cricle` as the `node_marker`.")
    end

    #get point to pixel scale
    xlims = ax.xaxis.attributes.limits[]
    xrange = ax.xaxis.attributes.endpoints[]
    ylims = ax.yaxis.attributes.limits[]
    yrange = ax.yaxis.attributes.endpoints[]
    dx = xlims[2] - xlims[1]
    dxpx = xrange[2][1] - xrange[1][1]
    dy = ylims[2] - ylims[1]
    dypx = yrange[2][2] - yrange[1][2]
    dpx = [dx/dxpx, dy/dypx]

    #node attr
    node_pos = gp.node_pos[]
    node_size = gp.node_size[]
    node_rad = node_size ./ 2 #radius

    #arrow attr
    arrow_pos = gp.arrow_pos[]
    arrow_size = gp.arrow_size[]
    arrow_rot = gp.arrow_rot[]
    arrow_rad = arrow_size ./ 2 #radius

    #scene and projection
    sc = Makie.parent_scene(gp)
    to_px(point) = project(sc, point)

    #edge and graph
    edge_paths = gp.edge_paths[]
    g = gp.graph[]

    for (i,e) in enumerate(edges(g))
        #update arrow rotation
        e_tan = tangent(edge_paths[i], t) #tangent at destination node
        e_tan_px = to_px(e_tan) - to_px(Point2(0,0)) #project to pixels
        θ = atan(e_tan_px[2], e_tan_px[1]) #tangent angle
        arrow_rot.rotation[i] = θ

        #update arrow position
        j = dst(e)
        d = node_rad[1] + arrow_rad[1] #distance between center of node and center of arrow
        p0 = node_pos[j]
        p1 = p0 .- d * [cos(θ),sin(θ)] .* dpx
        arrow_pos[i] = p1
    end
    
    gp.arrow_pos[] = gp.arrow_pos[]
    gp.arrow_rot[] = gp.arrow_rot[]

    return nothing
end