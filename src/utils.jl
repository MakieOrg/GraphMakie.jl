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
    @assert p.marker[] == '➤'
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
    scale_factor(marker)

Get base size (scaling) in pixels for `marker`.
"""
function scale_factor(marker)
    if marker in [Circle, Rect, Arrow]
        r = 1
    else #based of Makie.default_marker_map(), all of the markers have scale_factor = 0.75
        r = 0.75
    end

    return r
end

"""
    distance_between_markers(marker1, size1, marker2, size2, θ)

Calculate distance between 2 markers at an angle θ.
TODO: Implement for noncircular marker1.
"""
function distance_between_markers(marker1, size1, marker2, size2, θ)
    #NOTE: If marker1 is Circle or :circle, θ has no impact on the distance.
    marker1_scale = scale_factor(marker1)
    marker2_scale = scale_factor(marker2)
    d = marker1_scale*size1/2 + marker2_scale*size2/2

    return d
end