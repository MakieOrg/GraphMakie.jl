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
    ps = filter(p -> p isa Makie.Text && p[1][] == labels, gp.plots)
    if isempty(ps)
        return nothing
    elseif length(ps) == 1
        return ps[1]
    else
        error("Could not determine plot $ps")
    end
end

"""
    getattr(o::Observable, idx)

If observable wraps an AbstractVector or AbstractDict return
the value at idx. If dict has no key idx rerturn nothing.
Else return the one and only element.
"""
function getattr(o::Observable, idx)
    if o[] isa AbstractVector && !isa(o[], Point)
        return o[][idx]
    elseif o[] isa AbstractDict
        return get(o[], idx, nothing)
    else
        return o[]
    end
end

"""
    Pointf0(p::Point{N, T})

Convert Point{N, T} or NTuple{N, T} to Point{N, Float32}.
"""
Pointf0(p::Union{Point{N,T}, NTuple{N,T}}) where {N,T} = Point{N, Float32}(p)
Pointf0(p::Vararg{T,N}) where {N,T} = Point{N, Float32}(p)

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
    return Point2f0(x/norm, y/norm)
end

function plot_controlpoints!(ax::Axis, gp::GraphPlot)
    ep = get_edge_plot(gp)
    ep.plots[1] isa BezierSegments || return
    ep = ep.plots[1]
    paths = ep[:paths][]

    for (i, p) in enumerate(paths)
        p isa Line && continue
        color = getattr(gp.edge_color, i)
        for (j, c) in enumerate(p.commands)
            if c isa CurveTo
                segs = [p.commands[j-1].p, c.c1, c.p, c.c2]
                linesegments!(ax, segs; color, linestyle=:dot)
                scatter!(ax, [c.c1, c.c2]; color)
            end
        end
    end
end
