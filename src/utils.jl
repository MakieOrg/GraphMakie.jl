export get_edge_plot, get_arrow_plot, get_node_plot, get_nlabel_plot, get_elabel_plot

function get_edge_plot(gp::GraphPlot)
    p = gp.plots[1]
    @assert p isa BezierSegments
    return p
end

function get_arrow_plot(gp::GraphPlot)
    p = gp.plots[2]
    @assert p isa Scatter
    @assert p.marker[] == '➤'
    return p
end

function get_node_plot(gp::GraphPlot)
    p = gp.plots[3]
    @assert p isa Scatter
    @assert p[1][] == gp[:node_pos][]
    return p
end

get_nlabel_plot(gp::GraphPlot) = _get_label_plot(gp, gp.nlabels[])
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
the value at idx. Else return the one and only element.
"""
function getattr(o::Observable, idx)
    if o[] isa AbstractVector && !isa(o[], Point) || o[] isa AbstractDict
        return o[][idx]
    else
        return o[]
    end
end
