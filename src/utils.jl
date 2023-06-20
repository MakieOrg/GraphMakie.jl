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
    getedgekeys(gr::G, edgedat::D) where {G<:AbstractGraph, K<:AbstractEdge, D<:AbstractDict{K}, IsDirected{G}}

Return enumeration of edges for directed graph
"""
@traitfn function getedgekeys(gr::G, edgedat::D) where {G<:AbstractGraph, K<:AbstractEdge, D<:AbstractDict{K}; IsDirected{G}}
    return edges(gr)
end

"""
    getedgekeys(gr::G, edgedat::D) where {G<:AbstractGraph, K<:AbstractEdge, D<:AbstractDict{K}, IsDirected{G}}

Return enumeration of edges for undirected graph such that the user's keys are used

# Extended help
Wraps the `edges()` method such that the edges are referenced as the user defined them in the dictionary.
"""
@traitfn function getedgekeys(gr::G, edgedat::D) where {G<:AbstractGraph, K<:AbstractEdge, D<:AbstractDict{K}; !IsDirected{G}}
    Iterators.map(e -> reverse(e) ∈ keys(edgedat) ? reverse(e) : e , edges(gr))
end

"""
    getedgekeys(gr::AbstractGraph, <:AbstractDict{AbstractEdge})

Return enumeration of edge indices
"""
getedgekeys(gr::AbstractGraph, _) = 1:ne(gr)

"""
    getattr(o::Observable, idx, default=nothing)

If observable wraps an AbstractVector or AbstractDict return
the value at idx. If dict has no key idx returns default.
Else return the one and only element.
"""
getattr(o::Observable, idx, default=nothing) = getattr(o[], idx, default)

"""
    getattr(x, idx, default=nothing)

If `x` wraps an AbstractVector or AbstractDict return
the value at idx. If dict has no key idx return default.
Else return the one and only element.
"""
function getattr(x, idx, default=nothing)
    if x isa AbstractVector && !isa(x, Point)
        return x[idx]
    elseif x isa DefaultDict || x isa DefaultOrderedDict
        return getindex(x, idx)
    elseif x isa AbstractDict
        return get(x, idx, default)
    else
        return x === nothing ? default : x
    end
end

"""
    prep_vertex_attributes(oattr::Observable, ograph::Observable{<:AbstractGraph}, odefault::Observable)

Prepare the vertex attributes to be forwarded to the internal recipes.
If the attribute is a `Vector` or single value forward it as is (or the `odefault` value if isnothing).
If it is an `AbstractDict` expand it to a `Vector` using `indices`.
"""
function prep_vertex_attributes(oattr::Observable, ograph::Observable{<:AbstractGraph}, odefault::Observable=Observable(nothing))
    @lift begin
        if issingleattribute($oattr)
            isnothing($oattr) ? $odefault : $oattr
        elseif $oattr isa AbstractVector
            $oattr
        else
            [getattr($oattr, i, $odefault) for i in vertices($ograph)]
        end
    end
end

"""
    prep_edge_attributes(oattr::Observable, ograph::Observable{<:AbstractGraph}, odefault::Observable)

Prepare the edge attributes to be forwarded to the internal recipes.
If the attribute is a `Vector` or single value forward it as is (or the `odefault` value if isnothing).
If it is an `AbstractDict` expand it to a `Vector` using `indices`.
"""
function prep_edge_attributes(oattr::Observable, ograph::Observable{<:AbstractGraph}, odefault::Observable=Observable(nothing))
    @lift begin
        if issingleattribute($oattr)
            isnothing($oattr) ? $odefault : $oattr
        elseif $oattr isa AbstractVector
            $oattr
        else
            [getattr($oattr, i, $odefault) for i in getedgekeys($ograph, $oattr)]
        end
    end
end

"""
    issingleattribute(x)

Return `true` if `x` represents a single attribute value
"""
issingleattribute(x) = isa(x, Point) || (!isa(x, AbstractVector) && !isa(x, AbstractDict))

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
scale_factor(marker) = 1 #for 1x1 base sizes (Circle, Rect, Arrow)
function scale_factor(marker::Char)
    if marker == '➤'
        d = 0.675
    else
        d = 0.705 #set to the same value as :circle, but is really dependent on the Char
    end

    return d
end
function scale_factor(marker::Symbol)
    size_factor = 0.75 #Makie.default_marker_map() has all markers scaled by 0.75
    if marker == :circle #BezierCircle
        r = 0.47
    elseif marker in [:rect, :diamond, :vline, :hline] #BezierSquare
        rmarker = 0.95*sqrt(pi)/2/2
        r = sqrt(2*rmarker^2) #pithagoras to get radius of circle that circumscribes marker
    elseif marker in [:utriangle, :dtriangle, :ltriangle, :rtriangle] #Bezier Triangles
        r = 0.97/2
    elseif marker in [:star4, :star5, :star6, :star8] #Bezier Stars
        r = 0.6
    else #Bezier Crosses/Xs and Ngons
        r = 0.5
    end

    return 2*r*size_factor #get shape diameter
end

"""
    distance_between_markers(marker1, size1, marker2, size2)

Calculate distance between 2 markers.
TODO: Implement for noncircular marker1.
      (will require angle for line joining the 2 markers).
"""
function distance_between_markers(marker1, size1, marker2, size2)
    marker1_scale = scale_factor(marker1)
    marker2_scale = scale_factor(marker2)
    d = marker1_scale*size1/2 + marker2_scale*size2/2

    return d
end

"""
    point_near_dst(edge_path, p0::PT, d, to_px) where {PT}

Find point near destination node along `edge_path` a 
distance `d` pixels along the tangent line.
"""
function point_near_dst(edge_path, p0::PT, d, to_px) where {PT}
    pt = tangent(edge_path, 1) #edge tangent at dst node
    r = to_px(pt) - to_px(PT(0)) #direction vector in pixels
    scale_px = 1 ./ (to_px(PT(1)) - to_px(PT(0)))
    p1 = p0 - d*normalize(r)*scale_px

    return p1
end
