export get_edge_plot, get_arrow_plot, get_node_plot, get_nlabel_plot, get_elabel_plot, compute_auto_label_aligns

"Get the `EdgePlot` subplot from a `GraphPlot`."
get_edge_plot(gp::GraphPlot) = gp.edge_plot[]

"Get the scatter plot of the arrow heads from a `GraphPlot`."
get_arrow_plot(gp::GraphPlot) = gp.arrow_plot[]

"Get the scatter plot of the nodes from a `GraphPlot`."
get_node_plot(gp::GraphPlot) = gp.node_plot[]

"Get the text plot of the node labels from a `GraphPlot`."
get_nlabel_plot(gp::GraphPlot) = haskey(gp.attributes, :nlabels_plot) ? gp[:nlabels_plot][] : nothing
"Get the text plot of the edge labels from a `GraphPlot`."
get_elabel_plot(gp::GraphPlot) = haskey(gp.attributes, :elabels_plot) ? gp[:elabels_plot][] : nothing

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
getattr(o::Union{Observable,Makie.Computed}, idx, default=nothing) = getattr(o[], idx, default)

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
    prep_vertex_attributes(attr, graph::AbstractGraph, default_value)

Prepare the vertex attributes to be forwarded to the internal recipes.
If the attribute is a `Vector` or single value forward it as is (or the `default_value` if isnothing).
If it is an `AbstractDict` expand it to a `Vector` using vertex indices.
"""
function prep_vertex_attributes(attr, graph::AbstractGraph, default_value=nothing)
    if issingleattribute(attr)
        isnothing(attr) ? default_value : attr
    elseif attr isa AbstractVector
        attr
    else
        [getattr(attr, i, default_value) for i in vertices(graph)]
    end
end

"""
    prep_edge_attributes(attr, graph::AbstractGraph, default_value)

Prepare the edge attributes to be forwarded to the internal recipes.
If the attribute is a `Vector` or single value forward it as is (or the `default_value` if isnothing).
If it is an `AbstractDict` expand it to a `Vector` using edge indices.
"""
function prep_edge_attributes(attr, graph::AbstractGraph, default_value=nothing)
    if issingleattribute(attr)
        isnothing(attr) ? default_value : attr
    elseif attr isa AbstractVector
        attr
    else
        [getattr(attr, i, default_value) for i in getedgekeys(graph, attr)]
    end
end

"""
    issingleattribute(x)

Return `true` if `x` represents a single attribute value
"""
issingleattribute(x) = isa(x, Point) || (!isa(x, AbstractVector) && !isa(x, AbstractDict))

"""
    to_pointf32(p::Point{N, T})

Convert Point{N, T} or NTuple{N, T} to Point{N, Float32}.
"""
to_pointf32(p::Union{Point{N,T}, NTuple{N,T}}) where {N,T} = Point{N, Float32}(p)
to_pointf32(p::StaticVector{N, T}) where {N,T} = Point{N, Float32}(p)
to_pointf32(p::Vararg{T,N}) where {N,T} = Point{N, Float32}(p)
to_pointf32(p::Vector{T}) where {T} = Point{length(p), Float32}(p)

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

`size1` / `size2` may be scalars or `(width, height)` / `Vec2` marker sizes;
non-scalar sizes use `maximum` (larger axis / circumscribed extent).
"""
function distance_between_markers(marker1, size1, marker2, size2)
    marker1_scale = scale_factor(marker1)
    marker2_scale = scale_factor(marker2)
    d = marker1_scale * maximum(size1) / 2 +
        marker2_scale * maximum(size2) / 2

    return d
end

"""
    point_near_offset(edge_path, p0::PT, d, to_px, offset) where {PT}

Find point near the `offset ∈ [0, 1]` along `edge_path` a 
distance `d` pixels along the tangent line.
"""
function point_near_offset(edge_path, p0::PT, d, to_px, offset) where {PT}
    pt = tangent(edge_path, offset) #edge tangent along path
    r = to_px(pt) - to_px(PT(0)) #direction vector in pixels
    scale_px = 1 ./ (to_px(PT(1)) - to_px(PT(0)))
    p1 = p0 - d*normalize(r)*scale_px

    return p1
end

"""
    compute_auto_label_aligns(g::AbstractGraph, node_positions::AbstractVector)

Compute optimal label alignment for each node to avoid edge overlaps.

For each node, finds the largest angular gap between incident edges and places
the label in the middle of that gap.

# Arguments
- `g`: An `AbstractGraph` from Graphs.jl
- `node_positions`: Vector of node positions (Point2f or similar)

# Returns
- Vector of `(:horizontal, :vertical)` tuples for `nlabels_align`

# Notes
- For isolated nodes (no incident edges), returns `(:right, :bottom)` as default
- Works with both directed and undirected graphs
- Accounts for both incoming and outgoing edges
"""
function compute_auto_label_aligns(g::AbstractGraph, node_positions::AbstractVector)
    n = nv(g)
    aligns = Vector{Tuple{Symbol, Symbol}}(undef, n)

    for node in 1:n
        node_pos = node_positions[node]
        edge_angles = Float64[]

        # One angle per distinct neighbour (Graphs.jl all_neighbors); skip self-loops.
        for nbr in all_neighbors(g, node)
            nbr == node && continue
            nbr_pos = node_positions[nbr]
            dx = Float64(nbr_pos[1] - node_pos[1])
            dy = Float64(nbr_pos[2] - node_pos[2])
            push!(edge_angles, atan(dy, dx))
        end

        # Default alignment for isolated nodes
        if isempty(edge_angles)
            aligns[node] = (:right, :bottom)
            continue
        end
        
        # Find the largest angular gap between adjacent edges
        # Normalize angles (in radians) to [0, 2π] and sort
        angles_norm = sort([mod(a + 2π, 2π) for a in edge_angles])
        
        # Find gaps between adjacent angles (including wrap-around gap)
        best_gap = 0.0
        best_midpoint = -π/4  # Default: Southeast (bottom-right)
        
        for i in 1:length(angles_norm)
            next_i = i == length(angles_norm) ? 1 : i + 1
            # Gap from angles_norm[i] to angles_norm[next_i]
            if next_i == 1
                # Wrap-around gap
                gap = (2π - angles_norm[i]) + angles_norm[1]
                midpoint = angles_norm[i] + gap / 2
                if midpoint > π
                    midpoint -= 2π
                end
            else
                gap = angles_norm[next_i] - angles_norm[i]
                midpoint = angles_norm[i] + gap / 2
            end
            
            if gap > best_gap
                best_gap = gap
                best_midpoint = midpoint
            end
        end
        
        # Convert best_midpoint angle to alignment
        # Normalize to [0, 2π]
        angle = mod(best_midpoint + 2π, 2π)
        
        # Map angle to alignment (8 directions)
        # Each sector spans π/4 (45°), centered on cardinal/ordinal directions
        # Sector boundaries: 0=E, π/4=NE, π/2=N, 3π/4=NW, π=W, 5π/4=SW, 3π/2=S, 7π/4=SE
        #
        # IMPORTANT: Alignment semantics in Makie:
        # - (:left, :center) means label's LEFT edge at anchor → label extends RIGHT (East)
        # - (:right, :center) means label's RIGHT edge at anchor → label extends LEFT (West)
        # - (:center, :bottom) means label's BOTTOM at anchor → label extends UP (North)
        # - (:center, :top) means label's TOP at anchor → label extends DOWN (South)
        aligns[node] = if angle < π/8 || angle >= 15π/8
            (:left, :center)       # Label to East (right of node)
        elseif angle < 3π/8
            (:left, :bottom)       # Label to Northeast
        elseif angle < 5π/8
            (:center, :bottom)     # Label to North (above node)
        elseif angle < 7π/8
            (:right, :bottom)      # Label to Northwest
        elseif angle < 9π/8
            (:right, :center)      # Label to West (left of node)
        elseif angle < 11π/8
            (:right, :top)         # Label to Southwest
        elseif angle < 13π/8
            (:center, :top)        # Label to South (below node)
        else
            (:left, :top)          # Label to Southeast
        end
    end
    
    return aligns
end
