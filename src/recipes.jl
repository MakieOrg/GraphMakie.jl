using LinearAlgebra: normalize, ⋅, norm
using NetworkLayout: AbstractLayout, dim
export GraphPlot, graphplot, graphplot!, Arrow

const Arrow = Makie.Polygon(Point2f.([(-0.5,-0.5),(0.5,0),(-0.5,0.5),(-0.25,0)]))

"""
    graphplot(graph::AbstractGraph)
    graphplot!(ax, graph::AbstractGraph)

Creates a plot of the network `graph`. Consists of multiple steps:
- Layout the nodes: see `layout` attribute. The node position is accessible from outside
  the plot object `p` as an observable using `p[:node_pos]`.
- plot edges as `edgeplot`-plot
- if `arrow_show` plot arrowheads as `scatter`-plot
- plot nodes as `scatter`-plot
- if `nlabels!=nothing` plot node labels as `text`-plot
- if `elabels!=nothing` plot edge labels as `text`-plot

The main attributes for the subplots are exposed as attributes for `graphplot`.
Additional attributes for the `scatter`, `edgeplot` and `text` plots can be provided
as a named tuples to `node_attr`, `edge_attr`, `nlabels_attr` and `elabels_attr`.

Most of the arguments can be either given as a vector of length of the
edges/nodes or as a single value. One might run into errors when changing the
underlying graph and therefore changing the number of Edges/Nodes.

## Attributes
### Main attributes
- `layout=Spring()`: function `AbstractGraph->Vector{Point}` or `Vector{Point}` that determines the base layout.  Can also be any network layout from [NetworkLayout.jl](https://github.com/JuliaGraphs/NetworkLayout.jl), like `Spring`, `Stress`, `Spectral`, etc.
- `node_color=automatic`:
  Defaults to `scatter_theme.color` in absence of `ilabels`.
- `node_size=automatic`:
  Defaults to `scatter_theme.markersize` in absence of `ilabels`. Otherwise choses node size based on `ilabels` size.
- `node_marker=automatic`:
  Defaults to `scatter_theme.marker` in absence of `ilabels`.
- `node_strokewidth=automatic`
  Defaults to `scatter_theme.strokewidth` in absence of `ilabels`.
- `node_attr=(;)`: List of kw arguments which gets passed to the `scatter` command
- `edge_color=lineseg_theme.color`: Color for edges.
- `edge_width=lineseg_theme.linewidth`: Pass a vector with 2 width per edge to
  get pointy edges.
- `edge_linestyle=:solid`: Linestyle of edges. Can also be vector or dict for per-edge styling.
  When using different linestyles for different edges, GraphMakie
  creates separate line plots for each edge rather than combining them into one plot, which may reduce
  performance for graphs with many edges. For optimal performance with large graphs, use homogeneous
  linestyles.
- `edge_attr=(;)`: List of kw arguments which gets passed to the underlying `lines` command used for plotting edges.
- `arrow_show=Makie.automatic`: `Bool`, indicate edge directions with arrowheads?
  Defaults to `Graphs.is_directed(graph)`.
- `arrow_marker='➤'`
- `arrow_size=scatter_theme.markersize`: Size of arrowheads.
- `arrow_shift=0.5`: Shift arrow position from source (0) to dest (1) node.
  If `arrow_shift=:end`, the arrowhead will be placed on the surface of the destination node
  (assuming the destination node is circular).
- `arrow_attr=(;)`: List of kw arguments which gets passed to the `scatter` command

### Node labels
The position of each label is determined by the node position plus an offset in
data space.

- `nlabels=nothing`: `Vector{String}` with label for each node
- `nlabels_align=(:left, :bottom)`: Anchor of text field.
- `nlabels_distance=0.0`: Pixel distance from node in direction of align.
- `nlabels_color=labels_theme.color`
- `nlabels_offset=nothing`: `Point` or `Vector{Point}` (in data space)
- `nlabels_fontsize=labels_theme.fontsize`
- `nlabels_attr=(;)`: List of kw arguments which gets passed to the `text` command

### Inner node labels
Put labels inside the marker. If labels are provided, change default attributes to
`node_marker=Circle`, `node_strokewidth=1` and `node_color=:gray80`.
The `node_size` will match size of the `ilabels`.

- `ilabels=nothing`: `Vector` with label for each node
- `ilabels_color=labels_theme.color`
- `ilabels_fontsize=labels_theme.fontsize`
- `ilabels_attr=(;)`: List of kw arguments which gets passed to the `text` command

### Edge labels
The base position of each label is determined by `src + shift*(dst-src)`. The
additional `distance` parameter is given in pixels and shifts the text away from
the edge.

- `elabels=nothing`: `Vector{String}` with label for each edge
- `elabels_align=(:center, :center)`: Anchor of text field.
- `elabels_side = :left`: Side of the edge to put the edge label text
- `elabels_distance=Makie.automatic`: Pixel distance of anchor to edge. The direction is decided based on `elabels_side`
- `elabels_shift=0.5`: Position between src and dst of edge.
- `elabels_rotation=Makie.automatic`: Angle of text per label. If `nothing` this will be
  determined by the edge angle. If `automatic` it will also point upwards making it easy to read.
- `elabels_offset=nothing`: Additional offset in data space
- `elabels_color=labels_theme.color`
- `elabels_fontsize=labels_theme.fontsize`
- `elabels_attr=(;)`: List of kw arguments which gets passed to the `text` command

Self edges / loops:

- `selfedge_size=Makie.automatic()`: Size of selfloop (dict/vector possible).
- `selfedge_direction=Makie.automatic()`: Direction of center of the selfloop as `Point2` (dict/vector possible).
- `selfedge_width=Makie.automatic()`: Opening of selfloop in rad (dict/vector possible).
- Note: If valid waypoints are provided for selfloops, the selfedge attributes above will be ignored.

High level interface for curvy edges:
- `force_straight_edges=false`: If `true`, ignore all curvy edge attributes and draw all edges as straight lines.

- `curve_distance=0.1`:

    Specify a distance of the (now curved) line to the straight line *in data
    space*. Can be single value, array or dict. User provided `tangents` or
    `waypoints` will overrule this property.

- `curve_distance_usage=Makie.automatic()`:

    If `Makie.automatic()`, only plot double edges in a curvy way. Other options
    are `true` and `false`.

Tangents interface for curvy edges:

- `tangents=nothing`:

    Specify a pair of tangent vectors per edge (for src and dst). If `nothing`
    (or edge idx not in dict) draw a straight line.

- `tfactor=0.6`:

    Factor is used to calculate the bezier waypoints from the (normalized) tangents.
    Higher factor means bigger radius. Can be tuple per edge to specify different
    factor for src and dst.

- Note: Tangents are ignored on selfloops if no waypoints are provided.

Waypoints along edges:
- `waypoints=nothing`

    Specify waypoints for edges. This parameter should be given as a vector or
    dict. Waypoints will be crossed using natural cubic splines. The waypoints may
    or may not include the src/dst positions.

- `waypoint_radius=nothing`

    If the attribute `waypoint_radius` is `nothing` or `:spline` the waypoints will
    be crossed using natural cubic spline interpolation. If number (dict/vector
    possible), the waypoints won't be reached, instead they will be connected with
    straight lines which bend in the given radius around the waypoints.
"""
@recipe(GraphPlot, graph) do scene
    # TODO: figure out this whole theme business
    scatter_theme = default_theme(scene, Scatter)
    lineseg_theme = default_theme(scene, LineSegments)
    labels_theme = default_theme(scene, Makie.Text)
    Attributes(
        layout = Spring(),
        # node attributes (Scatter)
        node_color = automatic,
        node_size = automatic,
        node_marker = automatic,
        node_strokewidth = automatic,
        node_attr = (;),
        # edge attributes (LineSegements)
        edge_color = lineseg_theme.color,
        edge_width = lineseg_theme.linewidth,
        edge_linestyle = :solid,
        edge_attr = (;),
        # arrow attributes (Scatter)
        arrow_show = automatic,
        arrow_marker = '➤',
        arrow_size = scatter_theme.markersize,
        arrow_shift = 0.5,
        arrow_attr = (;),
        # node label attributes (Text)
        nlabels = nothing,
        nlabels_align = (:left, :bottom),
        nlabels_distance = 0.0,
        nlabels_color = labels_theme.color,
        nlabels_offset = nothing,
        nlabels_fontsize = labels_theme.fontsize,
        nlabels_attr = (;),
        # inner node labels
        ilabels = nothing,
        ilabels_color = labels_theme.color,
        ilabels_fontsize = labels_theme.fontsize,
        ilabels_attr = (;),
        # edge label attributes (Text)
        elabels = nothing,
        elabels_align = (:center, :center),
        elabels_side = :left,
        elabels_distance = automatic,
        elabels_shift = 0.5,
        elabels_rotation = automatic,
        elabels_offset = nothing,
        elabels_color = labels_theme.color,
        elabels_fontsize = labels_theme.fontsize,
        elabels_attr = (;),
        # self edge attributes
        force_straight_edges = false,
        selfedge_size = automatic,
        selfedge_direction = automatic,
        selfedge_width = automatic,
        curve_distance = 0.1,
        curve_distance_usage = automatic,
        tangents=nothing,
        tfactor=0.6,
        waypoints=nothing,
        waypoint_radius=nothing,
    )
end

function Makie.plot!(gp::GraphPlot)
    dfth = default_theme(gp.parent, GraphPlot)

    # create initial vertex positions, will be updated on changes to graph or layout
    # make node_position-Observable available as named attribute from the outside
    map!(gp.attributes, [:layout, :graph], :node_pos) do layout, graph
        if layout isa AbstractVector
            if length(layout) != nv(graph)
                throw(ArgumentError("The length of the layout vector does not match the number of nodes in the graph!"))
            else
                to_pointf32.(layout)
            end
        else
            [to_pointf32(p) for p in layout(graph)]
        end
    end

    sc = Makie.parent_scene(gp)
    add_input!(gp.attributes, :viewport, sc.viewport)
    add_input!(gp.attributes, :projectionview, sc.camera.projectionview)

    # function which projects the point in px space
    map!(gp.attributes, [:viewport, :projectionview], :to_px) do pxa, pv
        # project should transform to 2d point in px space
        (point) -> project(sc, point)
    end
    # get angle in px space from path p at point t
    map!(gp.attributes, :to_px, :to_angle) do tpx_func
        (path, p0, t) -> begin
            # TODO: maybe shorter tangent? For some perspectives this might give wrong angles in 3d
            p1 = p0 + tangent(path, t)
            any(isnan, p1) && return 0.0  # lines with zero lengths might lead to NaN tangents
            tpx = tpx_func(p1) - tpx_func(p0)
            atan(tpx[2], tpx[1])
        end
    end

    # plot inside labels
    scatter_theme = default_theme(sc, Scatter)

    if gp[:ilabels][] !== nothing
        map!(gp.attributes, :ilabels, :ilabels_text) do ilabels
            string.(ilabels)
        end

        ilabels_plot = text!(gp, gp[:node_pos];
            text=gp[:ilabels_text],
            align=(:center, :center),
            color=gp.ilabels_color,
            fontsize=gp.ilabels_fontsize,
            gp.ilabels_attr[]...)
        add_constant!(gp.attributes, :ilabels_plot, ilabels_plot) #make plotobj accessible

        # only shift very litte to mess less with 3d plots
        translate!(ilabels_plot, 0f32, 0f32, nextfloat(0f32))

        map!(gp.attributes, [:ilabels_plot, :ilabels_text, :ilabels_fontsize, :node_size], :node_size_m) do ilp, txt, ilabels_fontsize, node_size
            bbs = Makie.fast_string_boundingboxes(ilp)
            map(enumerate(bbs)) do (i, bb)
                _ns = getattr(node_size, i)
                if _ns == automatic
                    norm(bb.widths) + 0.1 * ilabels_fontsize
                else
                    _ns
                end
            end
        end
    else
        map!(gp.attributes, :node_size, :node_size_m) do node_size
            node_size === automatic ? scatter_theme.markersize : node_size
        end
    end

    map!(gp.attributes, [:node_color, :ilabels], :node_color_m) do node_color, ilabels
        if node_color === automatic
            ilabels !== nothing ? :gray80 : scatter_theme.color
        else
            node_color
        end
    end

    map!(gp.attributes, [:node_marker, :ilabels], :node_marker_m) do node_marker, ilabels
        if node_marker === automatic
            ilabels !== nothing ? Circle : scatter_theme.marker
        else
            node_marker
        end
    end

    map!(gp.attributes, [:node_strokewidth, :ilabels], :node_strokewidth_m) do node_strokewidth, ilabels
        if node_strokewidth === automatic
            ilabels !== nothing ? 1.0 : scatter_theme.strokewidth
        else
            node_strokewidth
        end
    end

    # compute initial edge paths; will be adjusted later if arrow_shift = :end
    # create array of paths triggered by node_pos changes
    # in case of a graph change the node_position will change anyway
    map!(gp.attributes, [:node_pos, :selfedge_size, :selfedge_direction, :selfedge_width, :curve_distance_usage, :curve_distance, :graph], :init_edge_paths) do pos, s, d, w, cdu, cd, g
        find_edge_paths(g, gp.attributes, pos)
    end

    # plot arrow heads
    map!(gp.attributes, [:init_edge_paths, :to_px, :arrow_shift, :node_marker_m, :node_size_m, :arrow_size, :graph], :arrow_shift_m) do paths, tpx, shift, nmarker, nsize, asize, g
        update_arrow_shift(g, gp, paths, tpx, nmarker, nsize, shift)
    end
    map!(gp.attributes, [:init_edge_paths, :arrow_shift_m, :node_pos], :arrow_pos) do paths, shift_m, np
        if !isempty(paths)
            broadcast(interpolate, paths, shift_m)
        else # if no edges return (empty) vector of points, broadcast yields Vector{Any} which can't be plotted
            Vector{eltype(np)}()
        end
    end
    map!(gp.attributes, [:init_edge_paths, :to_angle, :arrow_pos, :arrow_shift_m], :arrow_rot) do paths, tangle, apos, shift_m
        if !isempty(paths)
            Billboard(broadcast(tangle, paths, apos, shift_m))
        else
            Billboard(Float32[])
        end
    end

    # update edge paths to line up with arrow heads if arrow_shift = :end
    map!(gp.attributes, [:init_edge_paths, :arrow_pos, :arrow_shift], :edge_paths) do paths, apos, shift
        map(paths, apos, eachindex(apos)) do ep, ap, i
            if getattr(shift, i) == :end
                adjust_endpoint(ep, ap)
            else
                ep
            end
        end
    end

    # prepare edge plot attributes (makes them vectors of length ne(g) or single elements)
    map!(gp.attributes, [:edge_color, :graph], :edgeplot_color) do color, graph
        prep_edge_attributes(color, graph, dfth.edge_color[])
    end

    map!(gp.attributes, [:edge_width, :graph], :edgeplot_linewidth) do width, graph
        prep_edge_attributes(width, graph, dfth.edge_width[])
    end

    map!(gp.attributes, [:edge_linestyle, :graph], :edgeplot_linestyle) do style, graph
        prep_edge_attributes(style, graph, dfth.edge_linestyle[])
    end

    # actually plot edges
    edge_plot = edgeplot!(gp, gp[:edge_paths];
        color=gp[:edgeplot_color],
        linewidth=gp[:edgeplot_linewidth],
        linestyle=gp[:edgeplot_linestyle],
        gp.edge_attr[]...)
    add_constant!(gp.attributes, :edge_plot, edge_plot) #make plotobj accessible

    # arrow plots
    map!(gp.attributes, [:arrow_show, :graph], :arrow_show_m) do arrow_show, g
        arrow_show === automatic ? Graphs.is_directed(g) : arrow_show
    end

    # prepare arrow plot attributes
    map!(gp.attributes, [:arrow_marker, :graph], :arrowplot_marker) do marker, graph
        prep_edge_attributes(marker, graph, dfth.arrow_marker[])
    end

    map!(gp.attributes, [:arrow_size, :graph], :arrowplot_markersize) do size, graph
        prep_edge_attributes(size, graph, dfth.arrow_size[])
    end

    map!(gp.attributes, [:edge_color, :graph], :arrowplot_color) do color, graph
        prep_edge_attributes(color, graph, dfth.edge_color[])
    end

    arrow_plot = scatter!(gp,
        gp[:arrow_pos];
        marker = gp[:arrowplot_marker],
        markersize = gp[:arrowplot_markersize],
        color = gp[:arrowplot_color],
        rotation = gp[:arrow_rot],
        strokewidth = 0.0,
        markerspace = :pixel,
        visible = gp[:arrow_show_m],
        gp.arrow_attr[]...)
    add_constant!(gp.attributes, :arrow_plot, arrow_plot) #make plotobj accessible


    # prepare node plot attributes
    map!(gp.attributes, [:node_color_m, :graph], :nodeplot_color) do color, graph
        prep_vertex_attributes(color, graph, scatter_theme.color)
    end

    map!(gp.attributes, [:node_marker_m, :graph], :nodeplot_marker) do marker, graph
        prep_vertex_attributes(marker, graph, scatter_theme.marker)
    end

    map!(gp.attributes, [:node_size_m, :graph], :nodeplot_markersize) do size, graph
        prep_vertex_attributes(size, graph, scatter_theme.markersize)
    end

    map!(gp.attributes, [:node_strokewidth_m, :graph], :nodeplot_strokewidth) do width, graph
        prep_vertex_attributes(width, graph, scatter_theme.strokewidth)
    end

    vertex_plot = scatter!(gp, gp[:node_pos];
        color=gp[:nodeplot_color],
        marker=gp[:nodeplot_marker],
        markersize=gp[:nodeplot_markersize],
        strokewidth=gp[:nodeplot_strokewidth],
        gp[:node_attr][]...)
    add_constant!(gp.attributes, :node_plot, vertex_plot) #make plotobj accessible

    # plot node labels
    if gp.nlabels[] !== nothing
        map!(gp.attributes, [:node_pos, :nlabels_offset], :nlabels_positions) do np, offset
            if offset != nothing
                np .+ offset
            else
                copy(np)
            end
        end

        map!(gp.attributes, [:nlabels_align, :nlabels_distance], :nlabels_offset_processed) do align, distance
            if align isa Vector
                distance .* align_to_dir.(align)
            else
                distance .* align_to_dir(align)
            end
        end

        # prepare node labels attributes
        map!(gp.attributes, [:nlabels, :graph], :nlabels_text_processed) do labels, graph
            prep_vertex_attributes(labels, graph, "")
        end

        map!(gp.attributes, [:nlabels_align, :graph], :nlabels_align_processed) do align, graph
            prep_vertex_attributes(align, graph, dfth.nlabels_align[])
        end

        map!(gp.attributes, [:nlabels_color, :graph], :nlabels_color_processed) do color, graph
            prep_vertex_attributes(color, graph, dfth.nlabels_color[])
        end

        map!(gp.attributes, [:nlabels_fontsize, :graph], :nlabels_fontsize_processed) do fontsize, graph
            prep_vertex_attributes(fontsize, graph, dfth.nlabels_fontsize[])
        end

        nlabels_plot = text!(gp, gp[:nlabels_positions];
            text=gp[:nlabels_text_processed],
            align=gp[:nlabels_align_processed],
            color=gp[:nlabels_color_processed],
            offset=gp[:nlabels_offset_processed],
            fontsize=gp[:nlabels_fontsize_processed],
            gp.nlabels_attr[]...)
        add_constant!(gp.attributes, :nlabels_plot, nlabels_plot) #make plotobj accessible
    end

    # plot edge labels
    if gp.elabels[] !== nothing
        # positions: center point between nodes + offset + distance*normal + shift*edge direction
        map!(gp.attributes, [:edge_paths, :elabels_shift, :elabels_offset], :elabels_positions) do paths, shift, eloffset
            pos = broadcast(interpolate, paths, shift)

            if eloffset !== nothing
                pos .= pos .+ eloffset
            end
            pos
        end

        # rotations based on the edge_vec_px and opposite argument
        map!(gp.attributes, [:elabels_rotation, :to_angle, :elabels_positions, :edge_paths, :elabels_shift, :graph], :elabels_rotation_computed) do elabrots, tangle, pos, paths, shift, g
            rot = broadcast(tangle, paths, pos, shift)
            for i in 1:ne(g)
                valrot = getattr(elabrots, i, nothing)
                if valrot isa Real
                    # fix rotation to a single angle
                    rot[i] = valrot
                elseif valrot == automatic
                    # point the labels up
                    if (rot[i] > π/2 || rot[i] < - π/2)
                        rot[i] += π
                    end
                end
            end
            return rot
        end

        # calculate the offset in pixels in normal direction to the edge
        map!(gp.attributes, [:elabels_positions, :to_px, :elabels_distance, :elabels_side, :edge_paths, :elabels_shift, :graph, :elabels_fontsize, :edge_width], :elabels_offsets) do pos, tpx, dist, side, paths, shift, g, fontsize, edge_width
            tangent_px = broadcast(paths, pos, shift) do path, p0, t
                p1 = p0 + tangent(path, t)
                tpx(p1) - tpx(p0)
            end

            offsets = map(p -> Point(-p.data[2], p.data[1])/norm(p), tangent_px)
            offsets .= elabels_distance_offset(g, gp.attributes) .* offsets
        end

        # prepare edge labels attributes
        map!(gp.attributes, [:elabels, :graph], :elabels_text_processed) do labels, graph
            prep_edge_attributes(labels, graph, "")
        end

        map!(gp.attributes, [:elabels_align, :graph], :elabels_align_processed) do align, graph
            prep_edge_attributes(align, graph, dfth.elabels_align[])
        end

        map!(gp.attributes, [:elabels_color, :graph], :elabels_color_processed) do color, graph
            prep_edge_attributes(color, graph, dfth.elabels_color[])
        end

        map!(gp.attributes, [:elabels_fontsize, :graph], :elabels_fontsize_processed) do fontsize, graph
            prep_edge_attributes(fontsize, graph, dfth.elabels_fontsize[])
        end

        elabels_plot = text!(gp, gp[:elabels_positions];
            text=gp[:elabels_text_processed],
            rotation=gp[:elabels_rotation_computed],
            offset=gp[:elabels_offsets],
            align=gp[:elabels_align_processed],
            color=gp[:elabels_color_processed],
            fontsize=gp[:elabels_fontsize_processed],
            gp.elabels_attr[]...)
        add_constant!(gp.attributes, :elabels_plot, elabels_plot) #make plotobj accessible
    end

    return gp
end

"""
    elabels_distance_offset(g, attrs)

Returns the elabels_distance taking into consideration elabels_side
"""
function elabels_distance_offset(g, attrs)
    offs = zeros(ne(g))
    for i in 1:ne(g)
        attrval = getattr(attrs.elabels_distance, i, automatic)
        attrvalside = getattr(attrs.elabels_side, i, :left)
        if attrval isa Real
            if attrvalside == :left
                offs[i] = attrval
            elseif attrvalside == :right
                offs[i] = -attrval
            elseif attrvalside == :center
                offs[i] = zero(attrval)
            end
        elseif attrval == automatic
            offval = (getattr(attrs.elabels_fontsize, i) + getattr(attrs.edge_width, i))/2
            if attrvalside == :left
                offs[i] = offval
            elseif attrvalside == :right
                offs[i] = -offval
            elseif attrvalside == :center
                offs[i] = zero(offval)
            end
        end
    end
    return offs
end

"""
    find_edge_paths(g, attr, pos::AbstractVector{PT}) where {PT}

Returns an `AbstractPath` for each edge in the graph. Returns a vector of
paths. If `attr.force_straight_edges` is `true`, the paths will be just plain lines
"""
function find_edge_paths(g, attr, pos::AbstractVector{PT}) where {PT}
    # for straight_lines: return vector of Line rather than vector of AbstractPath
    if attr.force_straight_edges[]
        return map(edges(g)) do e
            p1, p2 = pos[src(e)], pos[dst(e)]
            Path(p1, p2)
        end
    end

    paths = Vector{AbstractPath{PT}}(undef, ne(g))
    for (i, e) in enumerate(edges(g))
        p1, p2 = pos[src(e)], pos[dst(e)]

        tangents = getattr(attr.tangents, i)
        tfactor = getattr(attr.tfactor, i)
        waypoints::Vector{PT} = getattr(attr.waypoints, i, PT[])
        if !isnothing(waypoints) && !isempty(waypoints) #remove p1 and p2 from waypoints if these are given
            waypoints[begin] == p1 && popfirst!(waypoints)
            waypoints[end] == p2 && pop!(waypoints)
        end

        cdu = getattr(attr.curve_distance_usage, i)
        if cdu === true
            curve_distance = getattr(attr.curve_distance, i, 0.0)
        elseif cdu === false
            curve_distance = 0.0
        elseif cdu === automatic
            if is_directed(g) && has_edge(g, dst(e), src(e))
                curve_distance = getattr(attr.curve_distance, i, 0.0)
            else
                curve_distance = 0.0
            end
        end

        if !isnothing(waypoints) && !isempty(waypoints) #there are waypoints
            radius = getattr(attr.waypoint_radius, i, nothing)
            if radius === nothing || radius === :spline
                paths[i] = Path(p1, waypoints..., p2; tangents, tfactor)
            elseif radius isa Real
                paths[i] = Path(radius, p1, waypoints..., p2)
            else
                throw(ArgumentError("Invalid radius $radius for edge $i!"))
            end
        elseif src(e) == dst(e) # selfedge
            size = getattr(attr.selfedge_size, i)
            direction = getattr(attr.selfedge_direction, i)
            width = getattr(attr.selfedge_width, i)
            paths[i] = selfedge_path(g, pos, src(e), size, direction, width)
        elseif !isnothing(tangents)
            paths[i] = Path(p1, p2; tangents, tfactor)
        elseif PT<:Point2 && !iszero(curve_distance)
            paths[i] = curved_path(p1, p2, curve_distance)
        else # straight line
            paths[i] = Path(p1, p2)
        end
    end

    return paths
end

"""
    selfedge_path(g, pos, v, size, direction, width)

Return a BezierPath for a selfedge.
"""
function selfedge_path(g, pos::AbstractVector{<:Point2}, v, size, direction, width)
    vp = pos[v]
    # get the vectors to all the neighbors
    ndirs = [pos[n] - vp for n in all_neighbors(g, v) if n != v]

    # angle and maximum width of loop
    γ, Δ = 0.0, 0.0

    if direction === automatic && !isempty(ndirs)
        angles = SVector{length(ndirs)}(atan(p[2], p[1]) for p in ndirs)
        angles = sort(angles)

        for i in 1:length(angles)
            α = angles[i]
            β = get(angles, i+1, 2π + angles[1])
            if β-α > Δ
                Δ = β-α
                γ = (β+α) / 2
            end
        end

        # set width of selfloop
        Δ = min(.7*Δ, π/2)
    elseif direction === automatic && isempty(ndirs)
        γ = π/2
        Δ = π/2
    else
        @assert direction isa Point2 "Direction of selfedge should be 2 dim vector ($direction)"
        γ = atan(direction[2], direction[1])
        Δ = π/2
    end

    if width !== automatic
        Δ = width
    end

    # the size (max distance to v) of loop
    # if there are no neighbors set to 0.5. Else half dist to nearest neighbor
    if size === automatic
        size = isempty(ndirs) ? 0.5 : minimum(norm.(ndirs)) * 0.5
    end

    # the actual length of the tagent vectors, magic number from `CurveTo`
    l = Float32( size/(cos(Δ/2) * 2*0.375) )
    t1 = vp + l * Point2f(cos(γ-Δ/2), sin(γ-Δ/2))
    t2 = vp + l * Point2f(cos(γ+Δ/2), sin(γ+Δ/2))

    return BezierPath([MoveTo(vp),
                       CurveTo(t1, t2, vp)])
end

function selfedge_path(g, pos::AbstractVector{<:Point3}, v, size, direction, width)
    error("Self edges in 3D not yet supported")
end

"""
    curved_path(p1, p2, curve_distance)

Return a BezierPath for a curved edge (not selfedge).
"""
function curved_path(p1::PT, p2::PT, curve_distance) where {PT}
    d = curve_distance
    s = norm(p2 - p1)
    γ = 2*atan(2 * d/s)
    a = (p2 - p1)/s * (4*d^2 + s^2)/(3s)

    m = @SMatrix[cos(γ) -sin(γ); sin(γ) cos(γ)]
    c1 = PT(p1 + m*a)
    c2 = PT(p2 - transpose(m)*a)

    return BezierPath([MoveTo(p1), CurveTo(c1, c2, p2)])
end

@recipe EdgePlot (paths,) begin
    Makie.documented_attributes(Lines)...
end

function Makie.plot!(p::EdgePlot)
    alllines = eltype(p[:paths][]) <: Line

    map!(p.attributes, :paths, [:points, :ranges]) do paths
        PT = ptype(eltype(paths))
        points = PT[]
        ranges = UnitRange{Int}[]
        for path in paths
            disc = discretize(path)
            pstart = length(points) + 1
            append!(points, disc)
            push!(points, PT(NaN)) # add NaN to separate segments
            pstop = pstart+length(disc)
            push!(ranges, pstart:pstop)
        end
        (points, ranges)
    end

    # if the user specified different linestyles, we need to fall back to plotting n `lines` rather than plotting
    split_edgeplots = !(p[:linestyle][] isa Union{Nothing,Symbol,Linestyle})
    add_constant!(p.attributes, :split_edgeplots, split_edgeplots)

    if !split_edgeplots
        # expand color and linewidth attributes (for curved edges with multiple segments)
        map!(_expand_args, p.attributes, [:color, :ranges], :color_expanded)
        map!(_expand_args, p.attributes, [:linewidth, :ranges], :linewidth_expanded)

        lines!(p, p.attributes, p[:points];
            color=p[:color_expanded],
            linewidth=p[:linewidth_expanded],
        )
    else
        # manually find colorrange
        if p[:colorrange][] === automatic && p[:color][] isa Union{Number, AbstractVector}
            minc = Inf
            maxc = -Inf
            for c in p[:color][]
                if c isa Number
                    minc = min(minc, c)
                    maxc = max(maxc, c)
                elseif c isa AbstractVector{<:Number} || c isa Tuple{<:Number}
                    minc = min(minc, minimum(c))
                    maxc = max(maxc, maximum(c))
                end
            end
            if minc != Inf || maxc != -Inf
                p[:colorrange][] = Float32.((minc, maxc))
            end
        end

        for i in 1:length(p[:paths][])
            thispoints = Symbol(:points,i)
            map!(p.attributes, [:points, :ranges], thispoints) do points, ranges
                view(points, ranges[i][1:end-1])
            end

            lines!(p, p.attributes, p[thispoints];
                color=_split_arg!(p.attributes, :color, i),
                linewidth=_split_arg!(p.attributes, :linewidth, i),
                linestyle=_split_arg!(p.attributes, :linestyle, i),
            )
        end
    end

    return p
end

function _expand_args(args::Union{AbstractVector, AbstractDict}, ranges)
    N_paths = length(ranges)
    N_points = N_paths > 0 ? ranges[end][end] : 0
    allstraight = N_paths*3 == N_points
    if args isa AbstractVector && length(args) != N_paths
        throw(ArgumentError("The length of the args vector $args does not match the number of edges!"))
    end

    elT = eltype(args) <:Tuple ? eltype(eltype(args)) : eltype(args)
    elT = elT <: Integer ? Float32 : elT # convert integers to floats for interpolation
    expanded = Vector{elT}(undef, N_points)
    for i in 1:N_paths
        attr = getattr(args, i)
        if attr isa Union{Tuple, AbstractVector}
            if length(attr) == length(ranges[i]) - 1
                expanded[ranges[i][1:end-1]] .= attr
                expanded[ranges[i][end]] = attr[end] # last point is always the same
            elseif eltype(attr) <: Number && length(attr) == 2 # interpolate between numbers
                expanded[ranges[i][1:end-1]].= range(attr[1], attr[2], length=length(ranges[i])-1)
                expanded[ranges[i][end]] = attr[end] # last point is always the same
            else
                throw(ArgumentError("Don't know how to map $(attr) to the $(length(ranges[i])) points in this edge."))
            end
        else
            expanded[ranges[i]] .= args[i]
        end
    end
    expanded
end
_expand_args(arg, ranges) = arg
function _split_arg!(cg::Makie.ComputeGraph, name, i)
    splitname = Symbol(name, i)
    map!(cg, [name, Symbol(:points, i)], splitname) do prop, pointsi
        attr = getattr(prop, i)
        # interpolate numeric values for intermediate points
        if attr isa Union{Tuple,AbstractVector} && eltype(attr) <: Number && length(attr) == 2
            attr = range(attr[1], attr[2], length=length(pointsi))
        end
        attr
    end
    cg[splitname]
end

"""
    update_arrow_shift(g, gp, edge_paths::Vector{<:AbstractPath{PT}}, to_px) where {PT}

Checks `arrow_shift` attr so that `arrow_shift = :end` gets transformed so that the arrowhead for that edge
lands on the surface of the destination node.
"""
function update_arrow_shift(g, gp, edge_paths::Vector{<:AbstractPath{PT}}, to_px, node_markers, node_sizes, shift) where {PT}
    arrow_shift = Vector{Float32}(undef, ne(g))

    for (i,e) in enumerate(edges(g))
        t = getattr(shift, i, 0.5)
        if t === :end
            j = dst(e)
            p0 = getattr(gp.node_pos, j)
            node_marker = getattr(node_markers, j)
            node_size = getattr(node_sizes, j)
            arrow_marker = getattr(gp.arrow_marker, i)
            arrow_size = getattr(gp.arrow_size, i)
            d = distance_between_markers(node_marker, node_size, arrow_marker, arrow_size)
            p1 = point_near_dst(edge_paths[i], p0, d, to_px)
            t = inverse_interpolate(edge_paths[i], p1)
            if isnan(t)
                @warn """
                    Shifting arrowheads to destination nodes failed.
                    This can happen when the markers are inadequately scaled (e.g., when zooming out too far).
                    Arrow shift has been reset to 0.5.
                """
                t = 0.5
            end
        end
        arrow_shift[i] = t
    end

    return arrow_shift
end

function update_arrow_shift(g, gp, edge_paths::Vector{<:AbstractPath{<:Point3}}, to_px, shift)
    arrow_shift = Vector{Float32}(undef, ne(g))

    for (i,e) in enumerate(edges(g))
        t = getattr(shift, i, 0.5)
        if t === :end #not supported because to_px does not give pixels in 3D space (would need to map 3D coordinates to pixels...?)
            error("`arrow_shift = :end` not supported for 3D plots.")
        end
        arrow_shift[i] = t
    end

    return arrow_shift
end

function Makie.preferred_axis_type(plot::Plot{GraphMakie.graphplot})
    if haskey(plot.kw, :layout)
        layout = plot.kw[:layout]
        dim = _dimensionality(layout, plot[1][])
        dim == 3 && return LScene
        dim == 2 && return Axis
    end
    Axis
end

_dimensionality(obs::Observable, g) = _dimensionality(obs[], g)
_dimensionality(layout::AbstractLayout, _) = dim(layout)
_dimensionality(layout::AbstractArray, _) = length(first(layout))
_dimensionality(layout, g) = length(first(layout(g)))
