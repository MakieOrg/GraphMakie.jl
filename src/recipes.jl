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
- `layout=Spring()`: function `AbstractGraph->Vector{Point}` or `Vector{Point}` determines the base layout
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
- `edge_attr=(;)`: List of kw arguments which gets passed to the `linesegments` command
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

### Curvy edges & self edges/loops
- `edge_plottype=Makie.automatic()`: Either `automatic`, `:linesegments` or
  `:beziersegments`. `:beziersegments` are much slower for big graphs!

Self edges / loops: 

- `selfedge_size=Makie.automatic()`: Size of selfloop (dict/vector possible).
- `selfedge_direction=Makie.automatic()`: Direction of center of the selfloop as `Point2` (dict/vector possible).
- `selfedge_width=Makie.automatic()`: Opening of selfloop in rad (dict/vector possible).
- Note: If valid waypoints are provided for selfloops, the selfedge attributes above will be ignored.

High level interface for curvy edges:

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
        edge_plottype = automatic,
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
    graph = gp[:graph]
    
    dfth = default_theme(gp.parent, GraphPlot)

    # create initial vertex positions, will be updated on changes to graph or layout
    # make node_position-Observable available as named attribute from the outside
    gp[:node_pos] = @lift if $(gp.layout) isa AbstractVector
        if length($(gp.layout)) != nv($graph)
            throw(ArgumentError("The length of the layout vector does not match the number of nodes in the graph!"))
        else
            Pointf.($(gp.layout))
        end
    else
        [Pointf(p) for p in ($(gp.layout))($graph)]
    end

    sc = Makie.parent_scene(gp)

    # function which projects the point in px space
    to_px = lift(gp, sc.viewport, sc.camera.projectionview) do pxa, pv
        # project should transform to 2d point in px space
        (point) -> project(sc, point)
    end
    # get angle in px space from path p at point t
    to_angle = @lift (path, p0, t) -> begin
        # TODO: maybe shorter tangent? For some perspectives this might give wrong angles in 3d
        p1 = p0 + tangent(path, t)
        tpx = $to_px(p1) - $to_px(p0)
        atan(tpx[2], tpx[1])
    end

    node_pos = gp[:node_pos]

    # plot inside labels
    scatter_theme = default_theme(sc, Scatter)

    if gp[:ilabels][] !== nothing
        positions = node_pos
        
        ilabels_plot = text!(gp, positions;
            text=@lift(string.($(gp.ilabels))),
            align=(:center, :center),
            color=gp.ilabels_color,
            fontsize=gp.ilabels_fontsize,
            gp.ilabels_attr...)

        translate!(ilabels_plot, 0, 0, 1)

        node_size_m = lift(ilabels_plot.plots[1][1], gp.ilabels_fontsize, gp.node_size) do glyphcollections, ilabels_fontsize, node_size
            map(enumerate(glyphcollections)) do (i, gc)
                _ns = getattr(node_size, i)
                if _ns == automatic
                    rect = Rect2f(boundingbox(gc, Quaternion((1,0,0,0))))
                    norm(rect.widths) + 0.1 * ilabels_fontsize
                else
                    _ns
                end
            end
        end
    else
        node_size_m = @lift $(gp.node_size) === automatic ? scatter_theme.markersize[] : $(gp.node_size)
    end

    node_color_m = @lift if $(gp.node_color) === automatic
        gp.ilabels[] !== nothing ? :gray80 : scatter_theme.color[]
    else
        $(gp.node_color)
    end
    
    node_marker_m = @lift if $(gp.node_marker) === automatic
        gp.ilabels[] !== nothing ? Circle : scatter_theme.marker[]
    else
        $(gp.node_marker)
    end

    node_strokewidth_m = @lift if $(gp.node_strokewidth) === automatic
        gp.ilabels[] !== nothing ? 1.0 : scatter_theme.strokewidth[]
    else
        $(gp.node_strokewidth)
    end

    # create array of pathes triggered by node_pos changes
    # in case of a graph change the node_position will change anyway
    gp[:edge_paths] = lift(node_pos, gp.selfedge_size,
                      gp.selfedge_direction, gp.selfedge_width, gp.curve_distance_usage, gp.curve_distance) do pos, s, d, w, cdu, cd
        find_edge_paths(graph[], gp.attributes, pos)
    end
    edge_paths = gp[:edge_paths]

    # plot edges
    edge_plot = edgeplot!(gp, edge_paths;
        color=prep_edge_attributes(gp.edge_color, graph, dfth.edge_color),
        linewidth=prep_edge_attributes(gp.edge_width, graph, dfth.edge_width),
        gp.edge_attr...)

    # plot arrow heads
    arrow_shift_m = lift(edge_paths, to_px, gp.arrow_shift, node_marker_m, node_size_m, gp.arrow_size) do paths, tpx, shift, nmarker, nsize, asize
        update_arrow_shift(graph[], gp, paths, tpx, node_marker_m, node_size_m, shift)
    end
    arrow_pos = @lift if !isempty(edge_paths[])
        broadcast(interpolate, edge_paths[], $arrow_shift_m)
    else # if no edges return (empty) vector of points, broadcast yields Vector{Any} which can't be plotted
        Vector{eltype(node_pos[])}()
    end
    arrow_rot = @lift if !isempty(edge_paths[])
        Billboard(broadcast($to_angle, edge_paths[], $arrow_pos, $(arrow_shift_m)))
    else
        Billboard(Float32[])
    end
    arrow_show_m = @lift $(gp.arrow_show) === automatic ? Graphs.is_directed($graph) : $(gp.arrow_show)
    arrow_heads = scatter!(gp,
        arrow_pos;
        marker = prep_edge_attributes(gp.arrow_marker, graph, dfth.arrow_marker),
        markersize = prep_edge_attributes(gp.arrow_size, graph, dfth.arrow_size),
        color=prep_edge_attributes(gp.edge_color, graph, dfth.edge_color),
        rotations = arrow_rot,
        strokewidth = 0.0,
        markerspace = :pixel,
        visible = arrow_show_m,
        gp.arrow_attr...)

    # plot vertices
    vertex_plot = scatter!(gp, node_pos;
        color=prep_vertex_attributes(node_color_m, graph, scatter_theme.color),
        marker=prep_vertex_attributes(node_marker_m, graph, scatter_theme.marker),
        markersize=prep_vertex_attributes(node_size_m, graph, scatter_theme.markersize),
        strokewidth=prep_vertex_attributes(node_strokewidth_m, graph, scatter_theme.strokewidth),
        gp.node_attr...)

    # plot node labels
    if gp.nlabels[] !== nothing
        positions = @lift begin
            if $(gp.nlabels_offset) != nothing
                $node_pos .+ $(gp.nlabels_offset)
            else
                copy($node_pos)
            end
        end

        offset = @lift if $(gp.nlabels_align) isa Vector
            $(gp.nlabels_distance) .* align_to_dir.($(gp.nlabels_align))
        else
            $(gp.nlabels_distance) .* align_to_dir($(gp.nlabels_align))
        end

        nlabels_plot = text!(gp, positions;
            text=prep_vertex_attributes(gp.nlabels, graph, Observable("")),
            align=prep_vertex_attributes(gp.nlabels_align, graph, dfth.nlabels_align),
            color=prep_vertex_attributes(gp.nlabels_color, graph, dfth.nlabels_color),
            offset=offset,
            fontsize=prep_vertex_attributes(gp.nlabels_fontsize, graph, dfth.nlabels_fontsize),
            gp.nlabels_attr...)
    end

    # plot edge labels
    if gp.elabels[] !== nothing
        # positions: center point between nodes + offset + distance*normal + shift*edge direction
        positions = @lift begin
            pos = broadcast(interpolate, $edge_paths, $(gp.elabels_shift))

            if $(gp.elabels_offset) !== nothing
                pos .= pos .+ $(gp.elabels_offset)
            end
            pos
        end

        # rotations based on the edge_vec_px and opposite argument
        rotation = lift(gp.elabels_rotation, to_angle, positions) do elabrots, tangle, pos
            rot = broadcast(tangle, edge_paths[], pos, gp.elabels_shift[])
            for i in 1:ne(graph[])
                valrot = getattr(gp.elabels_rotation, i, nothing)
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
        offsets = lift(positions, to_px, gp.elabels_distance, gp.elabels_side) do pos, tpx, dist, side
            tangent_px = broadcast(edge_paths[], pos, gp.elabels_shift[]) do path, p0, t
                p1 = p0 + tangent(path, t)
                tpx(p1) - tpx(p0)
            end

            offsets = map(p -> Point(-p.data[2], p.data[1])/norm(p), tangent_px)
            offsets .= elabels_distance_offset(graph[], gp.attributes) .* offsets
        end
        elabels_plot = text!(gp, positions;
            text=prep_edge_attributes(gp.elabels, graph, Observable("")),
            rotation=rotation,
            offset=offsets,
            align=prep_edge_attributes(gp.elabels_align, graph, dfth.elabels_align),
            color=prep_edge_attributes(gp.elabels_color, graph, dfth.elabels_color),
            fontsize=prep_edge_attributes(gp.elabels_fontsize, graph, dfth.elabels_fontsize),
            gp.elabels_attr...)
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

Returns an `AbstractPath` for each edge in the graph. Based on the `edge_plotype` attribute
this returns either arbitrary bezier curves or just lines.
"""
function find_edge_paths(g, attr, pos::AbstractVector{PT}) where {PT}
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

    plottype = attr[:edge_plottype][]

    if plottype === :beziersegments
        # user explicitly specified beziersegments
        return paths
    elseif plottype === :linesegments
        # user specified linesegments but there are non-lines
        if !isempty(paths)
            return straighten.(paths)
        else
            return Line{PT}[]
        end
    elseif plottype === automatic
        ls = try
            first(attr.edge_attr.linestyle[])
        catch
            nothing
        end
        if !isnothing(ls) && !isa(ls, Number)
            throw(ArgumentError("Linestyles per edge are only supported when passing `edge_plottype=:beziersegments` to `graphplot` explicitly."))
        end

        # try to narrow down the typ, i.e. just `Line`s
        # update :plottype in attr so this is persistent even if the graph changes
        T = isempty(paths) ? Line{PT} : mapreduce(typeof, promote_type, paths)

        if T <: Line
            attr[:edge_plottype][] = :linesegments
            return convert(Vector{T}, paths)
        elseif ne(g) > 500
            attr[:edge_plottype][] = :linesegments
            @warn "Since there are a lot of edges ($(ne(g)) > 500), they will be drawn as straight lines "*
                "even though they contain curvy edges. If you really want to plot them as "*
                "bezier curves pass `edge_plottype=:beziersegments` explicitly. This will have "*
                "much worse performance!"
            return straighten.(paths)
        else
            attr[:edge_plottype][] = :beziersegments
            return paths
        end
    else
        throw(ArgumentError("Invalid argument for edge_plottype: $plottype"))
    end
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

"""
    edgeplot(paths::Vector{AbstractPath})
    edgeplot!(sc, paths::Vector{AbstractPath})

Recipe to draw the edges. Attribute `plottype` can be either.
If `eltype(pathes) <: Line` draw using `linesegments`. If `<:AbstractPath`
draw using bezier segments.
All attributes are passed down to the actual recipe.
"""
@recipe(EdgePlot, paths) do scene
    Attributes()
end

function Makie.plot!(p::EdgePlot)
    N = length(p[:paths][])

    if eltype(p[:paths][]) <: Line
        PT = ptype(eltype(p[:paths][]))
        segs = Observable(Vector{PT}(undef, 2*N))

        update_segments!(segs, p[:paths][]) # first call to initialize
        on(p[:paths]) do paths # update if pathes change
            update_segments!(segs, paths)
        end

        linesegments!(p, segs; p.attributes...)
    else
        beziersegments!(p, p[:paths]; p.attributes...)
    end

    return p
end

function update_segments!(segs, paths)
    N = length(paths)
    if length(segs[]) != 2*N
        resize!(segs[], 2*N)
    end
    for (i, p) in enumerate(paths)
        segs[][2*i - 1] = interpolate(p, 0.0)
        segs[][2*i]     = interpolate(p, 1.0)
    end
    notify(segs)
end


"""
    beziersegments(paths::Vector{AbstractPath})
    beziersegments!(sc, paths::Vector{AbstractPath})

Recipe to draw bezier pathes. Each path will be descritized and ploted with a
separate `lines` plot. Scalar attributes will be used for all subplots. If you
provide vector attributs of same length als the pathes the i-th `lines` subplot
will see the i-th attribute.

TODO: Beziersegments plot won't work if the number of pathes changes.
"""
@recipe(BezierSegments, paths) do scene
    Attributes(default_theme(scene, LineSegments)...)
end

function Makie.plot!(p::BezierSegments)
    N = length(p[:paths][])
    PT = ptype(eltype(p[:paths][]))
    attr = p.attributes

    if N == 0 # no edges, still need to plot some atomic element otherwise recipe does not work
        lines!(p, PT[])
        return p
    end

    # set colorange automaticially if needed
    if attr[:color][] isa AbstractArray{<:Number}
        numbers = attr[:color][]
        colorrange = get(attr, :colorrange, nothing) |> to_value

        if colorrange === Makie.automatic
            attr[:colorrange] = extrema(numbers)
        end
    end

    # array which holds observables for the discretized values
    # this is needed so the `lines!` subplots will update
    disc = [Observable{Vector{PT}}() for i in 1:N]

    update_discretized!(disc, p[:paths][]) # first call to initialize
    on(p[:paths]) do paths # update if pathes change
        update_discretized!(disc, paths)
    end

    # plot all the lines
    for i in 1:N
        # each subplot will pick its specic arguments if there are vectors
        specific_attr = Attributes()

        for (k, v) in attr
            if k === :color && v[] isa AbstractVector{<:Number} && length(v[]) == N
                colormap = attr[:colormap][] |> to_colormap
                colorrange = attr[:colorrange][]
                specific_attr[k] = @lift Makie.interpolated_getindex(colormap,
                                                                     Float64($v[i]),
                                                                     colorrange)
            elseif k === :linestyle && v[] isa AbstractVector && first(v[]) isa Number
                # linestyle may be vector of numbers to specify pattern, in that case don't split
                specific_attr[k] = v
            elseif v[] isa AbstractVector && length(v[]) == N
                specific_attr[k] = @lift $v[i]
            else
                specific_attr[k] = v
            end
        end
        lines!(p, disc[i]; specific_attr...)
    end

    return p
end

function update_discretized!(disc, pathes)
    for (i, p) in enumerate(pathes)
        disc[i][] = discretize(p)
    end
end

"""
    update_arrow_shift(g, gp, edge_paths::Vector{<:AbstractPath{PT}}, to_px) where {PT}

Checks `arrow_shift` attr so that `arrow_shift = :end` gets transformed so that the arrowhead for that edge
lands on the surface of the destination node.
"""
function update_arrow_shift(g, gp, edge_paths::Vector{<:AbstractPath{PT}}, to_px, node_markers, node_sizes, shift) where {PT}
    arrow_shift = Vector{Float32}(undef, ne(g))

    warn_nan = false
    warn_3d = false
    for (i,e) in enumerate(edges(g))
        t = getattr(shift, i, 0.5)
        if t === :end
            if PT <: Point2
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
                    warn_nan = true
                    t = 0.5
                end
            else
                warn_3d = true
                t = 0.5
            end
        end
        arrow_shift[i] = t
    end
    if warn_nan
        @warn """
        Shifting arrowheads to destination nodes failed.
        This can happen when the markers are inadequately scaled (e.g., when zooming out too far).
        Arrow shift has been reset to 0.5.
        """
    end
    if warn_3d
        @warn "`arrow_shift=:end` not implemented for 3d plot!"
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

_dimensionality(layout::AbstractLayout, _) = dim(layout)
_dimensionality(layout::AbstractArray, _) = length(first(layout))
_dimensionality(layout, g) = length(first(layout(g)))
