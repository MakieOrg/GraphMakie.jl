using LinearAlgebra: normalize, ⋅, norm
export GraphPlot, graphplot, graphplot!

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
- `layout=Spring()`: function `AbstractGraph->Vector{Point}` determines the base layout
- `node_color=scatter_theme.color`
- `node_size=scatter_theme.markersize`
- `node_marker=scatter_theme.marker`
- `node_attr=(;)`: List of kw arguments which gets passed to the `scatter` command
- `edge_color=lineseg_theme.color`: Color for edges.
- `edge_width=lineseg_theme.linewidth`: Pass a vector with 2 width per edge to
  get pointy edges.
- `edge_attr=(;)`: List of kw arguments which gets passed to the `linesegments` command
- `arrow_show=Makie.automatic`: `Bool`, indicate edge directions with arrowheads?
  Defaults to `Graphs.is_directed(graph)`.
- `arrow_size=scatter_theme.markersize`: Size of arrowheads.
- `arrow_shift=0.5`: Shift arrow position from source (0) to dest (1) node.
- `arrow_attr=(;)`: List of kw arguments which gets passed to the `scatter` command

### Node labels
The position of each label is determined by the node position plus an offset in
data space.

- `nlabels=nothing`: `Vector{String}` with label for each node
- `nlabels_align=(:left, :bottom)`: Anchor of text field.
- `nlabels_distance=0.0`: Pixel distance from node in direction of align.
- `nlabels_color=labels_theme.color`
- `nlabels_offset=nothing`: `Point` or `Vector{Point}` (in data space)
- `nlabels_textsize=labels_theme.textsize`
- `nlabels_attr=(;)`: List of kw arguments which gets passed to the `text` command

### Edge labels
The base position of each label is determinded by `src + shift*(dst-src)`. The
additional `distance` parameter is given in pixels and shifts the text away from
the edge.

- `elabels=nothing`: `Vector{String}` with label for each edge
- `elabels_align=(:center, :bottom)`: Anchor of text field.
- `elabels_distance=0.0`: Pixel distance of anchor to edge.
- `elabels_shift=0.5`: Position between src and dst of edge.
- `elabels_opposite=Int[]`: List of edge indices, for which the label should be
  displayed on the opposite side
- `elabels_rotation=nothing`: Angle of text per label. If `nothing` this will be
  determined by the edge angle!
- `elabels_offset=nothing`: Additional offset in data space
- `elabels_color=labels_theme.color`
- `elabels_textsize=labels_theme.textsize`
- `elabels_attr=(;)`: List of kw arguments which gets passed to the `text` command

### Curvy edges & self edges/loops
- `edge_plottype=Makie.automatic()`: Either `automatic`, `:linesegments` or
  `:beziersegments`. `:beziersegments` are much slower for big graphs!

Self edges / loops:

- `selfedge_size=Makie.automatic()`: Size of self-edge-loop (dict/vector possible).
- `selfedge_direction=Makie.automatic()`: Direction of self-edge-loop as `Point2` (dict/vector possible).
- `selfedge_width=Makie.automatic()`: Opening of selfloop in rad (dict/vector possible).

High level interface for curvy edges:

- `curve_distance=0.1`:

    Specify a distance of the (now curved) line to the straight line *in data
    space*. Can be single value, array or dict. User proivded `tangents` or
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

Waypoints along edges:
- `waypoints=nothing`

    Specify waypoints for edges. This parameter should be given as a vector or
    dict. Waypoints will be crossed using natural cubic splines. The waypoints may
    or may not include the src/dst positions.

- `waypoint_radius=nothing`: If number (dict/vector possible) bent lines within radius of waypoints.
"""
@recipe(GraphPlot, graph) do scene
    # TODO: figure out this whole theme business
    scatter_theme = default_theme(scene, Scatter)
    lineseg_theme = default_theme(scene, LineSegments)
    labels_theme = default_theme(scene, Makie.Text)
    Attributes(
        layout = Spring(),
        # node attributes (Scatter)
        node_color = scatter_theme.color,
        node_size = scatter_theme.markersize,
        node_marker = scatter_theme.marker,
        node_attr = (;),
        # edge attributes (LineSegements)
        edge_color = lineseg_theme.color,
        edge_width = lineseg_theme.linewidth,
        edge_attr = (;),
        # arrow attributes (Scatter)
        arrow_show = automatic,
        arrow_size = scatter_theme.markersize,
        arrow_shift = 0.5,
        arrow_attr = (;),
        # node label attributes (Text)
        nlabels = nothing,
        nlabels_align = (:left, :bottom),
        nlabels_distance = 0.0,
        nlabels_color = labels_theme.color,
        nlabels_offset = nothing,
        nlabels_textsize = labels_theme.textsize,
        nlabels_attr = (;),
        # edge label attributes (Text)
        elabels = nothing,
        elabels_align = (:center, :bottom),
        elabels_distance = 0.0,
        elabels_shift = 0.5,
        elabels_opposite = Int[],
        elabels_rotation = automatic,
        elabels_offset = nothing,
        elabels_color = labels_theme.color,
        elabels_textsize = labels_theme.textsize,
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

    # create initial vertex positions, will be updated on changes to graph or layout
    # make node_position-Observable available as named attribute from the outside
    gp[:node_pos] = @lift [Pointf(p) for p in ($(gp.layout))($graph)]

    node_pos = gp[:node_pos]

    # create array of pathes triggered by node_pos changes
    # in case of a graph change the node_position will change anyway
    edge_paths = lift(node_pos, gp.selfedge_size,
                      gp.selfedge_direction, gp.selfedge_width) do pos, s, d, w
        find_edge_paths(graph[], gp.attributes, pos)
    end

    sc = Makie.parent_scene(gp)

    # function which projects the point in px space
    to_px = lift(sc.px_area, sc.camera.projectionview) do pxa, pv
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

    # plot edges
    edge_plot = edgeplot!(gp, edge_paths;
                          plottype=gp.edge_plottype,
                          color=gp.edge_color,
                          linewidth=gp.edge_width,
                          gp.edge_attr...)

    # plott arrow heads
    arrow_pos = @lift if !isempty($edge_paths)
        broadcast(interpolate, $edge_paths, $(gp.arrow_shift))
    else # if no edges return (empty) vector of points, broadcast yields Vector{Any} which can't be plotted
        Vector{eltype(node_pos[])}()
    end
    arrow_rot = @lift Billboard(broadcast($to_angle, edge_paths[], $arrow_pos, gp.arrow_shift[]))
    arrow_show = @lift $(gp.arrow_show) === automatic ? Graphs.is_directed($graph) : $(gp.arrow_show)
    arrow_heads = scatter!(gp,
                           arrow_pos;
                           marker = '➤',
                           markersize = gp.arrow_size,
                           color = gp.edge_color,
                           rotations = arrow_rot,
                           strokewidth = 0.0,
                           markerspace = Pixel,
                           visible = arrow_show,
                           gp.arrow_attr...)

    # plot vertices
    vertex_plot = scatter!(gp, node_pos;
                           color=gp.node_color,
                           marker=gp.node_marker,
                           markersize=gp.node_size,
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

        nlabels_plot = text!(gp, gp.nlabels;
                             position=positions,
                             align=gp.nlabels_align,
                             color=gp.nlabels_color,
                             offset=offset,
                             textsize=gp.nlabels_textsize,
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
        rotation = @lift begin
            if $(gp.elabels_rotation) isa Real
                # fix rotation to a single angle
                rot = $(gp.elabels_rotation)
            elseif $(gp.elabels_rotation) == automatic
                #point the labels up
                rot = broadcast($to_angle, edge_paths[], $positions, gp.elabels_shift[])
                upsidedownlabels = (rot .> π/2) .| (rot .< -π/2)
                rot[upsidedownlabels] .+= π
                gp.elabels_align[] isa Vector || (gp.elabels_align[] = fill(gp.elabels_align[], ne(gp.graph[]) ))
                begin
                    local dtopbot = Dict(:top => :bottom, :bottom => :top)
                    gp.elabels_align[] = [(x,upsidedownlabels[i] ? dtopbot[y] : y) for (i,(x,y)) in enumerate(gp.elabels_align[])]
                end
            else
                # determine rotation for each position
                rot = broadcast($to_angle, edge_paths[], $positions, gp.elabels_shift[])
                for i in $(gp.elabels_opposite)
                    rot[i] += π
                end
                # if there are user provided rotation for some labels, use those
                if $(gp.elabels_rotation) isa Vector
                    for (i, α) in enumerate($(gp.elabels_rotation))
                        α !== nothing && (rot[i] = α)
                    end
                end
            end
            return rot
        end

        # calculate the offset in pixels in normal direction to the edge
        offsets = @lift begin
            tangent_px = broadcast(edge_paths[], $positions, gp.elabels_shift[]) do path, p0, t
                p1 = p0 + tangent(path, t)
                $to_px(p1) - $to_px(p0)
            end

            offsets = map(p -> Point(-p.data[2], p.data[1])/norm(p), tangent_px)
            offsets .= $(gp.elabels_distance) .* offsets
            for i in $(gp.elabels_opposite)
                offsets[i] = -1.0 * offsets[i] # flip offset if in opposite
            end
            offsets
        end

        elabels_plot = text!(gp, gp.elabels;
                             position=positions,
                             rotation=rotation,
                             offset=offsets,
                             align=gp.elabels_align,
                             color=gp.elabels_color,
                             textsize=gp.elabels_textsize,
                             gp.elabels_attr...)
    end

    return gp
end


"""
    find_edge_paths(g, attr, pos::AbstractVector{PT}) where {PT}

Returns an `AbstractPath` for each edge in the graph.
"""
function find_edge_paths(g, attr, pos::AbstractVector{PT}) where {PT}
    paths = Vector{AbstractPath{PT}}(undef, ne(g))

    ne(g) == 0 && return paths

    for (i, e) in enumerate(edges(g))
        if e.src == e.dst # selfedge
            size = getattr(attr.selfedge_size, i)
            direction = getattr(attr.selfedge_direction, i)
            width = getattr(attr.selfedge_width, i)
            paths[i] = selfedge_path(g, pos, e.src, size, direction, width)
        else # no selfedge
            p1, p2 = pos[e.src], pos[e.dst]
            tangents = getattr(attr.tangents, i)
            tfactor = getattr(attr.tfactor, i)
            waypoints::Vector{PT} = getattr(attr.waypoints, i, PT[])

            cdu = getattr(attr.curve_distance_usage, i)
            if cdu === true
                curve_distance = getattr(attr.curve_distance, i, 0.0)
            elseif cdu === false
                curve_distance = 0.0
            elseif cdu === automatic
                if is_directed(g) && has_edge(g, e.dst, e.src)
                    curve_distance = getattr(attr.curve_distance, i, 0.0)
                else
                    curve_distance = 0.0
                end
            end

            if !isnothing(waypoints) && !isempty(waypoints) #there are waypoints
                # the waypoints may already include the endpoints
                waypoints[begin] == p1 && popfirst!(waypoints)
                waypoints[end] == p2 && pop!(waypoints)

                radius = getattr(attr.waypoint_radius, i, nothing)

                if isempty(waypoints) || radius === nothing || radius === :spline
                    paths[i] = Path(p1, waypoints..., p2; tangents, tfactor)
                elseif radius isa Real
                    paths[i] = Path(radius, p1, waypoints..., p2)
                else
                    throw(ArgumentError("Invalid radius $radius for edge $i!"))
                end
            elseif !isnothing(tangents)
                paths[i] = Path(p1, p2; tangents, tfactor)
            elseif PT<:Point2 && !iszero(curve_distance)
                d = curve_distance
                s = norm(p2 - p1)
                γ = 2*atan(2 * d/s)
                a = (p2 - p1)/s * (4*d^2 + s^2)/(3s)

                m = @SMatrix[cos(γ) -sin(γ); sin(γ) cos(γ)]
                c1 = PT(p1 + m*a)
                c2 = PT(p2 - transpose(m)*a)

                commands = [MoveTo(p1), CurveTo(c1, c2, p2)]
                paths[i] = BezierPath(commands)
            else # straight line
                paths[i] = Path(p1, p2)
            end
        end
    end

    # try to narrow down the typ, i.e. just `Line`s
    T = mapreduce(typeof, promote_type, paths)
    return convert(Vector{T}, paths)
end

"""
    selfedge_path(g, pos, v, size, direction, width)

Return a Path for the
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
    edgeplot(paths::Vector{AbstractPath})
    edgeplot!(sc, paths::Vector{AbstractPath})

Recipe to draw the edges. Attribute `plottype` can be either

- `:linesegments`: Draw edges as `linesegments` (just straight lines)
- `:beziersegments`: Draw edges as `beziersegments` (separate `lines` plot for
  each edge. Much slower than `linesegments`!)
- `Makie.automatic()`: Only use `beziersegements` if there is at least one curvy line.
"""
@recipe(EdgePlot, paths) do scene
    Attributes(
        plottype = automatic
    )
end

function Makie.plot!(p::EdgePlot)
    N = length(p[:paths][])
    # remove plottype from attributes, otherwise attributs can't be passed on
    plottype = pop!(p.attributes, :plottype)[]
    if plottype === automatic
        justlines = all(isline, p[:paths][])

        if justlines
            plottype = :linesegments
        elseif N > 500
            @warn "Since there are a lot of edges ($N), they will be drawn as straight lines "*
                "even though they contain curvy edges. If you realy want to plot them as "*
                "bezier curves pass `edge_plottype=:beziersegments` explicitly. This will have "*
                "much worse performance!"
            plottype = :linesegments
        else
            plottype = :beziersegments
        end
    end

    if plottype === :linesegments
        PT = ptype(eltype(p[:paths][]))
        segs = Observable(Vector{PT}(undef, 2*N))

        update_segments!(segs, p[:paths][]) # first call to initialize
        on(p[:paths]) do paths # update if pathes change
            update_segments!(segs, paths)
        end
        linesegments!(p, segs; p.attributes...)
    elseif plottype === :beziersegments
        beziersegments!(p, p[:paths]; p.attributes...)
    else
        error("Unknown plottype $(repr(plottype)) for edge plot")
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
