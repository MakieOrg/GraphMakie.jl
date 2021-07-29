using LinearAlgebra: normalize, ⋅, norm
export GraphPlot, graphplot, graphplot!

"""
    graphplot(graph::AbstractGraph)
    graphplot!(ax, graph::AbstractGraph)

Creates a plot of the network `graph`. Consists of multiple steps:
- Layout the nodes: see `layout` attribute. The node position is accesible from outside
  the plot object `p` as an observable using `p[:node_positions]`.
- plot edges as `linesegments`-plot
- if `arrows_show` plot arrowheads as `scatter`-plot
- plot nodes as `scatter`-plot
- if `nlabels!=nothing` plot node labels as `text`-plot
- if `elabels!=nothing` plot edge labels as `text`-plot

The main attributes for the subplots are exposed as attributes for `graphplot`.
Additional attributes for the `scatter`, `linesegments` and `text` plots can be provided
as a named tuples to `node_attr`, `edge_attr`, `nlabels_attr` and `elabels_attr`.

Most of the arguments can be either given as a vector of length of the
edges/nodes or as a single value. One might run into errors when changing the
underlying graph and therefore changing the number of Edges/Nodes.

## Attributes
### Main attributes
- `layout=Spring()`: function `AbstractGraph->Vector{Point}` determines the base layout
- `node_color=scatter_theme.color`
- `node_size=scatter_theme.markersize`
- `node_marker=scatter_theme.marker`k
- `node_attr=(;)`: List of kw arguments which gets passed to the `scatter` command
- `edge_color=lineseg_theme.color`: Pass a vector with 2 colors per edge to get
  color gradients.
- `edge_width=lineseg_theme.linewidth`: Pass a vector with 2 width per edge to
  get pointy edges.
- `edge_attr=(;)`: List of kw arguments which gets passed to the `linesegments` command
- `arrow_show=Makie.automatic`: `Bool`, indicate edge directions with arrowheads?
  Defaults to `true` for `SimpleDiGraph` and `false` otherwise.
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
        elabels_rotation = nothing,
        elabels_offset = nothing,
        elabels_color = labels_theme.color,
        elabels_textsize = labels_theme.textsize,
        elabels_attr = (;),
    )
end

function Makie.plot!(gp::GraphPlot)
    graph = gp[:graph]

    # create initial vertex positions, will be updated on changes to graph or layout
    # make node_position-Observable available as named attribute from the outside
    gp[:node_positions] = @lift [toF32(Point(p)) for p in ($(gp.layout))($graph)]

    node_pos = gp[:node_positions]

    # create array of pathes triggered by node_pos changes
    # in case of a graph change the node_position will change anyway
    edge_paths = @lift find_edge_paths(graph[], $node_pos)

    sc = Makie.parent_scene(gp)

    # function which projects the point in px space
    to_px = lift(sc.px_area, sc.camera.projectionview) do pxa, pv
        # project should transform to 2d point in px space
        (point) -> project(sc, point)
    end
    # get angle in px space from path p at point t
    to_angle = @lift (path, t) -> begin
        p0 = interpolate(path, t)
        p1 = p0 + tangent(path, t)
        tpx = $to_px(p1) - $to_px(p0)
        atan(tpx[2], tpx[1])
    end

    # plot edges
    edge_plot = beziersegments!(gp, edge_paths;
                                color=gp.edge_color,
                                linewidth=gp.edge_width,
                                gp.edge_attr...)

    # plott arrow heads
    arrow_pos = @lift broadcast((p, t) -> interpolate(p, t),
                                $edge_paths, $(gp.arrow_shift))
    arrow_rot = @lift Billboard(broadcast((p, t) -> $to_angle(p, t),
                                          $edge_paths, $(gp.arrow_shift)))
    arrow_show = @lift $(gp.arrow_show) === automatic ? $graph isa SimpleDiGraph : $(gp.arrow_show)
    arrow_heads = scatter!(gp,
                           arrow_pos,
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
        # rotations based on the edge_vec_px and opposite argument
        rotation = @lift begin
            if $(gp.elabels_rotation) isa Real
                # fix rotation to a single angle
                rot = $(gp.elabels_rotation)
            else
                # determine rotation for each position
                rot = broadcast((p, t) -> $to_angle(p, t),
                                $edge_paths, $(gp.elabels_shift))

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

        # positions: center point between nodes + offset + distance*normal + shift*edge direction
        positions = @lift begin
            pos = broadcast((p, t) -> interpolate(p, t), $edge_paths, $(gp.elabels_shift))

            if $(gp.elabels_offset) !== nothing
                pos .= pos .+ $(gp.elabels_offset)
            end
            pos
        end

        # calculate the offset in pixels in normal direction to the edge
        offsets = @lift begin
            tangent_px = broadcast($edge_paths, $(gp.elabels_shift)) do path, t
                p0 = interpolate(path, t)
                p1 = p0 + tangent(path, t)
                tpx = $to_px(p1) - $to_px(p0)
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

toF32(p::Point{2, T}) where T = Point2f0(p)
toF32(p::Point{3, T}) where T = Point3f0(p)

function align_to_dir(align)
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
    return Point(x/norm, y/norm)
end

function find_edge_paths(g, pos::AbstractVector{PT}) where {PT}
    paths = Vector{BezierPath{PT}}(undef, ne(g))
    for (i, e) in enumerate(edges(g))
        if e.src == e.dst
            paths[i] = selfedge_path(g, pos, e.src)
        else
            paths[i] = BezierPath(pos[e.src], pos[e.dst])
        end
    end
    return paths
end

function selfedge_path(g, pos::AbstractVector{Point2f0}, v)
    vp = pos[v]
    # get the vectors to all the neighbors
    ndirs = [pos[n] - vp for n in neighbors(g, v) if n != v]

    angles = SVector{length(ndirs)}(atan(p[2], p[1]) for p in ndirs)
    angles = sort(angles)

    # search biggest gap between neighbor edges
    γ, Δ = 0.0, 0.0
    for i in 1:length(angles)
        α = angles[i]
        β = get(angles, i+1, 2π + angles[1])
        if β-α > Δ
            Δ = (β-α) / 2
            γ = (β+α) / 2
        end
    end

    # set width of selfloop
    Δ = min(.7*Δ, π/4)

    # the size (max dis. to v) of loop
    size = minimum(norm.(ndirs)) * 0.75
    # the actual length of the tagent vectors, magic number from `CurveTo`
    l = Float32( size/cos(Δ) * 1/0.75 )
    t1 = vp + l * Point2f0(cos(γ-Δ), sin(γ-Δ))
    t2 = vp + l * Point2f0(cos(γ+Δ), sin(γ+Δ))

    return BezierPath([MoveTo(vp),
                       CurveTo(t1, t2, vp)])
end

function selfedge_path(g, pos::AbstractVector{Point3f0}, v)
    error("Self edges in 3D not yet supported")
end

"""
    beziersegments(paths::Vector{BezierPath})
    beziersegments!(sc, paths::Vector{BezierPath})

Recipe to draw bezier pathes. Each path will be descritized and ploted with a
separate `lines` plot. Scalar attributes will be used for all subplots. If you
provide vector attributs of same length als the pathes the i-th `lines` subplot
will see the i-th attribute.
"""
@recipe(BezierSegments, paths) do scene
    Attributes(default_theme(scene, Lines)...)
end

function Makie.plot!(p::BezierSegments)
    N = length(p[:paths][])
    PT = ptype(eltype(p[:paths][]))
    attr = p.attributes

    # every argument which is array of length N will be put into the individual attr.
    # and the i-th lines plot gehts the i-th attribute
    individual_attr = Attributes()
    for (k, v) in attr
        if k == :color && v[] isa AbstractArray{<:Number}
            # in case of colors as numbers transfom those
            # https://github.com/JuliaPlots/Makie.jl/blob/e59f4f0785a190146623ff235eeeacc56b04de6c/CairoMakie/src/utils.jl#L100-L113
            numbers = attr[:color][]

            colormap = get(attr, :colormap, nothing) |> to_value |> to_colormap
            colorrange = get(attr, :colorrange, nothing) |> to_value

            if colorrange === Makie.automatic
                attr[:colorrange] = colorrange = extrema(numbers)
            end

            individual_attr[k] = Makie.interpolated_getindex.(Ref(colormap),
                                                              Float64.(numbers), # ints don't work in interpolated_getindex
                                                              Ref(colorrange))
        elseif v[] isa AbstractArray && length(v[]) == N
            individual_attr[k] = v
        end
    end

    # array which holds observables for the discretized values
    # this is needed so the `lines!` subplots will update
    disc = [Observable{Vector{PT}}() for i in 1:N]
    function update_discretized!(disc, pathes)
        for (i, p) in enumerate(pathes)
            disc[i][] = discretize(p)
        end
    end
    update_discretized!(disc, p[:paths][]) # first call to initialize
    on(p[:paths]) do paths # update if pathes change
        update_discretized!(disc, paths)
    end

    # plot all the lines
    for i in 1:N
        specific_attr = Attributes()
        for (k, v) in individual_attr
            specific_attr[k] = @lift $v[i]
        end
        lines!(p, disc[i]; attr..., specific_attr...)
    end

    return p
end
