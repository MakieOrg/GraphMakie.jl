"""
    graphplot(graph::AbstractGraph)

Creates a plot of the network `graph`. Consists of multiple steps:
- Layout the nodes: the `layout` attribute is has to be a function `f(adj_matrix)::pos`
  where `pos` is either an array of `Point2f0` or `(x, y)` tuples
- plot edges as `linesegments`-plot
- plot nodes as `scatter`-plot

The main attributes for the subplots are exposed as attributes for `graphplot`.
Additional attributes for the `scatter` or `linesegments` plots can be provided
as a named tuples to `node_attr` and `edge_attr`.

## Attributes
$(ATTRIBUTES)
"""
@recipe(GraphPlot, graph) do scene
    scatter_theme = default_theme(scene, Scatter)
    lineseg_theme = default_theme(scene, LineSegments)
    Attributes(
        layout = NetworkLayout.Spring.layout,
        node_color = scatter_theme.color,
        node_size = scatter_theme.markersize,
        node_marker = scatter_theme.marker,
        node_attr = (;),
        edge_color = lineseg_theme.color,
        edge_width = lineseg_theme.linewidth,
        edge_attr = (;),
    )
end

function AbstractPlotting.plot!(gp::GraphPlot)
    graph = gp[:graph]

    # create initial vertex positions, will be updated on changes to graph or layout
    node_positions = @lift begin
        A = adjacency_matrix($graph)
        [Point2f0(p) for p in ($(gp.layout))(A)]
    end

    # create views into the node_positions, will be updated node_position changes
    # in case of a graph change the node_position will change anyway
    edge_segments = @lift begin
        indices = vec([getfield(e, s) for s in (:src, :dst), e in edges(graph[])])
        ($node_positions)[indices]
    end

    # in case the edge_with is same as the number of edges
    # create a new observable which doubles the values for compat with line segments
    if length(gp.edge_width[]) == ne(graph[])
        lineseg_width = @lift repeat($(gp.edge_width), inner=2)
    else
        lineseg_width = gp.edge_width
    end

    # plot edges
    edge_plot = linesegments!(gp, edge_segments;
                              color=gp.edge_color,
                              linewidth=lineseg_width,
                              gp.edge_attr...)

    # plot vertices
    vertex_plot = scatter!(gp, node_positions;
                           color=gp.node_color,
                           marker=gp.node_marker,
                           markersize=gp.node_size,
                           gp.node_attr...)
    return gp
end
