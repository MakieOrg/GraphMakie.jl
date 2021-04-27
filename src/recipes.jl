@recipe(GraphPlot, graph) do scene
    Attributes(
        layout = NetworkLayout.Spring.layout,
        nodecolor = :red,
        nodesize = 20,
        marker = :circle,
        edgecolor = :black,
        edgewidth = 2.5,
    )
end

function AbstractPlotting.plot!(gp::GraphPlot)
    graph = gp[:graph]

    # create initial vertex positions
    node_points = @lift begin
        A = adjacency_matrix($graph)
        ($(gp.layout))(A) # array of Point2f0
    end

    edge_segments = @lift begin
        indices = vec([getfield(e, s) for s in (:src, :dst), e in edges($graph)])
        view($node_points, indices)
    end

    # plot edges
    edge_plot = linesegments!(gp, edge_segments,
                              color=gp.edgecolor,
                              linewidth=gp.edgewidth)

    # plot vertices
    vertex_plot = scatter!(gp, node_points,
                           color=gp.nodecolor,
                           marker=gp.marker,
                           markersize=gp.nodesize)
    return gp
end
