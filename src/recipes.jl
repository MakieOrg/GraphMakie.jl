
function AbstractPlotting.plot!(graph_plot::Combined{Any, S} where S<:Tuple{<: AbstractGraph})
    # create initial vertex positions
    graph = graph_plot[1]

    layout_node = map(graph) do graph
        A = adjacency_matrix(graph)
        return [Point2f0(p[1], p[2]) for p in NetworkLayout.Spring.layout(A)]
    end

    # plot vertices
    vertex_plot = scatter!(graph_plot, layout_node)
        # plot edges
    lines = map(vertex_plot[1]) do x
        indices = vec([getfield(e, src) for src in (:src, :dst), e in edges(graph[])])
        return view(x, indices)
    end
    linesegments!(graph_plot, lines)
    # Plot vertices over lines
    reverse!(graph_plot.plots)
    return graph_plot
end
