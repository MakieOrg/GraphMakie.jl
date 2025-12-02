#=
# Visualising Makie's `ComputeGraph`
A `ComputeGraph` is Makies internal representation of the transformations that a
plot performs on its input data to reach a representation that can be handed to
backend for plotting.

In this example, we draw the topology of two `ComputeGraph`. Being able to visualise
the interdependence of some calculations can sometimes help in debugging Makie recipes.
=#

using CairoMakie
CairoMakie.activate!(; type="svg") #hide
using GraphMakie
using Graphs
using NetworkLayout

# ## Converting the `ComputeGraph`

# First, we convert the `ComputeGraph`, which does not implement the `Graph` interface
# to a bipartite `SimpleDiGraph` as well as some edge and node data

function convert_to_simple_graph(compute_graph)
    g = SimpleDiGraph()
    labels = Any[]
    for (k, v) in compute_graph.inputs
        add_vertex!(g)
        push!(labels, (:input, k))
    end
    for (k, v) in compute_graph.outputs
        k in keys(compute_graph.inputs) && continue
        add_vertex!(g)
        push!(labels, (:output, k))
    end
    found_calculations = Any[]
    ## every output is the result of a computation, so no recursion needed.
    for (k, i) in compute_graph.outputs
        calc = i.parent
        if calc isa Makie.ComputePipeline.ComputeEdge
            if !(calc in found_calculations)
                push!(found_calculations, calc)
                calc_id = (:calculation, calc.outputs[1].name)
                add_vertex!(g)
                push!(labels, calc_id)
                for j in calc.inputs
                    src_id = findfirst(labels) do label
                        return label in [(:input, j.name), (:output, j.name)]
                    end
                    add_edge!(g, src_id, nv(g))
                end
                for j in calc.outputs
                    dst_id = findfirst(labels) do label
                        return label in [(:input, j.name), (:output, j.name), (:calculation, j.name)]
                    end
                    add_edge!(g, nv(g), dst_id)
                end
            end
        end
    end
    return (g, labels)
end
nothing #hide

# we define a function which does the conversion and plotting for us

function plot_compute_graph!(ax, compute_graph; kwargs...)
    g, labels = convert_to_simple_graph(compute_graph)
    node_colors = map(labels) do i
        if i[1] == :calculation
            :black
        elseif i[1] == :input
            :darkorange
        else
            :red
        end
    end
    plot = graphplot!(ax, g; node_color=node_colors,
                      edge_color=[labels[e.src][1] == :calculation ? :orange : :darkgrey for e in edges(g)],
                      nlabels=[i[1] == :calculation ? "" : string(i[2]) for i in labels],
                      nlabels_align=(:center, :bottom), kwargs...)
    return plot
end
nothing #hide

# ## Simple ComputeGraph

# First, we demonstrate a simple `ComputeGraph`, constructed from two inputs which map onto two outputs

graph = Makie.ComputeGraph()
Makie.add_input!(graph, :input1, 1)
Makie.add_input!((key, value) -> Float32(value), graph, :input2, 2)

Makie.register_computation!(graph, [:input1, :input2], [:output]) do inputs, changed, cached
    input1, input2 = inputs
    return (input1[] + input2[],)
end

Makie.register_computation!(graph, [:input1, :output], [:output2]) do inputs, changed, cached
    input1, output = inputs
    return (input1[]^2 * output,)
end

# plotting the simple compute graph

f = Figure()
ax = Axis(f[1, 1])
g_plot = plot_compute_graph!(ax, graph; arrow_size=20, arrow_shift=0.9)
f

# ## `ComputeGraph` of a the `graphplot!` recipe

# plotting the compute graph of the above `graphplot!` plot

f = Figure()
ax = Axis(f[1, 1])
plot_compute_graph!(ax, g_plot.attributes; nlabels_fontsize=9, layout=Spring(; C=8.0, iterations=2000))
f