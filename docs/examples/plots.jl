#=
# Plotting Graphs with `GraphMakie.jl`
Plotting your first `AbstractGraph` from [`LightGraphs.jl`](https://juliagraphs.org/LightGraphs.jl/latest/)
is as simple as
=#
using CairoMakie
CairoMakie.activate!(type="png") # hide
set_theme!(resolution=(800, 400)) #hide
AbstractPlotting.inline!(true) # hide
using GraphMakie
using LightGraphs
import Random; Random.seed!(2) # hide

g = wheel_graph(10)
f, ax, p = graphplot(g)
hidedecorations!(ax); hidespines!(ax)
ax.aspect = DataAspect()
f # hide

#=
The `graphplot` command is a recipe which wraps several steps
- layout the graph in 2D space using a layout function,
- create a `scatter` plot for the nodes and
- create a `linesegments` plot for the edges.

The default layout is `NetworkLayout.Spring.layout` from
[`NetworkLayout.jl`](https://github.com/JuliaGraphs/NetworkLayout.jl). The
layout attribute can be any function which takes the adjacency matrix of the
graph an returns a list of `(x,y)` tuples or `Point2f0` objects.

Besides that there are some common attributes which are forwarded to the
underlying plot commands. See [`graphplot`](@ref).
=#
using GraphMakie.NetworkLayout

g = SimpleGraph(5)
add_edge!(g, 1, 2); add_edge!(g, 2, 4);
add_edge!(g, 4, 3); add_edge!(g, 3, 2);
add_edge!(g, 2, 5); add_edge!(g, 5, 4);
add_edge!(g, 4, 1); add_edge!(g, 1, 5);

## define some edge colors
edgecolors = [:black for i in 1:ne(g)]
edgecolors[4] = edgecolors[7] = :red

f, ax, p = graphplot(g, layout=NetworkLayout.Circular.layout,
                     node_color=[:black, :red, :red, :red, :black],
                     edge_color=edgecolors)

hidedecorations!(ax); hidespines!(ax)
ax.aspect = DataAspect()
f #hide

#=
We can interactively change the attributes as usual with Makie.
=#

fixed_layout(_) = [(0,0), (0,1), (0.5, 1.5), (1,1), (1,0)]
## set new layout
p.layout = fixed_layout; autolimits!(ax)
## change edge width & color
p.edge_width = 5.0
p.edge_color[][3] = :green;
p.edge_color = p.edge_color[] # trigger observable
f #hide

#=
## Nodelabels
=#
Random.seed!(2)
g = wheel_graph(10)

colors = [:black for i in 1:nv(g)]
colors[1] = :red

f, ax, p = graphplot(g,
                     nlabels=repr.(1:nv(g)),
                     nlabels_color=colors,
                     nlabels_align=(:center,:center))
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect()
f # hide

# This is not to nice, lets change the offsets based on the node_positions

offsets = 0.10 * (p[:node_positions][] .- p[:node_positions][][1])
offsets[1] = Point2f0(0, 0.3)
p.nlabels_offset[] = offsets
autolimits!(ax)
f # hide
