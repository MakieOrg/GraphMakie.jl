#=
# Plotting Graphs with `GraphMakie.jl`
## The `graphplot` Command
Plotting your first `AbstractGraph` from [`LightGraphs.jl`](https://juliagraphs.org/LightGraphs.jl/latest/)
is as simple as
=#
using CairoMakie
CairoMakie.activate!(type="png") # hide
set_theme!(resolution=(800, 400)) #hide
CairoMakie.inline!(true) # hide
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
## Adding Node Labels
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

# This is not very nice, lets change the offsets based on the `node_positions`

offsets = 0.15 * (p[:node_positions][] .- p[:node_positions][][1])
offsets[1] = Point2f0(0, 0.3)
p.nlabels_offset[] = offsets
autolimits!(ax)
f # hide

# ## Adding Edge Labels
Random.seed!(42)
g = barabasi_albert(6, 2)

labels =  repr.(1:ne(g))

f, ax, p = graphplot(g, elabels=labels,
                     elabels_color=[:black for i in 1:ne(g)],
                     edge_color=[:black for i in 1:ne(g)])
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect()
f # hide

#=
The position of the edge labels is determined by several plot arguments.
Each label is placed in the middle of the edge and rotated to match the edge rotation.

Note: Since the `text` is displayed in the `screen` system, this rotations only really works
for `DataAspect()`! See [the Makie docs](https://makie.juliaplots.org/stable/plotting_functions/text.html).

The rotaion for each label can be overwritten with the `elabels_rotation` argument.
=#
p.elabels_rotation[] = Vector{Union{Nothing, Float64}}(nothing, ne(g))
p.elabels_rotation[][5] = 0.0 # set absolute rotation angle for label 5
p.elabels_rotation[] = p.elabels_rotation[]
nothing #hide

#=
The position of each label can be modified using different arguments.
  - `elabels_opposite` is a vector if edge indices which tells the plot to display labels on the oppisite site.
  - `elabels_offset` will be added to the middle of the edge in axis coordinates
  - `elabels_distance` increses/decreases the normal distance to the edge
  - `elabels_shift` shifts the label along the edge
  - `elabels_align` tells the `text` plot where to place the text in relation to the position
=#
p.elabels_opposite[] = [4,7]

p.elabels_offset[] = [Point2f0(0.0, 0.0) for i in 1:ne(g)]
p.elabels_offset[][5] = Point2f0(-0.4,0)
p.elabels_offset[] = p.elabels_offset[]

p.elabels_shift[] = [0.5 for i in 1:ne(g)]
p.elabels_shift[][4] = 0.7
p.elabels_shift[][3] = 0.6
p.elabels_shift[] = p.elabels_shift[]

p.elabels_distance[] = zeros(ne(g))
p.elabels_distance[][8] = -.3
p.elabels_distance[] = p.elabels_distance[]

f # hide

# Is it a bird?
p.edge_width[] = 3.0
p.elabels_color[] = [:green, :red, :green, :green, :red, :goldenrod1, :green, :goldenrod1]
p.edge_color[] = [:green, :black, :green, :green, :black, :goldenrod1, :green, :goldenrod1]
p.elabels[] =["left", "There\n should\n be an add-\n itional node\n for the wings!", "right", "leg", "O", "beak", "leg", "weird"]
xlims!(ax, (-1.5,2.5))
f #hide
