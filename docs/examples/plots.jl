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
- layout the graph in space using a layout function,
- create a `scatter` plot for the nodes and
- create a `linesegments` plot for the edges.

The default layout is `Spring()` from
[`NetworkLayout.jl`](https://github.com/JuliaGraphs/NetworkLayout.jl). The
layout attribute can be any function which takes an `AbstractGraph` and returns
a list of `Point{dim,Ptype}` (see [`GeometryBasics.jl`](https://github.com/JuliaGeometry/GeometryBasics.jl)
objects where `dim` determines the dimensionality of the plot.

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

f, ax, p = graphplot(g, layout=Shell(),
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
All possible arguments are described in the docs of the [`graphplot`](@ref) function.

Basicially, each label is placed in the middle of the edge and rotated to match the edge rotation.
The rotaion for each label can be overwritten with the `elabels_rotation` argument.
=#
p.elabels_rotation[] = Vector{Union{Nothing, Float64}}(nothing, ne(g))
p.elabels_rotation[][5] = 0.0 # set absolute rotation angle for label 5
p.elabels_rotation[] = p.elabels_rotation[]
nothing #hide

#=
One can shift the label along the edge with the `elabels_shift` argument and determine the distance
in pixels using the `elabels_distance` argument.
=#
p.elabels_opposite[] = [1,2,8,6]

p.elabels_offset[] = [Point2f0(0.0, 0.0) for i in 1:ne(g)]
p.elabels_offset[][5] = Point2f0(-0.4,0)
p.elabels_offset[] = p.elabels_offset[]

p.elabels_shift[] = [0.5 for i in 1:ne(g)]
p.elabels_shift[][1] = 0.6
p.elabels_shift[][7] = 0.4
p.elabels_shift[] = p.elabels_shift[]

p.elabels_distance[] = zeros(ne(g))
p.elabels_distance[][8] = 15
p.elabels_distance[] = p.elabels_distance[]

f # hide

#=
## Indicate Edge Direction

It is possible to put arrows on the edges using the `arrow_show` parameter. This parameter
is `true` for `SimpleDiGraph` by default. The position and size of each arrowhead can be
change using the `arrow_shift` and `arrow_size` parameters.
=#
g = wheel_digraph(10)
arrow_size = [10+i for i in 1:ne(g)]
arrow_shift = range(0.1, 0.8, length=ne(g))
f, ax, p = graphplot(g; arrow_size, arrow_shift)
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect()
f # hide

#=
## Plot Graphs in 3D
If the layout returns points in 3 dimensions, the plot will be in 3D. However this is a bit
experimental. Feel free to file an issue if there are any problems.
=#
set_theme!(resolution=(800, 800)) #hide
g = smallgraph(:cubical)
elabels_shift = [0.5 for i in 1:ne(g)]
elabels_shift[[2,7,8,9]] .= 0.3
elabels_shift[10] = 0.25
graphplot(g; layout=Spring(dim=3, seed=5),
          elabels="Edge ".*repr.(1:ne(g)),
          elabels_textsize=12,
          elabels_opposite=[3,5,7,8,12],
          elabels_shift,
          elabels_distance=3,
          arrow_show=true,
          arrow_shift=0.9,
          arrow_size=15)

#=
Using [`JSServe.jl`](https://github.com/SimonDanisch/JSServe.jl) and [`WGLMakie.jl`](https://github.com/JuliaPlots/WGLMakie.jl)
we can also add some interactivity:
=#
using JSServe
Page(exportable=true, offline=true)
#
using WGLMakie
WGLMakie.activate!()
set_theme!(resolution=(800, 600))
g = smallgraph(:dodecahedral)
graphplot(g, layout=Spring(dim=3), node_size=100)
