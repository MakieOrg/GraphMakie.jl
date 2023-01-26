#=
# Images for Reference Testing

This is a more or less uncurated list of figures which are used for reference
testing. They might be interesting but they are most probably not.
=#
# ## Test different `nlabels_align` modes
using Graphs, GraphMakie, CairoMakie, NetworkLayout
CairoMakie.activate!(type="png") # hide
set_theme!(resolution=(400, 400)) #hide
g = SimpleGraph(9)
nlabels_align = [(:right, :bottom),
                 (:center, :bottom),
                 (:left, :bottom),
                 (:right, :center),
                 (:center, :center),
                 (:left, :center),
                 (:right, :top),
                 (:center, :top),
                 (:left, :top)]

nlabels = repr.(nlabels_align)
nlabels_fontsize = 10

fig, ax, p = graphplot(g; layout=SquareGrid(), nlabels, nlabels_align, nlabels_fontsize)
graphplot!(g; layout=SquareGrid(), nlabels, nlabels_align, nlabels_distance=30,
           nlabels_color=:red, nlabels_fontsize)
hidedecorations!(ax); xlims!(-2,4); ylims!(-4,2)
@save_reference fig

# change the align
pop!(ax.scene.plots) # remove red plot
p[:nlabels_distance] = 10
p[:nlabels_align][] = [(:left, :bottom) for i in 1:nv(g)]
p[:nlabels][] = ["↙" for i in 1:nv(g)]
p[:nlabels_offset][] = Point2(0.1,0.2)
@save_reference fig

# variable offsets
p[:nlabels_distance] = 0
p[:nlabels_align][] = [(:center, :center) for i in 1:nv(g)]
p[:nlabels][] = ["×" for i in 1:nv(g)]
p[:nlabels_color][] = :red
p[:nlabels_offset][] = [Point2(.1*cos(-2π/9*i),.1*sin(-2π/9*i)) for i in 1:nv(g)]
@save_reference fig

# ## Edge label placement
g = path_graph(4)
elabels = ["a" for i in 1:ne(g)]
elabels_align = (:center, :center)
fig, ax, p = graphplot(g; layout=SquareGrid(), elabels, elabels_align,
                       elabels_color=:red)
hidedecorations!(ax)
@save_reference fig

# change the positioning of the edge labels
p[:elabels][] = repr.(edges(g))
p[:elabels_shift][] = [0.25, 0.5, 0.75]
p[:elabels_rotation][] = [π/8, 0, -π/8]
p[:elabels_offset][] = Point2(0.05,0.05)
p[:elabels_fontsize][] = 10
autolimits!(ax)
@save_reference fig

# ## Changes of node positions
# edge and node labels follow graph movement
g = complete_digraph(3)
elabels = repr.(edges(g))
nlabels = repr.(1:nv(g))
fig, ax, p = graphplot(g; elabels, nlabels, elabels_fontsize=10)
@save_reference fig
#
limits!(ax, ax.finallimits[]) # freeze the limits
p[:node_pos][] = Point2f.([(1., -.5), (-1.,0.), (-1.,-1.)])
@save_reference fig

# ## Change of the underlying graph
gn = Observable(SimpleDiGraph(3))
fig, ax, p = graphplot(gn)
hidedecorations!(ax)
@save_reference fig

add_edge!(gn[], 1, 2)
add_edge!(gn[], 1, 3)
add_edge!(gn[], 2, 3)
notify(gn)
autolimits!(ax)
@save_reference fig

# add another edge
add_edge!(gn[], 2, 1)
notify(gn)
@save_reference fig

# ## Different combinations of edge width and edgelabel distance
g = path_graph(3)
layout(y) = _ -> Point2f.([(0,-y),(1,-y),(2,-y)])
elabels = ["Edge 1", "Edge 2"]
node_color = :red

fig, ax, p = graphplot(g; layout=layout(0), elabels, node_color)
graphplot!(g; layout=layout(1), elabels, node_color, edge_width=10, elabels_fontsize=25)
graphplot!(g; layout=layout(2), elabels, node_color, edge_width=25, elabels_fontsize=15)
graphplot!(g; layout=layout(3), elabels, node_color, edge_width=25, elabels_fontsize=[15,25])
graphplot!(g; layout=layout(4), elabels, node_color, edge_width=[10,25], elabels_fontsize=25)
autolimits!(ax); hidedecorations!(ax); hidespines!(ax); ylims!(-5,1)
@save_reference fig

# ## Different linestyles per edge
fig = Figure()
graphplot(fig[1,1],
          DiGraph([Edge(1 => 2), Edge(2 => 3)]),
          edge_attr = (; linestyle = [:dot, :dash]),
          edge_plottype = :beziersegments,
          )

graphplot(fig[1,2],
          DiGraph([Edge(1 => 2), Edge(2 => 1)]),
          edge_attr = (; linestyle = [:dot, :dash]),
          edge_plottype = :beziersegments,
          )

graphplot(fig[2,1],
          DiGraph([Edge(1 => 2), Edge(2 => 3), Edge(3=>4), Edge(4=>1)]),
          edge_attr = (; linestyle = [0.5, 1.0, 1.5, 2.5]),
          edge_plottype = :beziersegments,
          )
@save_reference fig

# ##self loop with waypoints
g1 = SimpleDiGraph(1)
add_edge!(g1, 1, 1) #add self loop
fig, ax, p = graphplot(g1, layout = _ -> [(0,0)], waypoints = [[(1,-1),(1,1),(-1,1),(-1,-1)]])
@save_reference fig

# ##shift arrows to nodes
fig,ax,p=graphplot(SimpleDiGraph(ones(2,2)),node_size=50,arrow_size=20,curve_distance=0.5)
move_arrows_to_nodes!(ax, p; t=0.99)
@save_reference fig