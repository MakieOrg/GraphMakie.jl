#=
# Images for Reference Testing

This is a more or less uncurated list of figures which are used for reference
testing. They might be interesting but they are most probably not.
=#
using Graphs, GraphMakie, CairoMakie, NetworkLayout
CairoMakie.activate!(type="png") # hide
set_theme!(resolution=(400, 400)) #hide
CairoMakie.inline!(true) # hide
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
nlabels_textsize = 10

fig, ax, p = graphplot(g; layout=SquareGrid(), nlabels, nlabels_align, nlabels_textsize)
graphplot!(g; layout=SquareGrid(), nlabels, nlabels_align, nlabels_distance=30,
           nlabels_color=:red, nlabels_textsize)
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
p[:nlabels][] = ["𐄂" for i in 1:nv(g)]
p[:nlabels_color][] = :red
p[:nlabels_offset][] = [Point2(.1*cos(-2π/9*i),.1*sin(-2π/9*i)) for i in 1:nv(g)]
@save_reference fig

# draw some edge labels
g = path_graph(4)
elabels = ["𐄂" for i in 1:ne(g)]
elabels_align = (:center, :center)
fig, ax, p = graphplot(g; layout=SquareGrid(), elabels, elabels_align,
                       elabels_color=:red, elabels_rotation=nothing)
hidedecorations!(ax)
@save_reference fig

# change the positioning of the edge labels
p[:elabels][] = repr.(edges(g))
p[:elabels_shift][] = [0.25, 0.5, 0.75]
p[:elabels_rotation][] = [π/8, 0, -π/8]
p[:elabels_offset][] = Point2(0.1,0.1)
p[:elabels_textsize][] = 10
autolimits!(ax)
@save_reference fig

# edge and node labels follow graph movement
g = complete_digraph(3)
elabels = repr.(edges(g))
nlabels = repr.(1:nv(g))
fig, ax, p = graphplot(g; elabels, nlabels, elabels_textsize=10, elabels_rotation=nothing)
@save_reference fig
p[:node_pos][] = Point2f.([(1., -.5), (-1.,0.), (-1.,-1.)])
@save_reference fig

# dynamical changes of the underlying graph
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
