#=
# Images for Reference Testing
This is a more or less uncurated list of figures which are used for reference testing.

## Node label positioning
=#
using Graphs, GraphMakie, CairoMakie, NetworkLayout
CairoMakie.activate!(type="png") # hide
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

fig, ax, p = graphplot(g; layout=SquareGrid(), nlabels, nlabels_align)
graphplot!(g; layout=SquareGrid(), nlabels, nlabels_align, nlabels_distance=60, nlabels_color=:red)
hidedecorations!(ax); xlims!(-2,4); ylims!(-4,2)
@save_reference fig
