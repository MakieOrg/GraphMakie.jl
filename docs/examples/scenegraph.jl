#=
# Scene graph in Makie
In Makie, plots and scenes are stored as a tree.  Scenes can hold child Scenes and Plots, and plots can hold other plots.

In this example, we create a simple plot (`streamplot`) and render the scene graph for it.
=#

using CairoMakie
CairoMakie.activate!(type="svg") #hide
using GraphMakie
using Graphs

# ## Extracting the scene graph

# First, we extract the scene graph by "walking" down the tree.  This uses multiple dispatch to dispatch based on scenes and plots.
# Scenes and Plots can both hold plots, and the leaf nodes can be either Scenes with no plots (empty Scenes) or atomic plots, i.e., 
# plots which can be rendered directly by the backend.

# This function simply initializes the graph and labels, and begins the traversal.
function walk_tree(scene)
    g = SimpleDiGraph()
    labels = Any[]
    walk_tree!(g, labels, scene)
    return (g, labels)
end
nothing #hide

# Now, we can walk down the Scene tree.  Scenes can have child Scenes as well as child Plots, 
# but in terms of semantic order we walk down the Scene tree before looking at the Scene's attached
# plots.
function walk_tree!(g, labels, scene::Scene)
    add_vertex!(g)
    top_vertex = vertices(g)[end]

    push!(labels, label_str(scene))

    for child_scene in scene.children
        child = walk_tree!(g, labels, child_scene)
        add_edge!(g, top_vertex, child)
    end

    for child_plot in scene.plots
        child = walk_tree!(g, labels, child_plot)
        add_edge!(g, top_vertex, child)
    end

    return top_vertex
end

function walk_tree!(g, labels, plot)
    add_vertex!(g)
    top_vertex = vertices(g)[end]

    push!(labels, label_str(plot))

    for child_plot in plot.plots
        child = walk_tree!(g, labels, child_plot)
        add_edge!(g, top_vertex, child)
    end

    return top_vertex
end

## This is a utility function for the label, to avoid
## the cruft that comes from excessive type printing.
label_str(::Scene) = "Scene"
label_str(::Plot{F}) where {F} = string(F) # get only the plot func, not the argument type
nothing #hide


# ## Creating the plot

# This is a simple streamplot in an LScene, which has the simplest axis (Axis3 is more complex!)

fig, ax, plt = streamplot(-2..2, -2..2; axis = (type = LScene,),) do x::Point2
    Point2(x[2], 4x[1])
end

# Let's walk down the tree with our previous `walk_tree` function:

newg, newl = walk_tree(fig.scene)
## This is for convenience later:
nlabels_align = [(:left, :center) for v in vertices(newg)]
nothing #hide
# We start out by plotting the graph itself.
f, a, p = graphplot(
    newg; 
    layout=GraphMakie.Buchheim(),
    nlabels=newl,
    nlabels_distance=10,
    nlabels_fontsize=30,
    nlabels_align,
    tangents=((0,-1),(0,-1)),
    figure = (; size = (900, 600)),
    axis = (limits = (-2.5, 2, -16, 2),)
)
hidedecorations!(a); hidespines!(a)
# Now, we add some improvements to the labels and positions (this is fairly minor):
_nlabels = deepcopy(newl)
_nlabels[1] = "Scene (root)"
_nlabels[2] = "Scene (lscene.blockscene)"
_nlabels[3] = "Scene (LScene)"
p.nlabels[] = _nlabels
fig

for v in vertices(newg)
    if isempty(inneighbors(newg, v)) # root
        nlabels_align[v] = (:center,:bottom)
    elseif isempty(outneighbors(newg, v)) #leaf
        nlabels_align[v] = (:center,:top)
    else
        self = p[:node_pos][][v]
        parent = p[:node_pos][][inneighbors(newg, v)[1]]
        if self[1] < parent[1] # left branch
            nlabels_align[v] = (:right,:center)
        end
    end
end
p.nlabels_align = nlabels_align
nothing #hide
# ### Final figure
@save_reference f
