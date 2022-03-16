#=
# Add interactions to your graph plot

In this example you will see, how to register interactions with your graph plot.
This tutorial will make use of the more basic [Interaction Interface](@ref).
If you just want to move nodes check out the [Predefined Interactions](@ref). The
implementation of those is quit similar to what is shown in this tutorial.

We star with a simple wheel graph again. This time we use arrays for some attributes
because we want to change them later in the interactions for individual nodes/edges.
=#
using CairoMakie
CairoMakie.activate!(type="png") # hide
set_theme!(resolution=(800, 400)) #hide
CairoMakie.inline!(true) # hide
using GraphMakie
using Graphs
using CairoMakie.Colors

g = wheel_graph(10)
f, ax, p = graphplot(g,
                     edge_width = [2.0 for i in 1:ne(g)],
                     edge_color = [colorant"gray" for i in 1:ne(g)],
                     node_size = [10 for i in 1:nv(g)],
                     node_color = [colorant"red" for i in 1:nv(g)])
hidedecorations!(ax); hidespines!(ax)
ax.aspect = DataAspect()
@save_reference f # hide

# Later on we want to enable drag interactions, therefore we disable the default
# `:rectanglezoom` interaction
deregister_interaction!(ax, :rectanglezoom)
nothing #hide

# ## Hover interactions
# At first, let's add some hover interaction for our nodes using the [`NodeHoverHandler`](@ref)
# constructor. We need to define a action function with the signature `fun(state, idx, event, axis)`.
# We use the action to make the nodes bigger on hover events.
function node_hover_action(state, idx, event, axis)
    @info idx #hide
    p.node_size[][idx] = state ? 20 : 10
    p.node_size[] = p.node_size[] # trigger observable
end
nhover = NodeHoverHandler(node_hover_action)
register_interaction!(ax, :nhover, nhover)

function set_cursor!(p) #hide
    direction = Point2f(-0.1, 0.2) #hide
    arrows!([p-direction], [direction], linewidth=3, arrowsize=20, lengthscale=0.8) #hide
end #hide
nodepos = copy(p[:node_pos][]) #hide
set_cursor!(nodepos[5] + Point2f(0.05, 0)) #hide
p.node_size[][5] = 20; p.node_size[] = p.node_size[] #hide
@save_reference f #hide

# Please run the script locally with `GLMakie.jl` if you want to play with the Graph 🙂
# The edge hover interaction is quite similar:

pop!(ax.scene.plots) #hide
p.node_size[][5] = 10; p.node_size[] = p.node_size[] #hide
function edge_hover_action(state, idx, event, axis)
    @info idx #hide
    p.edge_width[][idx]= state ? 5.0 : 2.0
    p.edge_width[] = p.edge_width[] # trigger observable
end
ehover = EdgeHoverHandler(edge_hover_action)
register_interaction!(ax, :ehover, ehover)

set_cursor!((nodepos[4]+nodepos[1])/2) #hide
p.edge_width[][3] = 5.0; p.edge_width[] = p.edge_width[] #hide
@save_reference f #hide

# ## Click interactions
# In a similar fashion we might change the color of nodes and lines by click.
function node_click_action(idx, args...)
    p.node_color[][idx] = rand(RGB)
    p.node_color[] = p.node_color[]
end
nclick = NodeClickHandler(node_click_action)
register_interaction!(ax, :nclick, nclick)

function edge_click_action(idx, args...)
    p.edge_color[][idx] = rand(RGB)
    p.edge_color[] = p.edge_color[]
end
eclick = EdgeClickHandler(edge_click_action)
register_interaction!(ax, :eclick, eclick)

p.edge_color[][3] = colorant"blue"; p.edge_color[] = p.edge_color[] #hide
p.node_color[][7] = colorant"yellow" #hide
p.node_color[][2] = colorant"brown" #hide
p.node_color[][9] = colorant"pink" #hide
p.node_color[][6] = colorant"green" #hide
p.node_color[] = p.node_color[] #hide
@save_reference f #hide

# ## Drag interactions
pop!(ax.scene.plots) #hide
p.edge_width[][3] = 2.0; p.edge_width[] = p.edge_width[] #hide
function node_drag_action(state, idx, event, axis)
    p[:node_pos][][idx] = event.data
    p[:node_pos][] = p[:node_pos][]
end
ndrag = NodeDragHandler(node_drag_action)
register_interaction!(ax, :ndrag, ndrag)

p[:node_pos][][1] = nodepos[1] + Point2f(1.0,0.5) #hide
p[:node_pos][] = p[:node_pos][] #hide
set_cursor!(p[:node_pos][][1] + Point2f(0.05, 0)) #hide
p.node_size[][1] = 20; p.node_size[] = p.node_size[] #hide
@save_reference f # hide

# The last example is not as straight forward. By dragging an edge we want to
# change the positions of both attached nodes. Therefore we need some more state
# inside the action. We can achieve this with a callable struct.
pop!(ax.scene.plots) #hide
p.node_size[][1] = 10; p.node_size[] = p.node_size[] #hide
mutable struct EdgeDragAction
    init::Union{Nothing, Point2f} # save click position
    src::Union{Nothing, Point2f}  # save src vertex position
    dst::Union{Nothing, Point2f}  # save dst vertex position
    EdgeDragAction() = new(nothing, nothing, nothing)
end
function (action::EdgeDragAction)(state, idx, event, axis)
    edge = collect(edges(g))[idx]
    if state == true
        if action.src===action.dst===action.init===nothing
            action.init = event.data
            action.src = p[:node_pos][][edge.src]
            action.dst = p[:node_pos][][edge.dst]
        end
        offset = event.data - action.init
        p[:node_pos][][edge.src] = action.src + offset
        p[:node_pos][][edge.dst] = action.dst + offset
        p[:node_pos][] = p[:node_pos][] # trigger change
    elseif state == false
        action.src = action.dst = action.init =  nothing
    end
end
edrag = EdgeDragHandler(EdgeDragAction())
register_interaction!(ax, :edrag, edrag)

p[:node_pos][][9] = nodepos[9] + Point2f(0.9,1.0) #hide
p[:node_pos][][10] = nodepos[10] + Point2f(0.9,1.0) #hide
p[:node_pos][] = p[:node_pos][] #hide
pm = (p[:node_pos][][9] + p[:node_pos][][10])/2 #hide
set_cursor!(pm) #hide
p.edge_width[][18] = 5.0; p.edge_width[] = p.edge_width[] #hide
@save_reference f # hide
