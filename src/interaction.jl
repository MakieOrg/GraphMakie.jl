using Makie: ScenePlot
import Makie.MakieLayout: registration_setup!, process_interaction

export NodeHoverHandler, EdgeHoverHandler
export NodeClickHandler, EdgeClickHandler
export NodeDragHandler, EdgeDragHandler

"""
    convert_selection(element, idx) = (element, idx)
    convert_selection(element::LineSegments, idx) = (element, Int(idx/2))

In case of `LineSegements` the reported `idx` by `pick` is allways two
times the edge index (since LineSegments have 2 points per Line).
"""
convert_selection(element, idx) = (element, idx)
convert_selection(element::LineSegments, idx) = (element, Int(idx/2))

"""
    abstract type GraphInteraction

This abstract type supertypes other `GraphInteractions`. Some `GraphInteractions` need
references to the scatter and linesegments plot of the `GraphPlot`. By overloding
`set_nodeplot!(T<:GraphInteraction)` and `set_edgeplot!` those references will be set
during `registration_setup!`.
"""
abstract type GraphInteraction end;
set_nodeplot!(::GraphInteraction, plot) = nothing
set_edgeplot!(::GraphInteraction, plot) = nothing

function registration_setup!(parent, inter::GraphInteraction)
    @assert parent isa Axis "GraphInteraction has to be registered to an Axis!"
    # in case of multiple graph plots in one axis this won't work
    gplots = filter(p->p isa GraphPlot, parent.scene.plots)
    @assert length(gplots)==1 "There has to be exactly one GraphPlot in Axis!"

    # in case of multiple plots of type Scatter/Linesegments in GraphPlot this won't work
    scatter = filter(p->p isa Scatter, gplots[1].plots)
    @assert length(scatter)==1 "There has to be exactly one Scatter-Plot in GraphPlot!"
    set_nodeplot!(inter, scatter[1])

    lines = filter(p->p isa LineSegments, gplots[1].plots)
    @assert length(lines)==1 "There has to be exactly one LineSegments-Plot in GraphPlot!"
    set_edgeplot!(inter, lines[1])
end

"""
    mutable struct HoverHandler{P<:ScenePlot, F} <: GraphInteraction

Object to handle hovers on `plot::P`.
"""
mutable struct HoverHandler{P<:ScenePlot, F} <: GraphInteraction
    idx::Union{Nothing, Int}
    plot::Union{Nothing, P}
    fun::F
end
set_nodeplot!(h::HoverHandler{Scatter}, plot) = h.plot = plot
set_edgeplot!(h::HoverHandler{LineSegments}, plot) = h.plot = plot

"""
    NodeHoverHandler(fun)

Initializes `HoverHandler` for nodes. Calls function

    fun(hoverstate, idx, event, axis)

with `hoverstate=true` on hover and `false` at the end of hover. `idx` is the node index.

# Example
```
julia> g = wheel_digraph(10)
julia> f, ax, p = graphplot(g, node_size = [20 for i in 1:nv(g)])
julia> function action(state, idx, event, axis)
           p.node_size[][idx] = state ? 40 : 20
           p.node_size[] = p.node_size[] #trigger observable
       end
julia> register_interaction!(ax, :nodehover, NodeHoverHandler(action))
```
"""
NodeHoverHandler(fun::F) where F = HoverHandler{Scatter, F}(nothing, nothing, fun)

"""
    EdgeHoverHandler(fun)

Initializes `HoverHandler` for edges. Calls function

    fun(hoverstate, idx, event, axis)

with `hoverstate=true` on hover and `false` at the end of hover. `idx` is the edge index.

# Example
```
julia> g = wheel_digraph(10)
julia> f, ax, p = graphplot(g, edge_width = [3.0 for i in 1:ne(g)])
julia> function action(state, idx, event, axis)
           p.edge_width[][idx] = state ? 6.0 : 3.0
           p.edge_width[] = p.edge_width[] #trigger observable
       end
julia> register_interaction!(ax, :edgehover, EdgeHoverHandler(action))
```
"""
EdgeHoverHandler(fun::F) where F = HoverHandler{LineSegments, F}(nothing, nothing, fun)

function process_interaction(handler::HoverHandler, event::MouseEvent, axis)
    if event.type === MouseEventTypes.over
        (element, idx) = convert_selection(mouse_selection(axis.scene)...)
        if element == handler.plot
            if handler.idx === nothing
                handler.idx = idx
                handler.fun(true, handler.idx, event, axis)
            end
        else
            if handler.idx !== nothing
                handler.fun(false, handler.idx, event, axis)
                handler.idx = nothing
            end
        end
    end
end


"""
    mutable struct DragHandler{P<:ScenePlot, F} <: GraphInteraction

Object to handle left mous drags on `plot::P`.
"""
mutable struct DragHandler{P<:ScenePlot, F} <: GraphInteraction
    idx::Union{Nothing, Int}
    plot::Union{Nothing, P}
    fun::F
end
set_nodeplot!(h::DragHandler{Scatter}, plot) = h.plot = plot
set_edgeplot!(h::DragHandler{LineSegments}, plot) = h.plot = plot

"""
    NodeDragHandler(fun)

Initializes `DragHandler` for Nodes. Calls function

    fun(dragstate, idx, event, axis)

where `dragstate=true` during the drag and `false` at the end of the drag,
the last time `fun` is triggered. `idx` is the node index.

# Example
```
julia> g = wheel_digraph(10)
julia> f, ax, p = graphplot(g, node_size=20)
julia> deregister_interaction!(ax, :rectanglezoom)
julia> function action(state, idx, event, axis)
           p[:node_positions][][idx] = event.data
           p[:node_positions][] = p[:node_positions][]
       end
julia> register_interaction!(ax, :nodedrag, NodeDragHandler(action))
```
"""
NodeDragHandler(fun::F) where F = DragHandler{Scatter, F}(nothing, nothing, fun)

"""
    EdgeDragHandler(fun)

Initializes `DragHandler` for Edges. Calls function

    fun(dragstate, idx, event, axis)

where `dragstate=true` during the drag and `false` at the end of the drag,
the last time `fun` is triggered. `idx` is the edge index.

# Example
```
julia> g = wheel_digraph(10)
julia> f, ax, p = graphplot(g, edge_width=3)
julia> deregister_interaction!(ax, :rectanglezoom)
julia> mutable struct EdgeDragAction
           init::Union{Nothing, Point2f0} # save click position
           src::Union{Nothing, Point2f0}  # save src vertex position
           dst::Union{Nothing, Point2f0}  # save dst vertex position
           EdgeDragAction() = new(nothing, nothing, nothing)
       end
julia> function (action::EdgeDragAction)(state, idx, event, axis)
           edge = collect(edges(g))[idx]
           if state == true
               if action.src===action.dst===action.init===nothing
                   action.init = event.data
                   action.src = p[:node_positions][][edge.src]
                   action.dst = p[:node_positions][][edge.dst]
               end
               offset = event.data - action.init
               p[:node_positions][][edge.src] = action.src + offset
               p[:node_positions][][edge.dst] = action.dst + offset
               p[:node_positions][] = p[:node_positions][] # trigger change
           elseif state == false
               action.src = action.dst = action.init =  nothing
           end
       end
julia> handler = EdgeDragHandler(EdgeDragAction())
julia> register_interaction!(ax, :edgedrag, handler)
```
"""
EdgeDragHandler(fun::F) where F = DragHandler{LineSegments, F}(nothing, nothing, fun)

function process_interaction(handler::DragHandler, event::MouseEvent, axis)
    if handler.idx === nothing # not in drag state
        if event.type === MouseEventTypes.leftdragstart && handler.idx === nothing
            # TODO: idealy this would take the position of the most recent leftdown event!
            # however i am not quite sure how to use mouse_selection on px or data points
            (element, idx) = convert_selection(mouse_selection(axis.scene)...)
            if element === handler.plot
                handler.idx = idx
            end
        end
    elseif handler.idx !== nothing # drag state
        if event.type === MouseEventTypes.leftdrag
            handler.fun(true, handler.idx, event, axis)
        elseif event.type === MouseEventTypes.leftdragstop
            handler.fun(false, handler.idx, event, axis)
            handler.idx = nothing
        end
    end
end


"""
    mutable struct ClickHandler{P<:ScenePlot, F} <: GraphInteraction

Object to handle left mouse clicks on `plot::P`.
"""
mutable struct ClickHandler{P<:ScenePlot, F} <: GraphInteraction
    plot::Union{Nothing, P}
    fun::F
end
set_nodeplot!(h::ClickHandler{Scatter}, plot) = h.plot = plot
set_edgeplot!(h::ClickHandler{LineSegments}, plot) = h.plot = plot

"""
    NodeClickHandler(fun)

Initializes `ClickHandler` for nodes. Calls function

    fun(idx, event, axis)

on left-click events where `idx` is the node index.

# Example
```
julia> using Makie.Colors
julia> g = wheel_digraph(10)
julia> f, ax, p = graphplot(g, node_size=30, node_color=[colorant"red" for i in 1:nv(g)])
julia> function action(idx, event, axis)
           p.node_color[][idx] = rand(RGB)
           p.node_color[] = p.node_color[]
       end
julia> register_interaction!(ax, :nodeclick, NodeClickHandler(action))
```
"""
NodeClickHandler(fun::F) where F = ClickHandler{Scatter, F}(nothing, fun)

"""
    EdgeClickHandler(fun)

Initializes `ClickHandler` for edges. Calls function

    fun(idx, event, axis)

on left-click events where `idx` is the edge index.

# Example
```
julia> using Makie.Colors
julia> g = wheel_digraph(10)
julia> f, ax, p = graphplot(g, edge_width=4, edge_color=[colorant"black" for i in 1:ne(g)])
julia> function action(idx, event, axis)
           p.edge_color[][idx] = rand(RGB)
           p.edge_color[] = p.edge_color[]
       end
julia> register_interaction!(ax, :edgeclick, EdgeClickHandler(action))
```
"""
EdgeClickHandler(fun::F) where F = ClickHandler{LineSegments, F}(nothing, fun)

function process_interaction(handler::ClickHandler, event::MouseEvent, axis)
    if event.type === MouseEventTypes.leftclick
        (element, idx) = convert_selection(mouse_selection(axis.scene)...)
        if element == handler.plot
            handler.fun(idx, event, axis)
        end
    end
end
