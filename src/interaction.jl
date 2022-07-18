using Makie: ScenePlot
import Makie: registration_setup!, process_interaction

export NodeHoverHandler, EdgeHoverHandler
export NodeHoverHighlight, EdgeHoverHighlight
export NodeClickHandler, EdgeClickHandler
export NodeDragHandler, EdgeDragHandler
export NodeDrag, EdgeDrag

"""
    convert_selection(element, idx)

`mouse_selection` returns the basic plot type. In case of `Lines` check if it is
part of a `BezierSegments` and convert selection accordingly. In case of `LineSegments`
check if it is part of a `EdgePlot` and convert idx.
"""
convert_selection(element, idx) = (element, idx)
function convert_selection(element::LineSegments, idx)
    if element.parent isa EdgePlot
        return (element.parent, Int(idx / 2))
    end
    return (element, idx)
end
function convert_selection(element::Lines, idx)
    if element.parent isa BezierSegments && element.parent.parent isa EdgePlot
        bezier = element.parent
        idx = findfirst(isequal(element), bezier.plots)
        return (bezier.parent, idx)
    end
    return (element, idx)
end

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
    # in case of multiple graph plots in one axis this won't work
    gplots = filter(p -> p isa GraphPlot, parent.scene.plots)
    if length(gplots) !== 1
        @warn "There are multiple GraphPlots, register interaction to first!"
    end
    gplot = gplots[1]
    set_nodeplot!(inter, get_node_plot(gplot))
    set_edgeplot!(inter, get_edge_plot(gplot))
end

####
#### Hover Interaction
####

"""
    mutable struct HoverHandler{P<:ScenePlot, F} <: GraphInteraction

Object to handle hovers on `plot::P`. Calls `fun` on hover.
"""
mutable struct HoverHandler{P<:ScenePlot,F} <: GraphInteraction
    idx::Union{Nothing,Int}
    plot::Union{Nothing,P}
    fun::F
end
set_nodeplot!(h::HoverHandler{Scatter}, plot) = h.plot = plot
set_edgeplot!(h::HoverHandler{EdgePlot}, plot) = h.plot = plot

function process_interaction(handler::HoverHandler, event::MouseEvent, axis)
    if event.type === MouseEventTypes.over
        (element, idx) = convert_selection(mouse_selection(axis.scene)...)
        if element == handler.plot
            if handler.idx === nothing
                handler.idx = idx
                ret = handler.fun(true, handler.idx, event, axis)
                return ret isa Bool ? ret : false
            end
        else
            if handler.idx !== nothing
                ret = handler.fun(false, handler.idx, event, axis)
                handler.idx = nothing
                return ret isa Bool ? ret : false
            end
        end
    end
    return false
end

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
NodeHoverHandler(fun::F) where {F} = HoverHandler{Scatter,F}(nothing, nothing, fun)

"""
    NodeHoverHeighlight(p::GraphPlot, factor=2)

Magnifies the `node_size` of node under cursor by `factor`.

# Example
```
julia> g = wheel_graph(10)
julia> f, ax, p = graphplot(g, node_size = [20 for i in 1:nv(g)])
julia> register_interaction!(ax, :nodehover, NodeHoverHighlight(p))
```
"""
function NodeHoverHighlight(p::GraphPlot, factor=2)
    @assert p.node_size[] isa Vector{<:Real} "`node_size` object needs to be an Vector{<:Real} for this interaction to work!"
    action = (state, idx, _, _) -> begin
        old = p.node_size[][idx]
        p.node_size[][idx] = state ? old * factor : old / factor
        p.node_size[] = p.node_size[] #trigger observable
    end
    return NodeHoverHandler(action)
end

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
EdgeHoverHandler(fun::F) where {F} = HoverHandler{EdgePlot,F}(nothing, nothing, fun)

"""
    EdgeHoverHeighlight(p::GraphPlot, factor=2)

Magnifies the `edge_width` of edge under cursor by `factor`.
If `arrow_size isa Vector{<:Real}` it also magnefies the arrow scatter.

# Example
```
julia> g = wheel_digraph(10)
julia> f, ax, p = graphplot(g, edge_width = [3 for i in 1:ne(g)],
                               arrow_size=[10 for i in 1:ne(g)])
julia> register_interaction!(ax, :nodehover, EdgeHoverHighlight(p))
```
"""
function EdgeHoverHighlight(p::GraphPlot, factor=2)
    @assert p.edge_width[] isa Vector{<:Real} "`edge_width` object needs to be an Vector{<:Real} for this interaction to work!"
    scale_arrows = p.arrow_size[] isa Vector{<:Real}

    if length(p.edge_width[]) == ne(p[:graph][])
        action = (state, idx, _, _) -> begin
            old = p.edge_width[][idx]
            p.edge_width[][idx] = state ? old * factor : old / factor
            p.edge_width[] = p.edge_width[] #trigger observable
            if scale_arrows
                old = p.arrow_size[][idx]
                p.arrow_size[][idx] = state ? old * factor : old / factor
                p.arrow_size[] = p.arrow_size[] #trigger observable
            end
        end
    elseif length(p.edge_width[]) == 2 * ne(p[:graph][])
        action = (state, idx, _, _) -> begin
            oldA = p.edge_width[][2 * idx - 1]
            oldB = p.edge_width[][2 * idx]
            p.edge_width[][2 * idx - 1] = state ? oldA * factor : oldA / factor
            p.edge_width[][2 * idx] = state ? oldB * factor : oldB / factor
            p.edge_width[] = p.edge_width[] #trigger observable
            if scale_arrows
                old = p.arrow_size[][idx]
                p.arrow_size[][idx] = state ? old * factor : old / factor
                p.arrow_size[] = p.arrow_size[] #trigger observable
            end
        end
    else
        error("Can not make sense of `length(p.edge_width)`")
    end

    return EdgeHoverHandler(action)
end

####
#### Drag Interaction
####

"""
    mutable struct DragHandler{P<:ScenePlot, F} <: GraphInteraction

Object to handle left mous drags on `plot::P`.
"""
mutable struct DragHandler{P<:ScenePlot,F} <: GraphInteraction
    dragstate::Bool
    idx::Int
    plot::Union{Nothing,P}
    fun::F
end
set_nodeplot!(h::DragHandler{Scatter}, plot) = h.plot = plot
set_edgeplot!(h::DragHandler{EdgePlot}, plot) = h.plot = plot

function process_interaction(handler::DragHandler, event::MouseEvent, axis)
    if handler.dragstate == false # not in drag state
        if event.type === MouseEventTypes.leftdown
            # on leftdown save idx if happens to be on right element/plot
            (element, idx) = convert_selection(mouse_selection(axis.scene)...)
            handler.idx = element === handler.plot ? idx : 0
        elseif event.type === MouseEventTypes.leftdragstart && handler.idx != 0
            # if idx!=0 the last leftdown was on right element!
            handler.dragstate = true
        end
    elseif handler.dragstate == true # drag state
        if event.type === MouseEventTypes.leftdrag
            ret = handler.fun(true, handler.idx, event, axis)
            return ret isa Bool ? ret : false
        elseif event.type === MouseEventTypes.leftdragstop
            ret = handler.fun(false, handler.idx, event, axis)
            handler.idx = 0
            handler.dragstate = false
            return ret isa Bool ? ret : false
        end
    end
    return false
end

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
           p[:node_pos][][idx] = event.data
           p[:node_pos][] = p[:node_pos][]
       end
julia> register_interaction!(ax, :nodedrag, NodeDragHandler(action))
```
"""
NodeDragHandler(fun::F) where {F} = DragHandler{Scatter,F}(false, 0, nothing, fun)

"""
    NodeDrag(p::GraphPlot)

Allows drag and drop of Nodes. Please deregister the `:rectanglezoom` interaction.

# Example
```
julia> g = wheel_graph(10)
julia> f, ax, p = graphplot(g, node_size = [10 for i in 1:nv(g)])
julia> deregister_interaction!(ax, :rectanglezoom)
julia> register_interaction!(ax, :nodehover, NodeHoverHighlight(p))
julia> register_interaction!(ax, :nodedrag, NodeDrag(p))
```
"""
function NodeDrag(p)
    action = (state, idx, event, _) -> begin
        p[:node_pos][][idx] = event.data
        p[:node_pos][] = p[:node_pos][]
    end
    return NodeDragHandler(action)
end

"""
    EdgeDragHandler(fun)

Initializes `DragHandler` for Edges. Calls function

    fun(dragstate, idx, event, axis)

where `dragstate=true` during the drag and `false` at the end of the drag,
the last time `fun` is triggered. `idx` is the edge index.

See [`EdgeDrag`](@ref) for a concrete implementation.
```
"""
EdgeDragHandler(fun::F) where {F} = DragHandler{EdgePlot,F}(false, 0, nothing, fun)

"""
    EdgeDrag(p::GraphPlot)

Allows drag and drop of Edges. Please deregister the `:rectanglezoom` interaction.

# Example
```
julia> g = wheel_graph(10)
julia> f, ax, p = graphplot(g, edge_width = [3 for i in 1:ne(g)])
julia> deregister_interaction!(ax, :rectanglezoom)
julia> register_interaction!(ax, :edgehover, EdgeHoverHighlight(p))
julia> register_interaction!(ax, :edgedrag, EdgeDrag(p))
```
"""
function EdgeDrag(p)
    action = EdgeDragAction(p)
    return EdgeDragHandler(action)
end

mutable struct EdgeDragAction{PT<:GraphPlot}
    p::PT
    init::Union{Nothing,Point2f} # save click position
    src::Union{Nothing,Point2f}  # save src vertex position
    dst::Union{Nothing,Point2f}  # save dst vertex position
    EdgeDragAction(p::T) where {T} = new{T}(p, nothing, nothing, nothing)
end

function (action::EdgeDragAction)(state, idx, event, _)
    edge = collect(edges(action.p[:graph][]))[idx]
    if state == true
        if action.src === action.dst === action.init === nothing
            action.init = event.data
            action.src = action.p[:node_pos][][edge.src]
            action.dst = action.p[:node_pos][][edge.dst]
        end
        offset = event.data - action.init
        action.p[:node_pos][][edge.src] = action.src + offset
        action.p[:node_pos][][edge.dst] = action.dst + offset
        action.p[:node_pos][] = action.p[:node_pos][] # trigger change
    elseif state == false
        action.src = action.dst = action.init = nothing
    end
end

####
#### Click Interaction
####

"""
    mutable struct ClickHandler{P<:ScenePlot, F} <: GraphInteraction

Object to handle left mouse clicks on `plot::P`.
"""
mutable struct ClickHandler{P<:ScenePlot,F} <: GraphInteraction
    plot::Union{Nothing,P}
    fun::F
end
set_nodeplot!(h::ClickHandler{Scatter}, plot) = h.plot = plot
set_edgeplot!(h::ClickHandler{EdgePlot}, plot) = h.plot = plot

function process_interaction(handler::ClickHandler, event::MouseEvent, axis)
    if event.type === MouseEventTypes.leftclick
        (element, idx) = convert_selection(mouse_selection(axis.scene)...)
        if element == handler.plot
            ret = handler.fun(idx, event, axis)
            return ret isa Bool ? ret : false
        end
    end
    return false
end

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
NodeClickHandler(fun::F) where {F} = ClickHandler{Scatter,F}(nothing, fun)

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
EdgeClickHandler(fun::F) where {F} = ClickHandler{EdgePlot,F}(nothing, fun)
