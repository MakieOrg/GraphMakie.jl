```@meta
CurrentModule = GraphMakie
```

# GraphMakie
This is the Documentation for [GraphMakie](https://github.com/MakieOrg/GraphMakie.jl).

This Package consists of two parts: a plot recipe for graphs types from
[Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) and some helper
functions to add interactions to those plots.

!!! note "about Graphs.jl and LightGraphs.jl"
    Starting from v0.3 `GraphMakie.jl` switches from `LightGraphs.jl` to `Graphs.jl` as the for the underlying graph package. See this [discourse post](https://discourse.julialang.org/t/lightgraphs-jl-transition/69526/17) for more information. If you want to use `LightGraphs.jl` please specifically `] add GraphMakie@0.2`!

There are also [plot examples](generated/plots.md) and [interaction examples](generated/interactions.md) pages.

## The `graphplot` Recipe
```@docs
graphplot
```

## Network Layouts
The layout algorithms are provided by [`NetworkLayout.jl`](https://github.com/JuliaGraphs/NetworkLayout.jl). See
the [docs](https://juliagraphs.org/NetworkLayout.jl/stable/) for a list of available layouts.

A layout has to be a function `f(g::AbstractGraph) -> pos::Vector{Point}`. You can also provide your
own layouts or use other packages like [`LayeredLayouts.jl`](https://github.com/oxinabox/LayeredLayouts.jl) for
DAG (see also the [Dependency Graph of a Package](@ref) example).
```julia
using LayeredLayouts
function mylayout(g::SimpleGraph)
   xs, ys, _ = solve_positions(Zarate(), g)
   return Point.(zip(xs, ys))
end
```

## Predefined Interactions
`GraphMakie.jl` provides some pre-built interactions to enable drag&drop of nodes and edges as well as highlight on hover.

To try them all use the following code in a `GLMakie` environment.
```julia
using GLMakie
using GraphMakie
using Graphs
g = wheel_graph(10)
f, ax, p = graphplot(g, edge_width=[3 for i in 1:ne(g)],
                     node_size=[10 for i in 1:nv(g)])

deregister_interaction!(ax, :rectanglezoom)
register_interaction!(ax, :nhover, NodeHoverHighlight(p))
register_interaction!(ax, :ehover, EdgeHoverHighlight(p))
register_interaction!(ax, :ndrag, NodeDrag(p))
register_interaction!(ax, :edrag, EdgeDrag(p))
```

```@docs
NodeHoverHighlight
EdgeHoverHighlight
NodeDrag
EdgeDrag
```

## Interaction Interface
`GraphMakie.jl` provides some helper functions to register interactions to your graph plot.
There are special interaction types for hovering, clicking and draging nodes and edges.
For more information on the axis interaction please consult the [`Makie.jl` docs](https://docs.makie.org/stable/examples/blocks/axis/index.html#custom_interactions).

The general idea is to create some handler type, provide some action function and register it
as an interaction with the axes.

### Click Interactions
```@docs
NodeClickHandler
EdgeClickHandler
```
### Hover Interactions
```@docs
NodeHoverHandler
EdgeHoverHandler
```

### Drag Interactions
```@docs
NodeDragHandler
EdgeDragHandler
```
