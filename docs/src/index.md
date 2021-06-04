```@meta
CurrentModule = GraphMakie
```

# GraphMakie
This is the Documentation for [GraphMakie](https://github.com/JuliaPlots/GraphMakie.jl).

This Package consists of two parts: a plot recipe for graphs types from
[LightGraphs.jl](https://juliagraphs.org/LightGraphs.jl/latest/) and some helper
functions to add interactions to those plots.

There are also [plot examples](generated/plots.md) and [interaction examples](generated/interactions.md) pages.

## The `graphplot` Recipe
```@docs
graphplot
```

## Predefined Interactions
`GraphMakie.jl` provides some pre-built interactions to enable drag&drop of nodes and edges as well as highlight on hover.

To try them all use the following code in a `GLMakie` environment.
```julia
using GLMakie
using GraphMakie
using LightGraphs
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
For more information on the axis interaction please consult the [`Makie.jl` docs](https://makie.juliaplots.org/dev/makielayout/axis.html#Custom-Interactions).

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
