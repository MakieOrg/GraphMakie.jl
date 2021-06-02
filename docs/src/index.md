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

## Interactions
`GraphMakie.jl` provides some helper functions to register interactions to your graph plot.
There are special Interaction types for hovering, clicking and draging nodes and edges.
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
