# GraphMakie

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](http://juliaplots.org/GraphMakie.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://juliaplots.org/GraphMakie.jl/dev/)
[![Build Status](https://github.com/JuliaPlots/GraphMakie.jl/workflows/CI/badge.svg)](https://github.com/JuliaPlots/GraphMakie.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaPlots/GraphMakie.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaPlots/GraphMakie.jl)

Plotting graphs¹... with [Makie](https://github.com/JuliaPlots/Makie.jl)!

This package is just a plotting recipe, you still need to use one of the Makie backends.

`GraphMakie` is in an early development stage and might break often. Any
contribution such as suggesting features, raising issues and opening PRs are
welcome!

Starting from v0.3 `GraphMakie.jl` switches from `LightGraphs.jl` to `Graphs.jl` for the graph representation. See this [discourse post](https://discourse.julialang.org/t/lightgraphs-jl-transition/69526/17) for more information. If you want to use `LightGraphs.jl` please specifically `] add GraphMakie@0.2`!

## Installation
``` julia
pkg> add GraphMakie
```

## Basic usage
```julia
using GLMakie
using GraphMakie
using Graphs
g = complete_graph(10)
graphplot(g)
```


----------------------------
¹the networky type with nodes and edges
