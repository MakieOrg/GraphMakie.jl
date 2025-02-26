# GraphMakie

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](http://graph.makie.org/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://graph.makie.org/dev/)
[![Build Status](https://github.com/MakieOrg/GraphMakie.jl/workflows/CI/badge.svg)](https://github.com/MakieOrg/GraphMakie.jl/actions)
[![Coverage](https://codecov.io/gh/MakieOrg/GraphMakie.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/MakieOrg/GraphMakie.jl)

Plotting graphs¹... with [Makie](https://github.com/MakieOrg/Makie.jl)!

This package is just a plotting recipe; you still need to use one of the Makie backends. `GraphMakie.jl` uses graph datatypes from [`Graphs.jl`](https://github.com/JuliaGraphs/Graphs.jl) and layout algorithms from [`NetworkLayout.jl`](https://github.com/JuliaGraphs/NetworkLayout.jl). Check out [our docs](http://graph.makie.org/stable) for detailed information on how to use this package.

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
