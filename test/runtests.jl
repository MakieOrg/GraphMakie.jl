using GraphMakie
using GraphMakie.LightGraphs
using GraphMakie.AbstractPlotting
using Test

@testset "GraphMakie.jl" begin
    g = wheel_digraph(10)
    f, ax, p = plot(g)
    f = plot_graph(g)
end
