using CairoMakie #md
CairoMakie.activate!(type="png") #md
CairoMakie.inline!(true) #md
using GLMakie #src
GLMakie.activate!() #src

set_theme!(resolution=(800, 600)) #hide
using GraphMakie
using LightGraphs
using GraphMakie.Makie.Colors
using LayeredLayouts
using PkgDeps

function depgraph(root)
    packages = [root]
    connections = Vector{Pair{Int,Int}}()

    for pkg in packages
        println("Check ", pkg)
        pkgidx = findfirst(isequal(pkg), packages)
        deps = direct_dependencies(pkg)

        for dep in keys(deps)
            idx = findfirst(isequal(dep), packages)
            if idx === nothing
                push!(packages, dep)
                idx = lastindex(packages)
            end
            push!(connections, idx => pkgidx)
        end
    end
    g = SimpleDiGraph(length(packages))
    t
    for c in connections
        add_edge!(g, c)
    end
    return (packages, g)
end

(packages, g) = depgraph("Revise")
N = length(packages)
xs, ys, paths = solve_positions(Zarate(), g)

# define layout as function adj_matrix -> Vector{Point}
# support for paths not yet implemented
lay = _ -> Point.(zip(xs,ys))

# manually tweak some of the lable aligns
align = [(:right, :bottom) for i in 1:N]
align[1] = (:left, :bottom)
align[3] = align[13] = (:left, :top)
align[6] = (:center, :bottom)
align[10] = (:right, :top)

# shift "CodeTracking" node in data space
offset = [Point2f0(0,0) for i in 1:N]
offset[6] = Point(0.0, 0.4)

f, ax, p = graphplot(g; layout=lay, arrow_size=15, edge_color=:gray,
                     nlabels=packages,
                     nlabels_align=align,
                     nlabels_offset=offset,
                     node_size=[9.0 for i in 1:N],
                     edge_width=[3 for i in 1:ne(g)])
ax.title = "Dep Graph of Revise.jl"
xlims!(ax, -0.5, 5.5)
hidedecorations!(ax); hidespines!(ax)
f #md

deregister_interaction!(ax, :rectanglezoom) #src
register_interaction!(ax, :nodehover, NodeHoverHighlight(p)) #src
register_interaction!(ax, :edgehover, EdgeHoverHighlight(p)) #src
register_interaction!(ax, :edrag, EdgeDrag(p)) #src
register_interaction!(ax, :ndrag, NodeDrag(p)) #src
