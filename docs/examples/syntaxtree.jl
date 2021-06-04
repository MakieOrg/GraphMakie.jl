#=
# AST of a Julia function definition

In this example we are going to plot an abstract syntax tree of a Julia function using the
Bucheim Layout from [`NetworkLayout.jl`](https://github.com/JuliaGraphs/NetworkLayout.jl).
=#
using CairoMakie
CairoMakie.activate!(type="png") #hide
CairoMakie.inline!(true) #hide
set_theme!(resolution=(800, 600)) #hide
using LightGraphs
using GraphMakie
using NetworkLayout
using CairoMakie

#=
The following code, which walks the AST and creats a `SimpleDiGraph` was taken and slightly
modified from [`TreeView.jl`](https://github.com/JuliaTeX/TreeView.jl). Thanks!
=#
function walk_tree(ex; show_call=true)
    g = SimpleDiGraph()
    labels = Any[]
    walk_tree!(g, labels, ex, show_call)
    return (g, labels)
end

function walk_tree!(g, labels, ex, show_call)
    add_vertex!(g)
    top_vertex = vertices(g)[end]

    where_start = 1  # which argument to start with

    if !(show_call) && ex.head == :call
        f = ex.args[1]   # the function name
        push!(labels, f)
        where_start = 2   # drop "call" from tree
    else
        push!(labels, ex.head)
    end

    for i in where_start:length(ex.args)
        if isa(ex.args[i], Expr)
            child = walk_tree!(g, labels, ex.args[i], show_call)
            add_edge!(g, top_vertex, child)
        elseif !isa(ex.args[i], LineNumberNode)
            add_vertex!(g)
            n = vertices(g)[end]
            add_edge!(g, top_vertex, n)
            push!(labels, ex.args[i])
        end
    end

    return top_vertex
end

#=
The expression we want to look at is the recursive definition of the Fibonacci sequence.
=#
expr = quote
    function fib(n)
        if n > 1
            return fib(n-1) + fib(n-2)
        else
            return n
        end
    end
end

g, labels = walk_tree(expr, show_call=true)
nlabels_align = [(:left, :bottom) for v in vertices(g)]
buch = Buchheim()
fig, ax, p = graphplot(g;  layout=buch,
                       nlabels=repr.(labels),
                       nlabels_distance=5,
                       nlabels_align)
hidedecorations!(ax); hidespines!(ax)
fig #hide

# This does not look nice yet! Lets tweak the `align` parameter of the nodes labels...
for v in vertices(g)
    if isempty(inneighbors(g, v)) # root
        nlabels_align[v] = (:center,:bottom)
    elseif isempty(outneighbors(g, v)) #leaf
        nlabels_align[v] = (:center,:top)
    else
        self = p[:node_positions][][v]
        parent = p[:node_positions][][inneighbors(g, v)[1]]
        if self[1] < parent[1] # left branch
            nlabels_align[v] = (:right,:bottom)
        end
    end
end
p.nlabels_align = nlabels_align
fig # hide
