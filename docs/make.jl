using GraphMakie
using Documenter
using Literate
using CairoMakie

# preload the deps from the examples to supress precompilation output in docs
using JSServe
using NetworkDynamics
using LayeredLayouts
using Graphs
using PkgDeps

DocMeta.setdocmeta!(GraphMakie, :DocTestSetup, :(using GraphMakie); recursive=true)

# generate examples
examples = [
    joinpath(@__DIR__, "examples", "plots.jl"),
    joinpath(@__DIR__, "examples", "interactions.jl"),
    joinpath(@__DIR__, "examples", "depgraph.jl"),
    joinpath(@__DIR__, "examples", "truss.jl"),
    joinpath(@__DIR__, "examples", "syntaxtree.jl"),
]
OUTPUT = joinpath(@__DIR__, "src", "generated")
isdir(OUTPUT) && rm(OUTPUT, recursive=true)
mkpath(OUTPUT)

for ex in examples
    Literate.markdown(ex, OUTPUT)
end

makedocs(; modules=[GraphMakie], authors="Simon Danisch, Hans Würfel",
         repo="https://github.com/JuliaPlots/GraphMakie.jl/blob/{commit}{path}#{line}",
         sitename="GraphMakie.jl",
         format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true",
                                canonical="https://JuliaPlots.github.io/GraphMakie.jl", assets=String[]),
         pages=["Home" => "index.md",
                "Examples" => [
                    "Feature Walkthrough" => "generated/plots.md",
                    "Interaction Examples" => "generated/interactions.md",
                    "Dependency Graph" => "generated/depgraph.md",
                    "Stress on Truss" => "generated/truss.md",
                    "Julia AST" => "generated/syntaxtree.md",
                ]
                ])

# if gh_pages branch gets to big, check out
# https://juliadocs.github.io/Documenter.jl/stable/man/hosting/#gh-pages-Branch
deploydocs(;repo="github.com/JuliaPlots/GraphMakie.jl",
           push_preview=true)
