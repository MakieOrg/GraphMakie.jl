using GraphMakie
using Documenter
using Literate
using CairoMakie

# preload the deps from the examples to supress precompilation output in docs
using JSServe
using NetworkDynamics
using LayeredLayouts
using Graphs
using RegistryInstances

DocMeta.setdocmeta!(GraphMakie, :DocTestSetup, :(using GraphMakie); recursive=true)

# generate examples
example_dir = joinpath(@__DIR__, "examples")
outdir = joinpath(@__DIR__, "src", "generated")
isdir(outdir) && rm(outdir, recursive=true)
mkpath(outdir)

for example in filter(contains(r".jl$"), readdir(example_dir, join=true))
    Literate.markdown(example, outdir;
                      preprocess=c->replace(c, "@save_reference " => ""))
end

makedocs(; modules=[GraphMakie], authors="Simon Danisch, Hans Würfel",
         repo="https://github.com/MakieOrg/GraphMakie.jl/blob/{commit}{path}#{line}",
         sitename="GraphMakie.jl",
         format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true",
                                canonical="https://graph.makie.org", assets=String[],
                                size_threshold_ignore=["generated/plots.md"]),
         pages=["Home" => "index.md",
                "Examples" => [
                    "Feature Walkthrough" => "generated/plots.md",
                    "Interaction Examples" => "generated/interactions.md",
                    "Dependency Graph" => "generated/depgraph.md",
                    "Stress on Truss" => "generated/truss.md",
                    "Julia AST" => "generated/syntaxtree.md",
                    "Decision Tree" => "generated/decisiontree.md",
                    "Reference Tests" => "generated/reftests.md",
                ],
                "🔗 Layouts (`NetworkLayout.jl`)" => "networklayout_forward.md",
                ],
         warnonly=[:missing_docs])

# if gh_pages branch gets to big, check out
# https://juliadocs.github.io/Documenter.jl/stable/man/hosting/#gh-pages-Branch
deploydocs(;repo="github.com/MakieOrg/GraphMakie.jl",
           push_preview=true)
