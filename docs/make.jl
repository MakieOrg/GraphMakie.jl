using GraphMakie
using Documenter
using Literate
using CairoMakie

DocMeta.setdocmeta!(GraphMakie, :DocTestSetup, :(using GraphMakie); recursive=true)

# generate examples
examples = [
    joinpath(@__DIR__, "examples", "plots.jl"),
]
OUTPUT = joinpath(@__DIR__, "src", "generated")
isdir(OUTPUT) && rm(OUTPUT, recursive=true)
mkpath(OUTPUT)

for ex in examples
    Literate.markdown(ex, OUTPUT)
end

makedocs(; modules=[GraphMakie], authors="Simon Danisch",
         repo="https://github.com/JuliaPlots/GraphMakie.jl/blob/{commit}{path}#{line}",
         sitename="GraphMakie.jl",
         format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true",
                                canonical="https://JuliaPlots.github.io/GraphMakie.jl", assets=String[]),
         pages=["Home" => "index.md",
                "Plot Examples" => "generated/plots.md"])

deploydocs(;repo="github.com/JuliaPlots/GraphMakie.jl",
           push_preview=true)
