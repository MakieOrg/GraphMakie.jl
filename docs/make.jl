using GraphMakie
using Documenter

DocMeta.setdocmeta!(GraphMakie, :DocTestSetup, :(using GraphMakie); recursive=true)

makedocs(;
    modules=[GraphMakie],
    authors="Simon Danisch",
    repo="https://github.com/SimonDanisch/GraphMakie.jl/blob/{commit}{path}#{line}",
    sitename="GraphMakie.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://SimonDanisch.github.io/GraphMakie.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SimonDanisch/GraphMakie.jl",
)
