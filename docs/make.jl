using PerManova
using Documenter

DocMeta.setdocmeta!(PerManova, :DocTestSetup, :(using PerManova); recursive=true)

makedocs(;
    modules=[PerManova],
    authors="Arthur Newbury",
    repo="https://github.com/EvoArt/PerManova.jl/blob/{commit}{path}#{line}",
    sitename="PerManova.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EvoArt.github.io/PerManova.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EvoArt/PerManova.jl",
)
