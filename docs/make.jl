using PERMANOVA
using Documenter

DocMeta.setdocmeta!(PERMANOVA, :DocTestSetup, :(using PERMANOVA); recursive=true)

makedocs(;
    modules=[PERMANOVA],
    authors="Arthur Newbury",
    repo="https://github.com/EvoArt/PERMANOVA.jl/blob/{commit}{path}#{line}",
    sitename="PERMANOVA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EvoArt.github.io/PERMANOVA.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EvoArt/PERMANOVA.jl",
)
