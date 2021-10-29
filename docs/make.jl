using PerMANOVA
using Documenter

DocMeta.setdocmeta!(PerMANOVA, :DocTestSetup, :(using PerMANOVA); recursive=true)

makedocs(;
    modules=[PerMANOVA],
    authors="Arthur Newbury",
    repo="https://github.com/EvoArt/PerMANOVA.jl/blob/{commit}{path}#{line}",
    sitename="PerMANOVA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EvoArt.github.io/PerMANOVA.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EvoArt/PerMANOVA.jl",
)
