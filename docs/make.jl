using hPFMD
using Documenter

DocMeta.setdocmeta!(hPFMD, :DocTestSetup, :(using hPFMD); recursive=true)

makedocs(;
    modules=[hPFMD],
    authors="zwu <w415146142@gmail.com> and contributors",
    repo="https://github.com/Chenghao Wu/hPFMD.jl/blob/{commit}{path}#{line}",
    sitename="hPFMD.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Chenghao Wu.github.io/hPFMD.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Chenghao Wu/hPFMD.jl",
    devbranch="main",
)
