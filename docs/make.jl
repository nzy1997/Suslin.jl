using SuslinStability
using Documenter

DocMeta.setdocmeta!(SuslinStability, :DocTestSetup, :(using SuslinStability); recursive=true)

makedocs(;
    modules=[SuslinStability],
    authors="nzy1997",
    sitename="SuslinStability.jl",
    format=Documenter.HTML(;
        canonical="https://nzy1997.github.io/SuslinStability.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/nzy1997/SuslinStability.jl",
    devbranch="main",
)
