using Suslin
using Documenter

DocMeta.setdocmeta!(Suslin, :DocTestSetup, :(using Suslin); recursive=true)

makedocs(;
    modules=[Suslin],
    authors="nzy1997",
    sitename="Suslin.jl",
    format=Documenter.HTML(;
        canonical="https://nzy1997.github.io/Suslin.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/nzy1997/Suslin.jl",
    devbranch="main",
)
