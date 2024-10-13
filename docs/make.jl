using ECCO
using Documenter

DocMeta.setdocmeta!(ECCO, :DocTestSetup, :(using ECCO); recursive=true)

makedocs(;
    modules=[ECCO],
    authors="gaelforget <gforget@mit.edu> and contributors",
    sitename="ECCO.jl",
    format=Documenter.HTML(;
        canonical="https://gaelforget.github.io/ECCO.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/gaelforget/ECCO.jl",
    devbranch="main",
)
