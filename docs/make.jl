using ECCO
using Documenter

DocMeta.setdocmeta!(ECCO, :DocTestSetup, :(using ECCO); recursive=true)

makedocs(;
    modules=[ECCO,ECCO.toy_problems,ECCO.Zygote_examples,
		ECCO.glacier_model,ECCO.Lorenz_models],
    authors="gaelforget <gforget@mit.edu> and contributors",
    sitename="ECCO.jl",
    format=Documenter.HTML(;
        canonical="https://gaelforget.github.io/ECCO.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API" => "API.md",
        "context" => "context.md",
    ],
    warnonly = [:cross_references,:missing_docs],
)

deploydocs(;
    repo="github.com/gaelforget/ECCO.jl",
    devbranch="main",
)
