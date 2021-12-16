using Mueller
using Documenter

DocMeta.setdocmeta!(Mueller, :DocTestSetup, :(using Mueller); recursive=true)

makedocs(;
    modules=[Mueller],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/mileslucas/Mueller.jl/blob/{commit}{path}#{line}",
    sitename="Mueller.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mileslucas.github.io/Mueller.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mileslucas/Mueller.jl",
    devbranch="main",
)
