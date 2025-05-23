using Mueller
using Documenter
using Documenter.Remotes: GitHub

DocMeta.setdocmeta!(Mueller, :DocTestSetup, :(using Mueller); recursive=true)

makedocs(;
    modules=[Mueller],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo=GitHub("JuliaPhysics/Mueller.jl"),
    sitename="Mueller.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaphysics.github.io/Mueller.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API/Reference" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/JuliaPhysics/Mueller.jl",
    devbranch="main",
    push_preview=true,
)
