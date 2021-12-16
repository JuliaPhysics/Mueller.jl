using Mueller
using Documenter

DocMeta.setdocmeta!(Mueller, :DocTestSetup, :(using Mueller); recursive=true)

makedocs(;
    modules=[Mueller],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/JuliaPhysics/Mueller.jl/blob/{commit}{path}#{line}",
    sitename="Mueller.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaphysics.github.io/Mueller.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API/Reference" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/JuliaPhysics/Mueller.jl",
    devbranch="main",
)
