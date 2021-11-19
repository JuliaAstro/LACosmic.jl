using LACosmic
using Documenter

DocMeta.setdocmeta!(LACosmic, :DocTestSetup, :(using LACosmic); recursive=true)

makedocs(;
    modules=[LACosmic],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/JuliaAstro/LACosmic.jl/blob/{commit}{path}#{line}",
    sitename="LACosmic.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaastro.github.io/LACosmic.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API/Reference" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/JuliaAstro/LACosmic.jl",
    devbranch="main",
    push_preview=true
)
