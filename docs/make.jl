using LACosmic
using Documenter
using Documenter.Remotes: GitHub

DocMeta.setdocmeta!(LACosmic, :DocTestSetup, :(using LACosmic); recursive=true)

include("pages.jl")

makedocs(;
    modules=[LACosmic],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo=GitHub("JuliaAstro/LACosmic.jl"),
    sitename="LACosmic.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaastro.github.io/LACosmic.jl",
        assets=String[],
    ),
    pages=pages,
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/JuliaAstro/LACosmic.jl",
    devbranch="main",
    push_preview=true
)
