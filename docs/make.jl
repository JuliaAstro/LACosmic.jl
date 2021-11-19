using LACosmic
using Documenter

DocMeta.setdocmeta!(LACosmic, :DocTestSetup, :(using LACosmic); recursive=true)

makedocs(;
    modules=[LACosmic],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/mileslucas/LACosmic.jl/blob/{commit}{path}#{line}",
    sitename="LACosmic.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mileslucas.com/LACosmic.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mileslucas/LACosmic.jl",
    devbranch="main",
    push_preview=true
)
