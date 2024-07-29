using DoViP_Lib
using Documenter

DocMeta.setdocmeta!(DoViP_Lib, :DocTestSetup, :(using DoViP_Lib); recursive=true)

makedocs(;
    modules=[DoViP_Lib],
    authors="Cristina Moraru",
    sitename="DoViP_Lib.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
