using FrequencyDomainDecomposition
using Documenter

makedocs(;
    modules=[FrequencyDomainDecomposition],
    authors="hzgzh",
    repo="https://github.com/hzgzh/FrequencyDomainDecomposition.jl/blob/{commit}{path}#L{line}",
    sitename="FrequencyDomainDecomposition.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://hzgzh.github.io/FrequencyDomainDecomposition.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/hzgzh/FrequencyDomainDecomposition.jl",
)
