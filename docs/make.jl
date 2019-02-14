using Documenter, PerronFrobenius
ENV["GKSwstype"] = "100"
push!(LOAD_PATH,"../src/")

PAGES = [
    "Overview" => "index.md",
    "Estimators" =>
        ["Overview" => "estimators/estimators.md",
        "Grid estimator" => "estimators/direct.md"
        ]
]

makedocs(
    modules = [PerronFrobenius],
    format = :markdown,
    sitename = "PerronFrobenius.jl",
    authors = "Kristian Agasøster Haaga",
    pages = PAGES
    # Use clean URLs, unless built as a "local" build
    #html_prettyurls = !("local" in ARGS),
    #html_canonical = "https://kahaaga.github.io/TimeseriesSurrogates.jl/latest/"
)

# deploydocs(
#     repo   = "github.com/kahaaga/PerronFrobenius.jl.git",
#     julia  = "0.6",
#     target = "build",
#     deps = nothing,
#     make = nothing,
#     osname = "linux"
# )
