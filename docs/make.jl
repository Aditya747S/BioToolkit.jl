using Documenter
using BioToolkit

makedocs(
    sitename = "BioToolkit.jl",
    # format = Documenter.HTML(),
    format = Documenter.HTML(
        edit_link = nothing,
        assets = ["assets/custom.css"],
        size_threshold = 2048_000,
        size_threshold_warn = 1000_000,
    ),
    modules = [BioToolkit],
    authors = "Aditya Sharma",
    checkdocs = :none,
    
    pages = [
        "Home" => "index.md",
        # "Guide" => "guide.md",
        "API Reference" => "api.md",  
        # "Alignment" => "alignment.md",
        "Benchmarks" => "benchmark_report.md"
    ]
)

if get(ENV, "CI", "") == "true" || get(ENV, "GITHUB_ACTIONS", "") == "true"
    deploydocs(
        repo = "github.com/your-username/BioToolkit.jl.git",
        devbranch = "main"
    )
end