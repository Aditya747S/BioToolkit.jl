using Documenter
using BioToolkit

makedocs(
    sitename = "BioToolkit.jl",
    # format = Documenter.HTML(),
    format = Documenter.HTML(
        edit_link = nothing,
        assets = ["assets/custom.css"],
    ),
    modules = [BioToolkit],
    authors = "Aditya Sharma",
    
    pages = [
        "Home" => "index.md",
        # "Guide" => "guide.md",
        "API Reference" => "api.md",  
        # "Alignment" => "alignment.md",
        "Benchmarks" => "benchmark_report.md"
    ]
)

deploydocs(
    repo = "github.com/your-username/BioToolkit.jl.git",
    devbranch = "main"
)