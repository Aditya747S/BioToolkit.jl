using Pkg

Pkg.activate(joinpath(@__DIR__, "docs"))
Pkg.instantiate()
include(joinpath(@__DIR__, "docs", "make.jl"))
