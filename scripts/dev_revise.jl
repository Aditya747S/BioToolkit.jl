using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Revise
using BioToolkit

println("Revise is active. Edit source files and keep this Julia session running.")
