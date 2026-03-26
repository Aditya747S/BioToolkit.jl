using BioToolkit
using Distributions, Random

println("Testing Coalescent Simulation...")
n = 5
ne = 1000
mu = 0.001
seq_len = 50

tree, haps = BioToolkit.simulate_coalescent(n; ne=ne, mu=mu, seq_len=seq_len)

println("Tree: ", tree)
println("Haplotypes: ", length(haps))
println("Success!")
