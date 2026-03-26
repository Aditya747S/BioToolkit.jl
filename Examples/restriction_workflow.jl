using BioToolkit

sequence = "TTTGAATTCGGATCCAAAGGCCATGACGTC"

println("Restriction workflow")
println("  enzymes in catalog: ", length(BioToolkit.restriction_enzyme_names()))
println("  EcoRI: ", BioToolkit.restriction_enzyme("EcoRI"))
println("  BamHI: ", BioToolkit.restriction_enzyme("BamHI"))

hits = BioToolkit.find_restriction_sites(sequence, "EcoRI")
println("  EcoRI sites: ", hits)
println("  digest fragments: ", BioToolkit.digest_sequence(sequence, ["EcoRI", "BamHI"]))
println("  digest map keys: ", collect(keys(BioToolkit.restriction_digest_map(sequence, ["EcoRI", "BamHI", "SmaI"]))))
