using BioToolkit

function print_preview(label, text; max_lines=8)
    lines = split(chomp(text), '\n')
    println(label)
    for line in lines[1:min(end, max_lines)]
        println(line)
    end
    if length(lines) > max_lines
        println("  ...")
    end
end

println("Entrez workflow")

try
    pubmed = BioToolkit.entrez_pubmed_search("BRCA1[Title] AND Homo sapiens[Organism]"; retmax=3)
    println("  PubMed hits: ", pubmed.count)
    println("  PubMed ids: ", pubmed.ids)

    if !isempty(pubmed.ids)
        summary = BioToolkit.entrez_summary("pubmed", pubmed.ids[1:min(2, length(pubmed.ids))])
        println("  PubMed summary keys: ", collect(keys(summary)))
    end

    taxonomy = BioToolkit.entrez_taxonomy_search("Homo sapiens"; retmax=3)
    println("  Taxonomy ids: ", taxonomy.ids)

    nuccore = BioToolkit.entrez_nucleotide_search("BRCA1 Homo sapiens[Organism]"; retmax=1)
    println("  Nucleotide ids: ", nuccore.ids)

    if !isempty(nuccore.ids)
        fasta = BioToolkit.entrez_fetch_fasta("nuccore", nuccore.ids[1:1])
        print_preview("  FASTA preview:", fasta)

        genbank = BioToolkit.entrez_fetch_genbank("nuccore", nuccore.ids[1:1])
        print_preview("  GenBank preview:", genbank)
    end
catch err
    println("  Live Entrez call failed: ", err)
    println("  This usually means the network is unavailable or NCBI rate-limited the request.")
end
