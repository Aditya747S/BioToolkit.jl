using Test
using Plots
using BioToolkit

@testset "Root entrypoint" begin
	@test isdefined(BioToolkit, :gc_content)
	@test isdefined(BioToolkit, :window_coverage)
	@test isdefined(BioToolkit, :write_fastq)
	@test isdefined(BioToolkit, :SeqRecordLite)
	@test isdefined(BioToolkit, :pairwise_align)
	@test isdefined(BioToolkit, :motif_counts)
	@test isdefined(BioToolkit, :parse_gff_record)
	@test isdefined(BioToolkit, :annotate_genbank_record)
	@test isdefined(BioToolkit, :read_genbank)
	@test isdefined(BioToolkit, :Restriction)
	@test isdefined(BioToolkit, :Entrez)
	@test isdefined(BioToolkit, :Medline)
	@test isdefined(BioToolkit, :Compass)
	@test isdefined(BioToolkit, :Proteomics)
	@test isdefined(BioToolkit, :Metabolomics)
	@test isdefined(BioToolkit, :SystemsBio)
	@test isdefined(BioToolkit, :GWAS)
	@test isdefined(BioToolkit, :gwas_linear_scan)
	@test isdefined(BioToolkit, :manhattan_plot)
	@test isdefined(BioToolkit, :KEGG)
	@test isdefined(BioToolkit, :Pathway)
	@test isdefined(BioToolkit, :SCOP)
	@test isdefined(BioToolkit, :CATH)
end
