# Install DESeq2
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("pasilla")

library(DESeq2)

# --- 1. Generate Sparse Synthetic Data ---
set.seed(42)
n_genes <- 1000
n_samples <- 10

# High variance/sparsity (size=0.1 creates many zeros)
counts <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 0.1), nrow = n_genes)
rownames(counts) <- paste0("gene_", 1:n_genes)
colnames(counts) <- paste0("sample_", 1:n_samples)

# Sample Metadata
coldata <- data.frame(
  row.names = colnames(counts),
  condition = factor(rep(c("A", "B"), each = 5))
)

# Add Differential Expression
counts[1:50, coldata$condition == "B"] <- counts[1:50, coldata$condition == "B"] * 4

# Save for Julia
write.csv(counts, "counts_matrix(2).csv")
write.csv(coldata, "col_data.csv")
cat("Saved data for Julia.\n")

# --- 2. Run DESeq2 with 'poscounts' ---
cat("Running R DESeq2 (poscounts mode)...\n")
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

# CRITICAL FIX: Use 'poscounts' to handle the zeros
dds <- estimateSizeFactors(dds, type = "poscounts") 
dds <- DESeq(dds) # Continue with the rest of the pipeline

# Get Results
res <- results(dds, contrast = c("condition", "B", "A"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Save R Results
write.csv(res_df, "r_results.csv", row.names = FALSE)
cat("Saved R results.\n")
print(summary(res))
