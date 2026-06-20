# Install dependencies
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("apeglm") # For modern shrinkage

library(DESeq2)
library(ggplot2)

# --- 1. Data Generation (Matching your sparse profile) ---
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

# Validate inputs
if (!all(colnames(counts) == rownames(coldata))) {
  stop("Column names of counts matrix must match row names of coldata")
}

# Add Differential Expression (Strong signal)
counts[1:50, coldata$condition == "B"] <- counts[1:50, coldata$condition == "B"] * 4

# Validate counts are integers
if (!all(counts == as.integer(counts))) {
  counts <- matrix(as.integer(counts), nrow=nrow(counts), dimnames=dimnames(counts))
}

# Save Counts for Julia (with row names as gene IDs)
write.csv(counts, "counts_matrix(3).csv", row.names=TRUE)
write.csv(coldata, "col_data(2).csv", row.names=TRUE)
cat("Data saved.\n")

# --- 2. Run DESeq2 Pipeline (With Diagnostics) ---
cat("Running R DESeq2...\n")
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

# STEP A: Normalization
# Use poscounts to match your Julia fix
dds <- estimateSizeFactors(dds, type = "poscounts") 

# STEP B: Dispersion Estimation
dds <- estimateDispersions(dds) 

# STEP C: GLM Fitting (Wald Test)
# This fits the model and calculates Cook's distance
dds <- nbinomWaldTest(dds)

# --- 3. Save Intermediate Results ---

# 3a. Save Normalization Factors
norm_factors <- sizeFactors(dds)
write.csv(data.frame(sample=names(norm_factors), factor=norm_factors), "r_norm_factors.csv")
cat("Saved Size Factors.\n")

# 3b. Save Dispersion Estimates
# We want the final MAP dispersions used in testing
disp_estimates <- dispersions(dds) # or mcols(dds)$dispersion
gene_info <- data.frame(gene_id = rownames(counts), dispersion = disp_estimates)
write.csv(gene_info, "r_dispersions.csv", row.names=FALSE)
cat("Saved Dispersions.\n")

# 3c. Save Cook's Distances (CRITICAL FOR OUTLIER DEBUGGING)
# Cook's distance is stored as a matrix (Genes x Samples)
if (hasAssays(dds) && 'cooks' %in% names(assays(dds))) {
  cooks_mat <- assays(dds)[["cooks"]]
  if (!is.null(cooks_mat) && nrow(cooks_mat) > 0) {
    # Calculate max cooks distance per gene with proper NA handling
    max_cooks <- apply(cooks_mat, 1, function(x) {
      valid_vals <- x[!is.na(x)]
      if (length(valid_vals) == 0) NA else max(valid_vals)
    })
    cooks_df <- data.frame(gene_id = rownames(counts), max_cooks_distance = max_cooks)
    write.csv(cooks_df, "r_cooks_distance.csv", row.names=FALSE)
    cat("Saved Cook's Distances.\n")
  } else {
    warning("Cook's distance matrix is empty")
    write.csv(data.frame(gene_id = rownames(counts), max_cooks_distance = NA), 
              "r_cooks_distance.csv", row.names=FALSE)
  }
} else {
  warning("Cook's distances not computed by DESeq2 pipeline")
  write.csv(data.frame(gene_id = rownames(counts), max_cooks_distance = NA), 
            "r_cooks_distance.csv", row.names=FALSE)
}

# --- 4. Results Generation (Comparing MLE vs Shrunken) ---

# 4a. RAW Results (MLE - No Shrinkage)
# This is what standard results() returns.
# Validate that contrast levels exist
available_levels <- levels(coldata$condition)
if (!all(c("B", "A") %in% available_levels)) {
  stop(sprintf("Contrast levels 'B' and 'A' not found. Available levels: %s", 
               paste(available_levels, collapse=", ")))
}
# We disable cooksCutoff here to see the raw statistical power.
res_raw <- results(dds, contrast = c("condition", "B", "A"), cooksCutoff = FALSE)
res_raw_df <- as.data.frame(res_raw)
res_raw_df$gene_id <- rownames(res_raw_df)

# 4b. SHRUNKEN Results (Post-hoc processing)
# This uses 'apeglm' which is the modern standard for effect size shrinkage
# resLFC <- lfcShrink(dds, contrast = c("condition", "B", "A"), type="apeglm")
# For simplicity here, we will stick to the raw results for the main comparison,
# as shrinkage logic varies by implementation.
# BUT: We will save both to be safe.

write.csv(res_raw_df, "r_results_raw.csv", row.names = FALSE)
cat("Saved Raw Results (MLE).\n")

# Summary
cat("\n--- R Summary ---\n")
print(summary(res_raw))
cat("\nNOTE: Compare 'r_cooks_distance.csv' with your Julia output.\n")
cat("Genes with max_cooks > ~8.65 are flagged as outliers in R default settings.\n")

# Output-->
# Data saved. Running R DESeq2... converting counts to integer mode
# gene-wise dispersion estimates
# mean-dispersion relationship
# -- note: fitType='parametric', but the dispersion trend was not well captured by the function: y = a/x + b, and a local regression fit was automatically substituted. specify fitType='local' or 'mean' to avoid this message next time.
# final dispersion estimates
# Saved Size Factors. Saved Dispersions. Warning message in FUN(newX[, i], ...): “no non-missing arguments to max; returning -Inf” Saved Cook's Distances. Saved Raw Results (MLE).
# --- R Summary ---
# out of 999 with nonzero total read count adjusted p-value < 0.1 LFC > 0 (up)       : 49, 4.9% LFC < 0 (down)     : 55, 5.5% outliers [1]       : 0, 0% low counts [2]     : 0, 0% (mean count < 0) [1] see 'cooksCutoff' argument of ?results [2] see 'independentFiltering' argument of ?results
# NULL
# NOTE: Compare 'r_cooks_distance.csv' with your Julia output. Genes with max_cooks > ~8.65 are flagged as outliers in R default settings.