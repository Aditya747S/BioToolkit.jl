# Install packages
if (!require("Matrix")) install.packages("Matrix")
if (!require("ggplot2")) install.packages("ggplot2")
library(Matrix)
library(ggplot2)

set.seed(42)

# 1. Generate Data
n_genes <- 1000
n_cells <- 500
latent_dim <- 4

# (Your generation logic - simplified for reproducibility)
gene_loadings <- matrix(0, nrow=n_genes, ncol=latent_dim)
gene_groups <- sample(1:latent_dim, n_genes, replace=TRUE)
for (i in 1:n_genes) {
  gene_loadings[i, gene_groups[i]] <- 2.5 + runif(1)
  gene_loadings[i, -gene_groups[i]] <- 0.1 * rnorm(latent_dim-1)
}
cell_scores <- matrix(rnorm(latent_dim * n_cells), nrow=latent_dim)
signal <- gene_loadings %*% cell_scores / sqrt(latent_dim)
baseline <- 35.0 + 4.0 * rnorm(n_genes)
noisy_counts <- round(baseline + 14.0 * signal + 0.35 * matrix(rnorm(n_genes * n_cells), nrow=n_genes))
noisy_counts[noisy_counts < 2] <- 0
counts <- Matrix(noisy_counts, sparse=TRUE)

# 2. Save Data for Julia
write.csv(as.matrix(counts), "counts_matrix.csv")
cat("Saved counts_matrix.csv for Julia.\n")

# 3. Run R Benchmark
normalize_counts <- function(counts) {
  cs <- colSums(counts)
  cs[cs==0] <- 1
  log1p(t(t(counts) / cs) * 10000)
}

procrustes_align <- function(source, target) {
  sm <- colMeans(source); tm <- colMeans(target)
  s_c <- sweep(source, 2, sm); t_c <- sweep(target, 2, tm)
  svd_res <- svd(crossprod(s_c, t_c))
  rot <- svd_res$u %*% t(svd_res$v)
  sc <- sum(svd_res$d) / max(sum(s_c^2), .Machine$double.eps)
  list(aligned=sc * (s_c %*% rot) + matrix(tm, nrow=nrow(source), ncol=ncol(source), byrow=T), 
       rot=rot, scale=sc, sm=sm, tm=tm)
}

n_comp <- 30
ref_idx <- 1:400
query_idx <- 401:500

full_norm <- normalize_counts(counts)
full_pca <- prcomp(t(full_norm), center=T, scale.=F, rank.=n_comp)
full_truth <- full_pca$x

counts_ref <- counts[, ref_idx]
counts_query <- counts[, query_idx]

ref_norm <- normalize_counts(counts_ref)
ref_pca <- prcomp(t(ref_norm), center=T, scale.=F, rank.=n_comp)

# Projection
proj_query <- scale(t(normalize_counts(counts_query)), ref_pca$center, F) %*% ref_pca$rotation
proj_ref <- scale(t(ref_norm), ref_pca$center, F) %*% ref_pca$rotation

# Alignment
align_res <- procrustes_align(proj_ref, full_truth[ref_idx, ])
aligned_query <- align_res$scale * (sweep(proj_query, 2, align_res$sm) %*% align_res$rot) + matrix(align_res$tm, nrow=100, ncol=30, byrow=T)

# Report
cat("R RMSE:", sqrt(mean((aligned_query - full_truth[query_idx, ])^2)), "\n")
cat("R Cor PC1:", cor(aligned_query[,1], full_truth[query_idx, 1]), "\n")


# Output -->
# Saved counts_matrix.csv for Julia.
# R RMSE: 0.1181703 
# R Cor PC1: 1 