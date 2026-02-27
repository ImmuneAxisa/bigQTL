
# =========================================================================
# HELPER FUNCTION: Add SNPs as covariates
# =========================================================================

add_snps_to_covariates <- function(bigsnp, snp_names, covar_df,
                                   ind.row = NULL) {
  "
  Extract SNPs from bigSNP object and add them as columns to covariate data frame.
  Only extracts genotypes for the individuals specified by ind.row.

  Args:
    bigsnp:    bigSNP object with $genotypes (FBM) and $map (SNP metadata)
    snp_names: Character vector of SNP names to extract
    covar_df:  Data frame of covariates (samples x covariates).
               Must have the same number of rows as length(ind.row).
    ind.row:   Integer vector of row indices to extract. If NULL, all rows.

  Returns:
    Data frame with original covariates + SNP columns
  "
  
  library(tidyverse)
  
  if (is.null(ind.row)) {
    ind.row <- bigstatsr::rows_along(bigsnp$genotypes)
  }
  
  # Match SNP names to indices in map
  snp_indices <- match(snp_names, bigsnp$map$marker.ID)
  
  if (any(is.na(snp_indices))) {
    missing <- snp_names[is.na(snp_indices)]
    stop(sprintf("Could not find %d SNPs in bigSNP map: %s",
                 length(missing), paste(head(missing, 3), collapse = ", ")))
  }
  
  # Extract SNP genotypes for kept individuals only and convert to tibble
  snp_genos <- bigsnp$genotypes[ind.row, snp_indices, drop = FALSE] %>%
    as_tibble() %>%
    setNames(snp_names)
  
  # Bind to covariates
  covar_combined <- bind_cols(covar_df, snp_genos)
  
  return(covar_combined)
}


# =========================================================================
# HELPER FUNCTION: Get cis SNPs for a gene on the fly
# =========================================================================

get_cis_snps <- function(bigsnp, gene_chr, gene_start, gene_end, cis_window = 1e6) {
  "
  Get SNP indices and names that are in cis with a given gene.

  Args:
    bigsnp: bigSNP object with $map containing chromosome and physical.pos
    gene_chr: Chromosome of the gene
    gene_start: Start position of gene
    gene_end: End position of gene
    cis_window: Padding around gene coordinates (default 1e6)

  Returns:
    List with elements: indices (SNP column indices), names (SNP names)
  "
  
  # Define cis window
  window_start <- gene_start - cis_window
  window_end <- gene_end + cis_window
  
  # Filter SNPs on same chromosome within window
  snp_map <- bigsnp$map
  
  cis_mask <- (snp_map$chromosome == gene_chr) &
    (snp_map$physical.pos >= window_start) &
    (snp_map$physical.pos <= window_end)
  
  if (!any(cis_mask)) {
    return(list(indices = integer(0), names = character(0)))
  }
  
  cis_indices <- which(cis_mask)
  cis_names <- snp_map$marker.ID[cis_mask]
  
  return(list(indices = cis_indices, names = cis_names))
}


# =========================================================================
# HELPER FUNCTION: Get feature indices from bigFeatures
# =========================================================================

get_feature_indices <- function(bigfeatures, feature_names) {
  "
  Match feature names to indices in bigFeatures object.

  Args:
    bigfeatures: bigFeatures object
    feature_names: Character vector of feature names

  Returns:
    Integer vector of column indices (features are columns)
  "
  
  indices <- match(feature_names, bigfeatures$colData$feature_name)
  
  if (any(is.na(indices))) {
    missing <- feature_names[is.na(indices)]
    stop(sprintf("Could not find %d features: %s",
                 length(missing), paste(head(missing, 3), collapse = ", ")))
  }
  
  return(indices)
}


# =========================================================================
# HELPER FUNCTION: Genotype PCA using bigsnpr
# =========================================================================

compute_geno_pcs <- function(bigsnp,
                             keep_ids = NULL,
                             n_pcs = 5,
                             ncores = 1) {
  "
  Run PCA on genotype data using bigsnpr, optionally subsetting individuals.

  Uses snp_autoSVD which performs PCA with automatic removal of
  long-range LD regions and outlier variants (robust to LD structure).

  Args:
    bigsnp:   bigSNP object (as returned by snp_attach / snp_readBed)
    keep_ids: Character vector of individual IDs to include. Must match
              values in bigsnp$fam$sample.ID. If NULL (default), all
              individuals are used.
    n_pcs:    Number of principal components to retain (default 5)
    ncores:   Number of cores for parallel computation (default 1)

  Returns:
    Matrix of genotype PCs (n_kept_samples x n_pcs), with columns named
    genoPC1, genoPC2, ..., genoPCn.
    Row names are set to the sample IDs of the kept individuals.
  "
  
  library(bigsnpr)
  
  message("Computing genotype PCs with snp_autoSVD...")
  
  # Convert IDs to row indices
  if (!is.null(keep_ids)) {
    ind.row <- which(bigsnp$fam$sample.ID %in% keep_ids)
    if (length(ind.row) == 0) {
      stop("None of the provided keep_ids matched bigsnp$fam$sample.ID")
    }
    n_missing <- length(keep_ids) - length(ind.row)
    if (n_missing > 0) {
      warning(sprintf("%d IDs in keep_ids were not found in bigsnp$fam$sample.ID", n_missing))
    }
    message(sprintf("  Subsetting to %d of %d individuals",
                    length(ind.row), nrow(bigsnp$fam)))
  } else {
    ind.row <- rows_along(bigsnp$genotypes)
    message(sprintf("  Using all %d individuals", length(ind.row)))
  }
  
  svd_result <- snp_autoSVD(
    G         = bigsnp$genotypes,
    infos.chr = bigsnp$map$chromosome,
    infos.pos = bigsnp$map$physical.pos,
    ind.row   = ind.row,
    k         = n_pcs,
    ncores    = ncores
  )
  
  # Extract sample scores (U * D), equivalent to flashpca $vectors
  geno_pcs <- predict(svd_result)
  
  geno_pcs <- geno_pcs[, seq_len(n_pcs), drop = FALSE]
  colnames(geno_pcs) <- paste0("genoPC", seq_len(n_pcs))
  rownames(geno_pcs) <- bigsnp$fam$sample.ID[ind.row]
  
  message(sprintf("  Retained %d genotype PCs for %d samples", n_pcs, nrow(geno_pcs)))
  
  return(geno_pcs)
}


# =========================================================================
# HELPER FUNCTION: Feature (expression/phenotype) PCA using prcomp
# =========================================================================

compute_feature_pcs <- function(feature_matrix,
                                keep_ids = NULL,
                                n_pcs = 5,
                                n_top_features = 5000,
                                mean_threshold = 6,
                                exclude_pattern = "^chr[XY]",
                                exclude_names = NULL) {
  "
  Run PCA on feature (e.g. expression) data using prcomp, optionally
  subsetting to specific samples.

  Filters features by mean expression, excludes sex-chromosome and
  user-specified features, selects top variable features, then runs
  prcomp on mean-centered data.

  Args:
    feature_matrix: Numeric matrix of features (features x samples).
                    Column names must be sample identifiers.
    keep_ids:       Character vector of sample IDs to include. Must match
                    column names of feature_matrix. If NULL (default), all
                    samples are used.
    n_pcs:          Number of principal components to retain (default 5)
    n_top_features: Number of top-variable features to use for PCA (default 5000)
    mean_threshold: Minimum row mean to include a feature (default 6)
    exclude_pattern: Regex pattern for feature names to exclude,
                     e.g. sex chromosomes (default '^chr[XY]')
    exclude_names:  Optional character vector of feature names to exclude
                    (e.g. ribosomal, mitochondrial genes). Default NULL.

  Returns:
    Matrix of feature PCs (n_kept_samples x n_pcs), with columns named
    featurePC1, featurePC2, ..., featurePCn.
    Row names are set to the sample IDs of the kept samples.
  "
  
  library(MatrixGenerics)
  
  message("Computing feature PCs with prcomp...")
  
  # Subset samples if requested
  if (!is.null(keep_ids)) {
    sample_mask <- colnames(feature_matrix) %in% keep_ids
    if (!any(sample_mask)) {
      stop("None of the provided keep_ids matched column names of feature_matrix")
    }
    n_missing <- length(keep_ids) - sum(sample_mask)
    if (n_missing > 0) {
      warning(sprintf("%d IDs in keep_ids were not found in feature_matrix columns", n_missing))
    }
    feature_matrix <- feature_matrix[, sample_mask, drop = FALSE]
    message(sprintf("  Subsetting to %d samples", ncol(feature_matrix)))
  } else {
    message(sprintf("  Using all %d samples", ncol(feature_matrix)))
  }
  
  # Compute per-feature statistics on the (possibly subsetted) matrix
  feat_means <- rowMeans(feature_matrix)
  feat_sds <- rowSds(feature_matrix)
  feature_names <- rownames(feature_matrix)
  
  # Build logical mask
  keep <- (feat_means > mean_threshold)
  
  if (!is.null(exclude_pattern) && nchar(exclude_pattern) > 0) {
    keep <- keep & !grepl(exclude_pattern, feature_names)
  }
  
  if (!is.null(exclude_names) && length(exclude_names) > 0) {
    keep <- keep & !(feature_names %in% exclude_names)
  }
  
  if (sum(keep) == 0) {
    stop("No features remain after filtering. Check mean_threshold, exclude_pattern, and exclude_names.")
  }
  
  # Among passing features, select top N by SD
  candidate_sds <- feat_sds[keep]
  candidate_names <- feature_names[keep]
  
  n_select <- min(n_top_features, length(candidate_names))
  top_idx <- order(candidate_sds, decreasing = TRUE)[seq_len(n_select)]
  selected_features <- candidate_names[top_idx]
  
  message(sprintf("  Selected %d features for PCA (from %d candidates passing filters)",
                  n_select, length(candidate_names)))
  
  # Run PCA: transpose (samples x features), mean-center, then prcomp
  pca_result <- feature_matrix[selected_features, , drop = FALSE] %>%
    t() %>%
    scale(scale = FALSE) %>%
    prcomp()
  
  # Extract scores
  feature_pcs <- pca_result$x[, seq_len(min(n_pcs, ncol(pca_result$x))), drop = FALSE]
  colnames(feature_pcs) <- paste0("featurePC", seq_len(ncol(feature_pcs)))
  
  message(sprintf("  Retained %d feature PCs for %d samples", ncol(feature_pcs), nrow(feature_pcs)))
  
  return(feature_pcs)
}


# =========================================================================
# HELPER FUNCTION: Rank Inverse Normal Transform (RINT)
# =========================================================================

apply_rint <- function(feature_matrix) {
  "
  Apply Rank Inverse Normal Transformation (RINT) to each feature (row).

  Uses RNOmni::RankNorm which maps ranks to quantiles of a standard
  normal distribution. This is applied per feature (row-wise) to
  remove distributional artefacts before association testing.

  Args:
    feature_matrix: Numeric matrix (features x samples)

  Returns:
    Numeric matrix of same dimensions, with each row RINT-transformed
  "
  
  library(RNOmni)
  
  message("Applying RINT to features...")
  
  rint_matrix <- t(apply(feature_matrix, 1, RankNorm))
  
  message(sprintf("  Transformed %d features across %d samples",
                  nrow(rint_matrix), ncol(rint_matrix)))
  
  return(rint_matrix)
}