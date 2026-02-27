
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