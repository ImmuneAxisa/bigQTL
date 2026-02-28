
# =========================================================================
# HELPER FUNCTION: Add SNPs as covariates
# =========================================================================

#' Add SNPs as covariates to a data frame
#'
#' Extracts SNP genotypes from a bigSNP object and appends them as columns
#' to a covariate data frame. Row names of \code{covar_df} are used to match
#' samples against \code{bigsnp$fam$sample.ID}.
#'
#' @param bigsnp bigSNP object with \code{$genotypes} (FBM) and \code{$map}
#'   (data frame with column \code{marker.ID})
#' @param snp_names Character vector of SNP names to extract
#' @param covar_df Data frame of covariates (samples x covariates). Row names
#'   must be sample IDs present in \code{bigsnp$fam$sample.ID}.
#'
#' @return Data frame with original covariates plus one column per SNP
#' @export
add_snps_to_covariates <- function(bigsnp, snp_names, covar_df) {
  
  # Derive row indices into genotype FBM from covar_df row names
  ind.row <- match(rownames(covar_df), bigsnp$fam$sample.ID)
  if (any(is.na(ind.row))) {
    missing <- rownames(covar_df)[is.na(ind.row)]
    stop(sprintf("%d sample ID(s) from covar_df not found in bigsnp$fam$sample.ID: %s",
                 length(missing), paste(head(missing, 3), collapse = ", ")))
  }
  
  # Match SNP names to indices in map
  snp_indices <- match(snp_names, bigsnp$map$marker.ID)
  
  if (any(is.na(snp_indices))) {
    missing <- snp_names[is.na(snp_indices)]
    stop(sprintf("Could not find %d SNPs in bigSNP map: %s",
                 length(missing), paste(head(missing, 3), collapse = ", ")))
  }
  
  # Extract SNP genotypes for kept individuals only and convert to data frame
  snp_genos <- as.data.frame(bigsnp$genotypes[ind.row, snp_indices, drop = FALSE])
  names(snp_genos) <- snp_names
  
  # Bind to covariates
  covar_combined <- cbind(covar_df, snp_genos)
  
  return(covar_combined)
}


# =========================================================================
# HELPER FUNCTION: Get cis SNPs for a gene on the fly
# =========================================================================

#' Get cis SNPs for a gene
#'
#' Returns the column indices and names of SNPs that fall within the cis
#' window of a given gene.
#'
#' @param bigsnp bigSNP object with \code{$map} containing \code{chromosome}
#'   and \code{physical.pos}
#' @param gene_chr Chromosome of the gene (character or integer)
#' @param gene_start Start position of gene (numeric)
#' @param gene_end End position of gene (numeric)
#' @param cis_window Padding around gene coordinates in base pairs (default 1e6)
#'
#' @return Named list with elements: \code{indices} (integer vector of SNP
#'   column indices) and \code{names} (character vector of SNP names)
#' @export
get_cis_snps <- function(bigsnp, gene_chr, gene_start, gene_end, cis_window = 1e6) {
  
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

#' Match feature names to indices in a bigFeatures object
#'
#' @param bigfeatures bigFeatures object with \code{$colData} containing
#'   a \code{feature_name} column
#' @param feature_names Character vector of feature names to look up
#'
#' @return Integer vector of column indices into the feature FBM
#' @export
get_feature_indices <- function(bigfeatures, feature_names) {
  
  indices <- match(feature_names, bigfeatures$colData$feature_name)
  
  if (any(is.na(indices))) {
    missing <- feature_names[is.na(indices)]
    stop(sprintf("Could not find %d features: %s",
                 length(missing), paste(head(missing, 3), collapse = ", ")))
  }
  
  return(indices)
}
