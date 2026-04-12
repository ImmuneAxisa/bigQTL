
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
# HELPER FUNCTION: Get cis SNPs for a phenotype on the fly
# =========================================================================

#' Get cis SNPs for a phenotype
#'
#' Returns the column indices and names of SNPs that fall within the cis
#' window of a given phenotype.
#'
#' @param bigsnp bigSNP object with \code{$map} containing \code{chromosome}
#'   and \code{physical.pos}
#' @param pheno_chr Chromosome of the phenotype (character or integer)
#' @param pheno_start Start position of phenotype (numeric)
#' @param pheno_end End position of phenotype (numeric)
#' @param cis_window Padding around phenotype coordinates in base pairs (default 1e6)
#'
#' @return Named list with elements: \code{indices} (integer vector of SNP
#'   column indices) and \code{names} (character vector of SNP names)
#' @export
get_cis_snps <- function(bigsnp, pheno_chr, pheno_start, pheno_end, cis_window = 1e6) {
  
  # Define cis window
  window_start <- pheno_start - cis_window
  window_end <- pheno_end + cis_window
  
  # Filter SNPs on same chromosome within window
  snp_map <- bigsnp$map
  
  cis_mask <- (snp_map$chromosome == pheno_chr) &
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
# HELPER FUNCTION: Get pheno indices from bigPheno
# =========================================================================

#' Match phenotype names to indices in a bigPheno object
#'
#' @param bigpheno bigPheno object with \code{$colData} containing
#'   a \code{pheno_name} column
#' @param pheno_names Character vector of phenotype names to look up
#'
#' @return Integer vector of column indices into the phenotype FBM
#' @export
get_pheno_indices <- function(bigpheno, pheno_names) {
  
  indices <- match(pheno_names, bigpheno$colData$pheno_name)
  
  if (any(is.na(indices))) {
    missing <- pheno_names[is.na(indices)]
    stop(sprintf("Could not find %d phenotypes: %s",
                 length(missing), paste(head(missing, 3), collapse = ", ")))
  }
  
  return(indices)
}


# =========================================================================
# HELPER FUNCTION: Get SNP indices from bigSNP
# =========================================================================

#' Match SNP names to indices in a bigSNP object
#'
#' @param bigsnp bigSNP object with \code{$map} containing
#'   a \code{marker.ID} column
#' @param snp_names Character vector of SNP names to look up
#'
#' @return Integer vector of column indices into the genotype FBM
#' @export
get_snp_indices <- function(bigsnp, snp_names) {
  
  indices <- match(snp_names, bigsnp$map$marker.ID)
  
  if (any(is.na(indices))) {
    missing <- snp_names[is.na(indices)]
    stop(sprintf("Could not find %d SNPs in bigSNP map: %s",
                 length(missing), paste(head(missing, 3), collapse = ", ")))
  }
  
  return(indices)
}
