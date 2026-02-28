
# =========================================================================
# HELPER FUNCTION: Genotype PCA using bigsnpr
# =========================================================================

#' Compute genotype principal components
#'
#' Runs PCA on genotype data using bigsnpr, optionally subsetting individuals.
#' Uses \code{bigsnpr::snp_autoSVD()} which performs PCA with automatic removal
#' of long-range LD regions and outlier variants (robust to LD structure).
#'
#' @param bigsnp bigSNP object with \code{$genotypes} (FBM), \code{$fam}
#'   (data frame with column \code{sample.ID}), and \code{$map} (data frame
#'   with columns \code{chromosome} and \code{physical.pos})
#' @param keep_ids Character vector of individual IDs to include. Must match
#'   values in \code{bigsnp$fam$sample.ID}. If NULL (default), all individuals
#'   are used.
#' @param n_pcs Number of principal components to retain (default 5)
#' @param ncores Number of cores for parallel computation (default 1)
#'
#' @return Matrix of genotype PCs (n_kept_samples x n_pcs), with columns named
#'   genoPC1, genoPC2, ..., genoPCn. Row names are the sample IDs of the kept
#'   individuals.
#' @export
compute_geno_pcs <- function(bigsnp, keep_ids = NULL, n_pcs = 5, ncores = 1) {

  message("Computing genotype PCs with snp_autoSVD...")

  if (!is.null(keep_ids)) {
    ind.row <- which(bigsnp$fam$sample.ID %in% keep_ids)
    if (length(ind.row) == 0) stop("None of the provided keep_ids matched bigsnp$fam$sample.ID")
    n_missing <- length(keep_ids) - length(ind.row)
    if (n_missing > 0) warning(sprintf("%d IDs in keep_ids were not found in bigsnp$fam$sample.ID", n_missing))
    message(sprintf("  Subsetting to %d of %d individuals", length(ind.row), nrow(bigsnp$fam)))
  } else {
    ind.row <- bigstatsr::rows_along(bigsnp$genotypes)
    message(sprintf("  Using all %d individuals", length(ind.row)))
  }

  svd_result <- bigsnpr::snp_autoSVD(
    G = bigsnp$genotypes, infos.chr = bigsnp$map$chromosome,
    infos.pos = bigsnp$map$physical.pos, ind.row = ind.row,
    k = n_pcs, ncores = ncores
  )

  geno_pcs <- stats::predict(svd_result)
  geno_pcs <- geno_pcs[, seq_len(n_pcs), drop = FALSE]
  colnames(geno_pcs) <- paste0("genoPC", seq_len(n_pcs))
  rownames(geno_pcs) <- bigsnp$fam$sample.ID[ind.row]
  message(sprintf("  Retained %d genotype PCs for %d samples", n_pcs, nrow(geno_pcs)))
  return(geno_pcs)
}


# =========================================================================
# HELPER FUNCTION: Feature PCA using bigstatsr (operates on FBM directly)
# =========================================================================

#' Compute feature principal components on file-backed matrix
#'
#' Runs PCA on feature (e.g. expression) data using bigstatsr, operating
#' directly on the FBM without loading data into memory. Uses
#' \code{bigstatsr::big_colstats()} to compute per-column variance on disk,
#' selects the top \code{n_top_features} by variance, then runs
#' \code{bigstatsr::big_randomSVD()} for PCA.
#'
#' @param bigfeatures bigFeatures object with \code{$features} (FBM,
#'   samples x features) and \code{$rowData} (data frame with column
#'   \code{sample_name})
#' @param keep_ids Character vector of sample IDs to include. Must match
#'   values in \code{bigfeatures$rowData$sample_name}. If NULL (default),
#'   all samples are used.
#' @param n_pcs Number of principal components to retain (default 5)
#' @param n_top_features Number of top-variance features to use for PCA
#'   (default 5000)
#' @param exclude_col_idx Integer vector of column indices to exclude from
#'   PCA (e.g. sex-chromosome features). Default NULL.
#' @param ncores Number of cores for parallel computation (default 1)
#'
#' @return Matrix of feature PCs (n_kept_samples x n_pcs), with columns named
#'   featurePC1, featurePC2, ..., featurePCn. Row names are the sample IDs
#'   of the kept samples.
#' @export
compute_feature_pcs <- function(bigfeatures, keep_ids = NULL, n_pcs = 5,
                                n_top_features = 5000, exclude_col_idx = NULL, ncores = 1) {

  message("Computing feature PCs with big_randomSVD...")

  fbm <- bigfeatures$features

  if (!is.null(keep_ids)) {
    ind.row <- which(bigfeatures$rowData$sample_name %in% keep_ids)
    if (length(ind.row) == 0) stop("None of the provided keep_ids matched bigfeatures$rowData$sample_name")
    n_missing <- length(keep_ids) - length(ind.row)
    if (n_missing > 0) warning(sprintf("%d IDs in keep_ids were not found in bigfeatures$rowData$sample_name", n_missing))
    message(sprintf("  Subsetting to %d of %d samples", length(ind.row), nrow(bigfeatures$rowData)))
  } else {
    ind.row <- bigstatsr::rows_along(fbm)
    message(sprintf("  Using all %d samples", length(ind.row)))
  }

  all_cols <- bigstatsr::cols_along(fbm)
  if (!is.null(exclude_col_idx)) {
    candidate_cols <- setdiff(all_cols, exclude_col_idx)
    message(sprintf("  Excluded %d columns, %d candidates remain", length(exclude_col_idx), length(candidate_cols)))
  } else {
    candidate_cols <- all_cols
  }

  col_stats <- bigstatsr::big_colstats(fbm, ind.row = ind.row, ind.col = candidate_cols, ncores = ncores)
  n_select <- min(n_top_features, length(candidate_cols))
  top_idx <- order(col_stats$var, decreasing = TRUE)[seq_len(n_select)]
  ind.col <- candidate_cols[top_idx]
  message(sprintf("  Selected %d features by variance for PCA", n_select))

  svd_result <- bigstatsr::big_randomSVD(
    X = fbm, fun.scaling = bigstatsr::big_scale(center = TRUE, scale = TRUE),
    ind.row = ind.row, ind.col = ind.col, k = n_pcs, ncores = ncores
  )

  feature_pcs <- stats::predict(svd_result)
  feature_pcs <- feature_pcs[, seq_len(min(n_pcs, ncol(feature_pcs))), drop = FALSE]
  colnames(feature_pcs) <- paste0("featurePC", seq_len(ncol(feature_pcs)))
  rownames(feature_pcs) <- bigfeatures$rowData$sample_name[ind.row]
  message(sprintf("  Retained %d feature PCs for %d samples", ncol(feature_pcs), nrow(feature_pcs)))
  return(feature_pcs)
}


# =========================================================================
# HELPER FUNCTION: Rank Inverse Normal Transform (RINT)
# =========================================================================

#' Rank Inverse Normal Transformation (Blom transform)
#'
#' Maps the ranks of \code{x} to quantiles of the standard normal distribution.
#' This is an inline replacement for \code{RNOmni::RankNorm()} with no external
#' dependencies.
#'
#' @param x Numeric vector to transform
#' @param k Blom offset constant (default 0.375). Controls continuity correction
#'   to avoid boundary quantiles; 0.375 is the recommended Blom constant.
#'
#' @return Numeric vector of same length as \code{x}, with values mapped to the
#'   standard normal distribution via rank-based quantile transformation.
#' @export
rint <- function(x, k = 0.375) {
  n <- length(x)
  r <- rank(x, ties.method = "average")
  stats::qnorm((r - k) / (n - 2 * k + 1))
}
