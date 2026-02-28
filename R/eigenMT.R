# =========================================================================
# eigenMT: R implementation of eigenMT multiple testing correction
#
# Implements the eigenMT method (Davis et al., Am J Hum Genet 2016,
# https://doi.org/10.1016/j.ajhg.2016.01.029) for cis-QTL studies,
# rewritten in R and integrated with the bigstatsr/bigsnpr ecosystem.
#
# Original eigenMT repository: https://github.com/joed3/eigenMT
#
# Uses bigstatsr::big_cor() for C++-backed correlation computation
# directly on the file-backed genotype matrix (FBM), avoiding
# materialisation of genotype submatrices into R memory.
# =========================================================================


# =========================================================================
# HELPER FUNCTION: Ledoit-Wolf shrinkage (OAS formula)
# =========================================================================

#' Apply Ledoit-Wolf shrinkage to a correlation matrix
#'
#' Applies Ledoit-Wolf shrinkage to a pre-computed sample correlation matrix
#' using the OAS formula (Chen et al. 2010). Shrinkage target is the
#' identity matrix.
#'
#' @param R Square correlation matrix (p x p)
#' @param n Sample size used to compute R
#'
#' @return Shrunk correlation matrix of same dimensions as R
#' @export
lw_shrink_cor <- function(R, n) {

  p <- ncol(R)
  if (p == 1) return(matrix(1, 1, 1))

  # Degenerate case: all SNPs in perfect LD
  if (all(abs(R - 1) < 1e-10)) {
    return(matrix(1, p, p))
  }

  # OAS formula for identity target on correlation scale
  trR2 <- sum(R^2)
  trR  <- p

  numerator   <- (1 - 2 / p) * trR2 + trR^2
  denominator <- (n + 1 - 2 / p) * (trR2 - trR^2 / p)

  if (abs(denominator) < 1e-15) return(R)

  alpha <- max(0, min(1, numerator / denominator))

  shrunk <- (1 - alpha) * R
  diag(shrunk) <- 1

  return(shrunk)
}


# =========================================================================
# HELPER FUNCTION: Count eigenvalues for variance threshold
# =========================================================================

#' Count eigenvalues needed to explain a variance threshold
#'
#' Counts the minimum number of top eigenvalues needed to explain at least
#' \code{var_thresh} fraction of the total variance.
#'
#' @param eigenvalues Numeric vector of eigenvalues
#' @param var_thresh Variance fraction threshold (e.g. 0.99)
#'
#' @return Integer count of eigenvalues needed to reach the threshold
#' @export
count_eigenvalues <- function(eigenvalues, var_thresh) {

  eigenvalues <- sort(eigenvalues, decreasing = TRUE)
  total_var <- sum(eigenvalues)
  if (total_var <= 0) return(length(eigenvalues))

  target <- total_var * var_thresh
  running_sum <- 0
  counter <- 0L

  for (ev in eigenvalues) {
    running_sum <- running_sum + ev
    counter <- counter + 1L
    if (running_sum >= target) break
  }
  return(counter)
}


# =========================================================================
# CONVENIENCE FUNCTION: Apply eigenMT Bonferroni correction
# =========================================================================

#' Apply eigenMT Bonferroni correction
#'
#' Applies eigenMT Bonferroni correction to a nominal p-value.
#'
#' @param pvalue Numeric scalar or vector of nominal p-values
#' @param m_eff Effective number of independent tests (integer)
#'
#' @return Corrected p-value(s) capped at 1: \code{pmin(pvalue * m_eff, 1)}
#' @export
eigenMT_correct <- function(pvalue, m_eff) {

  pmin(pvalue * m_eff, 1)
}


# =========================================================================
# MAIN FUNCTION: Per-gene eigenMT M_eff computation
# =========================================================================

#' Compute eigenMT effective number of tests for a single gene
#'
#' Computes the effective number of independent tests (M_eff) for a single
#' gene using the eigenMT method. Splits cis-SNPs into disjoint windows of
#' size \code{window}, computes the shrunk correlation matrix in each window
#' via \code{bigstatsr::big_cor()} + Ledoit-Wolf shrinkage, eigendecomposes,
#' and counts eigenvalues needed to explain \code{var_thresh} of total
#' variance. Sums M_eff across windows.
#'
#' @param bigsnp bigSNP object with \code{$genotypes} (FBM) and \code{$map}
#' @param snp_indices Integer vector of column indices into the genotype FBM
#'   for cis-SNPs of this gene
#' @param ind.row Integer vector of row indices (individuals) to use
#' @param var_thresh Variance fraction threshold for eigenvalue counting
#'   (default 0.99)
#' @param window Maximum SNP window size for LD block computation (default 200)
#'
#' @return Integer M_eff: effective number of independent tests for this gene
#' @export
eigenMT_gene <- function(bigsnp, snp_indices, ind.row,
                         var_thresh = 0.99, window = 200) {

  M <- length(snp_indices)
  if (M == 0) return(0L)
  if (M == 1) return(1L)

  snp_indices <- sort(snp_indices)
  m_eff <- 0L
  start <- 1L

  while (start <= M) {
    stop <- min(start + window - 1L, M)
    win_size <- stop - start + 1L

    if (win_size == 1L) {
      m_eff <- m_eff + 1L
      break
    }

    win_indices <- snp_indices[start:stop]

    # Compute correlation matrix directly on the FBM (C++ backed)
    raw_cor <- bigstatsr::big_cor(
      X       = bigsnp$genotypes,
      ind.row = ind.row,
      ind.col = win_indices
    )[]

    # Apply Ledoit-Wolf shrinkage
    shrunk_cor <- lw_shrink_cor(raw_cor, n = length(ind.row))

    # Eigenvalues
    eigs <- eigen(shrunk_cor, symmetric = TRUE, only.values = TRUE)$values
    eigs[eigs < 0] <- 0

    m_eff <- m_eff + count_eigenvalues(eigs, var_thresh)
    start <- stop + 1L
  }

  return(m_eff)
}


# =========================================================================
# BATCH FUNCTION: Compute M_eff for all genes
# =========================================================================

#' Compute eigenMT M_eff for all genes in batch
#'
#' Computes eigenMT M_eff for each gene in \code{features_coord}. Uses
#' \code{get_cis_snps()} to identify cis-SNPs per gene, then calls
#' \code{eigenMT_gene()} for each. Parallelises over genes via
#' \code{parallel::mclapply}.
#'
#' @param bigsnp bigSNP object with \code{$genotypes} (FBM), \code{$fam},
#'   and \code{$map}
#' @param features_coord Data frame with columns: feature_name, chromosome,
#'   start, end
#' @param ind.row Integer vector of row indices (individuals) into the
#'   genotype FBM
#' @param cis_window Padding around gene coordinates for cis-SNP lookup
#'   (default 1e6)
#' @param var_thresh Variance fraction threshold for eigenvalue counting
#'   (default 0.99)
#' @param eigenmt_window SNP window size for LD block computation (default 200)
#' @param ncores Number of cores for parallel computation (default 1)
#'
#' @return Data frame with columns: feature_name, n_cis_snps, m_eff
#' @export
eigenMT_batch <- function(bigsnp, features_coord, ind.row,
                          cis_window = 1e6, var_thresh = 0.99,
                          eigenmt_window = 200, ncores = 1) {

  results <- parallel::mclapply(seq_len(nrow(features_coord)), function(i) {
    gene_row <- features_coord[i, ]

    cis_result <- get_cis_snps(
      bigsnp     = bigsnp,
      gene_chr   = gene_row$chromosome,
      gene_start = gene_row$start,
      gene_end   = gene_row$end,
      cis_window = cis_window
    )

    snp_indices <- cis_result$indices
    n_cis <- length(snp_indices)

    m_eff <- eigenMT_gene(
      bigsnp      = bigsnp,
      snp_indices = snp_indices,
      ind.row     = ind.row,
      var_thresh  = var_thresh,
      window      = eigenmt_window
    )

    data.frame(
      feature_name = gene_row$feature_name,
      n_cis_snps   = n_cis,
      m_eff        = m_eff,
      stringsAsFactors = FALSE
    )
  }, mc.cores = ncores)

  do.call(rbind, results)
}
