# =========================================================================
# eigenMT: R implementation of eigenMT multiple testing correction
#
# Implements the eigenMT method (Davis et al., Am J Hum Genet 2016,
# https://doi.org/10.1016/j.ajhg.2016.01.029) for cis-QTL studies,
# rewritten in R and integrated with the bigstatsr/bigsnpr ecosystem.
#
# Original eigenMT repository: https://github.com/joed3/eigenMT
# R port with nlshrink: https://github.com/ImmuneAxisa/eigenMT
#
# Two shrinkage methods are supported:
#   "basic"   - Ledoit-Wolf 2004 analytical shrinkage implemented in base R.
#               Applied to the biased sample covariance of the raw genotype
#               matrix (matching the approach of Python sklearn LedoitWolf).
#               Results are very close to Python and require no extra packages.
#   "nlshrink" - Requires the nlshrink package (in
#               Suggests). Gives results identical to the ImmuneAxisa/eigenMT
#               R port and to the Python sklearn implementation.
# =========================================================================


# =========================================================================
# HELPER FUNCTION: Ledoit-Wolf shrinkage on raw genotype matrix (basic)
# =========================================================================

#' Apply Ledoit-Wolf shrinkage to a raw genotype matrix
#'
#' Applies the analytical Ledoit-Wolf 2004 shrinkage estimator to the
#' biased sample covariance of a raw genotype matrix and returns the
#' corresponding shrunk correlation matrix. This is structurally equivalent
#' to Python \code{sklearn.covariance.LedoitWolf} and gives results very
#' close to the Python eigenMT without requiring any additional packages.
#'
#' @param X Numeric matrix of genotypes: \code{n} samples x \code{p} SNPs.
#'   Missing values should be imputed before calling this function.
#'
#' @return Shrunk correlation matrix of dimension \code{p x p}.
#' @export
lw_shrink_geno <- function(X) {

  n <- nrow(X)
  p <- ncol(X)

  if (p == 1L) return(matrix(1, 1, 1))

  # Center columns (sklearn LedoitWolf centres by default)
  X_c <- sweep(X, 2L, colMeans(X), "-")

  # Biased sample covariance (1/n denominator, matching sklearn)
  S <- crossprod(X_c) / n

  # Degenerate case: all genotypes are effectively constant (monomorphic
  # window). The covariance matrix is near-zero, so treat as perfect LD
  # and return the all-ones correlation matrix (m_eff = 1).
  if (max(abs(S)) < 1e-10) {
    return(matrix(1, p, p))
  }

  trace_S  <- sum(diag(S))
  trace_S2 <- sum(S * S)

  mu <- trace_S / p

  denominator <- (n + 2) * (trace_S2 - trace_S^2 / p)

  # denominator is zero only when S is a scaled identity, i.e., no excess
  # correlation signal; skip shrinkage in that case (rho = 0).
  if (abs(denominator) < 1e-15) {
    rho <- 0
  } else {
    rho <- min(1, max(0, ((n - 2) / n * trace_S2 + trace_S^2) / denominator))
  }

  Sigma_hat <- (1 - rho) * S
  diag(Sigma_hat) <- diag(Sigma_hat) + rho * mu

  # Convert to correlation matrix.
  # A near-zero standard deviation means a (near-)monomorphic SNP in this
  # window. Setting its sd to 1 yields an identity row/column — treating it
  # as an independent test, which is conservative and numerically safe.
  sd_vec <- sqrt(diag(Sigma_hat))
  sd_vec[sd_vec < sqrt(.Machine$double.eps)] <- 1
  shrunk_cor <- Sigma_hat / outer(sd_vec, sd_vec)
  diag(shrunk_cor) <- 1

  return(shrunk_cor)
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
#' size \code{window}, computes the shrunk correlation matrix in each window,
#' eigendecomposes, and counts eigenvalues needed to explain \code{var_thresh}
#' of total variance. Sums M_eff across windows.
#'
#' Two shrinkage methods are available via \code{shrinkage_method}:
#' \describe{
#'   \item{\code{"basic"}}{Analytical Ledoit-Wolf 2004 shrinkage applied to the
#'     biased sample covariance of the raw genotype matrix (no extra
#'     dependencies). Results are very close to the Python eigenMT.}
#'   \item{\code{"nlshrink"}}{Non-parametric linear shrinkage via
#'     \code{nlshrink::linshrink_cov()} (requires the \pkg{nlshrink} package).
#'     Gives results identical to the Python sklearn LedoitWolf implementation.}
#' }
#'
#' @param bigsnp bigSNP object with \code{$genotypes} (FBM) and \code{$map}
#' @param snp_indices Integer vector of column indices into the genotype FBM
#'   for cis-SNPs of this gene
#' @param ind.row Integer vector of row indices (individuals) to use
#' @param var_thresh Variance fraction threshold for eigenvalue counting
#'   (default 0.99)
#' @param window Maximum SNP window size for LD block computation (default 200)
#' @param shrinkage_method Character string specifying the shrinkage estimator:
#'   \code{"basic"} (default) or \code{"nlshrink"} (requires \pkg{nlshrink}).
#'
#' @return Integer M_eff: effective number of independent tests for this gene
#' @export
eigenMT_gene <- function(bigsnp, snp_indices, ind.row,
                         var_thresh = 0.99, window = 200,
                         shrinkage_method = c("basic", "nlshrink")) {

  shrinkage_method <- match.arg(shrinkage_method)

  if (shrinkage_method == "nlshrink" &&
      !requireNamespace("nlshrink", quietly = TRUE)) {
    stop(
      'shrinkage_method = "nlshrink" requires the nlshrink package. ',
      'Install it with: install.packages("nlshrink")'
    )
  }

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

    # Extract raw genotype matrix for this window: n_samples x win_size
    geno_mat <- bigsnp$genotypes[ind.row, win_indices, drop = FALSE]
    storage.mode(geno_mat) <- "double"

    # Compute shrunk correlation matrix
    shrunk_cor <- if (shrinkage_method == "nlshrink") {
      shrunk_cov <- nlshrink::linshrink_cov(geno_mat)
      cov2cor(shrunk_cov)
    } else {
      lw_shrink_geno(geno_mat)
    }

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

#' Compute eigenMT M_eff for all phenotypes in batch
#'
#' Computes eigenMT M_eff for each phenotype in \code{pheno_coord}. Uses
#' \code{get_cis_snps()} to identify cis-SNPs per phenotype, then calls
#' \code{eigenMT_gene()} for each. Parallelises over phenotypes via
#' \code{parallel::mclapply}.
#'
#' @param bigsnp bigSNP object with \code{$genotypes} (FBM), \code{$fam},
#'   and \code{$map}
#' @param pheno_coord Data frame with columns: pheno_name, chromosome,
#'   start, end
#' @param ind.row Integer vector of row indices (individuals) into the
#'   genotype FBM
#' @param cis_window Padding around phenotype coordinates for cis-SNP lookup
#'   (default 1e6)
#' @param var_thresh Variance fraction threshold for eigenvalue counting
#'   (default 0.99)
#' @param eigenmt_window SNP window size for LD block computation (default 200)
#' @param shrinkage_method Character string specifying the shrinkage estimator
#'   passed to \code{eigenMT_gene()}: \code{"basic"} (default) or
#'   \code{"nlshrink"} (requires \pkg{nlshrink}).
#' @param ncores Number of cores for parallel computation (default 1)
#'
#' @return Data frame with columns: pheno_name, n_cis_snps, m_eff
#' @export
eigenMT_batch <- function(bigsnp, pheno_coord, ind.row,
                          cis_window = 1e6, var_thresh = 0.99,
                          eigenmt_window = 200,
                          shrinkage_method = c("basic", "nlshrink"),
                          ncores = 1) {

  shrinkage_method <- match.arg(shrinkage_method)

  results <- parallel::mclapply(seq_len(nrow(pheno_coord)), function(i) {
    pheno_row <- pheno_coord[i, ]

    cis_result <- get_cis_snps(
      bigsnp      = bigsnp,
      pheno_chr   = pheno_row$chromosome,
      pheno_start = pheno_row$start,
      pheno_end   = pheno_row$end,
      cis_window  = cis_window
    )

    snp_indices <- cis_result$indices
    n_cis <- length(snp_indices)

    m_eff <- eigenMT_gene(
      bigsnp           = bigsnp,
      snp_indices      = snp_indices,
      ind.row          = ind.row,
      var_thresh       = var_thresh,
      window           = eigenmt_window,
      shrinkage_method = shrinkage_method
    )

    data.frame(
      pheno_name = pheno_row$pheno_name,
      n_cis_snps = n_cis,
      m_eff      = m_eff,
      stringsAsFactors = FALSE
    )
  }, mc.cores = ncores)

  do.call(rbind, results)
}
