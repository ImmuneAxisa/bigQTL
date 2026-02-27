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

lw_shrink_cor <- function(R, n) {
  "
  Apply Ledoit-Wolf shrinkage to a pre-computed sample correlation matrix
  using the OAS formula (Chen et al. 2010). Shrinkage target is the
  identity matrix.

  Args:
    R: Square correlation matrix (p x p)
    n: Sample size used to compute R

  Returns:
    Shrunk correlation matrix of same dimensions as R
  "

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

count_eigenvalues <- function(eigenvalues, var_thresh) {
  "
  Count the minimum number of top eigenvalues needed to explain at least
  var_thresh fraction of the total variance.

  Args:
    eigenvalues: Numeric vector of eigenvalues
    var_thresh:  Variance fraction threshold (e.g. 0.99)

  Returns:
    Integer count of eigenvalues needed to reach the threshold
  "

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

eigenMT_correct <- function(pvalue, m_eff) {
  "
  Apply eigenMT Bonferroni correction to a nominal p-value.

  Args:
    pvalue: Numeric scalar or vector of nominal p-values
    m_eff:  Effective number of independent tests (integer)

  Returns:
    Corrected p-value(s) capped at 1: pmin(pvalue * m_eff, 1)
  "

  pmin(pvalue * m_eff, 1)
}


# =========================================================================
# MAIN FUNCTION: Per-gene eigenMT M_eff computation
# =========================================================================

eigenMT_gene <- function(bigsnp, snp_indices, ind.row,
                         var_thresh = 0.99, window = 200) {
  "
  Compute the effective number of independent tests (M_eff) for a single
  gene using the eigenMT method.

  Splits cis-SNPs into disjoint windows of size `window`, computes the
  shrunk correlation matrix in each window via bigstatsr::big_cor() +
  Ledoit-Wolf shrinkage, eigendecomposes, and counts eigenvalues needed
  to explain var_thresh of total variance. Sums M_eff across windows.

  Args:
    bigsnp:      bigSNP object with $genotypes (FBM) and $map
    snp_indices: Integer vector of column indices into the genotype FBM
                 for cis-SNPs of this gene
    ind.row:     Integer vector of row indices (individuals) to use
    var_thresh:  Variance fraction threshold for eigenvalue counting
                 (default 0.99)
    window:      Maximum SNP window size for LD block computation
                 (default 200)

  Returns:
    Integer M_eff: effective number of independent tests for this gene
  "

  library(bigstatsr)

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
    raw_cor <- big_cor(
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

eigenMT_batch <- function(bigsnp, features_coord, ind.row,
                          cis_window = 1e6, var_thresh = 0.99,
                          eigenmt_window = 200, ncores = 1) {
  "
  Compute eigenMT M_eff for each gene in features_coord.

  Uses get_cis_snps() from helpers.R to identify cis-SNPs per gene,
  then calls eigenMT_gene() for each. Parallelises over genes via
  mclapply.

  Args:
    bigsnp:          bigSNP object with $genotypes (FBM), $fam, and $map
    features_coord:  Data frame with columns: feature_name, chromosome,
                     start, end
    ind.row:         Integer vector of row indices (individuals) into the
                     genotype FBM
    cis_window:      Padding around gene coordinates for cis-SNP lookup
                     (default 1e6)
    var_thresh:      Variance fraction threshold for eigenvalue counting
                     (default 0.99)
    eigenmt_window:  SNP window size for LD block computation (default 200)
    ncores:          Number of cores for parallel computation (default 1)

  Returns:
    Tibble with columns: feature_name, n_cis_snps, m_eff
  "

  library(bigstatsr)
  library(tidyverse)

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

    tibble(
      feature_name = gene_row$feature_name,
      n_cis_snps   = n_cis,
      m_eff        = m_eff
    )
  }, mc.cores = ncores)

  bind_rows(results)
}
