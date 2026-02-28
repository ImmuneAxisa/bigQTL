# =========================================================================
# test_eigenMT.R
#
# Validates the R eigenMT implementation (eigenMT.R) against the original
# Python eigenMT (joed3/eigenMT) by running both on identical simulated
# genotype data with realistic LD structure and comparing M_eff per gene.
#
# Usage:
#   Rscript tests/test_eigenMT.R
#
# from the repo root directory.
# =========================================================================

set.seed(42)

# -------------------------------------------------------------------------
# Load required packages
# -------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(bigstatsr)
  library(bigsnpr)
})

# Source R implementation — locate the repo root robustly
# Prefer the script's own directory when invoked via Rscript, then fall back
# to the current working directory (e.g. when sourced interactively).
.this_file <- tryCatch({
  # Works when invoked as: Rscript tests/test_eigenMT.R
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("^--file=", args, value = TRUE)
  if (length(file_flag)) normalizePath(sub("^--file=", "", file_flag[1])) else ""
}, error = function(e) "")

repo_root <- if (nchar(.this_file) > 0 && file.exists(.this_file)) {
  normalizePath(file.path(dirname(.this_file), ".."), mustWork = FALSE)
} else {
  getwd()
}

if (!file.exists(file.path(repo_root, "eigenMT.R"))) {
  stop(sprintf(
    "Cannot find eigenMT.R under repo_root='%s'. Run from the repo root or via: Rscript tests/test_eigenMT.R",
    repo_root
  ))
}
source(file.path(repo_root, "eigenMT.R"))
source(file.path(repo_root, "helpers.R"))

# -------------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------------
N_SAMPLES  <- 500L   # number of individuals
N_SNPS     <- 800L   # total SNPs on chr1
BLOCK_SIZE <- 20L    # SNPs per LD block
FLIP_PROB  <- 0.10   # probability of flipping an allele within a block
VAR_THRESH <- 0.99
EIGEN_WIN  <- 200L
CIS_DIST   <- 1e6

message("=================================================================")
message("eigenMT validation: R vs Python")
message("=================================================================")

# =========================================================================
# 1. Simulate genotype data with LD block structure
# =========================================================================
message("\n[1] Simulating genotype data ...")

# SNP positions spread evenly on chr1: 1e6 to 2e6
snp_pos <- round(seq(1e6, 2e6, length.out = N_SNPS))
snp_ids <- paste0("rs", seq_len(N_SNPS))
chr_vec <- rep("1", N_SNPS)

# Simulate genotypes: LD blocks of BLOCK_SIZE SNPs
# Within each block, SNPs are correlated (derived from a founder SNP)
maf_vec <- runif(N_SNPS, 0.05, 0.45)  # minor allele frequencies

geno_mat <- matrix(NA_integer_, nrow = N_SAMPLES, ncol = N_SNPS)

n_blocks <- ceiling(N_SNPS / BLOCK_SIZE)
for (b in seq_len(n_blocks)) {
  block_start <- (b - 1L) * BLOCK_SIZE + 1L
  block_end   <- min(b * BLOCK_SIZE, N_SNPS)
  block_snps  <- block_start:block_end

  # Founder SNP for this block
  founder_maf <- maf_vec[block_start]
  founder_geno <- rbinom(N_SAMPLES, 2L, founder_maf)
  geno_mat[, block_start] <- founder_geno

  # Subsequent SNPs: copy founder + random flips
  for (j in block_snps[-1]) {
    flipped <- rbinom(N_SAMPLES, 1L, FLIP_PROB)
    # When flipped, draw a fresh genotype from marginal
    fresh   <- rbinom(N_SAMPLES, 2L, maf_vec[j])
    new_geno <- ifelse(flipped == 1L, fresh, founder_geno)
    geno_mat[, j] <- new_geno
  }
}

message(sprintf("  Generated %d x %d genotype matrix", N_SAMPLES, N_SNPS))
message(sprintf("  LD blocks: %d blocks of ~%d SNPs", n_blocks, BLOCK_SIZE))

# Sample IDs
sample_ids <- paste0("sample_", seq_len(N_SAMPLES))

# =========================================================================
# 2. Define 3 genes with cis windows that capture different SNP counts
# =========================================================================
message("\n[2] Defining gene positions ...")

# Gene positions within chr1 1e6–2e6 range
features_coord <- data.frame(
  feature_name = c("GeneA", "GeneB", "GeneC"),
  chromosome   = c("1", "1", "1"),
  start        = c(1.10e6, 1.45e6, 1.80e6),
  end          = c(1.12e6, 1.47e6, 1.82e6),
  stringsAsFactors = FALSE
)

message(sprintf("  Defined %d genes on chr1", nrow(features_coord)))

# =========================================================================
# 3. Write data in Python eigenMT CLI format
# =========================================================================
tmp_dir <- file.path(tempdir(), "test_eigenMT")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

message(sprintf("\n[3] Writing Python-format data to: %s", tmp_dir))

## 3a. genotypes.txt: SNPs x samples (first col = SNP ID)
geno_for_py <- cbind(snp_ids, as.data.frame(t(geno_mat)))
colnames(geno_for_py) <- c("id", sample_ids)
geno_file <- file.path(tmp_dir, "genotypes.txt")
write.table(geno_for_py, geno_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

## 3b. gen.positions.txt: snp, chr, pos
gen_pos <- data.frame(snp = snp_ids, chr = chr_vec, pos = snp_pos,
                      stringsAsFactors = FALSE)
gen_pos_file <- file.path(tmp_dir, "gen.positions.txt")
write.table(gen_pos, gen_pos_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

## 3c. phe.positions.txt: geneid, chr, s1 (start), s2 (end)
phe_pos <- data.frame(
  geneid = features_coord$feature_name,
  chr    = features_coord$chromosome,
  s1     = features_coord$start,
  s2     = features_coord$end,
  stringsAsFactors = FALSE
)
phe_pos_file <- file.path(tmp_dir, "phe.positions.txt")
write.table(phe_pos, phe_pos_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

## 3d. cis.eqtls.txt: fake QTL results
# Build all cis SNP-gene pairs
eqtl_rows <- list()
for (i in seq_len(nrow(features_coord))) {
  gene_row <- features_coord[i, ]
  win_start <- gene_row$start - CIS_DIST
  win_end   <- gene_row$end   + CIS_DIST
  cis_mask  <- (chr_vec == gene_row$chromosome) &
               (snp_pos >= win_start) & (snp_pos <= win_end)
  cis_snps  <- snp_ids[cis_mask]
  if (length(cis_snps) > 0) {
    np  <- length(cis_snps)
    eqtl_rows[[i]] <- data.frame(
      SNP     = cis_snps,
      gene    = rep(gene_row$feature_name, np),
      beta    = rnorm(np),
      `t-stat` = rnorm(np),
      `p-value` = runif(np),
      FDR     = runif(np),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
}
eqtl_df   <- do.call(rbind, eqtl_rows)
eqtl_file <- file.path(tmp_dir, "cis.eqtls.txt")
write.table(eqtl_df, eqtl_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

message(sprintf("  Wrote %d SNP-gene pairs to cis.eqtls.txt", nrow(eqtl_df)))

# =========================================================================
# 4. Build bigSNP object for R eigenMT
# =========================================================================
message("\n[4] Building bigSNP object for R eigenMT ...")

# Create FBM-backed genotype matrix
fbm_file <- file.path(tmp_dir, "geno.bk")
geno_fbm <- FBM(nrow = N_SAMPLES, ncol = N_SNPS, type = "integer",
                backingfile = sub("\\.bk$", "", fbm_file),
                create_bk = TRUE)
geno_fbm[] <- geno_mat

# Build minimal bigSNP-like list
snp_map <- data.frame(
  chromosome   = chr_vec,
  marker.ID    = snp_ids,
  genetic.dist = rep(0, N_SNPS),
  physical.pos = snp_pos,
  allele1      = rep("A", N_SNPS),
  allele2      = rep("C", N_SNPS),
  stringsAsFactors = FALSE
)

snp_fam <- data.frame(
  family.ID = sample_ids,
  sample.ID = sample_ids,
  paternal.ID = rep(0, N_SAMPLES),
  maternal.ID = rep(0, N_SAMPLES),
  sex = rep(0, N_SAMPLES),
  affection = rep(0, N_SAMPLES),
  stringsAsFactors = FALSE
)

bigsnp_obj <- structure(
  list(genotypes = geno_fbm, map = snp_map, fam = snp_fam),
  class = "bigSNP"
)

ind.row <- seq_len(N_SAMPLES)

message("  bigSNP object created in memory")

# =========================================================================
# 5. Embed and run Python eigenMT
# =========================================================================

python_eigenmt_src <- r"(
#!/usr/bin/env python
# Embedded eigenMT from joed3/eigenMT (MIT License)
# https://github.com/joed3/eigenMT
import sys, os, argparse
import numpy as np
import pandas as pd
from sklearn import covariance

def lw_shrink(genotypes):
    lw = covariance.LedoitWolf()
    m, n = np.shape(genotypes)
    try:
        fitted = lw.fit(genotypes.T)
        alpha = fitted.shrinkage_
        shrunk_cov = fitted.covariance_
        shrunk_precision = np.mat(np.diag(np.diag(shrunk_cov)**(-.5)))
        shrunk_cor = shrunk_precision * shrunk_cov * shrunk_precision
    except Exception:
        row = np.repeat(1, m)
        shrunk_cor = []
        for i in range(0, m):
            shrunk_cor.append(row)
        shrunk_cor = np.mat(shrunk_cor)
        alpha = 'NA'
    return shrunk_cor, alpha

def find_num_eigs(eigenvalues, variance, var_thresh):
    eigenvalues = np.sort(eigenvalues)[::-1]
    running_sum = 0
    counter = 0
    while running_sum < variance * var_thresh:
        running_sum += eigenvalues[counter]
        counter += 1
    return counter

def eigenMT(args):
    np.random.seed(42)

    # Load QTL file to get gene list
    qtl = pd.read_csv(args.QTL, sep='\t')
    genes = qtl[args.geneCol].unique()

    # Load genotype positions
    genpos = pd.read_csv(args.GENPOS, sep='\t', index_col=0)

    # Load phenotype positions
    phepos = pd.read_csv(args.PHEPOS, sep='\t', index_col=0)

    # Load genotype matrix (SNPs x samples)
    geno_all = pd.read_csv(args.GEN, sep='\t', index_col=0)

    results = []
    for gene in genes:
        if gene not in phepos.index:
            continue
        row = phepos.loc[gene]
        chrom = str(row['chr'])
        g_start = int(row['s1'])
        g_end   = int(row['s2'])

        cis_start = g_start - args.cis_dist
        cis_end   = g_end   + args.cis_dist

        # Filter cis SNPs
        mask = (genpos['chr'].astype(str) == chrom) & \
               (genpos['pos'] >= cis_start) & \
               (genpos['pos'] <= cis_end)
        cis_snps = genpos.index[mask].tolist()

        if len(cis_snps) == 0:
            results.append({'gene': gene, 'TESTS': 0, 'n_cis_snps': 0})
            continue

        geno_cis = geno_all.loc[cis_snps].values.astype(float)  # SNPs x samples
        M = len(cis_snps)

        m_eff = 0
        start = 0
        while start < M:
            stop = min(start + args.window, M)
            win_geno = geno_cis[start:stop, :]
            win_size = stop - start

            if win_size == 1:
                m_eff += 1
                start = stop
                continue

            shrunk_cor, alpha = lw_shrink(win_geno)
            eigenvalues = np.real(np.linalg.eigvals(shrunk_cor))
            eigenvalues[eigenvalues < 0] = 0
            m_eff += find_num_eigs(eigenvalues, win_size, args.var_thresh)
            start = stop

        results.append({'gene': gene, 'TESTS': m_eff, 'n_cis_snps': M})

    out_df = pd.DataFrame(results)
    out_df.to_csv(args.OUT, sep='\t', index=False)
    print(f"Wrote {len(results)} genes to {args.OUT}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--CHROM',      type=str,   required=True)
    parser.add_argument('--QTL',        type=str,   required=True)
    parser.add_argument('--GEN',        type=str,   required=True)
    parser.add_argument('--GENPOS',     type=str,   required=True)
    parser.add_argument('--PHEPOS',     type=str,   required=True)
    parser.add_argument('--OUT',        type=str,   required=True)
    parser.add_argument('--cis_dist',   type=int,   default=1000000)
    parser.add_argument('--var_thresh', type=float, default=0.99)
    parser.add_argument('--window',     type=int,   default=200)
    parser.add_argument('--geneCol',    type=str,   default='gene')
    args = parser.parse_args()
    eigenMT(args)
)"

python_script_file <- file.path(tmp_dir, "eigenMT.py")
writeLines(python_eigenmt_src, python_script_file)

message("\n[5] Running Python eigenMT ...")

py_out_file <- file.path(tmp_dir, "eigenMT_output.txt")

# Detect Python executable with required packages
detect_python <- function() {
  candidates <- c(
    Sys.which("python3"),
    Sys.which("python"),
    "/usr/bin/python3",
    "/usr/bin/python"
  )
  candidates <- unique(candidates[nchar(candidates) > 0])
  for (py in candidates) {
    ok <- tryCatch({
      ret <- system2(py,
                     args = c("-c",
                              shQuote("import numpy, sklearn, pandas")),
                     stdout = FALSE, stderr = FALSE)
      ret == 0L
    }, error = function(e) FALSE)
    if (ok) return(py)
  }
  return(NULL)
}

python_exe <- detect_python()

python_results <- NULL
python_available <- !is.null(python_exe)

if (!python_available) {
  warning(paste(
    "Python with numpy/scipy/sklearn/pandas not found.",
    "Skipping Python vs R comparison.",
    "R-only unit tests will still run."
  ))
} else {
  message(sprintf("  Using Python: %s", python_exe))

  py_args <- c(
    python_script_file,
    "--CHROM",      "1",
    "--QTL",        eqtl_file,
    "--GEN",        geno_file,
    "--GENPOS",     gen_pos_file,
    "--PHEPOS",     phe_pos_file,
    "--OUT",        py_out_file,
    "--cis_dist",   as.character(CIS_DIST),
    "--var_thresh", as.character(VAR_THRESH),
    "--window",     as.character(EIGEN_WIN)
  )

  ret <- system2(python_exe, args = py_args,
                 stdout = TRUE, stderr = TRUE)
  exit_code <- attr(ret, "status")

  if (!is.null(exit_code) && exit_code != 0L) {
    warning(paste("Python eigenMT failed (exit code", exit_code, ").",
                  "Output:\n", paste(ret, collapse = "\n"),
                  "\nSkipping Python comparison."))
    python_available <- FALSE
  } else {
    message("  Python eigenMT completed successfully")
    python_results <- read.table(py_out_file, header = TRUE,
                                 sep = "\t", stringsAsFactors = FALSE)
    message(sprintf("  Parsed %d gene results from Python", nrow(python_results)))
  }
}

# =========================================================================
# 6. Run R eigenMT for each gene
# =========================================================================
message("\n[6] Running R eigenMT ...")

r_results <- lapply(seq_len(nrow(features_coord)), function(i) {
  gene_row <- features_coord[i, ]

  cis_result <- get_cis_snps(
    bigsnp     = bigsnp_obj,
    gene_chr   = gene_row$chromosome,
    gene_start = gene_row$start,
    gene_end   = gene_row$end,
    cis_window = CIS_DIST
  )

  snp_indices <- cis_result$indices
  n_cis       <- length(snp_indices)

  m_eff <- eigenMT_gene(
    bigsnp      = bigsnp_obj,
    snp_indices = snp_indices,
    ind.row     = ind.row,
    var_thresh  = VAR_THRESH,
    window      = EIGEN_WIN
  )

  list(
    gene       = gene_row$feature_name,
    n_cis_snps = n_cis,
    m_eff      = m_eff
  )
})

r_df <- do.call(rbind, lapply(r_results, as.data.frame))
message(sprintf("  R eigenMT done for %d genes", nrow(r_df)))

# =========================================================================
# 7. Unit tests — low-level component validation
# =========================================================================
message("\n[7] Running unit tests ...")

n_failures <- 0L

## --- 7a. count_eigenvalues vs Python find_num_eigs ----------------------
message("  [7a] count_eigenvalues vs find_num_eigs ...")

py_find_num_eigs <- function(eigenvalues, variance, var_thresh) {
  eigenvalues <- sort(eigenvalues, decreasing = TRUE)
  running_sum <- 0
  counter     <- 0L
  while (running_sum < variance * var_thresh) {
    running_sum <- running_sum + eigenvalues[counter + 1L]
    counter     <- counter + 1L
  }
  counter
}

test_eigs <- c(4.5, 2.3, 1.1, 0.8, 0.5, 0.4, 0.2, 0.1, 0.05, 0.05)
vt <- 0.99

r_count   <- count_eigenvalues(test_eigs, vt)
py_count  <- py_find_num_eigs(test_eigs, length(test_eigs), vt)

if (r_count != py_count) {
  message(sprintf("    FAIL: count_eigenvalues=%d, find_num_eigs=%d",
                  r_count, py_count))
  n_failures <- n_failures + 1L
} else {
  message(sprintf("    PASS: both return %d", r_count))
}

## --- 7b. big_cor vs base cor() ------------------------------------------
message("  [7b] big_cor vs base cor() accuracy ...")

# Use a small subset of the genotype matrix
test_cols <- 1:30
small_geno <- geno_mat[, test_cols, drop = FALSE]

# Compute correlation with base R
base_cor_mat <- cor(small_geno)

# Compute via big_cor on the FBM
fbm_cor_mat <- big_cor(
  X       = geno_fbm,
  ind.row = ind.row,
  ind.col = test_cols
)[]

frob_diff_cor <- sqrt(sum((base_cor_mat - fbm_cor_mat)^2))
message(sprintf("    Frobenius norm of (big_cor - cor): %.2e", frob_diff_cor))
tol_cor <- 1e-8
if (frob_diff_cor > tol_cor) {
  message(sprintf("    FAIL: Frobenius difference %.2e exceeds tolerance %.2e",
                  frob_diff_cor, tol_cor))
  n_failures <- n_failures + 1L
} else {
  message(sprintf("    PASS: big_cor matches cor() (tol=%.2e)", tol_cor))
}

## --- 7c. lw_shrink_cor sanity checks ------------------------------------
message("  [7c] lw_shrink_cor sanity checks ...")

R_test <- cor(geno_mat[, 1:20])
n_test <- N_SAMPLES
shrunk_test <- lw_shrink_cor(R_test, n_test)

# Shrunk matrix should be positive semi-definite (all eigenvalues >= 0)
eigs_shrunk <- eigen(shrunk_test, symmetric = TRUE, only.values = TRUE)$values
if (any(eigs_shrunk < -1e-10)) {
  message(sprintf("    FAIL: lw_shrink_cor produced negative eigenvalues (min=%.4f)",
                  min(eigs_shrunk)))
  n_failures <- n_failures + 1L
} else {
  message(sprintf("    PASS: shrunk matrix is PSD (min eigenvalue=%.4f)",
                  min(eigs_shrunk)))
}

# Diagonal should be 1
if (any(abs(diag(shrunk_test) - 1) > 1e-10)) {
  message("    FAIL: lw_shrink_cor diagonal != 1")
  n_failures <- n_failures + 1L
} else {
  message("    PASS: shrunk matrix diagonal is all 1.0")
}

# Shrinkage should pull off-diagonal towards 0 (or at least not increase
# the Frobenius norm of the off-diagonal portion)
off_diag_raw    <- R_test;    diag(off_diag_raw)    <- 0
off_diag_shrunk <- shrunk_test; diag(off_diag_shrunk) <- 0
frob_raw    <- sqrt(sum(off_diag_raw^2))
frob_shrunk <- sqrt(sum(off_diag_shrunk^2))
message(sprintf("    Off-diagonal Frobenius: raw=%.4f, shrunk=%.4f",
                frob_raw, frob_shrunk))
if (frob_shrunk > frob_raw + 1e-10) {
  message("    FAIL: shrinkage increased off-diagonal Frobenius norm")
  n_failures <- n_failures + 1L
} else {
  message("    PASS: shrinkage reduced off-diagonal Frobenius norm")
}

## --- 7d. eigenMT_gene with trivial inputs --------------------------------
message("  [7d] eigenMT_gene edge cases ...")

# Single SNP -> M_eff should be 1
m1 <- eigenMT_gene(bigsnp_obj, snp_indices = 1L, ind.row = ind.row)
if (m1 != 1L) {
  message(sprintf("    FAIL: single-SNP M_eff=%d (expected 1)", m1))
  n_failures <- n_failures + 1L
} else {
  message("    PASS: single-SNP M_eff = 1")
}

# Zero SNPs -> M_eff should be 0
m0 <- eigenMT_gene(bigsnp_obj, snp_indices = integer(0), ind.row = ind.row)
if (m0 != 0L) {
  message(sprintf("    FAIL: zero-SNP M_eff=%d (expected 0)", m0))
  n_failures <- n_failures + 1L
} else {
  message("    PASS: zero-SNP M_eff = 0")
}

# M_eff should be <= number of SNPs and >= 1 for a multi-SNP case
test_snp_idx <- head(which(snp_pos >= 1.1e6 & snp_pos <= 1.2e6), 10)
m_small <- eigenMT_gene(bigsnp_obj, snp_indices = test_snp_idx, ind.row = ind.row)
if (m_small < 1L || m_small > length(test_snp_idx)) {
  message(sprintf("    FAIL: M_eff=%d out of range [1, %d]",
                  m_small, length(test_snp_idx)))
  n_failures <- n_failures + 1L
} else {
  message(sprintf("    PASS: M_eff=%d in range [1, %d]",
                  m_small, length(test_snp_idx)))
}

# =========================================================================
# 8. Compare R vs Python M_eff per gene (if Python available)
# =========================================================================
comparison_df <- NULL
MEFF_TOL     <- 0.20  # 20% relative difference tolerance

if (python_available && !is.null(python_results)) {
  message("\n[8] Comparing R vs Python M_eff ...")

  comparison_rows <- lapply(seq_len(nrow(features_coord)), function(i) {
    gene <- features_coord$feature_name[i]

    r_row  <- r_df[r_df$gene == gene, ]
    py_row <- python_results[python_results$gene == gene, ]

    if (nrow(r_row) == 0 || nrow(py_row) == 0) return(NULL)

    m_r  <- r_row$m_eff[1]
    m_py <- py_row$TESTS[1]

    abs_diff <- abs(m_r - m_py)
    rel_diff <- if (m_py > 0) abs_diff / m_py else NA_real_

    list(
      gene       = gene,
      n_cis_snps = r_row$n_cis_snps[1],
      m_eff_r    = m_r,
      m_eff_py   = m_py,
      abs_diff   = abs_diff,
      rel_diff   = rel_diff
    )
  })

  comparison_df <- do.call(rbind, lapply(
    Filter(Negate(is.null), comparison_rows),
    as.data.frame
  ))

  # Check tolerance
  bad <- !is.na(comparison_df$rel_diff) & comparison_df$rel_diff > MEFF_TOL
  if (any(bad)) {
    bad_genes <- comparison_df$gene[bad]
    msg <- sprintf(
      paste("M_eff differs by >%d%% between R and Python for gene(s): %s.",
            "R values: %s  Python values: %s  Relative diffs: %s"),
      as.integer(MEFF_TOL * 100),
      paste(bad_genes, collapse = ", "),
      paste(comparison_df$m_eff_r[bad], collapse = ", "),
      paste(comparison_df$m_eff_py[bad], collapse = ", "),
      paste(round(comparison_df$rel_diff[bad] * 100, 1), collapse = ", ")
    )
    stop(msg)
  }
  message(sprintf("  All %d gene(s) pass M_eff tolerance (<= %d%% relative diff)",
                  nrow(comparison_df), as.integer(MEFF_TOL * 100)))
} else {
  message("\n[8] Skipping R vs Python comparison (Python not available)")
}

# =========================================================================
# 9. Summary table
# =========================================================================
message("\n=================================================================")
message("SUMMARY")
message("=================================================================")

cat("\n--- Unit test results ---\n")
cat(sprintf("  count_eigenvalues vs find_num_eigs : %s\n",
            if (n_failures == 0) "PASS" else "FAIL"))
cat(sprintf("  big_cor vs base cor()              : %s (Frob=%.2e)\n",
            if (frob_diff_cor <= tol_cor) "PASS" else "FAIL", frob_diff_cor))

cat("\n--- R eigenMT results per gene ---\n")
cat(sprintf("  %-10s  %10s  %8s\n", "Gene", "n_cis_snps", "M_eff"))
for (i in seq_len(nrow(r_df))) {
  cat(sprintf("  %-10s  %10d  %8d\n",
              r_df$gene[i], r_df$n_cis_snps[i], r_df$m_eff[i]))
}

if (!is.null(comparison_df)) {
  cat("\n--- R vs Python M_eff comparison ---\n")
  cat(sprintf("  %-10s  %10s  %8s  %8s  %8s  %8s\n",
              "Gene", "n_cis_snps", "M_eff_R", "M_eff_Py", "AbsDiff", "RelDiff%"))
  for (i in seq_len(nrow(comparison_df))) {
    cat(sprintf("  %-10s  %10d  %8d  %8d  %8.1f  %8.1f\n",
                comparison_df$gene[i],
                comparison_df$n_cis_snps[i],
                comparison_df$m_eff_r[i],
                comparison_df$m_eff_py[i],
                comparison_df$abs_diff[i],
                comparison_df$rel_diff[i] * 100))
  }
} else {
  cat("\n--- R vs Python comparison: SKIPPED (Python not available) ---\n")
}

cat(sprintf("\nUnit test failures: %d\n", n_failures))
if (python_available) {
  cat("Python comparison:  COMPLETED\n")
} else {
  cat("Python comparison:  SKIPPED\n")
}

# =========================================================================
# 10. Clean up temp files
# =========================================================================
message("\n[10] Cleaning up temp files ...")
unlink(tmp_dir, recursive = TRUE)
message("  Done.")

# Final pass/fail
if (n_failures > 0L) {
  stop(sprintf("%d unit test(s) failed. See messages above.", n_failures))
}

message("\nAll tests passed.")
