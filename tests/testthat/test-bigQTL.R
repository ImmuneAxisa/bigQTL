# bigQTL() had zero test coverage (see dev/CODE_REVIEW.md section 3): its
# PC-computation-and-augmentation logic (compute_pheno_pcs()/
# compute_geno_pcs() -> cbind() -> run_conditional_qtl()) is exactly the
# sample-alignment-sensitive code flagged in section 2.4, and it's the
# function most likely to be a new user's first call into the package.
#
# These tests exercise it end-to-end on `dummy_qtl_data`, a realistic
# fixture with 5 phenotypes / 3,2,1,0,0 known independent causal signals
# (dev/generate_dummy_qtl_data.R), and lock in the sample-ID realignment
# invariant with a genuinely permuted design_base.

# compute_pheno_pcs()/compute_geno_pcs() use randomized SVD algorithms; a
# fixed seed right before bigQTL() makes recovery reproducible (see
# dev/generate_dummy_qtl_data.R for how this seed was validated).
run_bigqtl_on_dummy_data <- function(sample_order = NULL, ...) {
  d <- load_dummy_qtl_data(sample_order = sample_order)
  out_dir <- tempfile("bigqtl_test_")
  set.seed(1)
  list(
    result = bigQTL(
      bigpheno    = d$bigpheno,
      bigsnp      = d$bigsnp,
      pheno_coord = d$pheno_coord,
      design_base = d$design_base,
      n_pheno_pcs = 3,
      n_geno_pcs  = 3,
      min_snps    = 10,
      output_dir  = out_dir,
      verbose     = FALSE,
      ...
    ),
    out_dir = out_dir,
    ground_truth = d$ground_truth
  )
}

# Number of independent signals recovered by stepwise conditioning for one
# phenotype's rows of the stepwise output. Step 0 (marginal) is always
# present even when it fails pval_threshold, so max(step) alone can't tell
# "0 signals" apart from "1 signal, no further conditioning step passed";
# any step >= 1 present unambiguously implies step 0 passed.
count_recovered_signals <- function(df_pheno, pval_threshold = 1e-3) {
  max_step <- max(df_pheno$step)
  if (max_step >= 1) return(as.integer(max_step) + 1L)
  step0 <- df_pheno[df_pheno$step == 0, ]
  if (min(step0$pvalue) < pval_threshold) 1L else 0L
}

test_that("bigQTL runs end-to-end and returns the expected structure", {
  out <- run_bigqtl_on_dummy_data()
  on.exit(unlink(out$out_dir, recursive = TRUE), add = TRUE)

  expect_type(out$result, "list")
  expect_named(out$result, c("stepwise", "allbutone", "stepwise_dir", "allbutone_dir"))
  expect_true(inherits(out$result$stepwise, "Dataset"))
  expect_true(inherits(out$result$allbutone, "Dataset"))
})

test_that("bigQTL recovers the correct number of independent signals per phenotype", {
  out <- run_bigqtl_on_dummy_data()
  on.exit(unlink(out$out_dir, recursive = TRUE), add = TRUE)

  df <- as.data.frame(out$result$stepwise)
  recovered <- vapply(names(out$ground_truth$n_signals), function(p) {
    count_recovered_signals(df[df$pheno == p, ])
  }, integer(1))

  expect_equal(unname(recovered), unname(out$ground_truth$n_signals))
})

test_that("bigQTL recovers the true causal SNP as lead of its own step", {
  out <- run_bigqtl_on_dummy_data()
  on.exit(unlink(out$out_dir, recursive = TRUE), add = TRUE)

  df <- as.data.frame(out$result$stepwise)
  causal <- out$ground_truth$causal

  for (i in seq_len(nrow(causal))) {
    pheno_rows <- df[df$pheno == causal$pheno_name[i], ]
    # The causal SNP must appear as the minimum-pvalue ("lead") SNP of
    # exactly one stepwise step (order across steps is not guaranteed).
    lead_snps <- vapply(split(pheno_rows, pheno_rows$step), function(tbl) {
      tbl$snp[which.min(tbl$pvalue)]
    }, character(1))
    expect_true(causal$snp[i] %in% lead_snps)
  }
})

test_that("bigQTL is invariant to sample-order permutation between design_base and bigsnp/bigpheno", {
  d0 <- load_dummy_qtl_data()
  permuted_order <- sample(d0$sample_ids)

  out_a <- run_bigqtl_on_dummy_data(sample_order = NULL)
  on.exit(unlink(out_a$out_dir, recursive = TRUE), add = TRUE)
  out_b <- run_bigqtl_on_dummy_data(sample_order = permuted_order)
  on.exit(unlink(out_b$out_dir, recursive = TRUE), add = TRUE)

  cols <- c("step", "conditioning_snps", "snp", "beta", "se", "t_stat", "pvalue", "fdr", "pheno")
  df_a <- as.data.frame(out_a$result$stepwise)[, cols]
  df_b <- as.data.frame(out_b$result$stepwise)[, cols]

  # A SNP being conditioned on is also re-tested as one of the cis-SNPs
  # (run_stepwise() doesn't exclude it), which produces a rank-deficient
  # design (the SNP is perfectly collinear with its own covariate column).
  # bigstatsr correctly reports se/t_stat/pvalue as NaN there, but the
  # coefficient itself is an arbitrary solution of a singular system --
  # it's allowed to differ with floating-point row order and isn't part of
  # the invariant this test checks, so drop those rows before comparing.
  df_a <- df_a[!is.na(df_a$pvalue), ]
  df_b <- df_b[!is.na(df_b$pvalue), ]

  ord_a <- order(df_a$pheno, df_a$step, df_a$snp)
  ord_b <- order(df_b$pheno, df_b$step, df_b$snp)

  expect_equal(df_a[ord_a, ], df_b[ord_b, ], ignore_attr = TRUE)
})

test_that("bigQTL with a permuted design_base still recovers the correct signals", {
  d0 <- load_dummy_qtl_data()
  permuted_order <- sample(d0$sample_ids)

  out <- run_bigqtl_on_dummy_data(sample_order = permuted_order)
  on.exit(unlink(out$out_dir, recursive = TRUE), add = TRUE)

  df <- as.data.frame(out$result$stepwise)
  recovered <- vapply(names(out$ground_truth$n_signals), function(p) {
    count_recovered_signals(df[df$pheno == p, ])
  }, integer(1))

  expect_equal(unname(recovered), unname(out$ground_truth$n_signals))
})
