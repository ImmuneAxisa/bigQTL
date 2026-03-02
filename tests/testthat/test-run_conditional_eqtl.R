# Helper: build minimal synthetic data for run_conditional_eqtl() tests.
# Returns a list with bigsnp, bigfeatures, features_coord, and design_base.
make_eqtl_test_data <- function(n_samples = 20, n_snps = 5, seed = 42) {
  set.seed(seed)

  sample_ids <- paste0("s", seq_len(n_samples))

  # Genotype FBM (samples x SNPs)
  geno_mat <- matrix(
    sample(0:2, n_samples * n_snps, replace = TRUE),
    nrow = n_samples, ncol = n_snps
  )
  geno_fbm <- bigstatsr::as_FBM(geno_mat)

  # Minimal bigsnp-like list
  bigsnp <- list(
    genotypes = geno_fbm,
    fam = data.frame(sample.ID = sample_ids, stringsAsFactors = FALSE),
    map = data.frame(
      chromosome   = rep("1", n_snps),
      marker.ID    = paste0("rs", seq_len(n_snps)),
      physical.pos = as.integer(seq(950000L, 1050000L, length.out = n_snps)),
      stringsAsFactors = FALSE
    )
  )

  # Expression matrix wrapped in bigFeatures (samples x 1 gene)
  expr_mat <- matrix(rnorm(n_samples), nrow = n_samples, ncol = 1)
  rownames(expr_mat) <- sample_ids
  colnames(expr_mat) <- "geneA"
  bigfeatures <- bigFeatures(expr_mat)

  # Gene coordinates (SNPs fall within cis window of 1 Mb)
  features_coord <- data.frame(
    feature_name = "geneA",
    chromosome   = "1",
    start        = 1000000L,
    end          = 1000100L,
    stringsAsFactors = FALSE
  )

  # Design matrix (one covariate) with sample IDs as row names
  design_base <- data.frame(
    covar1 = rnorm(n_samples),
    row.names = sample_ids,
    stringsAsFactors = FALSE
  )

  list(
    bigsnp        = bigsnp,
    bigfeatures   = bigfeatures,
    features_coord = features_coord,
    design_base   = design_base,
    sample_ids    = sample_ids
  )
}

# =========================================================================
# Happy path
# =========================================================================

test_that("run_conditional_eqtl runs without error with minimal valid inputs", {
  d <- make_eqtl_test_data()
  out_dir <- tempfile("eqtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  expect_no_error(
    run_conditional_eqtl(
      bigfeatures    = d$bigfeatures,
      bigsnp         = d$bigsnp,
      features_coord = d$features_coord,
      design_base    = d$design_base,
      output_dir     = out_dir
    )
  )
})

# =========================================================================
# Return structure
# =========================================================================

test_that("run_conditional_eqtl returns a list with expected named elements", {
  d <- make_eqtl_test_data()
  out_dir <- tempfile("eqtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- run_conditional_eqtl(
    bigfeatures    = d$bigfeatures,
    bigsnp         = d$bigsnp,
    features_coord = d$features_coord,
    design_base    = d$design_base,
    output_dir     = out_dir
  )

  expect_type(result, "list")
  expect_named(result, c("stepwise", "allbutone", "stepwise_dir", "allbutone_dir"))
  expect_type(result$stepwise_dir, "character")
  expect_type(result$allbutone_dir, "character")
  # stepwise must be an Arrow Dataset
  expect_true(inherits(result$stepwise, "Dataset"))
})

# =========================================================================
# Conditional disabling
# =========================================================================

test_that("run_conditional_eqtl with do_conditioning=FALSE returns valid structure", {
  d <- make_eqtl_test_data()
  out_dir <- tempfile("eqtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- run_conditional_eqtl(
    bigfeatures     = d$bigfeatures,
    bigsnp          = d$bigsnp,
    features_coord  = d$features_coord,
    design_base     = d$design_base,
    do_conditioning = FALSE,
    output_dir      = out_dir
  )

  expect_named(result, c("stepwise", "allbutone", "stepwise_dir", "allbutone_dir"))
  expect_true(inherits(result$stepwise, "Dataset"))
})

test_that("run_conditional_eqtl with do_allbutone=FALSE returns NULL allbutone", {
  d <- make_eqtl_test_data()
  out_dir <- tempfile("eqtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- run_conditional_eqtl(
    bigfeatures    = d$bigfeatures,
    bigsnp         = d$bigsnp,
    features_coord = d$features_coord,
    design_base    = d$design_base,
    do_allbutone   = FALSE,
    output_dir     = out_dir
  )

  expect_null(result$allbutone)
})

# =========================================================================
# Output directories
# =========================================================================

test_that("run_conditional_eqtl creates stepwise and allbutone output directories", {
  d <- make_eqtl_test_data()
  out_dir <- tempfile("eqtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- run_conditional_eqtl(
    bigfeatures    = d$bigfeatures,
    bigsnp         = d$bigsnp,
    features_coord = d$features_coord,
    design_base    = d$design_base,
    output_dir     = out_dir
  )

  expect_true(dir.exists(result$stepwise_dir))
  expect_true(dir.exists(result$allbutone_dir))
})

# =========================================================================
# Error cases
# =========================================================================

test_that("run_conditional_eqtl errors when design_base has default integer row names", {
  d <- make_eqtl_test_data()
  out_dir <- tempfile("eqtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  rownames(d$design_base) <- NULL  # resets to "1", "2", ...

  expect_error(
    run_conditional_eqtl(
      bigfeatures    = d$bigfeatures,
      bigsnp         = d$bigsnp,
      features_coord = d$features_coord,
      design_base    = d$design_base,
      output_dir     = out_dir
    ),
    "meaningful row names"
  )
})

test_that("run_conditional_eqtl errors when design_base sample IDs missing from bigsnp", {
  d <- make_eqtl_test_data()
  out_dir <- tempfile("eqtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  rownames(d$design_base) <- paste0("x", seq_len(nrow(d$design_base)))

  expect_error(
    run_conditional_eqtl(
      bigfeatures    = d$bigfeatures,
      bigsnp         = d$bigsnp,
      features_coord = d$features_coord,
      design_base    = d$design_base,
      output_dir     = out_dir
    ),
    "not found in bigsnp"
  )
})

test_that("run_conditional_eqtl errors when design_base sample IDs missing from bigfeatures", {
  d <- make_eqtl_test_data()
  out_dir <- tempfile("eqtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  d$bigfeatures$rowData$sample_name <- paste0("y", seq_len(nrow(d$design_base)))

  expect_error(
    run_conditional_eqtl(
      bigfeatures    = d$bigfeatures,
      bigsnp         = d$bigsnp,
      features_coord = d$features_coord,
      design_base    = d$design_base,
      output_dir     = out_dir
    ),
    "not found in bigfeatures"
  )
})

test_that("run_conditional_eqtl errors when features_coord is missing required columns", {
  d <- make_eqtl_test_data()
  out_dir <- tempfile("eqtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  d$features_coord$chromosome <- NULL  # drop required column

  expect_error(
    run_conditional_eqtl(
      bigfeatures    = d$bigfeatures,
      bigsnp         = d$bigsnp,
      features_coord = d$features_coord,
      design_base    = d$design_base,
      output_dir     = out_dir
    ),
    "features_coord must contain columns"
  )
})
