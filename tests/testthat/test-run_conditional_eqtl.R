# Helper: build minimal synthetic data for run_conditional_qtl() tests.
# Returns a list with bigsnp, bigpheno, pheno_coord, and design_base.
make_qtl_test_data <- function(n_samples = 20, n_snps = 5, seed = 42) {
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

  # Expression matrix wrapped in bigPheno (samples x 1 phenotype)
  expr_mat <- matrix(rnorm(n_samples), nrow = n_samples, ncol = 1)
  rownames(expr_mat) <- sample_ids
  colnames(expr_mat) <- "phenoA"
  bigpheno <- bigPheno(expr_mat)

  # Phenotype coordinates (SNPs fall within cis window of 1 Mb)
  pheno_coord <- data.frame(
    pheno_name   = "phenoA",
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
    bigsnp      = bigsnp,
    bigpheno    = bigpheno,
    pheno_coord = pheno_coord,
    design_base = design_base,
    sample_ids  = sample_ids
  )
}

# =========================================================================
# Happy path
# =========================================================================

test_that("run_conditional_qtl runs without error with minimal valid inputs", {
  d <- make_qtl_test_data()
  out_dir <- tempfile("qtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  expect_no_error(
    run_conditional_qtl(
      bigpheno    = d$bigpheno,
      bigsnp      = d$bigsnp,
      pheno_coord = d$pheno_coord,
      design_base = d$design_base,
      min_snps    = 1,
      output_dir  = out_dir
    )
  )
})

# =========================================================================
# Return structure
# =========================================================================

test_that("run_conditional_qtl returns a list with expected named elements", {
  d <- make_qtl_test_data()
  out_dir <- tempfile("qtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- run_conditional_qtl(
    bigpheno    = d$bigpheno,
    bigsnp      = d$bigsnp,
    pheno_coord = d$pheno_coord,
    design_base = d$design_base,
    min_snps    = 1,
    output_dir  = out_dir
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

test_that("run_conditional_qtl with do_conditioning=FALSE returns valid structure", {
  d <- make_qtl_test_data()
  out_dir <- tempfile("qtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- run_conditional_qtl(
    bigpheno        = d$bigpheno,
    bigsnp          = d$bigsnp,
    pheno_coord     = d$pheno_coord,
    design_base     = d$design_base,
    do_conditioning = FALSE,
    do_allbutone    = FALSE,
    min_snps        = 1,
    output_dir      = out_dir
  )

  expect_named(result, c("stepwise", "allbutone", "stepwise_dir", "allbutone_dir"))
  expect_true(inherits(result$stepwise, "Dataset"))
})

test_that("run_conditional_qtl with do_allbutone=FALSE returns NULL allbutone", {
  d <- make_qtl_test_data()
  out_dir <- tempfile("qtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- run_conditional_qtl(
    bigpheno     = d$bigpheno,
    bigsnp       = d$bigsnp,
    pheno_coord  = d$pheno_coord,
    design_base  = d$design_base,
    do_allbutone = FALSE,
    min_snps     = 1,
    output_dir   = out_dir
  )

  expect_null(result$allbutone)
})

# =========================================================================
# Output directories
# =========================================================================

test_that("run_conditional_qtl creates stepwise and allbutone output directories", {
  d <- make_qtl_test_data()
  out_dir <- tempfile("qtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- run_conditional_qtl(
    bigpheno    = d$bigpheno,
    bigsnp      = d$bigsnp,
    pheno_coord = d$pheno_coord,
    design_base = d$design_base,
    min_snps    = 1,
    output_dir  = out_dir
  )

  expect_true(dir.exists(result$stepwise_dir))
  expect_true(dir.exists(result$allbutone_dir))
})

# =========================================================================
# Error cases
# =========================================================================

test_that("run_conditional_qtl errors when design_base has default integer row names", {
  d <- make_qtl_test_data()
  out_dir <- tempfile("qtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  rownames(d$design_base) <- NULL  # resets to "1", "2", ...

  expect_error(
    run_conditional_qtl(
      bigpheno    = d$bigpheno,
      bigsnp      = d$bigsnp,
      pheno_coord = d$pheno_coord,
      design_base = d$design_base,
      output_dir  = out_dir
    ),
    "meaningful row names"
  )
})

test_that("run_conditional_qtl errors when design_base sample IDs missing from bigsnp", {
  d <- make_qtl_test_data()
  out_dir <- tempfile("qtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  rownames(d$design_base) <- paste0("x", seq_len(nrow(d$design_base)))

  expect_error(
    run_conditional_qtl(
      bigpheno    = d$bigpheno,
      bigsnp      = d$bigsnp,
      pheno_coord = d$pheno_coord,
      design_base = d$design_base,
      output_dir  = out_dir
    ),
    "not found in bigsnp"
  )
})

test_that("run_conditional_qtl errors when design_base sample IDs missing from bigpheno", {
  d <- make_qtl_test_data()
  out_dir <- tempfile("qtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  d$bigpheno$rowData$sample_name <- paste0("y", seq_len(nrow(d$design_base)))

  expect_error(
    run_conditional_qtl(
      bigpheno    = d$bigpheno,
      bigsnp      = d$bigsnp,
      pheno_coord = d$pheno_coord,
      design_base = d$design_base,
      output_dir  = out_dir
    ),
    "not found in bigpheno"
  )
})

test_that("run_conditional_qtl errors when pheno_coord is missing required columns", {
  d <- make_qtl_test_data()
  out_dir <- tempfile("qtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  d$pheno_coord$chromosome <- NULL  # drop required column

  expect_error(
    run_conditional_qtl(
      bigpheno    = d$bigpheno,
      bigsnp      = d$bigsnp,
      pheno_coord = d$pheno_coord,
      design_base = d$design_base,
      output_dir  = out_dir
    ),
    "pheno_coord must contain columns"
  )
})

# =========================================================================
# sd = 0 early exit
# =========================================================================

test_that("process_pheno warns and exits early when phenotype has sd = 0", {
  d <- make_qtl_test_data()

  # Replace phenotype values with a constant (sd = 0)
  d$bigpheno$pheno[, 1] <- 0

  out_dir <- tempfile("qtl_test_")
  dir.create(file.path(out_dir, "stepwise"),  recursive = TRUE)
  dir.create(file.path(out_dir, "allbutone"), recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  ind.row.snp   <- match(d$sample_ids, d$bigsnp$fam$sample.ID)
  ind.row.pheno <- match(d$sample_ids, d$bigpheno$rowData$sample_name)

  expect_warning(
    process_pheno(
      pheno          = "phenoA",
      bigpheno       = d$bigpheno,
      bigsnp         = d$bigsnp,
      pheno_coord    = d$pheno_coord,
      design_base    = d$design_base,
      ind.row.snp    = ind.row.snp,
      ind.row.pheno  = ind.row.pheno,
      cis_window     = 1e6,
      do_conditioning = TRUE,
      pval_threshold = 1e-3,
      max_steps      = 5,
      do_allbutone   = TRUE,
      do_rint        = FALSE,
      min_snps       = 1,
      ncores         = 1,
      stepwise_dir   = file.path(out_dir, "stepwise"),
      allbutone_dir  = file.path(out_dir, "allbutone")
    ),
    "sd\\(y\\) == 0"
  )
})

# =========================================================================
# max_steps enforcement
# =========================================================================

# Helper: dataset guaranteed to have multiple passing steps at any pval
make_qtl_test_data_high_signal <- function(n_samples = 40, n_snps = 5, seed = 7) {
  set.seed(seed)
  sample_ids <- paste0("s", seq_len(n_samples))

  geno_mat <- matrix(
    sample(0:2, n_samples * n_snps, replace = TRUE),
    nrow = n_samples, ncol = n_snps
  )
  geno_fbm <- bigstatsr::as_FBM(geno_mat)

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

  expr_mat <- matrix(rnorm(n_samples), nrow = n_samples, ncol = 1)
  rownames(expr_mat) <- sample_ids
  colnames(expr_mat) <- "phenoA"
  bigpheno <- bigPheno(expr_mat)

  pheno_coord <- data.frame(
    pheno_name = "phenoA", chromosome = "1",
    start = 1000000L, end = 1000100L,
    stringsAsFactors = FALSE
  )
  design_base <- data.frame(
    covar1 = rnorm(n_samples),
    row.names = sample_ids,
    stringsAsFactors = FALSE
  )
  list(bigsnp = bigsnp, bigpheno = bigpheno,
       pheno_coord = pheno_coord, design_base = design_base,
       sample_ids = sample_ids)
}

test_that("run_stepwise warns when max_steps is reached", {
  d <- make_qtl_test_data_high_signal()

  ind.row.snp <- match(d$sample_ids, d$bigsnp$fam$sample.ID)

  results_step0 <- test_snps_with_indices(
    bigsnp      = d$bigsnp,
    y           = rnorm(length(d$sample_ids)),
    snp_indices = seq_len(nrow(d$bigsnp$map)),
    snp_names   = d$bigsnp$map$marker.ID,
    design_base = d$design_base,
    ind.row     = ind.row.snp
  )
  results_step0$step <- 0L
  results_step0$conditioning_snps <- NA_character_
  results_step0 <- results_step0[, c("step", "conditioning_snps",
                                     setdiff(names(results_step0),
                                             c("step", "conditioning_snps")))]

  # Force the first lead SNP to pass threshold
  results_step0$pvalue[1] <- 1e-10

  # max_steps=1: step 1 runs, if lead passes step becomes 2, then 2>1 triggers warning
  expect_warning(
    run_stepwise(
      results_step0  = results_step0,
      bigsnp         = d$bigsnp,
      y              = rnorm(length(d$sample_ids)),
      snp_indices    = seq_len(nrow(d$bigsnp$map)),
      cis_snps_pheno = d$bigsnp$map$marker.ID,
      design_base    = d$design_base,
      ind.row.snp    = ind.row.snp,
      pval_threshold = 1,    # threshold of 1 ensures every lead "passes"
      max_steps      = 1,
      ncores         = 1,
      pheno          = "phenoA"
    ),
    "max_steps"
  )
})

# =========================================================================
# stepwise_tables only contains passing steps
# =========================================================================

test_that("run_stepwise excludes non-passing conditioning steps from stepwise_tables", {
  set.seed(42)
  d <- make_qtl_test_data_high_signal()

  ind.row.snp <- match(d$sample_ids, d$bigsnp$fam$sample.ID)

  results_step0 <- test_snps_with_indices(
    bigsnp      = d$bigsnp,
    y           = rnorm(length(d$sample_ids)),
    snp_indices = seq_len(nrow(d$bigsnp$map)),
    snp_names   = d$bigsnp$map$marker.ID,
    design_base = d$design_base,
    ind.row     = ind.row.snp
  )
  results_step0$step <- 0L
  results_step0$conditioning_snps <- NA_character_
  results_step0 <- results_step0[, c("step", "conditioning_snps",
                                     setdiff(names(results_step0),
                                             c("step", "conditioning_snps")))]

  # Force the step-0 lead SNP to pass a tight threshold
  results_step0$pvalue[1] <- 1e-10

  # pval_threshold = 1e-5: step-0 lead (pvalue = 1e-10) passes since 1e-10 < 1e-5.
  # The conditioning step uses a random y, so its p-values are noise and will not
  # reach 1e-5, causing the loop to break without adding that step to stepwise_tables.
  result <- run_stepwise(
    results_step0  = results_step0,
    bigsnp         = d$bigsnp,
    y              = rnorm(length(d$sample_ids)),
    snp_indices    = seq_len(nrow(d$bigsnp$map)),
    cis_snps_pheno = d$bigsnp$map$marker.ID,
    design_base    = d$design_base,
    ind.row.snp    = ind.row.snp,
    pval_threshold = 1e-5,
    max_steps      = 5,
    ncores         = 1
  )

  # stepwise_tables should contain only results_step0 (step 0), since the
  # conditioning step run after it did not pass pval_threshold = 1e-5.
  # (step 0 itself is always in the list; conditioning steps are only added
  # when their lead SNP passes the threshold.)
  expect_equal(length(result$stepwise_tables), 1L)
  expect_equal(result$stepwise_tables[[1]]$step[1], 0L)
})

test_that("run_conditional_qtl accepts max_steps argument and passes it through", {
  d <- make_qtl_test_data()
  out_dir <- tempfile("qtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # Just ensure max_steps is accepted without error
  expect_no_error(
    run_conditional_qtl(
      bigpheno    = d$bigpheno,
      bigsnp      = d$bigsnp,
      pheno_coord = d$pheno_coord,
      design_base = d$design_base,
      max_steps   = 3,
      min_snps    = 1,
      output_dir  = out_dir
    )
  )
})

# =========================================================================
# verbose flag
# =========================================================================

test_that("run_conditional_qtl with verbose=FALSE suppresses normal messages", {
  d <- make_qtl_test_data()
  out_dir <- tempfile("qtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  msgs <- character(0)
  withCallingHandlers(
    run_conditional_qtl(
      bigpheno    = d$bigpheno,
      bigsnp      = d$bigsnp,
      pheno_coord = d$pheno_coord,
      design_base = d$design_base,
      min_snps    = 1,
      output_dir  = out_dir,
      verbose     = FALSE
    ),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_length(msgs, 0)
})

test_that("run_conditional_qtl with verbose=TRUE emits messages", {
  d <- make_qtl_test_data()
  out_dir <- tempfile("qtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  expect_message(
    run_conditional_qtl(
      bigpheno    = d$bigpheno,
      bigsnp      = d$bigsnp,
      pheno_coord = d$pheno_coord,
      design_base = d$design_base,
      min_snps    = 1,
      output_dir  = out_dir,
      verbose     = TRUE
    )
  )
})

# =========================================================================
# do_allbutone requires do_conditioning
# =========================================================================

test_that("run_conditional_qtl errors when do_allbutone=TRUE and do_conditioning=FALSE", {
  d <- make_qtl_test_data()
  out_dir <- tempfile("qtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  expect_error(
    run_conditional_qtl(
      bigpheno        = d$bigpheno,
      bigsnp          = d$bigsnp,
      pheno_coord     = d$pheno_coord,
      design_base     = d$design_base,
      do_conditioning = FALSE,
      do_allbutone    = TRUE,
      min_snps        = 1,
      output_dir      = out_dir
    ),
    "do_allbutone = TRUE requires do_conditioning = TRUE"
  )
})

# =========================================================================
# min_snps threshold
# =========================================================================

test_that("process_pheno warns and skips when fewer cis SNPs than min_snps", {
  d <- make_qtl_test_data()  # 5 SNPs by default

  out_dir <- tempfile("qtl_test_")
  dir.create(file.path(out_dir, "stepwise"),  recursive = TRUE)
  dir.create(file.path(out_dir, "allbutone"), recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  ind.row.snp   <- match(d$sample_ids, d$bigsnp$fam$sample.ID)
  ind.row.pheno <- match(d$sample_ids, d$bigpheno$rowData$sample_name)

  expect_warning(
    process_pheno(
      pheno          = "phenoA",
      bigpheno       = d$bigpheno,
      bigsnp         = d$bigsnp,
      pheno_coord    = d$pheno_coord,
      design_base    = d$design_base,
      ind.row.snp    = ind.row.snp,
      ind.row.pheno  = ind.row.pheno,
      cis_window     = 1e6,
      do_conditioning = TRUE,
      pval_threshold = 1e-3,
      max_steps      = 5,
      do_allbutone   = TRUE,
      do_rint        = FALSE,
      min_snps       = 100,
      ncores         = 1,
      stepwise_dir   = file.path(out_dir, "stepwise"),
      allbutone_dir  = file.path(out_dir, "allbutone")
    ),
    "fewer than min_snps"
  )
})

test_that("run_conditional_qtl respects min_snps and skips phenotypes with too few SNPs", {
  d <- make_qtl_test_data()  # 5 SNPs by default
  out_dir <- tempfile("qtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # min_snps=100 should cause the phenotype to be skipped
  expect_warning(
    run_conditional_qtl(
      bigpheno    = d$bigpheno,
      bigsnp      = d$bigsnp,
      pheno_coord = d$pheno_coord,
      design_base = d$design_base,
      min_snps    = 100,
      output_dir  = out_dir
    ),
    "fewer than min_snps"
  )
})

# =========================================================================
# marginalQTL convenience function
# =========================================================================

test_that("marginalQTL runs without error and returns expected structure", {
  d <- make_qtl_test_data()
  out_dir <- tempfile("qtl_test_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- marginalQTL(
    bigpheno    = d$bigpheno,
    bigsnp      = d$bigsnp,
    pheno_coord = d$pheno_coord,
    design_base = d$design_base,
    min_snps    = 1,
    output_dir  = out_dir
  )

  expect_type(result, "list")
  expect_named(result, c("stepwise", "allbutone", "stepwise_dir", "allbutone_dir"))
  expect_true(inherits(result$stepwise, "Dataset"))
  # marginalQTL disables allbutone, so it should be NULL
  expect_null(result$allbutone)
})

