# run_allbutone() had zero test coverage (see dev/CODE_REVIEW.md section 3):
# every existing fixture used random, uncorrelated genotypes/phenotypes, so
# `length(conditioning_snps) > 1` (process_pheno's gate for calling
# run_allbutone(), bigQTL.R:380) was never true under test and the
# function's body never executed.
#
# It's also the function containing the fragile last-SNP-reuse optimization
# described in dev/CODE_REVIEW.md section 1: stepwise_tables[[length(...)]]
# is reused for the last independent SNP on the assumption that it was
# computed by conditioning on exactly conditioning_snps[-length(...)]. These
# tests exercise both the reuse branch and the freshly-computed branch on
# `dummy_qtl_data`'s pheno1 (3 real signals) and pheno2 (2 real signals),
# which genuinely trigger run_allbutone() the way real data would.

run_stepwise_for_pheno <- function(d, pheno_name, pval_threshold = 1e-3) {
  ind.row <- match(d$sample_ids, d$bigsnp$fam$sample.ID)
  y <- rint(d$bigpheno$pheno[ind.row, get_pheno_indices(d$bigpheno, pheno_name)])
  pr <- d$pheno_coord[d$pheno_coord$pheno_name == pheno_name, ]
  cis <- get_cis_snps(d$bigsnp, pheno_chr = pr$chromosome,
                      pheno_start = pr$start, pheno_end = pr$end, cis_window = 1e6)

  results_step0 <- test_snps_with_indices(
    bigsnp = d$bigsnp, y = y, snp_indices = cis$indices, snp_names = cis$names,
    design_base = d$design_base, ind.row = ind.row, ncores = 1
  )
  results_step0$step <- 0L
  results_step0$conditioning_snps <- NA_character_

  stepwise <- run_stepwise(
    results_step0 = results_step0, bigsnp = d$bigsnp, y = y,
    snp_indices = cis$indices, cis_snps_pheno = cis$names,
    design_base = d$design_base, ind.row.snp = ind.row,
    pval_threshold = pval_threshold, max_steps = 5, ncores = 1, pheno = pheno_name
  )

  list(y = y, cis = cis, ind.row = ind.row, stepwise = stepwise)
}

test_that("run_allbutone is exercised with >= 2 real conditioning SNPs (pheno1, 3 signals)", {
  d <- load_dummy_qtl_data()
  sw <- run_stepwise_for_pheno(d, "pheno1")

  expect_equal(length(sw$stepwise$conditioning_snps), 3L)

  allbutone_tables <- run_allbutone(
    conditioning_snps = sw$stepwise$conditioning_snps,
    stepwise_tables   = sw$stepwise$stepwise_tables,
    bigsnp            = d$bigsnp, y = sw$y,
    snp_indices       = sw$cis$indices, cis_snps_pheno = sw$cis$names,
    design_base       = d$design_base, ind.row.snp = sw$ind.row, ncores = 1
  )

  expect_equal(length(allbutone_tables), 3L)

  for (tbl in allbutone_tables) {
    expect_named(tbl, c("snp", "beta", "se", "t_stat", "pvalue", "fdr",
                        "indep", "conditioning_snps"))
  }

  # Each independent SNP must be recovered as a strong hit when everything
  # else is conditioned out.
  conditioning_snps <- sw$stepwise$conditioning_snps
  for (i in seq_along(conditioning_snps)) {
    tbl <- allbutone_tables[[i]]
    own_row <- tbl[tbl$snp == conditioning_snps[i], ]
    expect_equal(nrow(own_row), 1L)
    expect_lt(own_row$pvalue, 1e-3)
    expect_false(grepl(conditioning_snps[i], tbl$conditioning_snps[1], fixed = TRUE))
  }
})

test_that("run_allbutone's last-SNP reuse matches a freshly-computed all-but-one result", {
  d <- load_dummy_qtl_data()
  sw <- run_stepwise_for_pheno(d, "pheno1")
  conditioning_snps <- sw$stepwise$conditioning_snps
  n <- length(conditioning_snps)

  allbutone_tables <- run_allbutone(
    conditioning_snps = conditioning_snps,
    stepwise_tables   = sw$stepwise$stepwise_tables,
    bigsnp            = d$bigsnp, y = sw$y,
    snp_indices       = sw$cis$indices, cis_snps_pheno = sw$cis$names,
    design_base       = d$design_base, ind.row.snp = sw$ind.row, ncores = 1
  )

  reused <- allbutone_tables[[n]]

  fresh <- test_snps_with_indices(
    bigsnp = d$bigsnp, y = sw$y, snp_indices = sw$cis$indices,
    snp_names = sw$cis$names, design_base = d$design_base, ind.row = sw$ind.row,
    snp_conditioning = conditioning_snps[-n], ncores = 1
  )

  ord_reused <- order(reused$snp)
  ord_fresh  <- order(fresh$snp)

  expect_equal(reused$snp[ord_reused], fresh$snp[ord_fresh])
  expect_equal(reused$beta[ord_reused], fresh$beta[ord_fresh])
  expect_equal(reused$pvalue[ord_reused], fresh$pvalue[ord_fresh])
  expect_equal(reused$conditioning_snps[1], paste(conditioning_snps[-n], collapse = ";"))
})

test_that("run_allbutone is exercised with exactly 2 real conditioning SNPs (pheno2, 2 signals)", {
  d <- load_dummy_qtl_data()
  sw <- run_stepwise_for_pheno(d, "pheno2")

  expect_equal(length(sw$stepwise$conditioning_snps), 2L)

  allbutone_tables <- run_allbutone(
    conditioning_snps = sw$stepwise$conditioning_snps,
    stepwise_tables   = sw$stepwise$stepwise_tables,
    bigsnp            = d$bigsnp, y = sw$y,
    snp_indices       = sw$cis$indices, cis_snps_pheno = sw$cis$names,
    design_base       = d$design_base, ind.row.snp = sw$ind.row, ncores = 1
  )

  expect_equal(length(allbutone_tables), 2L)

  conditioning_snps <- sw$stepwise$conditioning_snps
  for (i in seq_along(conditioning_snps)) {
    tbl <- allbutone_tables[[i]]
    own_row <- tbl[tbl$snp == conditioning_snps[i], ]
    expect_equal(nrow(own_row), 1L)
    expect_lt(own_row$pvalue, 1e-3)
  }

  # The i=1 table is freshly computed, conditioned on conditioning_snps[-1].
  expect_equal(allbutone_tables[[1]]$conditioning_snps[1], conditioning_snps[2])
  # The i=2 table is the reused last-stepwise-step table, conditioned on
  # conditioning_snps[-2].
  expect_equal(allbutone_tables[[2]]$conditioning_snps[1], conditioning_snps[1])
})

test_that("run_allbutone is not triggered for phenotypes with < 2 conditioning SNPs (pheno3, pheno4, pheno5)", {
  d <- load_dummy_qtl_data()

  for (pheno_name in c("pheno3", "pheno4", "pheno5")) {
    sw <- run_stepwise_for_pheno(d, pheno_name)
    expect_lt(length(sw$stepwise$conditioning_snps), 2L)
  }
})
