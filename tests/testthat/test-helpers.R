test_that("get_feature_indices returns correct indices", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  colnames(m) <- c("GeneA", "GeneB", "GeneC", "GeneD")
  bf <- bigFeatures(m)

  idx <- get_feature_indices(bf, c("GeneC", "GeneA"))
  expect_equal(idx, c(3L, 1L))
})

test_that("get_feature_indices errors on missing feature", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  colnames(m) <- c("GeneA", "GeneB", "GeneC")
  bf <- bigFeatures(m)

  expect_error(get_feature_indices(bf, "GeneX"), "Could not find")
})

test_that("get_cis_snps returns empty lists when no SNPs in window", {
  snp_map <- data.frame(
    chromosome   = c("1", "1", "2"),
    marker.ID    = c("rs1", "rs2", "rs3"),
    physical.pos = c(1000L, 2000L, 5000L),
    stringsAsFactors = FALSE
  )
  # Minimal bigsnp-like list
  bigsnp <- list(map = snp_map)

  result <- get_cis_snps(bigsnp, gene_chr = "3",
                         gene_start = 1000, gene_end = 2000,
                         cis_window = 500)
  expect_length(result$indices, 0)
  expect_length(result$names, 0)
})

test_that("get_cis_snps finds SNPs within window on correct chromosome", {
  snp_map <- data.frame(
    chromosome   = c("1", "1", "1", "2"),
    marker.ID    = c("rs1", "rs2", "rs3", "rs4"),
    physical.pos = c(900L, 1500L, 3200L, 1500L),
    stringsAsFactors = FALSE
  )
  bigsnp <- list(map = snp_map)

  # Gene on chr1, 1000-2000, window=500 -> [500, 2500]
  result <- get_cis_snps(bigsnp, gene_chr = "1",
                         gene_start = 1000, gene_end = 2000,
                         cis_window = 500)
  expect_equal(result$indices, c(1L, 2L))
  expect_equal(result$names, c("rs1", "rs2"))
})

test_that("get_cis_snps window boundary is inclusive", {
  snp_map <- data.frame(
    chromosome   = "1",
    marker.ID    = "rs_exact",
    physical.pos = 500L,
    stringsAsFactors = FALSE
  )
  bigsnp <- list(map = snp_map)

  # Gene 1000-2000, window=500 -> lower bound = 500 exactly
  result <- get_cis_snps(bigsnp, gene_chr = "1",
                         gene_start = 1000, gene_end = 2000,
                         cis_window = 500)
  expect_equal(result$names, "rs_exact")
})
