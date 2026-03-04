test_that("bigPheno creates object with correct class", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(m) <- paste0("s", 1:3)
  colnames(m) <- paste0("g", 1:4)
  bp <- bigPheno(m)
  expect_s3_class(bp, "bigPheno")
})

test_that("bigPheno has expected list elements", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(m) <- paste0("s", 1:3)
  colnames(m) <- paste0("g", 1:4)
  bp <- bigPheno(m)
  expect_named(bp, c("pheno", "rowData", "colData"))
})

test_that("bigPheno rowData has sample_name column when not supplied", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  rownames(m) <- c("sampleA", "sampleB")
  colnames(m) <- c("geneX", "geneY", "geneZ")
  bp <- bigPheno(m)
  expect_true("sample_name" %in% colnames(bp$rowData))
  expect_equal(bp$rowData$sample_name, c("sampleA", "sampleB"))
})

test_that("bigPheno colData has pheno_name column when not supplied", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  colnames(m) <- c("geneX", "geneY", "geneZ")
  bp <- bigPheno(m)
  expect_true("pheno_name" %in% colnames(bp$colData))
  expect_equal(bp$colData$pheno_name, c("geneX", "geneY", "geneZ"))
})

test_that("bigPheno generates default sample names when rownames are NULL", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  colnames(m) <- c("g1", "g2", "g3")
  bp <- bigPheno(m)
  expect_equal(bp$rowData$sample_name, c("sample_1", "sample_2"))
})

test_that("bigPheno generates default pheno names when colnames are NULL", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  bp <- bigPheno(m)
  expect_equal(bp$colData$pheno_name, c("pheno_1", "pheno_2", "pheno_3"))
})

test_that("bigPheno errors on non-matrix input", {
  expect_error(bigPheno(data.frame(a = 1:3)), "matrix must be a numeric matrix")
})

test_that("bigPheno errors when rowData has wrong number of rows", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  bad_rowData <- data.frame(sample_name = c("s1", "s2", "s3"))
  expect_error(bigPheno(m, rowData = bad_rowData), "rowData has")
})

test_that("bigPheno errors when colData has wrong number of rows", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  bad_colData <- data.frame(pheno_name = c("g1", "g2"))
  expect_error(bigPheno(m, colData = bad_colData), "colData has")
})

test_that("bigPheno errors when colData lacks pheno_name column", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  bad_colData <- data.frame(name = c("g1", "g2", "g3"))
  expect_error(bigPheno(m, colData = bad_colData), "pheno_name")
})

