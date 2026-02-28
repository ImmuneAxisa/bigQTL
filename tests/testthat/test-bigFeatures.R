test_that("bigFeatures creates object with correct class", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(m) <- paste0("s", 1:3)
  colnames(m) <- paste0("g", 1:4)
  bf <- bigFeatures(m)
  expect_s3_class(bf, "bigFeatures")
})

test_that("bigFeatures has expected list elements", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(m) <- paste0("s", 1:3)
  colnames(m) <- paste0("g", 1:4)
  bf <- bigFeatures(m)
  expect_named(bf, c("features", "rowData", "colData"))
})

test_that("bigFeatures rowData has sample_name column when not supplied", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  rownames(m) <- c("sampleA", "sampleB")
  colnames(m) <- c("geneX", "geneY", "geneZ")
  bf <- bigFeatures(m)
  expect_true("sample_name" %in% colnames(bf$rowData))
  expect_equal(bf$rowData$sample_name, c("sampleA", "sampleB"))
})

test_that("bigFeatures colData has feature_name column when not supplied", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  colnames(m) <- c("geneX", "geneY", "geneZ")
  bf <- bigFeatures(m)
  expect_true("feature_name" %in% colnames(bf$colData))
  expect_equal(bf$colData$feature_name, c("geneX", "geneY", "geneZ"))
})

test_that("bigFeatures generates default sample names when rownames are NULL", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  colnames(m) <- c("g1", "g2", "g3")
  bf <- bigFeatures(m)
  expect_equal(bf$rowData$sample_name, c("sample_1", "sample_2"))
})

test_that("bigFeatures generates default feature names when colnames are NULL", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  bf <- bigFeatures(m)
  expect_equal(bf$colData$feature_name, c("feature_1", "feature_2", "feature_3"))
})

test_that("bigFeatures errors on non-matrix input", {
  expect_error(bigFeatures(data.frame(a = 1:3)), "matrix must be a numeric matrix")
})

test_that("bigFeatures errors when rowData has wrong number of rows", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  bad_rowData <- data.frame(sample_name = c("s1", "s2", "s3"))
  expect_error(bigFeatures(m, rowData = bad_rowData), "rowData has")
})

test_that("bigFeatures errors when colData has wrong number of rows", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  bad_colData <- data.frame(feature_name = c("g1", "g2"))
  expect_error(bigFeatures(m, colData = bad_colData), "colData has")
})

test_that("bigFeatures errors when colData lacks feature_name column", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  bad_colData <- data.frame(name = c("g1", "g2", "g3"))
  expect_error(bigFeatures(m, colData = bad_colData), "feature_name")
})
