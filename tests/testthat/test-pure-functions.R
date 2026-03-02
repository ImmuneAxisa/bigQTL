test_that("rint returns a vector of the same length", {
  x <- c(1.5, 3.2, 2.1, 0.4, 5.0)
  out <- rint(x)
  expect_length(out, length(x))
})

test_that("rint output is approximately standard normal", {
  set.seed(42)
  x <- rnorm(200)
  out <- rint(x)
  # Mean should be near 0, sd near 1
  expect_lt(abs(mean(out)), 0.1)
  expect_lt(abs(sd(out) - 1), 0.1)
})

test_that("rint preserves rank ordering", {
  x <- c(3, 1, 4, 1, 5, 9, 2, 6)
  out <- rint(x)
  expect_equal(rank(x, ties.method = "average"), rank(out, ties.method = "average"))
})

test_that("rint handles ties without error", {
  x <- c(1, 1, 2, 2, 3)
  expect_no_error(rint(x))
})

test_that("count_eigenvalues returns 1 for a single eigenvalue", {
  expect_equal(count_eigenvalues(5, 0.99), 1L)
})

test_that("count_eigenvalues returns all for zero total variance", {
  eigs <- c(0, 0, 0)
  expect_equal(count_eigenvalues(eigs, 0.99), length(eigs))
})

test_that("count_eigenvalues returns correct count at threshold", {
  eigs <- c(30, 10, 8, 5, 1, 0.5)
  expect_equal(count_eigenvalues(eigs, 0.99), 5L)
})

test_that("count_eigenvalues sorts descending before counting", {
  eigs_unsorted <- c(1, 0.5, 10, 5, 2)
  eigs_sorted   <- c(10, 5, 2, 1, 0.5)
  expect_equal(count_eigenvalues(eigs_unsorted, 0.99),
               count_eigenvalues(eigs_sorted, 0.99))
})

test_that("eigenMT_correct caps at 1", {
  expect_equal(eigenMT_correct(0.5, 5), 1)
  expect_equal(eigenMT_correct(0.1, 5), 0.5)
})

test_that("eigenMT_correct returns pmin(p * m_eff, 1)", {
  pvals  <- c(0.01, 0.05, 0.1, 0.5)
  m_eff  <- 10L
  result <- eigenMT_correct(pvals, m_eff)
  expect_equal(result, pmin(pvals * m_eff, 1))
})

test_that("lw_shrink_cor returns identity for 1x1 matrix", {
  R <- matrix(1, 1, 1)
  out <- lw_shrink_cor(R, n = 100)
  expect_equal(out, matrix(1, 1, 1))
})

test_that("lw_shrink_cor diagonal is all 1", {
  set.seed(42)
  geno <- matrix(sample(0:2, 50 * 20, replace = TRUE), 50, 20)
  R <- cor(geno)
  out <- lw_shrink_cor(R, n = 50)
  expect_true(all(abs(diag(out) - 1) < 1e-10))
})

test_that("lw_shrink_cor output is positive semi-definite", {
  set.seed(42)
  geno <- matrix(sample(0:2, 50 * 20, replace = TRUE), 50, 20)
  R <- cor(geno)
  out <- lw_shrink_cor(R, n = 50)
  eigs <- eigen(out, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-10))
})

test_that("lw_shrink_cor handles perfect LD (all-ones matrix)", {
  R <- matrix(1, 5, 5)
  out <- lw_shrink_cor(R, n = 100)
  expect_equal(out, matrix(1, 5, 5))
})
