context("twodlda")

# Helper: create a simple 2-class tensor dataset
# Class 1: random matrices with positive mean shift
# Class 2: random matrices with negative mean shift
make_twodlda_tensor <- function(I = 8, J = 6, n_per_class = 5, seed = 42) {
  set.seed(seed)
  K <- 2 * n_per_class
  arr <- array(0, dim = c(I, J, K))
  # class 1: mean = +1
  for (k in 1:n_per_class) {
    arr[, , k] <- matrix(rnorm(I * J, mean = 1, sd = 0.5), I, J)
  }
  # class 2: mean = -1
  for (k in (n_per_class + 1):K) {
    arr[, , k] <- matrix(rnorm(I * J, mean = -1, sd = 0.5), I, J)
  }
  labels <- rep(c(1, 2), each = n_per_class)
  list(tnsr = as.tensor(arr), labels = labels, I = I, J = J, K = K)
}

test_that("twodlda requires Tensor input", {
  expect_error(twodlda(array(1:24, c(3, 4, 2)), labels = c(1, 2)))
})

test_that("twodlda rejects zero tensor", {
  tnsr <- as.tensor(array(0, dim = c(3, 4, 2)))
  expect_error(twodlda(tnsr, labels = c(1, 2)), "Zero tensor")
})

test_that("twodlda requires labels", {
  tnsr <- rand_tensor(c(4, 5, 6))
  expect_error(twodlda(tnsr), "labels")
})

test_that("twodlda requires 3-mode tensor", {
  tnsr <- rand_tensor(c(3, 3, 3, 3))
  expect_error(twodlda(tnsr, labels = c(1, 2, 3, 4)))
})

test_that("twodlda requires labels length to match mode 3", {
  dat <- make_twodlda_tensor()
  expect_error(twodlda(dat$tnsr, labels = c(1, 2)), "labels")
})

test_that("twodlda requires at least 2 classes", {
  dat <- make_twodlda_tensor()
  expect_error(twodlda(dat$tnsr, labels = rep(1, dat$K)))
})

test_that("twodlda returns correct structure", {
  dat <- make_twodlda_tensor()
  res <- twodlda(dat$tnsr, labels = dat$labels, r_ranks = 3, c_ranks = 2)
  expect_true(is.list(res))
  expect_true(all(c("L", "R", "Z", "conv", "est", "norm_percent",
                     "fnorm_resid", "all_resids") %in% names(res)))
  # L is I x r_ranks
  expect_true(is.matrix(res$L))
  expect_equal(nrow(res$L), dat$I)
  expect_equal(ncol(res$L), 3)
  # R is J x c_ranks
  expect_true(is.matrix(res$R))
  expect_equal(nrow(res$R), dat$J)
  expect_equal(ncol(res$R), 2)
  # Z is list of K matrices
  expect_true(is.list(res$Z))
  expect_equal(length(res$Z), dat$K)
  expect_equal(dim(res$Z[[1]]), c(3, 2))
  # est is Tensor
  expect_true(is(res$est, "Tensor"))
  expect_equal(res$est@modes, dat$tnsr@modes)
  # scalars
  expect_true(is.numeric(res$norm_percent))
  expect_true(is.numeric(res$fnorm_resid))
  expect_true(is.logical(res$conv))
  expect_true(is.numeric(res$all_resids))
})

test_that("twodlda projections separate classes", {
  dat <- make_twodlda_tensor()
  res <- twodlda(dat$tnsr, labels = dat$labels, r_ranks = 2, c_ranks = 2,
                 max_iter = 50)
  # Project each slice and compute class centroids
  n_per_class <- dat$K / 2
  proj <- lapply(1:dat$K, function(k) {
    as.vector(t(res$L) %*% dat$tnsr@data[, , k] %*% res$R)
  })
  centroid1 <- Reduce("+", proj[1:n_per_class]) / n_per_class
  centroid2 <- Reduce("+", proj[(n_per_class + 1):dat$K]) / n_per_class
  # Centroids should be separated
  expect_gt(sqrt(sum((centroid1 - centroid2)^2)), 0.1)
})

test_that("twodlda norm_percent is consistent with fnorm_resid", {
  dat <- make_twodlda_tensor()
  res <- twodlda(dat$tnsr, labels = dat$labels, r_ranks = 3, c_ranks = 2)
  expected_pct <- (1 - res$fnorm_resid / fnorm(dat$tnsr)) * 100
  expect_equal(res$norm_percent, expected_pct, tolerance = 1e-10)
})
