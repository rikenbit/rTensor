context("rescal")

# Helper: construct a known low-rank RESCAL tensor
# X_k = A %*% R_k %*% t(A),  k = 1,...,K
make_rescal_tensor <- function(N = 8, R = 3, K = 4, seed = 42) {
  set.seed(seed)
  A <- qr.Q(qr(matrix(rnorm(N * R), N, R)))
  R_list <- lapply(seq_len(K), function(k) matrix(rnorm(R * R), R, R))
  arr <- array(0, dim = c(N, N, K))
  for (k in seq_len(K)) {
    arr[, , k] <- A %*% R_list[[k]] %*% t(A)
  }
  list(tnsr = as.tensor(arr), A = A, R = R_list, N = N, R_comp = R, K = K)
}

test_that("rescal requires Tensor input", {
  expect_error(rescal(array(1:24, c(3, 4, 2)), num_components = 2))
})

test_that("rescal rejects zero tensor", {
  tnsr <- as.tensor(array(0, dim = c(3, 3, 2)))
  expect_error(rescal(tnsr, num_components = 2), "Zero tensor")
})

test_that("rescal requires num_components", {
  tnsr <- rand_tensor(c(4, 4, 3))
  expect_error(rescal(tnsr), "num_components")
})

test_that("rescal requires 3-mode tensor", {
  tnsr <- rand_tensor(c(3, 3, 3, 3))
  expect_error(rescal(tnsr, num_components = 2))
})

test_that("rescal requires square frontal slices", {
  tnsr <- rand_tensor(c(4, 5, 3))
  expect_error(rescal(tnsr, num_components = 2))
})

test_that("rescal returns correct structure", {
  dat <- make_rescal_tensor()
  res <- rescal(dat$tnsr, num_components = dat$R_comp, max_iter = 50)
  expect_true(is.list(res))
  expect_true(all(c("A", "R", "conv", "est", "norm_percent",
                     "fnorm_resid", "all_resids") %in% names(res)))
  # A is N x R matrix
  expect_true(is.matrix(res$A))
  expect_equal(nrow(res$A), dat$N)
  expect_equal(ncol(res$A), dat$R_comp)
  # R is list of K dense matrices (R x R)
  expect_true(is.list(res$R))
  expect_equal(length(res$R), dat$K)
  expect_equal(dim(res$R[[1]]), c(dat$R_comp, dat$R_comp))
  # est is Tensor
  expect_true(is(res$est, "Tensor"))
  expect_equal(res$est@modes, dat$tnsr@modes)
  # scalars
  expect_true(is.numeric(res$norm_percent))
  expect_true(is.numeric(res$fnorm_resid))
  expect_true(is.logical(res$conv))
  expect_true(is.numeric(res$all_resids))
})

test_that("rescal recovers low-rank tensor accurately", {
  dat <- make_rescal_tensor()
  res <- rescal(dat$tnsr, num_components = dat$R_comp, max_iter = 100, tol = 1e-8)
  expect_gt(res$norm_percent, 99)
  expect_lt(res$fnorm_resid / fnorm(dat$tnsr), 0.01)
})

test_that("rescal norm_percent is consistent with fnorm_resid", {
  dat <- make_rescal_tensor()
  res <- rescal(dat$tnsr, num_components = dat$R_comp, max_iter = 50)
  expected_pct <- (1 - res$fnorm_resid / fnorm(dat$tnsr)) * 100
  expect_equal(res$norm_percent, expected_pct, tolerance = 1e-10)
})

test_that("rescal all_resids eventually decreases", {
  dat <- make_rescal_tensor()
  res <- rescal(dat$tnsr, num_components = dat$R_comp, max_iter = 50)
  resids <- res$all_resids
  if (length(resids) > 2) {
    expect_lt(tail(resids, 1), resids[min(3, length(resids))])
  }
})
