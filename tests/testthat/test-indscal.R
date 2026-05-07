context("indscal")

# Helper: construct a known low-rank INDSCAL tensor
# X_k = A %*% D_k %*% t(A),  k = 1,...,K
make_indscal_tensor <- function(I = 8, R = 3, K = 4, seed = 42) {
  set.seed(seed)
  A <- qr.Q(qr(matrix(rnorm(I * R), I, R)))  # orthogonal columns
  D_list <- lapply(seq_len(K), function(k) diag(runif(R, 0.5, 2.0)))
  arr <- array(0, dim = c(I, I, K))
  for (k in seq_len(K)) {
    arr[, , k] <- A %*% D_list[[k]] %*% t(A)
  }
  list(tnsr = as.tensor(arr), A = A, D = D_list, I = I, R = R, K = K)
}

test_that("indscal requires Tensor input", {
  expect_error(indscal(array(1:24, c(3, 4, 2)), num_components = 2))
})

test_that("indscal rejects zero tensor", {
  tnsr <- as.tensor(array(0, dim = c(3, 3, 2)))
  expect_error(indscal(tnsr, num_components = 2), "Zero tensor")
})

test_that("indscal requires num_components", {
  tnsr <- rand_tensor(c(4, 4, 3))
  expect_error(indscal(tnsr), "num_components")
})

test_that("indscal requires 3-mode tensor", {
  tnsr <- rand_tensor(c(3, 3, 3, 3))
  expect_error(indscal(tnsr, num_components = 2))
})

test_that("indscal requires square frontal slices", {
  tnsr <- rand_tensor(c(4, 5, 3))
  expect_error(indscal(tnsr, num_components = 2))
})

test_that("indscal returns correct structure", {
  dat <- make_indscal_tensor()
  res <- indscal(dat$tnsr, num_components = dat$R, max_iter = 50)
  expect(is.list(res), "result should be a list")
  expect(all(c("A", "D", "conv", "est", "norm_percent",
               "fnorm_resid", "all_resids") %in% names(res)),
         "result missing expected names")
  # A is I x R matrix
  expect(is.matrix(res$A), "A should be a matrix")
  expect_equal(nrow(res$A), dat$I)
  expect_equal(ncol(res$A), dat$R)
  # D is list of K diagonal matrices
  expect(is.list(res$D), "D should be a list")
  expect_equal(length(res$D), dat$K)
  expect_equal(dim(res$D[[1]]), c(dat$R, dat$R))
  # est is Tensor
  expect(is(res$est, "Tensor"), "est should be a Tensor")
  expect_equal(res$est@modes, dat$tnsr@modes)
  # scalars
  expect(is.numeric(res$norm_percent), "norm_percent should be numeric")
  expect(is.numeric(res$fnorm_resid), "fnorm_resid should be numeric")
  expect(is.logical(res$conv), "conv should be logical")
  expect(is.numeric(res$all_resids), "all_resids should be numeric")
})

test_that("indscal recovers low-rank tensor accurately", {
  dat <- make_indscal_tensor()
  set.seed(123)
  res <- indscal(dat$tnsr, num_components = dat$R, max_iter = 200, tol = 1e-8)
  # Should explain > 95% of norm for exact low-rank data
  expect_gt(res$norm_percent, 95)
  # Frobenius residual should be small
  expect_lt(res$fnorm_resid / fnorm(dat$tnsr), 0.05)
})

test_that("indscal norm_percent is consistent with fnorm_resid", {
  dat <- make_indscal_tensor()
  res <- indscal(dat$tnsr, num_components = dat$R, max_iter = 50)
  expected_pct <- (1 - res$fnorm_resid / fnorm(dat$tnsr)) * 100
  expect_equal(res$norm_percent, expected_pct, tolerance = 1e-10)
})

test_that("indscal all_resids eventually decreases", {
  dat <- make_indscal_tensor()
  res <- indscal(dat$tnsr, num_components = dat$R, max_iter = 50)
  resids <- res$all_resids
  if (length(resids) > 2) {
    # Final residual should be smaller than the residual after warmup
    expect_lt(tail(resids, 1), resids[min(3, length(resids))])
  }
})
