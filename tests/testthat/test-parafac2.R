context("parafac2")

# Helper: construct a known PARAFAC2 tensor
# X_k = H_k %*% B %*% diag(c_k) %*% t(C_common_cols)
# with H_k^T H_k = Phi (constant across k)
# Simplified: all slices same size, X_k = H_k B D_k C^T
make_parafac2_tensor <- function(I = 8, J = 6, K = 4, R = 3, seed = 42) {
  set.seed(seed)
  # Shared profile matrix
  B <- matrix(rnorm(R * R), R, R)
  # Shared C matrix (J x R)
  C <- qr.Q(qr(matrix(rnorm(J * R), J, R)))
  # Per-slice orthogonal H_k (I x R) satisfying H_k^T H_k = Phi
  Phi <- diag(R)  # simplest case
  H_list <- lapply(seq_len(K), function(k) {
    qr.Q(qr(matrix(rnorm(I * R), I, R)))
  })
  # Per-slice weights
  c_vec <- lapply(seq_len(K), function(k) runif(R, 0.5, 2.0))
  arr <- array(0, dim = c(I, J, K))
  for (k in seq_len(K)) {
    arr[, , k] <- H_list[[k]] %*% B %*% diag(c_vec[[k]]) %*% t(C)
  }
  list(tnsr = as.tensor(arr), H = H_list, B = B, C = C, c_vec = c_vec,
       I = I, J = J, K = K, R = R)
}

test_that("parafac2 requires Tensor input", {
  expect_error(parafac2(array(1:24, c(3, 4, 2)), num_components = 2))
})

test_that("parafac2 rejects zero tensor", {
  tnsr <- as.tensor(array(0, dim = c(3, 4, 2)))
  expect_error(parafac2(tnsr, num_components = 2), "Zero tensor")
})

test_that("parafac2 requires num_components", {
  tnsr <- rand_tensor(c(4, 5, 3))
  expect_error(parafac2(tnsr), "num_components")
})

test_that("parafac2 requires 3-mode tensor", {
  tnsr <- rand_tensor(c(3, 3, 3, 3))
  expect_error(parafac2(tnsr, num_components = 2))
})

test_that("parafac2 returns correct structure", {
  dat <- make_parafac2_tensor()
  res <- parafac2(dat$tnsr, num_components = dat$R, max_iter = 50)
  expect(is.list(res), "result should be a list")
  expect(all(c("H", "B", "C", "D", "conv", "est", "norm_percent",
               "fnorm_resid", "all_resids") %in% names(res)),
         "result missing expected names")
  # H is list of K matrices (I x R)
  expect(is.list(res$H), "H should be a list")
  expect_equal(length(res$H), dat$K)
  expect_equal(nrow(res$H[[1]]), dat$I)
  expect_equal(ncol(res$H[[1]]), dat$R)
  # B is R x R matrix
  expect(is.matrix(res$B), "B should be a matrix")
  expect_equal(dim(res$B), c(dat$R, dat$R))
  # C is J x R matrix
  expect(is.matrix(res$C), "C should be a matrix")
  expect_equal(nrow(res$C), dat$J)
  expect_equal(ncol(res$C), dat$R)
  # D is list of K diagonal matrices (R x R)
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

test_that("parafac2 H_k satisfy cross-product constraint", {
  dat <- make_parafac2_tensor()
  res <- parafac2(dat$tnsr, num_components = dat$R, max_iter = 100, tol = 1e-6)
  # H_k^T H_k should be approximately equal across all k
  cross_prods <- lapply(res$H, function(Hk) t(Hk) %*% Hk)
  for (k in 2:dat$K) {
    expect_equal(cross_prods[[1]], cross_prods[[k]], tolerance = 1e-4)
  }
})

test_that("parafac2 recovers low-rank tensor accurately", {
  dat <- make_parafac2_tensor()
  set.seed(123)
  res <- parafac2(dat$tnsr, num_components = dat$R, max_iter = 200, tol = 1e-8)
  expect_gt(res$norm_percent, 95)
  expect_lt(res$fnorm_resid / fnorm(dat$tnsr), 0.05)
})

test_that("parafac2 norm_percent is consistent with fnorm_resid", {
  dat <- make_parafac2_tensor()
  res <- parafac2(dat$tnsr, num_components = dat$R, max_iter = 50)
  expected_pct <- (1 - res$fnorm_resid / fnorm(dat$tnsr)) * 100
  expect_equal(res$norm_percent, expected_pct, tolerance = 1e-10)
})

test_that("parafac2 all_resids eventually decreases", {
  dat <- make_parafac2_tensor()
  res <- parafac2(dat$tnsr, num_components = dat$R, max_iter = 50)
  resids <- res$all_resids
  if (length(resids) > 2) {
    expect_lt(tail(resids, 1), resids[min(3, length(resids))])
  }
})
