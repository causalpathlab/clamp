test_that("weighted SER with null weights", {
  set.seed(12345)
  p <- 10
  n <- 1000
  X <- matrix(rnorm(n * p), nrow = n)
  y <- 0.5*X[,1] + 0.8*X[,3] + 0.2*X[,7] + rnorm(n)

  out = compute_colstats(X,center = T, scale = T)
  attr(X,"scaled:center") = out$scaled_center
  attr(X,"scaled:scale") = out$scaled_scale
  attr(X, "d") = out$d

  res <- weighted_single_effect_regression(y=y, X=X, prior_varB=1)
  # res <- single_effect_regression(y=y, X=X, V=1)  ## internal function

  # expect_equal(res$logBF_model, 323.6682, tolerance = 0.001)
})

test_that("weighted SER with a unit weight matrix", {
  set.seed(12345)
  p <- 10
  n <- 1000
  X <- matrix(rnorm(n * p), nrow = n)
  y <- 0.5*X[,1] + 0.8*X[,3] + 0.2*X[,7] + rnorm(n)
  W <- matrix(1, nrow = n, ncol = p)

  out = compute_colstats(X,center = T, scale = T)
  attr(X,"scaled:center") = out$scaled_center
  attr(X,"scaled:scale") = out$scaled_scale
  attr(X, "d") = out$d

  res <- weighted_single_effect_regression(y=y, X=X, W=W, prior_varB=1)

  # expect_equal(res$logBF_model, 323.6682, tolerance = 0.001)
})

