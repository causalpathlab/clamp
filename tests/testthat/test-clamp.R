test_that("clamp", {
  set.seed(12345)
  p <- 10
  n <- 1000
  X <- matrix(rnorm(n * p), nrow = n)
  y <- 0.5*X[,1] + 0.8*X[,3] + 0.2*X[,7] + rnorm(n)

  res <- clamp(X, y)
  # summary.clamp(res)
  # This output should be equivalent to
  # res_susie <- susieR::susie(X, y)

  # expect_equal(sum(res$pip), 6.816112, tolerance = 1e-4)
  # expect_equal(res$elbo[res$niter], -1415.958, tolerance = 1e-3)

})
