test_that("A linear regression model on clamp() with W=NULL", {
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

test_that("A Poisson regression model on clamp() with W=NULL", {
  source.dir <- "~/Dropbox/paper-causal.susie/codes/clamp/R"
  file.sources <- list.files(source.dir, full.names = T, pattern="*.R")
  sapply(file.sources, source, .GlobalEnv)

  set.seed(123456789)
  nn <- 1000
  pp <- 5

  X <- matrix(rnorm(nn * pp), ncol = pp)
  effect_idx <- 2
  bb <- 2
  # data.frame(variable = paste0("X", effect_idx), effect_size = bb)
  h2 <- 1

  Eta <- sqrt(h2) * (X[ ,effect_idx, drop = F] %*% as.matrix(bb)) +
    sqrt(1-h2) * rnorm(nn)
  y1 <- rpois(nn, exp(Eta))

  # res_cl <- clamp(X, y1, family = "poisson", standardize = F,
  #                 estimate_residual_variance = F, estimate_prior_variance = F)
  # res_cl <- clamp(X, y1, family = "poisson", standardize = F,
  #                 estimate_residual_variance = F)
  # res_cl <- clamp(X, y1, family = "poisson", standardize = F,
  #                 estimate_residual_variance = F, robust_method = "huber")

  W <- matrix(2, nrow = nn, ncol = pp)

  res_cl <- clamp(X, y1, W = W,
                  family = "poisson", standardize = F,
                  estimate_residual_variance = F, robust_method = "huber")
})
