test_that("update each effect of the weighted SER (without W)", {
  set.seed(12345)
  p <- 10
  n <- 1000
  X <- matrix(rnorm(n * p), nrow = n)
  y <- 0.5*X[,1] + 0.8*X[,3] + 0.2*X[,7] + rnorm(n)

  # X <- scale(X)
  out = compute_colstats(X,center = TRUE,scale = TRUE)
  attr(X,"scaled:center") = out$scaled_center
  attr(X,"scaled:scale") = out$scaled_scale
  attr(X,"d") = out$d

  # Center and scale input.
  if (TRUE)
    y = y -  mean(y)

  s <- init_setup(n, p, maxL=10,
                  family="linear",
                  scaled_prior_variance=0.2,
                  residual_variance=1,
                  prior_inclusion_prob=NULL,
                  null_weight=0,
                  varY=as.numeric(var(y)),
                  standardize=TRUE)
  s <- init_finalize(s)

  res <- update_each_effect(X=X, y=y, s=s,
                            estimate_prior_method = "none")

  # expect_equal(sum(res$logBF), 486.4492, tolerance = 0.001)
})


test_that("update each effect of the weighted SER (with unified W)", {
  set.seed(12345)
  p <- 10
  n <- 1000
  X <- matrix(rnorm(n * p), nrow = n)
  y <- 0.5*X[,1] + 0.8*X[,3] + 0.2*X[,7] + rnorm(n)

  # X <- scale(X)
  out = compute_colstats(X,center = T, scale = T)
  attr(X,"scaled:center") = out$scaled_center
  attr(X,"scaled:scale") = out$scaled_scale
  attr(X,"d") = out$d

  # Center and scale input.
  if (TRUE)
    y_ = y -  mean(y)

  s <- init_setup(n, p, maxL=10,
                  family="linear",
                  scaled_prior_variance=0.2,
                  residual_variance=1,
                  prior_inclusion_prob=NULL,
                  null_weight=0,
                  varY=as.numeric(var(y)),
                  standardize=TRUE)
  s <- init_finalize(s)

  # weight matrix W
  W <- matrix(1, nrow = n, ncol = p)

  res <- update_each_effect(X=X, y=y_, s=s, W=W,
                            estimate_prior_method="none")

  # expect_equal(sum(res$logBF), 26.32327, tolerance = 0.1)  ##?!
})

##########
## test the equivalence of the `update_each_effect()` in `susieR`
