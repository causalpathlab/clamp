# test_that("compute the evidence of lower bound (ELBO)", {
#   set.seed(12345)
#   p <- 10
#   n <- 1000
#   X <- matrix(rnorm(n * p), nrow = n)
#   y <- 0.5*X[,1] + 0.8*X[,3] + 0.2*X[,7] + rnorm(n)
#
#   out = compute_colstats(X, center = TRUE, scale = TRUE)
#   attr(X,"scaled:center") = out$cm
#   attr(X,"scaled:scale") = out$csd
#   # attr(X,"d") = out$d
#
#   # Center and scale input.
#   if (TRUE)
#     y = y -  mean(y)
#
#   s <- init_setup(n, p, maxL=10,
#                   scaled_prior_variance=0.2,
#                   residual_variance=1,
#                   prior_inclusion_prob=NULL,
#                   as.numeric(var(y)),
#                   standardize=TRUE)
#   s <- init_finalize(s)
#
#   # update once
#   s <- update_each_effect_v2(X=X, y=y, s=s,
#                              robust_method="none")
#
# })
