## source functions
.libPaths("~/projects/Rlibs/")
library(glmnet)
library(susieR)
library(data.table)

`%&%` <- function(a, b) paste0(a, b)
if.needed <- function(.files, .code) {
  if(!all(file.exists(.files))){
    .code
  }
}

# getwd()
files_path <- "./R/"
r_file_names <- list.files(path = files_path)

for (i in 1 : length(r_file_names)) {
  source(files_path %&% r_file_names[i])
}

##################
# create filename of files that save the results.

# n: sample size
# p: number of variables
# h2: heritability
# seed: random seed
# .prefix: add identities such as date/"sim/real"
create_filename <- function(n, p, h2, seed, .prefix = NULL) {

  .filename <- "n" %&% n %&% "-p" %&% p %&% "-h" %&% (100*h2) %&%
    "-s" %&% seed %&% ".rds"

  if (!is.null(.prefix)) {
    .filename <- .prefix %&% "-" %&% .filename
  }
  return(.filename)
}


###################
# simulations

## read in arguments
args <- commandArgs(trailingOnly = TRUE)
ID <- as.numeric(args[1])
seed <- as.numeric(args[2])
nn <- as.numeric(args[3])
pp <- as.numeric(args[4])
h2 <- as.numeric(args[5])
n_effect_vars <- as.numeric(args[6])  ## number of effective variables

.result.dir <- "~/projects/causal-susie/results/" %&% Sys.Date() %&% "-synthetic/"
# .result.dir <- "./results/synthetic-" %&% Sys.Date() %&% "-" %&% model %&% "/"
if (!dir.exists(.result.dir)) {
  dir.create(.result.dir, recursive=TRUE, showWarnings=FALSE)
}

.filename <- create_filename(nn, pp, h2, seed, "ID"%&%ID)

if.needed(.result.dir %&% .filename, {

  # methods names
  ## name the columns (methods):
  methods_names <- c("cl.pip", "su_vn.pip", "su_ad.pip", "la.coef", "en.coef")

  expit <- function(eta) {
    return ( ifelse( eta > 0, 1 / (1+exp(-eta)), exp(eta) / (1+exp(eta)) ) )
  }

  ##################
  ## Data generation
  set.seed(seed)

  U <- rnorm(nn)  # confounder
  zeta <- ( (1:pp) - pp/2 ) / (pp+1)
  # zeta <- rnorm(pp, sd = 1)
  # hist(zeta)
  delta <- matrix(rnorm(nn*pp, 0, 0.1), nrow=nn)
  UU <- (outer(U, zeta) + delta)[, , drop=F]  # confounding effect on X.
  Pmat <- expit(UU)  # Prob matrix

  X <- sapply(1:pp, function(j) {rbinom(nn, 1, Pmat[,j])})
  X <- apply(X, 2, as.double)
  esp <- rnorm(nn)

  causal_vars <- sample.int(pp, size = 5)
  # coefs <- rnorm(length(causal_vars) + 1)
  coefs <- append(rep(1, times = length(causal_vars)), 1)  # all-one
  y <- sqrt(h2) * cbind(X[, causal_vars, drop=F], U) %*% as.matrix(coefs) + sqrt(1-h2) * esp

  # print(data.frame(variable = c(paste0("X", causal_vars), "U"),
  #                  effect_size = coefs))

  ##################
  ## Methods

  # Method 1: causal-susie
  ## step 1: estimate the propensity scores and construct the weight matrix
  Uhat <- rsvd(X, k=3)$u  ## 3 PCs as the substitute confounders.
  PS <- apply(X, 2,
              function(xcol) predict(glm(xcol ~ Uhat, family = binomial),
                                     type = "response"))
  Wmat <- ifelse(X == 1, 1/PS, 1/(1-PS))
  ## step 2: fit a causal-susie(?) model
  res_cl <- tryCatch(
    withCallingHandlers(
      clamp(X, y, Wmat, seed=seed, standardize=F, max_iter = 100),
      error = function(e) {write.to.log(sys.calls())},
      warning = function(w) {
        write.to.log(sys.calls())
        invokeRestart("muffleWarning")
        }),
    error = function(e) {
      message("clamp. seed: ", seed, "; error: ", e)
    }
    )

  # Method 2: susie (vanilla)
  res_su_vn <- susie(X, y)

  # Method 3: susie + regressing on residuals from (y ~ Uhat)
  yres <- residuals(lm(y ~ Uhat))
  res_su_ad <- susie(X, yres)

  # Method 4: lasso
  res_la <- cv.glmnet(X, y, family = "gaussian", alpha = 1)

  # Method 5: elastic net
  res_en <- cv.glmnet(X, y, family = "gaussian", alpha = 0.5)

  ##################
  ## Evaluation

  extract_pip <- function(.res) {
    if (is.null(.res$pip)) {
      # causal-susie raised an error
      pip <- rep(NA, times = pp)
    } else {
      pip <- .res$pip
      names(pip) <- names(.res$pip)
      # no intercept is included.
    }

    return(pip)
  }

  output <- list()

  # save setup
  output$setup <- list(
    seed = seed,
    n = nn,
    p = pp,
    h2 = h2,
    n_effect_vars = n_effect_vars,
    causal_vars = causal_vars
  )

  # save results
  output$selected <- data.table(
    extract_pip(res_cl),
    extract_pip(res_su_vn),
    extract_pip(res_su_ad),
    as.numeric(coef(res_la, s = "lambda.min"))[-1],
    as.numeric(coef(res_en, s = "lambda.min"))[-1]
  )
  colnames(output$selected) <- methods_names

  # saveRDS
  saveRDS(output, file = .result.dir %&% .filename)
})







