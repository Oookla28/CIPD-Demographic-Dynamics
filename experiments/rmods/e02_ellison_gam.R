#!/usr/bin/env Rscript
# E02 baseline: Ellison et al. (2026, Demography) Bayesian P-spline GAM + HMC.
# Faithful adaptation of their r/0_functions.r (bbase) and stan/mod_p.stan
# functional form (log-rate = P-spline in age, RW2-type penalty via prior on
# 2nd differences) to period age x parity ASFR matrices with missing cells.
#
# I/O (JSON):
#   in:  {"cases": [{"id":..., "age":[15..49], "y":[...], "y_na":[...bool]}, ...]}
#   out: {"cases": [{"id":..., "pred":[...], "q05":[...], "q95":[...]}, ...]}
#
# Usage: Rscript e02_ellison_gam.R input.json output.json [n_knots] [chains] [iter]

suppressPackageStartupMessages({
  library(rstan)
  library(jsonlite)
})
rstan_options(auto_write = TRUE)
options(mc.cores = 2)

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]; outfile <- args[2]
n_knots <- if (length(args) >= 3) as.integer(args[3]) else 8
n_chains <- if (length(args) >= 4) as.integer(args[4]) else 2
n_iter <- if (length(args) >= 5) as.integer(args[5]) else 800

## --- Ellison et al. (2026) P-spline basis (from their r/0_functions.r) ---
tpower <- function(x, t, p) (x - t)^p * (x > t)
bbase <- function(x, xl = min(x), xr = max(x), n.knots = 10, deg = 3) {
  nseg <- n.knots - 1
  dx <- (xr - xl) / nseg
  knots <- seq(xl - deg * dx, xr + deg * dx, len = n.knots + 2 * deg)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)
  B
}

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  vector[N] age;
  vector[N] y;
  matrix[N, K] B;
}
parameters {
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  beta[1] ~ normal(0, 5);
  beta[2] ~ normal(0, 5);
  for (k in 3:K)
    (beta[k] - 2 * beta[k-1] + beta[k-2]) ~ normal(0, 0.5);
  sigma ~ normal(0, 0.5);
  y ~ normal(B * beta, sigma);
}
generated quantities {
  vector[N] eta = B * beta;
}
"

sm <- stan_model(model_code = stan_code, model_name = "ellison_pspline")

inp <- fromJSON(infile, simplifyVector = FALSE)
ages_all <- unlist(inp$cases[[1]]$age)
res <- list()
for (cs in inp$cases) {
  id <- cs$id
  y <- unlist(cs$y)
  m <- !unlist(cs$y_na)        # observed mask
  if (sum(m) < 6) {            # too few observations: persistence fallback
    pred <- rep(NA_real_, length(y))
    res[[length(res) + 1]] <- list(id = id, pred = pred, q05 = pred, q95 = pred,
                                   method = "insufficient_data")
    next
  }
  B <- bbase(ages_all, n.knots = n_knots, deg = 3)
  sd_data <- list(N = sum(m), K = ncol(B), age = ages_all[m],
                  y = log(pmax(y[m], 1e-5)), B = B[m, , drop = FALSE])
  fit <- sampling(sm, data = sd_data, chains = n_chains, iter = n_iter,
                  refresh = 0, show_messages = FALSE)
  eta <- rstan::extract(fit, pars = "eta")$eta      # draws x N_obs
  Bfull <- B
  # posterior predictive for all ages: re-evaluate B%*%beta per draw
  beta_draws <- rstan::extract(fit, pars = "beta")$beta
  mu <- exp(beta_draws %*% t(Bfull))                # draws x 35
  pred <- apply(mu, 2, median)
  q05 <- apply(mu, 2, quantile, 0.05)
  q95 <- apply(mu, 2, quantile, 0.95)
  res[[length(res) + 1]] <- list(id = id, pred = as.numeric(pred),
                                 q05 = as.numeric(q05), q95 = as.numeric(q95),
                                 method = "ellison_gam_hmc")
}
write_json(list(cases = res), outfile, auto_unbox = TRUE, digits = 8)
cat("done:", length(res), "cases ->", outfile, "\n")
