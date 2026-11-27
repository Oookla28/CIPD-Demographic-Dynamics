#!/usr/bin/env Rscript
# Single event-study file -> HonestDiD relative-magnitudes sensitivity grid.
suppressPackageStartupMessages({ library(jsonlite); library(HonestDiD) })
args <- commandArgs(trailingOnly = TRUE)
es <- fromJSON(args[1])
betahat <- as.numeric(es$beta)
sigma <- matrix(unlist(es$vcov), nrow = length(betahat))
numPre <- as.integer(es$n_pre); numPost <- as.integer(es$n_post)
rows <- list()
for (M in c(0, 0.5, 1, 1.5, 2, 3)) {
  r <- tryCatch({
    o <- HonestDiD::createSensitivityResults_relativeMagnitudes(
      betahat = betahat, sigma = sigma, numPrePeriods = numPre,
      numPostPeriods = numPost, Mbarvec = M, method = "C-LF",
      gridPoints = 2000, grid.ub = 2, grid.lb = -2)
    data.frame(M = M, lb = o$lb[1], ub = o$ub[1], method = "FLCI_relativeMagnitudes")
  }, error = function(e) data.frame(M = M, lb = NA, ub = NA,
                                    method = paste("error:", conditionMessage(e))))
  rows[[length(rows) + 1]] <- r
}
write_json(list(sensitivity = do.call(rbind, rows)), args[2], auto_unbox = TRUE, digits = 10)
cat("ok\n")
