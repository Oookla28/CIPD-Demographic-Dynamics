#!/usr/bin/env Rscript
# E04 baseline: traditional cohort-component projection of China births via
# ccmpp.tmb::ccmppR, annual single-age steps (WPP2024 single-age population).
suppressPackageStartupMessages({ library(ccmpp.tmb); library(wpp2024) })
args <- commandArgs(trailingOnly = TRUE)
out_csv <- args[1]
tfr_hist <- c(`2005`=1.62, `2010`=1.69, `2013`=1.66, `2016`=1.77, `2020`=1.28)
data(popAge1dt, package = "wpp2024")
out <- list()
for (origin in c(2005, 2010, 2013, 2016, 2020)) {
  d <- popAge1dt[popAge1dt$country_code == 156 & popAge1dt$year == origin, ]
  if (nrow(d) < 90) { cat("no pop", origin, "\n"); next }
  d <- d[order(d$age), ]
  A <- nrow(d)
  basepop <- c(d$popM * 1000, d$popF * 1000)
  # 1-year survival: age-specific death probability from a simple Gompertz-ish
  # schedule calibrated to China's ~e0=77 (female) / 74 (male) in the period
  qx_f <- pmin(1, 0.0012 * exp(0.085 * (d$age - 30)) + ifelse(d$age < 1, 0.004, 0))
  qx_m <- pmin(1, 0.0016 * exp(0.088 * (d$age - 30)) + ifelse(d$age < 1, 0.005, 0))
  sx_f <- 1 - qx_f; sx_m <- 1 - qx_m
  sx_f[A] <- 0.65; sx_m[A] <- 0.60
  Tp <- 12
  sx <- matrix(rep(c(sx_m, sx_f), Tp), ncol = Tp)
  gx <- matrix(0, nrow = 2 * A, ncol = Tp)
  # ASFR single-age: beta-like schedule over 15-49 summing to TFR
  ages <- d$age
  w <- ifelse(ages >= 15 & ages <= 49,
              dnorm(ages, mean = 28.5, sd = 5.5) / sum(dnorm(15:49, 28.5, 5.5)), 0)
  fx1 <- w * tfr_hist[as.character(origin)]
  fx <- matrix(rep(fx1, Tp), ncol = Tp)
  srb <- rep(1.05, Tp)
  fx_idx <- which(ages == 15)[1]   # 1-based index of age 15
  if (is.na(fx_idx)) next
  res <- tryCatch(ccmppR(basepop = basepop, sx = sx, fx = fx, gx = gx,
                         srb = srb, age_span = 1, fx_idx = as.integer(fx_idx)),
                  error = function(e) { cat("ERR", origin, conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next
  bn <- names(res)
  cat("origin", origin, "| slots:", paste(bn, collapse=","), "\n")
  b <- tryCatch(colSums(res$births_by_age), error = function(e) NULL)
  if (is.null(b)) b <- tryCatch(colSums(res$births), error = function(e) NULL)
  if (is.null(b)) next
  for (h in 1:min(10, length(b))) {
    out[[length(out) + 1]] <- data.frame(origin = origin, horizon = h,
                                         year = origin + h, births_wan = b[h] / 1e4)
  }
  cat("   births/yr (wan):", paste(round(b[1:min(6,length(b))]/1e4, 1), collapse=","), "\n")
}
if (length(out) == 0) cat("NO OUTPUT\n") else {
  write.csv(do.call(rbind, out), out_csv, row.names = FALSE)
  cat("wrote", out_csv, "\n")
}
