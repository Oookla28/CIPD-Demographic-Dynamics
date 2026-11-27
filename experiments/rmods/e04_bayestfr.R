#!/usr/bin/env Rscript
# E04 baseline: bayesTFR projections at rolling origins; quantiles from the
# prediction object (5-year native resolution). Country indices taken from
# pred$tfr_matrix_reconstructed colnames (UN country codes).
suppressPackageStartupMessages({ library(bayesTFR) })
args <- commandArgs(trailingOnly = TRUE)
out_csv <- args[1]
origins <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else c(2005L, 2010L, 2016L)
n_iter <- if (length(args) >= 3) as.integer(args[3]) else 600L
data(UNlocations, package = "wpp2019")
iso_name <- c(JPN="Japan", KOR="Republic of Korea", DEUTNP="Germany", SWE="Sweden",
  NOR="Norway", FIN="Finland", DNK="Denmark", AUT="Austria", CZE="Czechia",
  NLD="Netherlands", HUN="Hungary", ITA="Italy", ESP="Spain",
  USA="United States of America", CAN="Canada", POL="Poland",
  GBRTENW="United Kingdom", CHE="Switzerland", TWN="China, Taiwan Province of China")
codes <- sapply(iso_name, function(nm) {
  i <- which(UNlocations$name == nm)
  if (length(i) == 0) NA_integer_ else UNlocations$country_code[i[1]] })
out <- list()
for (origin in origins) {
  cat("=== origin", origin, "\n")
  outdir <- tempfile()
  ok <- tryCatch({
    mcmc <- run.tfr.mcmc(nr.chains = 1, iter = n_iter, output.dir = outdir,
                         start.year = 1950, present.year = origin, wpp.year = 2019,
                         annual = FALSE, uncertainty = TRUE, verbose = FALSE)
    pred <- NULL
    for (ey in c(origin + 10, origin + 8)) {
      pred <- tryCatch(tfr.predict(mcmc, start.year = origin, end.year = ey,
                                   burnin = max(100, n_iter %/% 4), nr.traj = 200,
                                   verbose = FALSE, uncertainty = TRUE),
                       error = function(e) { cat("  predict end.year", ey, "ERR:", conditionMessage(e), "\n"); NULL })
      if (!is.null(pred)) break
    }
    if (is.null(pred)) stop("prediction failed")
    q <- pred$quantiles
    yrs <- as.integer(dimnames(q)[[3]])
    qlev <- as.numeric(dimnames(q)[[2]])
    tm <- pred$tfr_matrix_reconstructed
    cn <- colnames(tm)
    for (ctry in names(codes)) {
      cc <- codes[[ctry]]
      if (is.na(cc)) next
      idx <- which(cn == as.character(cc))
      if (length(idx) == 0) next
      for (h in 1:10) {
        yr <- origin + h
        yi <- which(yrs == yr)
        if (length(yi) == 0) next
        gq <- function(p) q[idx[1], which.min(abs(qlev - p)), yi[1]]
        out[[length(out) + 1]] <- data.frame(country = ctry, origin = origin,
                                             horizon = h, year = yr,
                                             median = gq(0.5), q05 = gq(0.05), q95 = gq(0.95))
      }
    }
    TRUE
  }, error = function(e) { cat("ERR", conditionMessage(e), "\n"); FALSE })
  unlink(outdir, recursive = TRUE)
}
if (length(out) == 0) cat("NO OUTPUT\n") else {
  res <- do.call(rbind, out)
  write.csv(res, out_csv, row.names = FALSE)
  cat("wrote", nrow(res), "rows ->", out_csv, "\n")
}
