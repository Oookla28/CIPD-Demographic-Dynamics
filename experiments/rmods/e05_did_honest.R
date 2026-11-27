#!/usr/bin/env Rscript
# E05 R component:
#  (1) HonestDiD Delta^SD(M) sensitivity for the two PyFixest event studies;
#  (2) bcallaway11/did att_gt on the dose-contrast design (2014 timing).
# Usage: Rscript e05_did_honest.R <panel_csv> <did_out.json> <honest_out.json>
suppressPackageStartupMessages({
  library(jsonlite)
  library(HonestDiD)
})

args <- commandArgs(trailingOnly = TRUE)
panel_csv <- args[1]
did_out <- args[2]
honest_out <- args[3]

## ---------------- (1) HonestDiD on exported event-study coefficients -------
es_files <- c(A_2014 = "results/logs/e05_event_study_A_2014.json",
              B_2016 = "results/logs/e05_event_study_B_2016.json")
sens_rows <- list()
for (nm in names(es_files)) {
  f <- es_files[[nm]]
  if (!file.exists(f)) next
  es <- fromJSON(f)
  betahat <- as.numeric(es$beta)
  sigma <- matrix(unlist(es$vcov), nrow = length(betahat))
  numPre <- as.integer(es$n_pre); numPost <- as.integer(es$n_post)
  for (M in c(0, 0.5, 1, 1.5, 2, 3)) {
    res <- tryCatch({
      out <- HonestDiD::createSensitivityResults_relativeMagnitudes(
        betahat = betahat, sigma = sigma,
        numPrePeriods = numPre, numPostPeriods = numPost,
        Mbarvec = M, method = "C-LF", gridPoints = 2000, grid.ub = 2, grid.lb = -2)
      data.frame(policy = nm, M = M, M_effective = M,
                 lb = out$lb[1], ub = out$ub[1],
                 method = "FLCI_relativeMagnitudes")
    }, error = function(e) {
      data.frame(policy = nm, M = M, M_effective = M, lb = NA, ub = NA,
                 method = paste("error:", conditionMessage(e)))
    })
    sens_rows[[length(sens_rows) + 1]] <- res
  }
  # smoothness Delta^SD with M calibrated to the pre-period second-difference scale
  pre_b <- betahat[1:numPre]
  if (length(pre_b) >= 3) {
    d2 <- abs(diff(pre_b, differences = 2))
    M_scale <- max(mean(d2), 1e-8)
  } else M_scale <- 1
  for (M_raw in c(0, 0.5, 1, 1.5, 2, 3)) {
    M <- M_raw * M_scale
    res <- tryCatch({
      out <- HonestDiD::createSensitivityResults(
        betahat = betahat, sigma = sigma,
        numPrePeriods = numPre, numPostPeriods = numPost,
        Mvec = M, method = "FLCI")
      data.frame(policy = nm, M = M_raw, M_effective = M, lb = out$lb[1], ub = out$ub[1],
                 method = "FLCI_smoothness")
    }, error = function(e) data.frame(policy = nm, M = M_raw, M_effective = M, lb = NA, ub = NA,
                                      method = paste("error:", conditionMessage(e))))
    sens_rows[[length(sens_rows) + 1]] <- res
  }
}
sens <- do.call(rbind, sens_rows)
write_json(list(sensitivity = sens), honest_out, auto_unbox = TRUE, digits = 10)
cat("HonestDiD rows:", nrow(sens), "\n")

## ---------------- (2) did::att_gt dose-contrast design ---------------------
suppressPackageStartupMessages(library(did))
d <- read.csv(panel_csv)
d$province <- as.numeric(factor(d$province))
out <- list()
run_one <- function(dat, intensity, policy_year, label) {
  dat$G <- ifelse(dat[[intensity]] > median(dat[[intensity]], na.rm = TRUE),
                  policy_year, 0)
  dat <- dat[!is.na(dat$log_births), ]
  res <- tryCatch({
    att <- did::att_gt(yname = "log_births", tname = "year", idname = "province",
                       gname = "G", data = dat, control_group = "nevertreated",
                       est_method = "dr", base_period = "universal",
                       print_details = FALSE, allow_unbalanced_panel = TRUE)
    agg <- did::aggte(att, type = "simple", na.rm = TRUE)
    list(label = label, att = agg$overall.att, se = agg$overall.se,
         crit_val = agg$crit.val.egt,
         ci_lo = agg$overall.att - agg$crit.val.egt * agg$overall.se,
         ci_hi = agg$overall.att + agg$crit.val.egt * agg$overall.se,
         group_time = as.data.frame(att$group), att_gt = as.numeric(att$att),
         se_gt = as.numeric(att$se))
  }, error = function(e) list(label = label, error = conditionMessage(e)))
  res
}
out[[1]] <- run_one(d, "D", 2014, "2014_dandui_timing_dose_contrast")
out[[2]] <- run_one(d, "E_onechild", 2016, "2016_quanmian_intensity_dose_contrast")
write_json(out, did_out, auto_unbox = TRUE, digits = 10)
for (o in out) {
  if (!is.null(o$error)) cat("did ERROR:", o$label, o$error, "\n")
  else cat("did OK:", o$label, "ATT=", round(o$att, 5), "SE=", round(o$se, 5), "\n")
}
