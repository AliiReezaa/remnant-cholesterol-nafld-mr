# Extended MR methods and diagnostics for the primary 52-variant instrument set (three NAFLD outcome datasets)
# Input : analysis/primary_clumped52/harmonised_<dataset>.tsv (written by the Python pipeline)
# Output: extended_analyses/R_out/*.csv and R_out/R_console_log.txt
base <- path.expand("~/Documents/02_PhD_Research/Fatty Liver")
setwd(file.path(base, "extended_analyses")); dir.create("R_out", showWarnings = FALSE)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
for (p in c("MendelianRandomization", "mr.raps")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
if (!requireNamespace("RadialMR", quietly = TRUE)) try(remotes::install_github("WSpiller/RadialMR", upgrade = "never"))
con_log <- file("R_out/R_console_log.txt", open = "wt"); sink(con_log, split = TRUE); sink(con_log, type = "message")
cat("Run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n"); cat(R.version.string, "\n")
suppressMessages({library(MendelianRandomization); library(mr.raps)})
have_radial <- requireNamespace("RadialMR", quietly = TRUE); if (have_radial) suppressMessages(library(RadialMR))
cat("MendelianRandomization", as.character(packageVersion("MendelianRandomization")), "| mr.raps", as.character(packageVersion("mr.raps")),
    "| RadialMR", if (have_radial) as.character(packageVersion("RadialMR")) else "NOT AVAILABLE", "\n")
sets <- c(Ghodsian = "GCST90091033", Fairfield = "GCST90054782", FinnGen = "finn-b-NAFLD")
est <- list(); diag <- list(); infl <- list()
tidy <- function(ds, method, b, se, p, n, extra = NA) data.frame(dataset = ds, method = method, nsnp = n, b = b, se = se, pval = p,
                                                              OR = exp(b), OR_lo = exp(b - 1.96 * se), OR_hi = exp(b + 1.96 * se), note = extra)
for (nm in names(sets)) {
  cat("\n=====", nm, "=====\n")
  d <- read.delim(file.path(base, "analysis", "primary_clumped52", paste0("harmonised_", sets[[nm]], ".tsv")))
  n <- nrow(d); cat("variants:", n, "\n")
  inp <- mr_input(bx = d$bx, bxse = d$sex, by = d$by, byse = d$sey, snps = d$SNP)
  set.seed(2026)
  r_ml <- mr_maxlik(inp);  est[[length(est) + 1]] <- tidy(nm, "Maximum likelihood", r_ml@Estimate, r_ml@StdError, r_ml@Pvalue, n)
  r_ls <- tryCatch(mr_lasso(inp), error = function(e) { cat("MR-Lasso failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(r_ls)) est[[length(est) + 1]] <- tidy(nm, "MR-Lasso", r_ls@Estimate, r_ls@StdError, r_ls@Pvalue, n, paste("valid instruments:", paste(r_ls@Valid, collapse = ",")))
  r_cm <- tryCatch(mr_conmix(inp), error = function(e) NULL)
  if (!is.null(r_cm)) est[[length(est) + 1]] <- tidy(nm, "Contamination mixture", r_cm@Estimate, (r_cm@CIUpper - r_cm@CILower) / (2 * 1.96), NA, n, paste0("CI ", round(exp(r_cm@CILower), 3), "-", round(exp(r_cm@CIUpper), 3), "; valid ", length(r_cm@Valid)))
  r_dv <- tryCatch(mr_divw(inp), error = function(e) NULL)
  if (!is.null(r_dv)) est[[length(est) + 1]] <- tidy(nm, "Debiased IVW", r_dv@Estimate, r_dv@StdError, r_dv@Pvalue, n)
  r_iv <- mr_ivw(inp, model = "random"); est[[length(est) + 1]] <- tidy(nm, "IVW (random effects, MendelianRandomization)", r_iv@Estimate, r_iv@StdError, r_iv@Pvalue, n)
  r_fx <- mr_ivw(inp, model = "fixed");  est[[length(est) + 1]] <- tidy(nm, "IVW (fixed effect)", r_fx@Estimate, r_fx@StdError, r_fx@Pvalue, n)
  for (loss in c("tukey", "huber")) {
    rr <- mr.raps(d$bx, d$by, d$sex, d$sey, over.dispersion = TRUE, loss.function = loss)
    est[[length(est) + 1]] <- tidy(nm, paste0("MR-RAPS (", loss, ", overdispersed)"), rr$beta.hat, rr$beta.se, 2 * pnorm(-abs(rr$beta.hat / rr$beta.se)), n)
  }
  # heterogeneity statistics
  Q <- r_iv@Heter.Stat[1]; Qp <- r_iv@Heter.Stat[2]; I2 <- max(0, (Q - (n - 1)) / Q)
  w <- 1 / d$sex^2; gbar <- sum(w * d$bx) / sum(w); Qgx <- sum(w * (d$bx - gbar)^2); I2gx <- (Qgx - (n - 1)) / Qgx
  r_eg <- mr_egger(inp)
  diag[[nm]] <- data.frame(dataset = nm, nsnp = n, Q_IVW = Q, Q_IVW_p = Qp, I2_IVW = I2, Q_Egger = r_eg@Heter.Stat[1], Q_Egger_p = r_eg@Heter.Stat[2],
                           I2_GX = I2gx, Egger_intercept = r_eg@Intercept, Egger_intercept_se = r_eg@StdError.Int, Egger_intercept_p = r_eg@Pvalue.Int,
                           phi_random_effects = r_iv@RSE^2)
  # F >= 30 and MAF >= 0.01 filters
  for (flt in c("F>=30", "MAF>=0.01")) {
    keep <- if (flt == "F>=30") d$F >= 30 else pmin(d$eaf_x, 1 - d$eaf_x) >= 0.01
    inp2 <- mr_input(bx = d$bx[keep], bxse = d$sex[keep], by = d$by[keep], byse = d$sey[keep], snps = d$SNP[keep])
    r2 <- mr_ivw(inp2, model = "random"); est[[length(est) + 1]] <- tidy(nm, paste("IVW, variants with", flt), r2@Estimate, r2@StdError, r2@Pvalue, sum(keep))
  }
  # influence diagnostics (Cook's distance, studentised residuals, hat values) from the weighted IVW regression
  fit <- lm(by ~ 0 + bx, weights = 1 / d$sey^2, data = d)
  infl[[nm]] <- data.frame(dataset = nm, SNP = d$SNP, cooks_d = cooks.distance(fit), student_resid = rstudent(fit), hat = hatvalues(fit),
                           cook_threshold_4_over_n = 4 / n)
  # radial MR
  if (have_radial) {
    rd <- format_radial(d$bx, d$by, d$sex, d$sey, d$SNP)
    rv <- ivw_radial(rd, alpha = 0.05, weights = 3, tol = 0.0001)
    cat("Radial IVW: estimate", rv$coef[1, 1], "SE", rv$coef[1, 2], "; outliers:", paste(rv$outliers$SNP, collapse = ", "), "\n")
    rq <- rv$qstatistic; cat("Radial Q:", rq, "\n")
    write.csv(rv$data, paste0("R_out/radial_data_", nm, ".csv"), row.names = FALSE)
    est[[length(est) + 1]] <- tidy(nm, "Radial IVW (modified second order weights)", rv$coef[1, 1], rv$coef[1, 2], rv$coef[1, 4], n,
                                   paste("radial outliers:", paste(rv$outliers$SNP, collapse = ";")))
  }
}
E <- do.call(rbind, est); D <- do.call(rbind, diag); I <- do.call(rbind, infl)
print(E[, c("dataset", "method", "nsnp", "OR", "OR_lo", "OR_hi", "pval", "note")], row.names = FALSE); print(D, row.names = FALSE)
write.csv(E, "R_out/extended_estimates.csv", row.names = FALSE); write.csv(D, "R_out/heterogeneity_I2.csv", row.names = FALSE)
write.csv(I, "R_out/influence_diagnostics.csv", row.names = FALSE)
cat("\n==== sessionInfo ====\n"); print(sessionInfo()); cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
sink(type = "message"); sink(); close(con_log)
cat("DONE - extended analyses written to extended_analyses/R_out\n")
