# Positive control, R arm: genetically predicted remnant cholesterol -> coronary artery disease (Nikpay 2015, GCST003116)
# Independent of the Python harmonisation: raw exposure + raw outcome tables are harmonised by TwoSampleMR itself.
base <- path.expand("~/Documents/02_PhD_Research/Fatty Liver/positive_control_CAD")
setwd(base); dir.create("R_out", showWarnings = FALSE)
for (p in c("remotes")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) remotes::install_github("MRCIEU/TwoSampleMR", upgrade = "never")
if (!requireNamespace("MRPRESSO", quietly = TRUE)) remotes::install_github("rondolab/MR-PRESSO", upgrade = "never")
library(TwoSampleMR); library(MRPRESSO)
con_log <- file("R_out/R_console_log.txt", open = "wt")   # everything printed to the console (stdout) and all messages/warnings
sink(con_log, split = TRUE); sink(con_log, type = "message")
cat("Run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n"); cat(R.version.string, "\n")
cat("MRPRESSO version:", as.character(packageVersion("MRPRESSO")), "\n")
cat("TwoSampleMR version:", as.character(packageVersion("TwoSampleMR")), "\n")
ex <- read.delim("inputs/rc_52_instruments_exposure.tsv")
ex$F <- (ex$beta / ex$standard_error)^2
exp_dat <- format_data(ex, type = "exposure", snp_col = "variant_id", beta_col = "beta", se_col = "standard_error",
                       effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                       eaf_col = "effect_allele_frequency", pval_col = "p_value", chr_col = "chr", pos_col = "pos")
exp_dat$exposure <- "Remnant cholesterol"; exp_dat$id.exposure <- "RC"
out <- read.delim("inputs/cad_outcome_for_52_instruments.tsv")
out_dat <- format_data(out, type = "outcome", snp_col = "rsid", beta_col = "beta", se_col = "se",
                       effect_allele_col = "ea", other_allele_col = "oa", eaf_col = "eaf", pval_col = "p")
out_dat$outcome <- "CAD (Nikpay 2015)"; out_dat$id.outcome <- "CAD"
dat <- harmonise_data(exp_dat, out_dat, action = 2)
cat("Harmonised SNPs kept for MR:", sum(dat$mr_keep), " of ", nrow(dat), "\n")
write.csv(dat, "R_out/harmonised_TwoSampleMR.csv", row.names = FALSE)
res <- mr(dat, method_list = c("mr_ivw_mre", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))
res$OR <- exp(res$b); res$OR_lo <- exp(res$b - 1.96 * res$se); res$OR_hi <- exp(res$b + 1.96 * res$se)
print(res[, c("method", "nsnp", "b", "se", "pval", "OR", "OR_lo", "OR_hi")])
write.csv(res, "R_out/mr_results_TwoSampleMR.csv", row.names = FALSE)
het <- mr_heterogeneity(dat); print(het); write.csv(het, "R_out/heterogeneity.csv", row.names = FALSE)
ple <- mr_pleiotropy_test(dat); print(ple); write.csv(ple, "R_out/egger_intercept.csv", row.names = FALSE)
loo <- mr_leaveoneout(dat); write.csv(loo, "R_out/leave_one_out.csv", row.names = FALSE)
d <- dat[dat$mr_keep, ]
set.seed(2026)
p <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure",
               OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = d, NbDistribution = 10000, SignifThreshold = 0.05)
print(p)  # full MR-PRESSO object
print(p$`Main MR results`); print(p$`MR-PRESSO results`$`Global Test`); print(p$`MR-PRESSO results`$`Outlier Test`)
saveRDS(p, "R_out/mrpresso_object.rds")
mr_ <- p$`Main MR results`
write.csv(mr_, "R_out/mrpresso_main.csv", row.names = FALSE)
cat("Distortion:\n"); print(p$`MR-PRESSO results`$`Distortion Test`)
cat("Outlier SNPs:", paste(d$SNP[p$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`], collapse = ", "), "\n")
cat("\n==== sessionInfo ====\n"); print(sessionInfo())
cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
sink(type = "message"); sink(); close(con_log)
cat("DONE - results in R_out\n")
