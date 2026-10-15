# MR-PRESSO re-run with the official R package (primary 52-variant instrument set)
# Input: harmonised_<outcome>.tsv files written by the Python pipeline
# Output: mrpresso_results.csv and mrpresso_log.txt in this folder
base <- path.expand("~/rc_nafld_mr/analysis")
setwd(file.path(base, "R"))
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
if (!requireNamespace("MRPRESSO", quietly = TRUE)) remotes::install_github("rondolab/MR-PRESSO", upgrade = "never")
library(MRPRESSO)
sets <- c(Ghodsian = "GCST90091033", Fairfield = "GCST90054782", FinnGen = "finn-b-NAFLD")
res <- list(); sink("mrpresso_log.txt", split = TRUE)
for (nm in names(sets)) {
  d <- read.delim(file.path(base, "primary_clumped52", paste0("harmonised_", sets[[nm]], ".tsv")))
  cat("\n=====", nm, "n variants =", nrow(d), "=====\n")
  set.seed(2026)
  p <- mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "sey", SdExposure = "sex",
                 OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = d,
                 NbDistribution = 10000, SignifThreshold = 0.05)
  print(p$`Main MR results`); print(p$`MR-PRESSO results`$`Global Test`)
  print(p$`MR-PRESSO results`$`Outlier Test`)
  mr <- p$`Main MR results`
  glob <- p$`MR-PRESSO results`$`Global Test`
  out <- p$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  res[[nm]] <- data.frame(dataset = nm, n = nrow(d),
    raw_beta = mr$`Causal Estimate`[1], raw_sd = mr$Sd[1], raw_p = mr$`P-value`[1],
    corrected_beta = mr$`Causal Estimate`[2], corrected_sd = mr$Sd[2], corrected_p = mr$`P-value`[2],
    global_RSSobs = glob$RSSobs, global_p = glob$Pvalue,
    outliers = if (is.null(out) || all(is.na(out))) "none" else paste(d$SNP[out], collapse = "; "))
}
sink()
write.csv(do.call(rbind, res), "mrpresso_results.csv", row.names = FALSE)
cat("Done. Results in mrpresso_results.csv\n")
