###############################################################################
#  GENETIC CORRELATION (LD SCORE REGRESSION): REMNANT CHOLESTEROL (RC) vs NAFLD
#  Study: two-sample Mendelian randomization of RC and NAFLD (Soleymani et al.)
#
#  PURPOSE
#  Run this file top to bottom in RStudio (Ctrl+Alt+R, or step with Ctrl+Enter).
#  It reproduces, from raw public downloads, every number of the genetic
#  correlation analysis: SNP heritability, rg of RC with two NAFLD GWAS and
#  with coronary artery disease (CAD, positive control), the external
#  benchmark with GLGC lipid GWAS, the sensitivity analyses, the bulk-genome
#  analysis, an independent cross-check and the instrument-composition check.
#
#  METHOD IN ONE PARAGRAPH
#  LD score regression (Bulik-Sullivan et al. 2015, Nat Genet 47:291 and 47:1236)
#  regresses chi2 (heritability) or z1*z2 (genetic covariance) of every HapMap3
#  SNP on its LD score. Slope = N*h2/M (or sqrt(N1N2)*rho_g/M), intercept =
#  inflation (or sample-overlap) term. rg = rho_g / sqrt(h2_1 * h2_2). Standard
#  errors come from a 200-block block jackknife.
#
#  IMPORTANT: THIS IS A LINE-BY-LINE R PORT OF THE OFFICIAL LDSC SOFTWARE
#  The reported numbers were produced with the official LDSC
#  code (github.com/bulik/ldsc, v2.0.x, Python 3 build). So that you can rerun
#  everything in R without Python, the functions in Sections 2 and 3
#  (munge_sumstats_R, ldsc_core, irwls, lstsq_jackknife, ldsc_hsq, ldsc_gencov,
#  ldsc_rg) re-implement the official munge_sumstats.py and
#  ldscore/{regressions,jackknife,irwls,sumstats}.py step by step. Section 13
#  prints the R results next to the official-software results for all 36 rg
#  comparisons and the 4 heritabilities: the largest difference is 0.00005
#  (the official numbers are rounded to 4 decimals), i.e. they agree.
#
#  REQUIREMENTS
#  R >= 4.0, packages data.table and R.utils (installed below if missing),
#  about 12 GB free disk (downloads ~3 GB, plus working files), at least 8 GB
#  RAM (tested on a 2-core machine with 7 GB), and internet for Section 1.
#  Run time: downloads (about 3 GB) plus roughly 15 minutes of computation.
#  Munged files are kept; set FORCE_REMUNGE <- TRUE before running to redo them.
#
#  HOW TO RUN
#  Open this file in RStudio, press Ctrl+Shift+Enter (source), or run it line by
#  line with Ctrl+Enter. To use another folder, set WORKDIR <- "/your/path"
#  in the console before running. Everything printed is also saved in logs/.
#  Results: results/ldsc_h2_R.csv, results/ldsc_rg_all_runs_R.csv,
#  results/instrument_composition_R.csv.
#
#  ANALYSIS DECISIONS
#   - Regression SNPs: HapMap3 SNPs, MHC excluded (1,187,349 SNPs).
#   - LD scores: 1000 Genomes Phase 3 European (S-LDSC reference, Zenodo 7768714).
#   - Primary run: free intercepts; h2 with the two-step estimator (cutoff 30,
#     LDSC default); rg with the one-step estimator (LDSC default).
#   - Case-control traits: N = cases + controls; liability-scale h2 uses
#     sample prevalence = cases / (cases + controls). Population prevalence:
#     NAFLD = sample prevalence (record-based diagnosis in a cohort), CAD = 0.05.
#     (rg itself does not depend on these prevalences.)
#   - No chi2 truncation in rg (LDSC default).
###############################################################################

## ---------------------------------------------------------------------------
## 0. SETUP
## ---------------------------------------------------------------------------
# 0.1 Choose the working folder (change this line if you like). Everything is
#     written inside it: downloads, munged files, results, logs.
if (!exists("FORCE_REMUNGE")) FORCE_REMUNGE <- FALSE          # TRUE = redo munging even if files exist
if (!exists("WORKDIR")) WORKDIR <- path.expand("~/RC_NAFLD_genetic_correlation")
dir.create(WORKDIR, showWarnings = FALSE, recursive = TRUE)
for (d in c("ref", "gwas", "munged", "results", "logs"))
  dir.create(file.path(WORKDIR, d), showWarnings = FALSE)
setwd(WORKDIR)

# 0.2 Packages
for (p in c("data.table", "R.utils"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({ library(data.table); library(R.utils) })
setDTthreads(0)                       # use all cores for reading big files
options(timeout = 24 * 3600)          # allow very long downloads
options(width = 200, digits = 7, scipen = 999)

# 0.3 Log everything that is printed to a file as well as to the console
LOGFILE <- file.path(WORKDIR, "logs", paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(LOGFILE, split = TRUE)
cat("Started:", format(Sys.time()), "\n"); print(sessionInfo())

## ---------------------------------------------------------------------------
## 1. DOWNLOAD THE PUBLIC DATA (skipped if a file is already present and its
##    MD5 checksum is correct)
## ---------------------------------------------------------------------------
md5_ok <- function(f, md5) unname(tools::md5sum(f)) == md5
get_file <- function(url, dest, md5) {
  if (file.exists(dest) && !md5_ok(dest, md5)) { message("Checksum mismatch, re-downloading ", basename(dest)); file.remove(dest) }
  if (!file.exists(dest)) {
    message("Downloading ", basename(dest), " ...")
    download.file(url, dest, mode = "wb", method = "libcurl")
  }
  if (!md5_ok(dest, md5)) stop("MD5 check FAILED for ", dest, ". Delete it and run this line again.")
  message("OK  ", basename(dest), "  (MD5 verified)")
  invisible(dest)
}

# 1.1 LD score reference files: S-LDSC reference set, Zenodo record 10.5281/zenodo.7768714
ZEN <- "https://zenodo.org/api/records/7768714/files/%s/content"
zen_files <- c("1000G_Phase3_weights_hm3_no_MHC.tgz" = "a98ac0f089ee285177544a3e6e721ca3",  # HapMap3 LD scores, no MHC
               "1000G_Phase3_ldscores.tgz"           = "e4352ccf778f296835d73985350a863b",  # LD scores summed over all 1000G SNPs (used in sensitivity B)
               "1000G_Phase3_frq.tgz"                = "ac29686ffd5b6378789857a522ebca77")  # alleles of the 1000G SNPs
zen_dirs <- c("1000G_Phase3_weights_hm3_no_MHC.tgz" = "1000G_Phase3_weights_hm3_no_MHC",
              "1000G_Phase3_ldscores.tgz"           = "LDscore",
              "1000G_Phase3_frq.tgz"                = "1000G_Phase3_frq")
for (f in names(zen_files)) {
  dest <- file.path("ref", f)
  get_file(sprintf(ZEN, f), dest, zen_files[[f]])
  if (!dir.exists(file.path("ref", zen_dirs[[f]]))) untar(dest, exdir = "ref")   # unpack once
}

# 1.2 GWAS summary statistics (all full genome-wide files)
EBI <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
gwas_dl <- list(
  # Remnant cholesterol, Richardson et al. 2022 (UK Biobank NMR, N = 115,082), GCST90092943
  RC        = list(url = paste0(EBI, "/GCST90092001-GCST90093000/GCST90092943/harmonised/35213538-GCST90092943-EFO_0010815.h.tsv.gz"),
                   file = "RC_GCST90092943.h.tsv.gz", md5 = "20717c9a039e3c8f04641fc927b86610"),
  # NAFLD, Ghodsian et al. 2021 (8,434 cases / 770,180 controls), GCST90091033
  Ghodsian  = list(url = paste0(EBI, "/GCST90091001-GCST90092000/GCST90091033/harmonised/34841290-GCST90091033-EFO_0003095.h.tsv.gz"),
                   file = "Ghodsian_GCST90091033.h.tsv.gz", md5 = "126496bbc3248f90319ff9406b1c76c3"),
  # NAFLD, Fairfield et al. 2022 (4,761 cases / 373,227 controls), GCST90054782. ORIGINAL author file (it contains the INFO column)
  Fairfield = list(url = paste0(EBI, "/GCST90054001-GCST90055000/GCST90054782/GCST90054782_buildGRCh37.tsv"),
                   file = "Fairfield_GCST90054782_buildGRCh37.tsv", md5 = "17881eddb45e2bd0c8f3cd617902d446"),
  # Coronary artery disease, Nikpay et al. 2015 CARDIoGRAMplusC4D (60,801 cases / 123,504 controls), GCST003116  (positive control)
  CAD       = list(url = paste0(EBI, "/GCST003001-GCST004000/GCST003116/harmonised/26343387-GCST003116-EFO_0000378.h.tsv.gz"),
                   file = "CAD_GCST003116.h.tsv.gz", md5 = "45b58f3ca315eca3d09cb2077cedc2d5"),
  # Global Lipids Genetics Consortium, Willer et al. 2013 (N ~ 188,000): external benchmark
  GLGC_TG   = list(url = paste0(EBI, "/GCST002001-GCST003000/GCST002216/jointGwasMc_TG.txt.gz"),
                   file = "GLGC_TG.txt.gz", md5 = "50e26af9107cca7b2e38fc7eae0ac767"),
  GLGC_LDL  = list(url = paste0(EBI, "/GCST002001-GCST003000/GCST002222/jointGwasMc_LDL%20(1).txt.gz"),
                   file = "GLGC_LDL.txt.gz", md5 = "c01937bd545acc888b119d3f8cb74d76")
)
for (nm in names(gwas_dl)) { g <- gwas_dl[[nm]]; get_file(g$url, file.path("gwas", g$file), g$md5) }

## ---------------------------------------------------------------------------
## 2. HAPMAP3 REGRESSION-SNP LIST WITH ALLELES (what munge_sumstats --merge-alleles uses)
##    SNP set = SNPs in the HapMap3 (MHC-excluded) LD-score weight files.
##    Alleles = A1/A2 of the same SNP in the 1000G Phase 3 EUR frequency files.
## ---------------------------------------------------------------------------
read_ld <- function(files) {
  # one LD-score file per chromosome; LDSC sorts by CHR, BP and drops duplicated rs numbers
  x <- rbindlist(lapply(files, fread, showProgress = FALSE))
  x <- x[order(CHR, BP)]
  unique(x, by = "SNP")[, .(SNP, CHR, BP, L2)]
}
w_files <- list.files("ref/1000G_Phase3_weights_hm3_no_MHC", pattern = "^weights\\.hm3_noMHC\\.[0-9]+\\.l2\\.ldscore\\.gz$", full.names = TRUE)
stopifnot(length(w_files) == 22)
LD_HM3 <- read_ld(w_files)                       # LD scores over HapMap3 SNPs (used as reference AND as regression weights)
M_HM3  <- nrow(LD_HM3)                           # number of SNPs the LD scores are summed over
cat("HapMap3 SNPs (MHC excluded):", M_HM3, "\n")          # expect 1,187,349

frq <- rbindlist(lapply(list.files("ref/1000G_Phase3_frq", pattern = "\\.frq$", full.names = TRUE), fread, showProgress = FALSE))
HM3 <- merge(LD_HM3[, .(SNP)], unique(frq[, .(SNP, A1, A2)], by = "SNP"), by = "SNP")
HM3[, MA := toupper(paste0(A1, A2))]             # 'merged alleles' string, e.g. "GA"
cat("HapMap3 SNPs with alleles:", nrow(HM3), "\n")       # expect 1,187,349
rm(frq); invisible(gc())

# Full-1000G LD scores (sensitivity analysis B): LD summed over ~5.96 million 1000G SNPs
l_files <- list.files("ref/LDscore", pattern = "^LDscore\\.[0-9]+\\.l2\\.ldscore\\.gz$", full.names = TRUE)
LD_FULL <- read_ld(l_files)
M_FULL  <- sum(sapply(list.files("ref/LDscore", pattern = "\\.l2\\.M_5_50$", full.names = TRUE), function(f) as.numeric(scan(f, quiet = TRUE))[1]))
cat("Full-1000G reference: SNPs with LD scores", nrow(LD_FULL), "; M (MAF>5%) =", M_FULL, "\n")   # expect M = 5,961,159

## ---------------------------------------------------------------------------
## 3. FUNCTIONS (R ports of the official LDSC code)
## ---------------------------------------------------------------------------

## 3.1 Allele logic, identical to ldscore/sumstats.py -------------------------
##     VALID  = strand-unambiguous SNP allele pairs (A/T and C/G pairs are excluded)
##     MATCH  = two SNPs have the same alleles (allowing strand flip and swap)
##     FLIP   = their alleles are swapped, so the sign of the z-score must be reversed
COMP <- c(A = "T", T = "A", C = "G", G = "C")
VALID_PAIRS <- c("AC", "AG", "CA", "CT", "GA", "GT", "TC", "TG")
.g  <- expand.grid(p = VALID_PAIRS, q = VALID_PAIRS, stringsAsFactors = FALSE)
.ch <- function(s, i) substr(s, i, i)
.p1 <- .ch(.g$p, 1); .p2 <- .ch(.g$p, 2); .q1 <- .ch(.g$q, 1); .q2 <- .ch(.g$q, 2)
.cq1 <- unname(COMP[.q1]); .cq2 <- unname(COMP[.q2])
.same   <- .p1 == .q1  & .p2 == .q2           # same strand, same order
.strand <- .p1 == .cq1 & .p2 == .cq2          # other strand, same order
.swap   <- .p1 == .q2  & .p2 == .q1           # same strand, swapped order
.swapst <- .p1 == .cq2 & .p2 == .cq1          # other strand, swapped order
.any <- .same | .strand | .swap | .swapst
MATCH_KEYS <- paste0(.g$p, .g$q)[.any]
FLIP_FLAG  <- (.swap | .swapst)[.any]
allele_check <- function(key) { i <- match(key, MATCH_KEYS); list(valid = !is.na(i), flip = FLIP_FLAG[i]) }

## 3.2 munge_sumstats_R: port of munge_sumstats.py ----------------------------
##     Input table columns: SNP, A1 (effect allele), A2, BETA, P, FRQ, optional INFO and N.
##     Steps and their order are those of the official program:
##       drop missing -> keep HapMap3 SNPs -> INFO >= 0.9 -> MAF > 0.01 -> 0 < P <= 1
##       -> unambiguous SNP alleles -> drop duplicated rs numbers -> N -> Z = sign(BETA) * sqrt(qchisq(P))
##       -> match alleles to the HapMap3 list -> Z written with 3 decimals (as the official program does)
munge_sumstats_R <- function(d, name, Nfixed = NULL, Ncas = NULL, Ncon = NULL, info_min = 0.9, maf_min = 0.01) {
  lg <- character(); L <- function(...) { m <- paste0(...); lg <<- c(lg, m); cat(m, "\n") }
  L("---- munge: ", name, " ----")
  L("Read ", format(nrow(d), big.mark = ","), " SNPs from the prepared table.")
  d <- copy(d)
  keepcols <- setdiff(names(d), "INFO")                       # INFO may be missing, everything else may not
  n0 <- nrow(d); d <- d[complete.cases(d[, ..keepcols])]; r_na <- n0 - nrow(d)
  n0 <- nrow(d); d <- d[SNP %chin% HM3$SNP];               r_merge <- n0 - nrow(d)
  r_info <- 0L
  if ("INFO" %in% names(d)) { n0 <- nrow(d); d <- d[!is.na(INFO) & INFO >= info_min]; r_info <- n0 - nrow(d) }
  n0 <- nrow(d); d <- d[FRQ >= 0 & FRQ <= 1 & pmin(FRQ, 1 - FRQ) > maf_min];            r_frq <- n0 - nrow(d)
  n0 <- nrow(d); d <- d[P > 0 & P <= 1];                                                   r_p   <- n0 - nrow(d)
  d[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]
  n0 <- nrow(d); d <- d[paste0(A1, A2) %chin% VALID_PAIRS];                                r_al  <- n0 - nrow(d)
  L("Removed ", format(r_merge, big.mark = ","), " SNPs not in the HapMap3 allele list.")
  L("Removed ", r_na, " SNPs with missing values.")
  L("Removed ", r_info, " SNPs with INFO <= ", info_min, ".")
  L("Removed ", r_frq, " SNPs with MAF <= ", maf_min, ".")
  L("Removed ", r_p, " SNPs with out-of-bounds p-values.")
  L("Removed ", r_al, " variants that were not SNPs or were strand-ambiguous.")
  L(format(nrow(d), big.mark = ","), " SNPs remain.")
  n0 <- nrow(d); d <- unique(d, by = "SNP")                       # keeps the first occurrence, like drop_duplicates
  L("Removed ", n0 - nrow(d), " SNPs with duplicated rs numbers (", format(nrow(d), big.mark = ","), " SNPs remain).")
  if ("N" %in% names(d)) {                                         # per-SNP sample size supplied (GLGC)
    n_min <- as.numeric(quantile(d$N, 0.9, type = 7)) / 1.5
    n0 <- nrow(d); d <- d[N >= n_min]
    L("Removed ", n0 - nrow(d), " SNPs with N < ", n_min, " (", format(nrow(d), big.mark = ","), " SNPs remain).")
  } else if (!is.null(Nfixed)) { d[, N := Nfixed]; L("Using N = ", Nfixed)
  } else { d[, N := Ncas + Ncon]; L("Using N_cas = ", Ncas, "; N_con = ", Ncon) }
  d[, Z := sqrt(qchisq(P, df = 1, lower.tail = FALSE))]          # P -> |Z|
  d[, Z := fifelse(BETA < 0, -Z, Z)]                              # sign from the effect size
  L("Median value of BETA was ", signif(median(d$BETA), 6), ", which seems sensible.")
  out <- merge(HM3[, .(SNP, MA)], d[, .(SNP, A1, A2, N, Z)], by = "SNP", all.x = TRUE, sort = FALSE)   # HapMap3 order
  chk <- allele_check(paste0(out$A1, out$A2, out$MA))
  has <- !is.na(out$A1); bad <- has & !chk$valid
  L("Removed ", sum(bad), " SNPs whose alleles did not match the HapMap3 list (", format(sum(has) - sum(bad), big.mark = ","), " SNPs remain).")
  out[bad, `:=`(A1 = NA_character_, A2 = NA_character_, N = NA_real_, Z = NA_real_)]
  chisq <- out$Z^2                                                 # summary statistics are taken before rounding, as in the official program
  out[, `:=`(Z = as.numeric(sprintf("%.3f", Z)), N = as.numeric(sprintf("%.3f", N)))]    # official files keep 3 decimals
  out[, MA := NULL]
  L("Writing summary statistics for ", format(nrow(out), big.mark = ","), " SNPs (", format(sum(!is.na(out$Z)), big.mark = ","), " with nonmissing beta).")
  L("Mean chi^2 = ", round(mean(chisq, na.rm = TRUE), 3))
  L("Lambda GC = ", round(median(chisq, na.rm = TRUE) / 0.4549, 3))
  L("Max chi^2 = ", round(max(chisq, na.rm = TRUE), 3))
  L(sum(chisq > 29, na.rm = TRUE), " Genome-wide significant SNPs (some may have been removed by filtering).")
  writeLines(lg, file.path("logs", paste0("munge_", name, ".log")))
  fwrite(out[, .(SNP, N, Z, A1, A2)], file.path("munged", paste0(name, ".sumstats.gz")), sep = "\t", na = "NA")
  setkey(out, NULL); out
}

## 3.3 Reading and preparing each GWAS (what the prepared-table step did) -------
##     Every file is reduced to: SNP, A1 (effect allele), A2, BETA, P, FRQ [, INFO, N]
##     P is recomputed from BETA/SE at full precision, because some catalogue files
##     store P with only 2 significant figures: P = 2*pnorm(-|BETA/SE|).
read_gwas <- function(name) {
  f <- file.path("gwas", gwas_dl[[name]]$file)
  if (name %in% c("RC", "Ghodsian", "CAD")) {
    d <- fread(f, select = c("variant_id", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error"),
               colClasses = c(variant_id = "character", effect_allele = "character", other_allele = "character"), showProgress = FALSE)
    setnames(d, c("SNP", "A1", "A2", "FRQ", "BETA", "SE"))
    d <- d[grepl("^rs[0-9]+$", SNP)]
  } else if (name == "Fairfield") {            # odds ratios with 95% CI: BETA = ln(OR), SE from the CI width
    d <- fread(f, select = c("variant_id", "effect_allele", "other_allele", "effect_allele_frequency", "odds_ratio", "ci_lower", "ci_upper", "info"),
               colClasses = c(variant_id = "character", effect_allele = "character", other_allele = "character"), showProgress = FALSE)
    setnames(d, c("SNP", "A1", "A2", "FRQ", "OR", "LL", "UL", "INFO"))
    d <- d[grepl("^rs[0-9]+$", SNP)]
    d <- d[OR > 0 & LL > 0 & UL > 0]
    d[, `:=`(BETA = log(OR), SE = (log(UL) - log(LL)) / (2 * 1.959964))]
    d[, c("OR", "LL", "UL") := NULL]
  } else {                                     # GLGC: SNP, alleles (lower case), beta, se, N, allele frequency
    d <- fread(f, select = c("rsid", "A1", "A2", "beta", "se", "N", "Freq.A1.1000G.EUR"), colClasses = c(rsid = "character", A1 = "character", A2 = "character"), showProgress = FALSE)
    setnames(d, c("SNP", "A1", "A2", "BETA", "SE", "N", "FRQ"))
    d <- d[complete.cases(d)]
  }
  d <- d[!is.na(BETA) & !is.na(SE) & !is.na(A1) & !is.na(A2) & SE > 0]
  d[, P := 2 * pnorm(-abs(BETA / SE))]
  d <- d[P > 0]
  d[, SE := NULL]
  setcolorder(d, intersect(c("SNP", "A1", "A2", "BETA", "P", "FRQ", "INFO", "N"), names(d)))
  d[]
}

## 3.5 Block jackknife (port of ldscore/jackknife.py) --------------------------
##     x: n x p design matrix, y: length-n response, sep: block boundaries (0 ... n).
jk_separators <- function(N, nb) floor(seq(0, N, length.out = nb + 1))      # np.floor(np.linspace(0, N, nb + 1))
jk_pseudovalues <- function(delete, est) {                                 # nb * est - (nb - 1) * delete
  nb <- nrow(delete); nb * matrix(est, nb, ncol(delete), byrow = TRUE) - (nb - 1) * delete
}
jk_summary <- function(pseudo) {                                           # jackknife estimate and covariance
  nb <- nrow(pseudo)
  list(est = colMeans(pseudo), cov = matrix(cov(pseudo) / nb, ncol(pseudo), ncol(pseudo)))
}
lstsq_jackknife <- function(x, y, sep) {
  nb <- length(sep) - 1; p <- ncol(x)
  xty <- matrix(0, nb, p); xtx <- array(0, c(nb, p, p))
  for (i in seq_len(nb)) {                                                 # block values of X'Y and X'X
    r <- (sep[i] + 1):sep[i + 1]; xb <- x[r, , drop = FALSE]
    xty[i, ] <- crossprod(xb, y[r]); xtx[i, , ] <- crossprod(xb)
  }
  xty_t <- colSums(xty); xtx_t <- apply(xtx, c(2, 3), sum)
  est <- solve(xtx_t, xty_t)                                               # whole-data estimate
  del <- matrix(0, nb, p)
  for (j in seq_len(nb)) del[j, ] <- solve(xtx_t - matrix(xtx[j, , ], p, p), xty_t - xty[j, ])   # leave-one-block-out
  s <- jk_summary(jk_pseudovalues(del, est))
  list(est = est, delete = del, cov = s$cov, se = sqrt(diag(s$cov)), separators = sep)
}

## 3.6 Iteratively re-weighted least squares (port of ldscore/irwls.py) --------
irwls_wls <- function(x, y, w) {                                           # weighted least squares, weights normalised to sum 1
  if (any(w <= 0)) stop("Weights must be > 0")
  w <- w / sum(w)
  qr.coef(qr(x * w, tol = 1e-12), y * w)
}
irwls <- function(x, y, update_func, nb, w, separators = NULL) {
  w <- sqrt(w)
  for (i in 1:2) w <- sqrt(update_func(irwls_wls(x, y, w)))                # two re-weighting passes, as in LDSC
  w <- w / sum(w)
  sep <- if (is.null(separators)) jk_separators(nrow(x), nb) else separators
  lstsq_jackknife(x * w, y * w, sep)
}

## 3.7 Regression weights ------------------------------------------------------
hsq_weights <- function(ld, w_ld, N, M, hsq, intercept = 1) {
  hsq <- min(max(hsq, 0), 1); ld <- pmax(ld, 1); w_ld <- pmax(w_ld, 1)
  c <- hsq * N / M
  (1 / (2 * (intercept + c * ld)^2)) * (1 / w_ld)
}
gencov_weights <- function(ld, w_ld, N1, N2, M, h1, h2, rho_g, int_gencov = 0, int_h1 = 1, int_h2 = 1) {
  h1 <- min(max(h1, 0), 1); h2 <- min(max(h2, 0), 1); rho_g <- min(max(rho_g, -1), 1)
  ld <- pmax(ld, 1); w_ld <- pmax(w_ld, 1)
  a <- N1 * h1 * ld / M + int_h1
  b <- N2 * h2 * ld / M + int_h2
  cc <- sqrt(N1 * N2) * rho_g * ld / M + int_gencov
  (1 / (a * b + cc^2)) * (1 / w_ld)
}

## 3.8 LD score regression core (port of LD_Score_Regression.__init__, 1 annotation)
##     y: chi^2 (h2) or z1*z2 (gencov);  ld: reference LD score;  w_ld: regression-SNP LD score;
##     N: sample size (sqrt(N1*N2) for gencov);  M: number of SNPs the LD score is summed over.
##     intercept = NULL -> free intercept; step1_ii -> two-step estimator (intercept from SNPs flagged TRUE).
ldsc_core <- function(y, ld, w_ld, N, M, nb, intercept = NULL, step1_ii = NULL, weight_fn, null_int) {
  Nbar <- mean(N); free <- is.null(intercept)
  i0 <- if (free) null_int else intercept
  tot_agg <- M * (mean(y) - i0) / mean(ld * N)                            # crude starting value of h2 (or gencov)
  initial_w <- weight_fn(ld, w_ld, NULL, tot_agg, i0)
  x <- ld * N / Nbar                                                       # keeps the condition number low
  if (free) { X <- cbind(x, 1); yp <- y } else { X <- matrix(x, ncol = 1); yp <- y - intercept }
  if (!is.null(step1_ii) && !free) stop("twostep is not compatible with a constrained intercept.")
  if (is.null(step1_ii)) {                                                 # one-step estimator
    upd <- function(a) weight_fn(ld, w_ld, NULL, M * a[1] / Nbar, if (free) a[2] else intercept)
    jk <- irwls(X, yp, upd, nb, initial_w)
  } else {                                                                 # two-step estimator
    ii <- step1_ii; n1 <- sum(ii)
    x1 <- X[ii, , drop = FALSE]; yp1 <- yp[ii]; w1 <- w_ld[ii]; iw1 <- initial_w[ii]
    # (as in LDSC, the step-1 weights use the N-scaled LD column of x1)
    upd1 <- function(a) weight_fn(x1[, 1], w1, ii, M * a[1] / Nbar, a[2])
    s1 <- irwls(x1, yp1, upd1, nb, iw1)
    s1_int <- s1$est[2]
    yp2 <- yp - s1_int; x2 <- X[, 1, drop = FALSE]
    upd2 <- function(a) weight_fn(ld, w_ld, NULL, M * a[1] / Nbar, s1_int)
    map0 <- which(ii) - 1                                                  # unmask the step-1 block boundaries
    sep <- c(0, map0[s1$separators[2:(length(s1$separators) - 1)] + 1], length(ii))
    s2 <- irwls(x2, yp2, upd2, nb, initial_w, separators = sep)
    cc <- sum(initial_w * x2) / sum(initial_w * x2^2)
    est <- c(s2$est, s1_int)
    del <- cbind(s2$delete[, 1] - cc * (s1$delete[, 2] - s1_int), s1$delete[, 2])
    sm <- jk_summary(jk_pseudovalues(del, est))
    jk <- list(est = est, delete = del, cov = sm$cov, se = sqrt(diag(sm$cov)))
  }
  coef <- jk$est[1] / Nbar
  coef_cov <- jk$cov[1, 1] / Nbar^2
  out <- list(coef = coef, tot = M * coef, tot_se = M * sqrt(coef_cov),
              tot_delete = jk$delete[, 1] * M / Nbar,
              intercept = if (free) jk$est[2] else intercept,
              intercept_se = if (free) jk$se[2] else NA_real_,
              nb = nb, n_snp = length(y))
  out
}

## 3.9 Heritability (class Hsq) --------------------------------------------------
ldsc_hsq <- function(chisq, ld, w_ld, N, M, nb, intercept = NULL, twostep = NULL) {
  wf <- function(ldv, wldv, ii, h, int) hsq_weights(ldv, wldv, if (is.null(ii)) N else N[ii], M, h, int)
  ii <- if (!is.null(twostep)) chisq < twostep else NULL
  r <- ldsc_core(chisq, ld, w_ld, N, M, nb, intercept, ii, wf, null_int = 1)
  r$mean_chisq <- mean(chisq); r$lambda_gc <- median(chisq) / 0.4549
  r$ratio <- if (is.null(intercept) && r$mean_chisq > 1) (r$intercept - 1) / (r$mean_chisq - 1) else NA_real_
  r$ratio_se <- if (is.null(intercept) && r$mean_chisq > 1) r$intercept_se / (r$mean_chisq - 1) else NA_real_
  r
}

## 3.10 Genetic covariance (class Gencov) ----------------------------------------
ldsc_gencov <- function(z1, z2, ld, w_ld, N1, N2, M, h1, h2, int1, int2, nb, intercept_gencov = NULL, twostep = NULL) {
  wf <- function(ldv, wldv, ii, rho, int) {
    n1 <- if (is.null(ii)) N1 else N1[ii]; n2 <- if (is.null(ii)) N2 else N2[ii]
    gencov_weights(ldv, wldv, n1, n2, M, h1, h2, rho, int, int1, int2)
  }
  ii <- if (!is.null(twostep)) (z1^2 < twostep & z2^2 < twostep) else NULL
  r <- ldsc_core(z1 * z2, ld, w_ld, sqrt(N1 * N2), M, nb, intercept_gencov, ii, wf, null_int = 0)
  r$mean_z1z2 <- mean(z1 * z2)
  r
}

## 3.11 Genetic correlation (class RG) -------------------------------------------
##      rg = gencov / sqrt(h2_1 * h2_2); SE by a ratio block jackknife (pseudovalues from delete values).
ldsc_rg <- function(z1, z2, ld, w_ld, N1, N2, M, nb = 200, twostep = NULL, int_h1 = NULL, int_h2 = NULL, int_gencov = NULL) {
  nb <- min(nb, length(z1))
  h1 <- ldsc_hsq(z1^2, ld, w_ld, N1, M, nb, int_h1, twostep)
  h2 <- ldsc_hsq(z2^2, ld, w_ld, N2, M, nb, int_h2, twostep)
  g  <- ldsc_gencov(z1, z2, ld, w_ld, N1, N2, M, h1$tot, h2$tot, h1$intercept, h2$intercept, nb, int_gencov, twostep)
  res <- list(h1 = h1, h2 = h2, g = g, rg = NA_real_, rg_se = NA_real_, z = NA_real_, p = NA_real_)
  if (h1$tot > 0 && h2$tot > 0) {
    rg <- g$tot / sqrt(h1$tot * h2$tot)
    den <- sqrt(h1$tot_delete * h2$tot_delete)
    pseudo <- nb * rg - (nb - 1) * g$tot_delete / den
    se <- sqrt(var(pseudo) / nb)
    res$rg <- rg; res$rg_se <- se; res$z <- rg / se; res$p <- pchisq((rg / se)^2, 1, lower.tail = FALSE)
  }
  res
}

## 3.12 Observed -> liability scale (h2_obs_to_liab, gencov_obs_to_liab) ----------
h2_obs_to_liab <- function(P, K) {                                         # multiplicative factor
  if (is.na(P) || is.na(K)) return(1)
  th <- qnorm(K, lower.tail = FALSE)
  K^2 * (1 - K)^2 / (P * (1 - P) * dnorm(th)^2)
}

## 3.13 Reading munged files, merging with LD scores and aligning alleles -------
read_munged <- function(name, dir = "munged") {
  d <- fread(file.path(dir, paste0(name, ".sumstats.gz")), showProgress = FALSE, na.strings = c("NA", "nan", "NaN", ""))
  d <- d[!is.na(Z) & !is.na(N)]
  d[, `:=`(N = as.numeric(N), Z = as.numeric(Z))]                         # integers would overflow in N1 * N2
  unique(d, by = "SNP")                                                    # drop duplicated rs numbers (first one kept)
}

## Prepare the regression data set exactly as LDSC does: SNPs in trait 1 (non-missing), merged with the reference
## LD scores (order of the LD-score file), then with the regression-SNP LD scores, then with trait 2.
## Returns a data.table with columns SNP, L2 (reference LD), W (regression LD), N1, Z1, N2, Z2 (aligned to trait 1 alleles).
merge_for_rg <- function(t1, t2, ref_ld, w_ld, verbose = TRUE, l1 = NULL, l2 = NULL) {
  a <- if (is.character(t1)) read_munged(t1) else t1
  b <- if (is.character(t2)) read_munged(t2) else t2
  i <- match(ref_ld$SNP, a$SNP); k <- !is.na(i)
  d <- data.table(SNP = ref_ld$SNP[k], L2 = ref_ld$L2[k], N1 = a$N[i[k]], Z1 = a$Z[i[k]], A1 = a$A1[i[k]], A2 = a$A2[i[k]])
  j <- match(d$SNP, w_ld$SNP); d <- d[!is.na(j)]; d[, W := w_ld$L2[j[!is.na(j)]]]
  i2 <- match(d$SNP, b$SNP); d <- d[!is.na(i2)]; i2 <- i2[!is.na(i2)]
  d[, `:=`(N2 = b$N[i2], Z2 = b$Z[i2], A1x = b$A1[i2], A2x = b$A2[i2])]
  d <- d[complete.cases(d)]
  chk <- allele_check(paste0(d$A1, d$A2, d$A1x, d$A2x))                   # same SNP alleles after strand/swap?
  d <- d[chk$valid]; flip <- chk$flip[chk$valid]
  d[, Z2 := fifelse(flip, -Z2, Z2)]                                        # swapped alleles -> reverse sign of z2
  d[, c("A1", "A2", "A1x", "A2x") := NULL]
  if (verbose) cat(sprintf("  %s x %s: %s SNPs after merging with LD scores and aligning alleles\n",
                           if (!is.null(l1)) l1 else t1, if (!is.null(l2)) l2 else t2, format(nrow(d), big.mark = ",")))
  d[]
}

## Run rg for one trait pair; one tidy row out. P/K = sample and population prevalence of trait 2 (NA for quantitative).
run_rg_pair <- function(t1, t2, ref_ld, w_ld, M, twostep = NULL, chisq_max = NULL, P2 = NA, K2 = NA, label = "", n1 = NULL, n2 = NULL) {
  d <- merge_for_rg(t1, t2, ref_ld, w_ld, l1 = n1, l2 = n2)
  if (!is.null(chisq_max)) d <- d[Z1^2 * Z2^2 < chisq_max^2]
  r <- ldsc_rg(d$Z1, d$Z2, d$L2, d$W, d$N1, d$N2, M, nb = 200, twostep = twostep)
  cl <- h2_obs_to_liab(P2, K2)
  data.table(run = label, trait1 = if (!is.null(n1)) n1 else t1, trait2 = if (!is.null(n2)) n2 else t2,
             n_snps = nrow(d), rg = r$rg, se = r$rg_se, z = r$z, p = r$p,
             h2_trait1 = r$h1$tot, h2_trait1_se = r$h1$tot_se, h2_trait2 = cl * r$h2$tot, h2_trait2_se = cl * r$h2$tot_se,
             h2_trait2_scale = if (is.na(P2)) "Observed" else "Liability",
             h2_int_trait1 = r$h1$intercept, h2_int_trait2 = r$h2$intercept, gcov_int = r$g$intercept, gcov_int_se = r$g$intercept_se)
}

## Heritability of one trait (LDSC --h2): two-step estimator, cutoff 30, free intercept (the LDSC default).
run_h2 <- function(name, ref_ld, w_ld, M, P = NA, K = NA) {
  a <- read_munged(name)
  i <- match(ref_ld$SNP, a$SNP); k <- !is.na(i)
  d <- data.table(SNP = ref_ld$SNP[k], L2 = ref_ld$L2[k], N = a$N[i[k]], Z = a$Z[i[k]])
  j <- match(d$SNP, w_ld$SNP); d <- d[!is.na(j)]; d[, W := w_ld$L2[j[!is.na(j)]]]
  r <- ldsc_hsq(d$Z^2, d$L2, d$W, d$N, M, min(200, nrow(d)), NULL, twostep = 30)
  cl <- h2_obs_to_liab(P, K)
  data.table(trait = name, scale = if (is.na(P)) "Observed" else "Liability", n_snps = nrow(d),
             h2 = cl * r$tot, h2_se = cl * r$tot_se, z_h2 = r$tot / r$tot_se,
             mean_chi2 = r$mean_chisq, lambda_gc = r$lambda_gc, intercept = r$intercept, intercept_se = r$intercept_se)
}

## ---------------------------------------------------------------------------
## 4. PREPARE AND MUNGE EVERY GWAS  (what 01_prep_sumstats.py and 02_munge.sh did)
##    Output: munged/<trait>.sumstats.gz  (SNP, N, Z, A1, A2 on the 1,187,349 HapMap3 SNPs)
##    Logs:   logs/munge_<trait>.log
## ---------------------------------------------------------------------------
cat("\n==== 4. MUNGING ====\n")
MUNGE <- list(                                  # sample sizes are those of the original studies
  RC        = list(Nfixed = 115082),
  Ghodsian  = list(Ncas = 8434, Ncon = 770180),
  Fairfield = list(Ncas = 4761, Ncon = 373227),
  CAD       = list(Ncas = 60801, Ncon = 123504),
  GLGC_TG   = list(),                           # per-SNP N is in the file
  GLGC_LDL  = list()
)
for (nm in names(MUNGE)) {
  if (file.exists(file.path("munged", paste0(nm, ".sumstats.gz"))) && !isTRUE(FORCE_REMUNGE)) { cat("munged/", nm, " exists, skipped (set FORCE_REMUNGE <- TRUE to redo)\n", sep = ""); next }
  d <- read_gwas(nm)
  invisible(munge_sumstats_R(d, nm, Nfixed = MUNGE[[nm]]$Nfixed, Ncas = MUNGE[[nm]]$Ncas, Ncon = MUNGE[[nm]]$Ncon))
  rm(d); invisible(gc())
}

## ---------------------------------------------------------------------------
## 5. PRIMARY ANALYSIS: SNP HERITABILITY (h2) OF EACH TRAIT
##    Two-step estimator (cutoff 30), free intercept, 200 jackknife blocks (LDSC --h2 default).
##    Case-control traits are reported on the liability scale (sample prevalence = cases / (cases + controls)).
## ---------------------------------------------------------------------------
cat("\n==== 5. SNP HERITABILITY ====\n")
PREV <- list(Ghodsian  = c(P = 8434  / (8434 + 770180),  K = 8434  / (8434 + 770180)),     # 0.010832
             Fairfield = c(P = 4761  / (4761 + 373227),  K = 4761  / (4761 + 373227)),     # 0.012595
             CAD       = c(P = 60801 / (60801 + 123504), K = 0.05))                         # 0.329888, population 0.05
getP <- function(nm) if (nm %in% names(PREV)) PREV[[nm]][["P"]] else NA
getK <- function(nm) if (nm %in% names(PREV)) PREV[[nm]][["K"]] else NA
H2 <- rbindlist(lapply(c("RC", "Ghodsian", "Fairfield", "CAD"), function(nm) run_h2(nm, LD_HM3, LD_HM3, M_HM3, getP(nm), getK(nm))))
print(H2[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
fwrite(H2, "results/ldsc_h2_R.csv")

## ---------------------------------------------------------------------------
## 6. PRIMARY ANALYSIS: GENETIC CORRELATION OF RC WITH NAFLD (two GWAS) AND WITH CAD
##    One-step estimator, free intercepts (LDSC --rg default), HapMap3 LD scores, M = 1,187,349.
## ---------------------------------------------------------------------------
cat("\n==== 6. PRIMARY rg ====\n")
rg_run <- function(t1, t2s, label, ref = LD_HM3, w = LD_HM3, M = M_HM3, prev = TRUE, ...) {   # prev = FALSE: h2 of trait 2 on the observed scale
  rbindlist(lapply(t2s, function(t2) run_rg_pair(t1, t2, ref, w, M, P2 = if (prev) getP(t2) else NA, K2 = if (prev) getK(t2) else NA, label = label, ...)))
}
RES <- list()
RES$primary   <- rg_run("RC", c("Ghodsian", "Fairfield", "CAD"), "rg_RC_vs_all")
RES$nafld_cc  <- rg_run("Ghodsian", "Fairfield", "rg_Ghodsian_vs_Fairfield")      # concordance of the two NAFLD GWAS
show <- function(x) print(x[, .(run, trait1, trait2, n_snps, rg = round(rg, 4), se = round(se, 4), z = round(z, 3), p = formatC(p, format = "e", digits = 2),
                                 h2_trait2 = round(h2_trait2, 4), h2_int_trait2 = round(h2_int_trait2, 4), gcov_int = round(gcov_int, 4), gcov_int_se = round(gcov_int_se, 4))])
show(RES$primary); show(RES$nafld_cc)

## ---------------------------------------------------------------------------
## 7. SENSITIVITY ANALYSES OF THE RC rg (RC is always trait 1)
##    A: no chi2 truncation (identical to the default, because rg has no default truncation)
##    B: reference LD scores from all ~6 million 1000G SNPs (M = 5,961,159), HapMap3 LD scores as weights
##    C: two-step estimator (cutoff 30) for the intercepts
## ---------------------------------------------------------------------------
cat("\n==== 7. SENSITIVITY ====\n")
OUT3 <- c("CAD", "Ghodsian", "Fairfield")
RES$sensA <- rg_run("RC", OUT3, "sens_A_chisqmax_off", chisq_max = 100000)
RES$sensB <- rg_run("RC", OUT3, "sens_B_full1000G_ref", ref = LD_FULL, w = LD_HM3, M = M_FULL)
RES$sensC <- rg_run("RC", OUT3, "sens_C_twostep", twostep = 30)
show(RES$sensA); show(RES$sensB); show(RES$sensC)

## ---------------------------------------------------------------------------
## 8. EXTERNAL BENCHMARK: GLGC LIPID GWAS (Willer 2013) vs CAD, RC and NAFLD
##    Literature expectation: rg(TG, CAD) and rg(LDL, CAD) clearly positive. A pipeline that gives that is working.
## ---------------------------------------------------------------------------
cat("\n==== 8. GLGC BENCHMARK ====\n")
RES$glgc_tg  <- rg_run("GLGC_TG",  c("CAD", "GLGC_LDL", "RC", "Ghodsian", "Fairfield"), "rg_GLGC_TG_vs_all", prev = FALSE)
RES$glgc_ldl <- rg_run("GLGC_LDL", c("CAD", "RC", "Ghodsian", "Fairfield"),             "rg_GLGC_LDL_vs_all", prev = FALSE)
show(RES$glgc_tg); show(RES$glgc_ldl)

## ---------------------------------------------------------------------------
## 9. ALIGNMENT SELF-TESTS: the same GWAS compared with itself after (i) swapping the two alleles and reversing the
##    sign of z, (ii) moving to the other strand, (iii) both, and (iv) reversing z only.
##    Expected rg = +1, +1, +1, -1. Proves that allele alignment and sign handling are right.
## ---------------------------------------------------------------------------
cat("\n==== 9. ALIGNMENT SELF-TESTS ====\n")
rc <- read_munged("RC")
cmp <- function(x) unname(COMP[x])
variants <- list(
  swapped    = copy(rc)[, `:=`(A1 = rc$A2, A2 = rc$A1, Z = -rc$Z)],
  strand     = copy(rc)[, `:=`(A1 = cmp(rc$A1), A2 = cmp(rc$A2))],
  swapstrand = copy(rc)[, `:=`(A1 = cmp(rc$A2), A2 = cmp(rc$A1), Z = -rc$Z)],
  negated    = copy(rc)[, Z := -rc$Z])
RES$selftest <- rbindlist(lapply(names(variants), function(v)
  run_rg_pair(rc, variants[[v]], LD_HM3, LD_HM3, M_HM3, label = "selftest", n1 = "RC", n2 = paste0("RC_", v))))
show(RES$selftest)
rm(rc, variants)

## ---------------------------------------------------------------------------
## 10. BULK-GENOME vs GENOME-WIDE-SIGNIFICANT LOCI
##     Is the near-zero rg(RC, CAD) driven by the large-effect loci or by the polygenic bulk?
##     "loci" = HapMap3 SNPs within 1 Mb of an RC genome-wide significant SNP (|z| > 5.45, P < 5e-8);
##     "bulk" = all other HapMap3 SNPs. The same SNP set is applied to every trait.
## ---------------------------------------------------------------------------
cat("\n==== 10. BULK vs LOCI ====\n")
rc <- read_munged("RC")
gws <- merge(rc[abs(Z) > 5.45, .(SNP)], LD_HM3[, .(SNP, CHR, BP)], by = "SNP")
near <- logical(nrow(LD_HM3))
for (ch in unique(gws$CHR)) {
  idx <- which(LD_HM3$CHR == ch); bp <- LD_HM3$BP[idx]; pos <- sort(gws$BP[gws$CHR == ch])
  k <- findInterval(bp, pos)                                               # nearest genome-wide-significant SNP, both sides
  lo <- abs(bp - pos[pmax(k, 1)]); hi <- abs(bp - pos[pmin(k + 1, length(pos))])
  near[idx] <- pmin(lo, hi) <= 1e6
}
cat("HapMap3 SNPs:", nrow(LD_HM3), "| within 1 Mb of an RC GWS SNP:", sum(near), "| bulk:", sum(!near), "| RC GWS SNPs:", nrow(gws), "\n")
SETS <- list(bulk = LD_HM3$SNP[!near], loci = LD_HM3$SNP[near])
sub <- function(nm, k) { d <- read_munged(nm); d[SNP %chin% SETS[[k]]] }
RES$bulkloci <- rbindlist(lapply(names(SETS), function(k) {
  r1 <- sub("RC", k); r2 <- sub("GLGC_LDL", k)
  a <- rbindlist(lapply(c("CAD", "GLGC_LDL", "GLGC_TG", "Ghodsian", "Fairfield"), function(t2)
    run_rg_pair(r1, sub(t2, k), LD_HM3, LD_HM3, M_HM3, label = paste0("rg_RC_vs_all_", k), n1 = paste0("RC_", k), n2 = paste0(t2, "_", k))))
  b <- rbindlist(lapply(c("CAD", "GLGC_TG"), function(t2)
    run_rg_pair(r2, sub(t2, k), LD_HM3, LD_HM3, M_HM3, label = paste0("rg_LDL_vs_CAD_TG_", k), n1 = paste0("GLGC_LDL_", k), n2 = paste0(t2, "_", k))))
  rbind(a, b)
}))
show(RES$bulkloci)
rm(rc)

## ---------------------------------------------------------------------------
## 11. INDEPENDENT CROSS-CHECK: a deliberately simple re-implementation of the LD score regression
##     (one weighted-regression pass, leave-one-chromosome-out jackknife). If this agrees with the full port above,
##     the numbers do not depend on the details of the IRWLS / 200-block machinery.
##       E[chi2_j]     = N * h2 * l_j / M + a        E[z1_j z2_j] = sqrt(N1 N2) * rho_g * l_j / M + b
## ---------------------------------------------------------------------------
cat("\n==== 11. INDEPENDENT CROSS-CHECK ====\n")
simple_wls <- function(y, x, w) { X <- cbind(1, x); solve(crossprod(X, w * X), crossprod(X, w * y)) }
simple_h2 <- function(y, x, N, M) {
  b <- simple_wls(y, x, 1 / x); hh <- max(b[2] * M / N, 0.001)
  w <- 1 / (x * (1 + N * hh * x / M)^2); list(b = simple_wls(y, x, w), w = w)
}
independent_rg <- function(a, b) {
  A <- read_munged(a); B <- read_munged(b)
  d <- merge_for_rg(A, B, LD_HM3, LD_HM3, verbose = FALSE)
  d <- merge(d, LD_HM3[, .(SNP, CHR)], by = "SNP", sort = FALSE)
  d <- d[Z1^2 < max(80, 0.001 * max(A$N)) & Z2^2 < max(80, 0.001 * max(B$N))]      # chi2 truncation used in this check only
  N1 <- mean(d$N1); N2 <- mean(d$N2); x <- d$L2; M <- M_HM3
  f1 <- simple_h2(d$Z1^2, x, N1, M); f2 <- simple_h2(d$Z2^2, x, N2, M)
  h21 <- f1$b[2] * M / N1; h22 <- f2$b[2] * M / N2
  y <- d$Z1 * d$Z2; w <- 1 / (x * sqrt((1 + N1 * h21 * x / M) * (1 + N2 * h22 * x / M)))
  bg <- simple_wls(y, x, w); gcov <- bg[2] * M / sqrt(N1 * N2)
  rgs <- sapply(1:22, function(ch) { m <- d$CHR != ch
    (simple_wls(y[m], x[m], w[m])[2] * M / sqrt(N1 * N2)) /
      sqrt(max(simple_wls(d$Z1[m]^2, x[m], f1$w[m])[2] * M / N1, 1e-6) * max(simple_wls(d$Z2[m]^2, x[m], f2$w[m])[2] * M / N2, 1e-6)) })
  rg <- gcov / sqrt(h21 * h22); pseudo <- 22 * rg - 21 * rgs
  data.table(pair = paste(a, "vs", b), n_snps = nrow(d), h2_1 = h21, h2_2 = h22, gcov = gcov, rg = rg, se = sd(pseudo) / sqrt(22),
             int1 = f1$b[1], int2 = f2$b[1], int_gencov = bg[1])
}
RES$independent <- rbindlist(list(independent_rg("RC", "RC"), independent_rg("RC", "CAD"), independent_rg("RC", "Ghodsian"),
                                  independent_rg("RC", "Fairfield"), independent_rg("Ghodsian", "Fairfield")))
print(RES$independent[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])

## ---------------------------------------------------------------------------
## 12. INSTRUMENT COMPOSITION: do the 48 MR instruments for RC mostly act through LDL cholesterol or triglycerides?
##     The 48 SNPs (RC-raising allele = ea, effect bx) are those of the primary MR analysis (clumped, harmonised with
##     the Ghodsian NAFLD GWAS). For each SNP the GLGC LDL and TG z-scores are aligned to the RC-raising allele.
## ---------------------------------------------------------------------------
cat("\n==== 12. INSTRUMENT COMPOSITION ====\n")
INSTR <- data.frame(
  SNP = c("rs72654473", "rs12151108", "rs157580", "rs11591147", "rs12740374", "rs693", "rs261290", "rs112875651", "rs58542926", "rs247616", "rs12916", "rs633695", "rs964184", "rs115478735", "rs77542162", "rs4299376", "rs102275", "rs4704834", "rs72848251", "rs1883711", "rs4665972", "rs3777411", "rs11206517", "rs2738447", "rs34042070", "rs926053", "rs77960347", "rs56325564", "rs472495", "rs12804228", "rs72631343", "rs553427", "rs11854242", "rs1169292", "rs838146", "rs6601299", "rs10260606", "rs2326077", "rs2737245", "rs653178", "rs72843670", "rs9616847", "rs4939883", "rs1500188", "rs35633876", "rs6602911", "rs72694391", "rs12620918"),
  ea  = c("C", "G", "G", "G", "G", "G", "T", "G", "C", "C", "T", "A", "G", "A", "A", "G", "T", "A", "G", "G", "T", "C", "T", "A", "C", "G", "A", "G", "G", "A", "C", "C", "C", "C", "C", "T", "G", "C", "G", "C", "T", "A", "T", "A", "G", "C", "T", "T"),
  oa  = c("A", "A", "A", "T", "T", "A", "C", "A", "T", "T", "C", "G", "C", "T", "G", "T", "C", "G", "A", "C", "C", "T", "G", "C", "G", "C", "G", "A", "T", "G", "G", "T", "T", "T", "T", "C", "C", "T", "T", "T", "C", "T", "C", "G", "T", "T", "C", "C"),
  bx  = c(0.26377, 0.152836, -0.0836374, 0.308882, 0.0961513, -0.0769739, 0.0759505, 0.068676, 0.118025, 0.0608862, -0.0580953, -0.0546603, 0.0710628, -0.0594489, -0.149245, 0.047141, 0.0453463, -0.0412619, -0.0484938, -0.108501, 0.0359412, 0.0483808, -0.0908334, -0.0328659, -0.0411867, 0.0501326, -0.134647, -0.0303526, -0.0314252, -0.0522524, 0.0429508, -0.0286973, 0.0317694, -0.030447, 0.0278811, -0.044935, -0.0337973, 0.027081, 0.0279349, -0.0247731, -0.0352589, -0.0250082, -0.0309854, 0.0240193, 0.0238608, -0.0245499, -0.0226707, -0.0295835),
  sex = c(0.00659909, 0.00630319, 0.00418106, 0.0156254, 0.00493378, 0.00410273, 0.0043224, 0.0042589, 0.00780678, 0.00438395, 0.00418643, 0.00452543, 0.00604505, 0.00530499, 0.0137929, 0.00439223, 0.00428851, 0.00427763, 0.00531677, 0.0120475, 0.00421805, 0.00573293, 0.011479, 0.00417665, 0.00528041, 0.00651738, 0.0179189, 0.00413217, 0.0042923, 0.00742254, 0.00611975, 0.00410145, 0.00462128, 0.00444764, 0.00411717, 0.0068214, 0.0052964, 0.00434692, 0.00458436, 0.00410821, 0.00592019, 0.00422751, 0.00532421, 0.00412981, 0.00411923, 0.00427046, 0.00411579, 0.00542051),
  stringsAsFactors = FALSE)
glgc_z <- function(nm) {
  d <- fread(file.path("gwas", gwas_dl[[nm]]$file), select = c("rsid", "A1", "A2", "beta", "se"), showProgress = FALSE,
             colClasses = c(rsid = "character", A1 = "character", A2 = "character"))
  setnames(d, c("SNP", "A1", "A2", "b", "se")); d <- d[SNP %chin% INSTR$SNP]; d[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]
  unique(d, by = "SNP")
}
inst <- copy(as.data.table(INSTR))
for (nm in c("GLGC_LDL", "GLGC_TG")) {
  g <- glgc_z(nm); i <- match(inst$SNP, g$SNP)
  a1 <- g$A1[i]; a2 <- g$A2[i]
  same <- !is.na(a1) & a1 == inst$ea & a2 == inst$oa; swap <- !is.na(a1) & a1 == inst$oa & a2 == inst$ea
  z <- g$b[i] / g$se[i]; bb <- g$b[i]
  z[!(same | swap)] <- NA; bb[!(same | swap)] <- NA                          # alleles that do not match (palindromic / strand) are left out
  z[swap & !is.na(z)] <- -z[swap & !is.na(z)]; bb[swap & !is.na(bb)] <- -bb[swap & !is.na(bb)]    # align to the RC-raising allele
  inst[[paste0(nm, "_z")]] <- z; inst[[paste0(nm, "_b")]] <- bb
}
ok <- complete.cases(inst[, .(GLGC_LDL_z, GLGC_TG_z)])
cat("Instruments with both GLGC LDL and TG available:", sum(ok), "of", nrow(inst), "\n")
cat("Genome-wide significant (|z| > 5.45) for LDL:", sum(abs(inst$GLGC_LDL_z) > 5.45, na.rm = TRUE),
    "| for TG:", sum(abs(inst$GLGC_TG_z) > 5.45, na.rm = TRUE), "\n")
cat("Instruments with |z_LDL| > |z_TG|:", sum(abs(inst$GLGC_LDL_z) > abs(inst$GLGC_TG_z), na.rm = TRUE), "of", sum(ok), "\n")
cat("Mean z for the RC-raising allele: LDL =", round(mean(inst$GLGC_LDL_z, na.rm = TRUE), 2), "| TG =", round(mean(inst$GLGC_TG_z, na.rm = TRUE), 2), "\n")
cat("Correlation of instrument effect on RC (bx) with effect on LDL:", round(cor(inst$bx[ok], inst$GLGC_LDL_b[ok]), 3),
    "| on TG:", round(cor(inst$bx[ok], inst$GLGC_TG_b[ok]), 3), "\n")
fwrite(inst, "results/instrument_composition_R.csv")

## ---------------------------------------------------------------------------
## 13. COMPARISON WITH THE OFFICIAL-SOFTWARE RESULTS
##     OFFICIAL_* below are the numbers produced by the official LDSC (v2.0.1, Python) for this study.
## ---------------------------------------------------------------------------
cat("\n==== 13. R vs OFFICIAL LDSC ====\n")
OFFICIAL_RG <- data.frame(
  run = c("rg_GLGC_LDL_vs_all", "rg_GLGC_LDL_vs_all", "rg_GLGC_LDL_vs_all", "rg_GLGC_LDL_vs_all", "rg_GLGC_TG_vs_all", "rg_GLGC_TG_vs_all", "rg_GLGC_TG_vs_all", "rg_GLGC_TG_vs_all", "rg_GLGC_TG_vs_all", "rg_Ghodsian_vs_Fairfield", "rg_LDL_vs_CAD_TG_bulk", "rg_LDL_vs_CAD_TG_bulk", "rg_LDL_vs_CAD_TG_loci", "rg_LDL_vs_CAD_TG_loci", "rg_RC_vs_all", "rg_RC_vs_all", "rg_RC_vs_all", "rg_RC_vs_all_bulk", "rg_RC_vs_all_bulk", "rg_RC_vs_all_bulk", "rg_RC_vs_all_bulk", "rg_RC_vs_all_bulk", "rg_RC_vs_all_loci", "rg_RC_vs_all_loci", "rg_RC_vs_all_loci", "rg_RC_vs_all_loci", "rg_RC_vs_all_loci", "sens_A_chisqmax_off", "sens_A_chisqmax_off", "sens_A_chisqmax_off", "sens_B_full1000G_ref", "sens_B_full1000G_ref", "sens_B_full1000G_ref", "sens_C_twostep", "sens_C_twostep", "sens_C_twostep"),
  trait1 = c("GLGC_LDL", "GLGC_LDL", "GLGC_LDL", "GLGC_LDL", "GLGC_TG", "GLGC_TG", "GLGC_TG", "GLGC_TG", "GLGC_TG", "Ghodsian", "GLGC_LDL_bulk", "GLGC_LDL_bulk", "GLGC_LDL_loci", "GLGC_LDL_loci", "RC", "RC", "RC", "RC_bulk", "RC_bulk", "RC_bulk", "RC_bulk", "RC_bulk", "RC_loci", "RC_loci", "RC_loci", "RC_loci", "RC_loci", "RC", "RC", "RC", "RC", "RC", "RC", "RC", "RC", "RC"),
  trait2 = c("CAD", "RC", "Ghodsian", "Fairfield", "CAD", "GLGC_LDL", "RC", "Ghodsian", "Fairfield", "Fairfield", "CAD_bulk", "GLGC_TG_bulk", "CAD_loci", "GLGC_TG_loci", "Ghodsian", "Fairfield", "CAD", "CAD_bulk", "GLGC_LDL_bulk", "GLGC_TG_bulk", "Ghodsian_bulk", "Fairfield_bulk", "CAD_loci", "GLGC_LDL_loci", "GLGC_TG_loci", "Ghodsian_loci", "Fairfield_loci", "CAD", "Ghodsian", "Fairfield", "CAD", "Ghodsian", "Fairfield", "CAD", "Ghodsian", "Fairfield"),
  rg = c(0.2239, 0.8442, 0.2425, 0.1584, 0.27, 0.3624, 0.4239, 0.6358, 0.3784, 1.0359, 0.1908, 0.3061, 0.5129, 0.4414, -0.1891, 0.0656, -0.0016, -0.0357, 0.8178, 0.2211, -0.2993, 0.0661, 0.2952, 0.9814, 1.1429, 0.1302, 0.0766, -0.0016, -0.1891, 0.0656, -0.0151, -0.5663, 0.0292, 0.0446, -0.1722, 0.0486),
  se = c(0.041, 0.0537, 0.1369, 0.0758, 0.0337, 0.0546, 0.1033, 0.1798, 0.0802, 0.1617, 0.0422, 0.0396, 0.1314, 0.1048, 0.1421, 0.0916, 0.0594, 0.0493, 0.0591, 0.0521, 0.1622, 0.0828, 0.1756, 0.0987, 0.8236, 0.3753, 0.3924, 0.0594, 0.1421, 0.0916, 0.0758, 0.5939, 0.1129, 0.0425, 0.1085, 0.0667),
  z = c(5.4541, 15.7191, 1.7708, 2.0906, 8.0155, 6.6427, 4.1027, 3.5372, 4.7196, 6.4063, 4.5228, 7.7274, 3.9044, 4.2134, -1.3302, 0.7157, -0.027, -0.7244, 13.8272, 4.246, -1.8451, 0.7979, 1.6811, 9.9391, 1.3877, 0.3469, 0.1953, -0.027, -1.3302, 0.7157, -0.1985, -0.9534, 0.2582, 1.0484, -1.5882, 0.729),
  p = c(4.9209e-08, 1.1196e-55, 0.0766, 0.0366, 1.0973e-15, 3.0807e-11, 4.0843e-05, 0.0004, 2.3626e-06, 1.4907e-10, 6.1017e-06, 1.0976e-14, 9.4459e-05, 2.5154e-05, 0.1835, 0.4742, 0.9785, 0.4688, 1.7471e-43, 2.1762e-05, 0.065, 0.4249, 0.0927, 2.8143e-23, 0.1652, 0.7287, 0.8451, 0.9785, 0.1835, 0.4742, 0.8427, 0.3404, 0.7962, 0.2945, 0.1122, 0.466),
  gcov_int = c(0.0186, 0.0325, -0.0061, -0.004, 0.0062, 0.1715, -0.0089, -0.008, -0.0039, 0.3245, 0.0113, 0.1934, 0.3074, -0.1866, 0.0035, -0.0078, 0.0144, 0.0005, -0.0047, 0.0077, 0.0063, -0.0054, 0.4171, 1.5903, -0.4817, -0.0786, -0.0854, 0.0144, 0.0035, -0.0078, 0.0154, 0.0038, -0.0046, 0.0055, 0.0044, -0.0069),
  gcov_int_se = c(0.0106, 0.0504, 0.0063, 0.006, 0.006, 0.0128, 0.0145, 0.0062, 0.007, 0.0061, 0.0049, 0.0058, 0.2949, 0.1665, 0.0051, 0.0057, 0.0084, 0.0047, 0.0067, 0.0067, 0.0047, 0.005, 0.2182, 1.6184, 0.3167, 0.081, 0.0835, 0.0084, 0.0051, 0.0057, 0.0078, 0.0047, 0.0054, 0.0047, 0.0045, 0.005),
  h2_int_trait2 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  stringsAsFactors = FALSE)
OFFICIAL_H2 <- data.frame(
  trait = c("CAD", "Fairfield", "Ghodsian", "RC"),
  h2 = c(0.0725, 0.1453, 0.0233, 0.1095),
  h2_se = c(0.0049, 0.0185, 0.0086, 0.0152),
  mean_chi2 = c(1.1243, 1.0788, 1.0309, 1.2448),
  lambda_gc = c(1.0466, 1.0649, 1.0285, 1.127),
  intercept = c(0.8811, 0.9964, 1.0063, 1.0251),
  intercept_se = c(0.0074, 0.0067, 0.007, 0.0095),
  stringsAsFactors = FALSE)
RG_R <- rbindlist(RES[c("primary", "nafld_cc", "sensA", "sensB", "sensC", "glgc_tg", "glgc_ldl", "bulkloci")], fill = TRUE)
CMP <- merge(as.data.table(OFFICIAL_RG), RG_R, by = c("run", "trait1", "trait2"), suffixes = c("_official", "_R"))
CMP[, `:=`(d_rg = rg_R - rg_official, d_se = se_R - se_official)]
print(CMP[, .(run, trait1, trait2, rg_official, rg_R = round(rg_R, 4), se_official, se_R = round(se_R, 4), d_rg = round(d_rg, 5), d_se = round(d_se, 5))])
cat("\nLargest absolute difference in rg:", signif(max(abs(CMP$d_rg)), 3), "| in SE:", signif(max(abs(CMP$d_se)), 3),
    "| comparisons:", nrow(CMP), "(official values are rounded to 4 decimals)\n")
H2C <- merge(as.data.table(OFFICIAL_H2), H2, by = "trait", suffixes = c("_official", "_R"))
print(H2C[, .(trait, h2_official, h2_R = round(h2_R, 4), se_official = h2_se_official, se_R = round(h2_se_R, 4),
              int_official = intercept_official, int_R = round(intercept_R, 4))])
fwrite(RG_R, "results/ldsc_rg_all_runs_R.csv")

## ---------------------------------------------------------------------------
## 14. THE NUMBERS FOR THE MANUSCRIPT
## ---------------------------------------------------------------------------
cat("\n==== 14. SUMMARY ====\n")
pr <- RES$primary
for (i in seq_len(nrow(pr)))
  cat(sprintf("rg(RC, %s) = %.3f (SE %.3f), z = %.2f, P = %.3f; h2 of %s on the liability scale = %.3f (SE %.3f)\n",
              pr$trait2[i], pr$rg[i], pr$se[i], pr$z[i], pr$p[i], pr$trait2[i], pr$h2_trait2[i], pr$h2_trait2_se[i]))
cat("Finished:", format(Sys.time()), "\n")
sink()
