# Genetic correlation (LD score regression): remnant cholesterol vs NAFLD

Exploratory analysis requested after the MR. Nothing here has been added to the manuscript.

Pipeline (run in this order, all logs in `logs/`):
- `00_make_hm3_alleles.py`  HapMap3 allele list (1,187,349 SNPs, MHC excluded) from the S-LDSC reference files (Zenodo 10.5281/zenodo.7768714, 1000G Phase 3 EUR).
- `01_prep_sumstats.py`  full GWAS Catalog summary statistics for RC (GCST90092943), Ghodsian (GCST90091033), Fairfield (GCST90054782, original author file for the INFO column) and CAD (GCST003116). MD5 of every download in `provenance/`.
- `02_munge.sh`  LDSC `munge_sumstats.py` (INFO > 0.9 where available, MAF > 0.01, HM3 alleles, strand-ambiguous SNPs removed).
- `03_run_ldsc.sh`  SNP heritability and rg (free intercepts, block jackknife, HM3 LD scores as reference and weights, M = 1,187,349).
- `04_independent_check.py`  independent numpy re-implementation of the regression (agrees with ldsc.py).
- `05_benchmark_glgc.sh`  external benchmark with GLGC (Willer 2013) TG and LDL against CAD.
- `06_sensitivity.sh`, `08_bulk_vs_loci.sh`  alternative LDSC settings and a bulk-genome analysis.
- `07_collect.py`  parses all LDSC logs into `results/*.csv`.
- `R/Genetic_correlation_RC_NAFLD_full_pipeline.R`  the whole analysis in one R script for RStudio (downloads with MD5 checks, munging, h2, rg, sensitivity, GLGC benchmark, bulk vs loci, self-tests, independent cross-check, instrument composition). It is an R port of the official LDSC code and reproduces all 36 rg values and 4 h2 values of the Python run (largest difference 0.00005, official values are rounded to 4 decimals). A console log and CSV outputs of one complete run are in `R/example_output/`.

Software: LDSC 2.0.1 (PyPI port of bulik/ldsc, Python 3.11, numpy 1.26.4, pandas 1.5.3, scipy 1.11.4). The port failed to read gzipped headers in `munge_sumstats.py`; the one-line fix is in `provenance/munge_sumstats_py3_gzip_patch.diff`.
Checks: self-comparison tests (same data with swapped alleles, other strand, negated) return rg = +1, +1, +1, -1; the independent implementation reproduces ldsc.py; GLGC TG and LDL give the expected positive rg with CAD.
Not run: FinnGen (894 cases in R5, too few for LDSC, and the file is not downloadable from this workspace).
