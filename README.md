# Remnant cholesterol and NAFLD: two-sample Mendelian randomization

Code, harmonised instrument data and results for the study of genetically predicted remnant cholesterol (RC) and non-alcoholic fatty liver disease (NAFLD).

Authors: Alireza Soleymani T., Mahdi Akbarzadeh, Danial Habibi (corresponding).

## Summary of result
Across three outcome GWAS (Ghodsian GCST90091033; Fairfield GCST90054782; FinnGen R5 finn-b-NAFLD), genetically predicted RC was not associated with NAFLD (IVW OR 1.05, 0.86 to 1.27; 1.00, 0.78 to 1.29; 0.85, 0.56 to 1.28). The study is reported as a null result.

## Data sources (not redistributed here)
- Exposure: Richardson et al. 2022 remnant cholesterol GWAS (Nightingale NMR, UK Biobank, European).
- Outcomes: GWAS Catalog GCST90091033 and GCST90054782; FinnGen R5 (finn-b-NAFLD, OpenGWAS).
- LD reference: 1000 Genomes European panel (see `code/dl_ld.sh`).

## Layout
- `code/` analysis scripts (Python) and `run_mrpresso.R` (R MR-PRESSO confirmation).
- `analysis/` harmonised data, MR results and leave-one-out output for the primary 52-variant set, the no-APOE sensitivity set and the original 48-variant set.
- `work/` clumped instruments and proxies used.
- `results/` summary table, Supplementary_Tables.xlsx, MR-PRESSO output, figures.

## Reproducing
Scripts assume `~/rc_nafld_mr/` with `inputs/` (raw GWAS files, not included), `work/` and `analysis/`. Run in order: s1_instruments, s1b_dedupe_pos, s2_compare, s3_extract, s3b_extract_proxy, s5_clump, s6_proxies, s6b_apply_proxies, s4_analysis (main MR), s7_reproduce, s8_figures, s9_pubfigs, s10_supp; validate_against_R.py and run_mrpresso.R are cross-checks.
`pip install -r requirements.txt`

## Licence
MIT (code). Underlying GWAS summary statistics remain under their original terms.
