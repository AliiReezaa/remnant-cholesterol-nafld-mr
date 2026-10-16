# Remnant cholesterol and NAFLD: two-sample Mendelian randomization

Code, harmonised instrument data and results for the study of genetically predicted remnant cholesterol (RC) and non-alcoholic fatty liver disease (NAFLD).

Authors: Alireza Soleymani T., Mahdi Akbarzadeh, Danial Habibi (corresponding).

## Summary of result
Across three outcome GWAS (Ghodsian GCST90091033; Fairfield GCST90054782; FinnGen R5 finn-b-NAFLD), genetically predicted RC was not associated with NAFLD (IVW OR 1.05, 0.86 to 1.27; 1.00, 0.78 to 1.29; 0.85, 0.56 to 1.28). The study is reported as a null result.

## Positive control (added after the main analysis)
To check instrument validity we tested the same 52 instruments against coronary artery disease (CARDIoGRAMplusC4D, Nikpay et al. 2015, GCST003116): IVW OR 1.58 (1.38 to 1.81) per 1-SD higher genetically predicted RC, P = 5e-11, consistent across Egger, weighted median, mode and outlier-corrected analyses. Everything is in `positive_control/` (scripts, Python and R console logs, three R runs, inputs, outputs, provenance with MD5, side-by-side Python vs R comparison). See `positive_control/README.md`. The 240 MB CAD file is not stored here (download link and MD5 in `positive_control/logs/data_download_and_provenance.md`).

## Data sources (not redistributed here)
- Exposure: Richardson et al. 2022 remnant cholesterol GWAS (Nightingale NMR, UK Biobank, European).
- Outcomes: GWAS Catalog GCST90091033 and GCST90054782; FinnGen R5 (finn-b-NAFLD, OpenGWAS).
- Positive control outcome: GWAS Catalog GCST003116 (Nikpay et al. 2015, coronary artery disease).
- LD reference: 1000 Genomes European panel (see `code/dl_ld.sh`).

## Manuscript
The manuscript and supplement (results/manuscript) are generated from the result files by `code/manuscript_build` (build_main.py, build_supp.py); figures by `code/s13_manuscript_figures.py` (results/figures_v2).

## Layout
- `code/` analysis scripts (Python) and `run_mrpresso.R` (R MR-PRESSO confirmation).
- `analysis/` harmonised data, MR results and leave-one-out output for the primary 52-variant set, the no-APOE sensitivity set and the original 48-variant set.
- `work/` clumped instruments and proxies used.
- `positive_control/` CAD positive control: scripts, logs (Python and R console), data provenance, inputs and outputs of all runs. `code/s11_extract_cad.py` and `code/s12_positive_control.py` are copies of its Python scripts.
- `results/` summary table, Supplementary_Tables.xlsx, MR-PRESSO output, figures.

## Reproducing
Scripts assume `~/rc_nafld_mr/` with `inputs/` (raw GWAS files, not included), `work/` and `analysis/`. Run in order: s1_instruments, s1b_dedupe_pos, s2_compare, s3_extract, s3b_extract_proxy, s5_clump, s6_proxies, s6b_apply_proxies, s4_analysis (main MR), s7_reproduce, s8_figures, s9_pubfigs, s10_supp; validate_against_R.py and run_mrpresso.R are cross-checks.
`pip install -r requirements.txt`

## Licence
MIT (code). Underlying GWAS summary statistics remain under their original terms.
