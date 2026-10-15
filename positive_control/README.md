# Positive control: genetically predicted remnant cholesterol and coronary artery disease (CAD)

Purpose: check that the 52 remnant cholesterol instruments used for the NAFLD analysis reproduce a known causal effect (remnant cholesterol raises coronary disease risk; Varbo 2013, Nordestgaard 2016 in the manuscript reference list). Done after the NAFLD analysis, not pre-registered. Reported in full whatever it showed.

## Result (final, R run 3, TwoSampleMR 0.6.22; Python agrees)
Odds ratio of CAD per 1-SD higher genetically predicted remnant cholesterol, 49 variants, Nikpay 2015 (60,801 cases, 123,504 controls)

| Method | OR (95% CI) | P |
|---|---|---|
| IVW (multiplicative random effects) | 1.58 (1.38 to 1.81) | 5.2e-11 |
| MR-Egger | 1.87 (1.46 to 2.40) | 9.4e-06 |
| Weighted median | 1.48 (1.32 to 1.65) | 5.6e-12 |
| Simple mode | 1.41 (1.09 to 1.82) | 0.012 |
| Weighted mode | 1.39 (1.22 to 1.60) | 1.7e-05 |

Cochran's Q 276.9 (48 df), P = 6.2e-34. Egger intercept -0.0115 (SE 0.0072), P = 0.12. MR-PRESSO (R): global test RSSobs 290.7, P < 1e-4; 6 outliers (rs115478735, rs12151108, rs12740374, rs1500188, rs3777411, rs653178); outlier-corrected OR 1.49 (1.35 to 1.65), P = 1.4e-9; distortion test P = 0.175. MR-PRESSO (Python): 5 outliers, outlier-corrected OR 1.54 (1.39 to 1.71). Steiger filter (Python) removed rs653178: IVW OR 1.60 (1.42 to 1.81). Leave-one-out IVW OR 1.51 to 1.62, all P < 1e-8. Mean F 150, 49-variant R2 6.4%.

## Folder guide
- `logs/data_download_and_provenance.md` data source, file size, MD5 check, substitutions, citations
- `logs/python_s11_extract_cad.log`, `logs/python_s12_positive_control.log` console logs of the two Python scripts
- `s11_extract_cad.py` extracts the instruments from the CAD summary statistics
- `s12_positive_control.py` Python MR (same pipeline and library `mr_lib.py` as the NAFLD analysis, in `../code/`)
- `run_positive_control_TwoSampleMR.R` final R script (run 3)
- `R_out/` outputs of run 3, including `R_console_log.txt` (everything printed to the console, with session info)
- `R_out_run1_first_run/` + `run_positive_control_TwoSampleMR_run1.R`: first R run (partial logging)
- `R_out_run2_logged_unseeded/` + `run_positive_control_TwoSampleMR_run2_unseeded.R`: second run with full logging, before the random seed was set for the bootstrap
- `python_out/` Python results (mr_results.tsv, harmonised data, leave-one-out, summary.json)
- `positive_control_python_vs_R.csv` side-by-side comparison
- `inputs/` all tables fed to both analyses and the rows pulled from the public CAD file
- `cad_substitutes_used.tsv` instruments that needed a substitute

## Run history (nothing hidden)
1. Run 1 (R): results as above, console output sunk to `R_log.txt` only (no warnings/session info).
2. Run 2 (R): logging expanded. Point estimates and MR-PRESSO identical to run 1. The bootstrap standard errors of the weighted median and mode estimators changed slightly because `mr()` was called without a random seed (e.g. simple mode P 0.015 in run 1 vs 0.012 in run 2).
3. Run 3 (R, final): `set.seed(2026)` added before `mr()`. IVW, Egger, heterogeneity, intercept and MR-PRESSO identical to runs 1 and 2. The values in this README come from run 3.
4. Python: the pipeline was run twice (original and a logged re-run); outputs are identical.
5. Differences between Python and R: weighted mode (Python 1.41 vs R 1.39; Python's mode is evaluated on a 512-point grid and its simple and weighted modes fell on the same grid point), bootstrap SEs for the median and modes, and MR-PRESSO (Python found 5 outliers, R 6, as in the NAFLD analysis the R package is taken as primary).

## Caveats
- The CAD GWAS is of mixed ancestry (mainly European, with South and East Asian participants); the exposure GWAS is European.
- 2 instruments absent by rsID were matched by identical chromosome:position and alleles (merged rsIDs), 2 by LD proxies, 2 were dropped, 1 dropped as an ambiguous palindrome (49 of 52 used).
- Heterogeneity is large (APOE, LPA-region and other lipid loci); the effect persists after outlier removal and with Egger and Steiger filtering.
- This is a validity check for the instruments, not a new causal claim.
- The 240 MB CAD file is not stored in the repo. Download link and MD5 are in `logs/data_download_and_provenance.md`.

## Citations
- Nikpay M, et al. A comprehensive 1,000 Genomes-based genome-wide association meta-analysis of coronary artery disease. Nat Genet 2015;47:1121-1130. doi:10.1038/ng.3396
- Richardson TG, et al. PLoS Biol 2022;20:e3001547. doi:10.1371/journal.pbio.3001547
- Sollis E, et al. The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource. Nucleic Acids Res 2023;51:D977-D985. doi:10.1093/nar/gkac1010
- Hemani G, et al. The MR-Base platform supports systematic causal inference across the human phenome. eLife 2018;7:e34408 (TwoSampleMR). doi:10.7554/eLife.34408
- Verbanck M, et al. Nat Genet 2018;50:693-698 (MR-PRESSO). doi:10.1038/s41588-018-0099-7
- Hemani G, Tilling K, Davey Smith G. PLoS Genet 2017;13:e1007081 (Steiger). doi:10.1371/journal.pgen.1007081
- Varbo A, et al. J Am Coll Cardiol 2013;61:427-436. doi:10.1016/j.jacc.2012.08.1026
- Auton A, et al. A global reference for human genetic variation. Nature 2015;526:68-74 (LD reference). doi:10.1038/nature15393
- R Core Team. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna; 2025 (R 4.4.3). https://www.R-project.org/
