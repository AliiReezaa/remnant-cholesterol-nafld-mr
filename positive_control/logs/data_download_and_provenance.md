# Data provenance: CAD outcome for the positive control

Study: Nikpay M, Goel A, Won HH, et al. A comprehensive 1,000 Genomes-based genome-wide association meta-analysis of coronary artery disease. Nat Genet 2015;47(10):1121-1130. doi:10.1038/ng.3396. PMID 26343387.
Resource: GWAS Catalog accession GCST003116 (CARDIoGRAMplusC4D 1000 Genomes-based GWAS), EFO_0000378 (coronary artery disease). 60,801 cases and 123,504 controls (counts from the study README at https://grasp.nhlbi.nih.gov/downloads/FullResults/2015/2015_Nikpey_CardioGramC4D/NikpeyREADME.txt), mixed ancestry (mainly European, South Asian and East Asian).
File downloaded (harmonised, genome build 37), 15 October 2026:
  https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003116/harmonised/26343387-GCST003116-EFO_0000378-build37.f.tsv.gz
Method: `curl -s -C - -o cad_nikpay2015_GCST003116_build37.f.tsv.gz <URL>`, repeated until complete (the transfer was interrupted several times and resumed).
Size on disk: 239,499,711 bytes. Rows: 9,455,778 variants (9,455,779 lines including header).
MD5 of the downloaded file: 0d296695219d8f4867bee625ac031ff0
MD5 published by EMBL-EBI in the same directory (md5sum.txt): 0d296695219d8f4867bee625ac031ff0  (match)
`gzip -t` passed.
Columns: variant_id, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, nstudy, odds_ratio, ci_lower, ci_upper. beta is the log odds ratio of CAD per copy of the effect allele.

The full file (240 MB) is not stored in the GitHub repository (file-size limit). The rows used are stored in inputs/:
  cad_source_rows_direct_hits.tsv            the 46 instrument rows found by rsID
  cad_source_rows_substitute_candidates.tsv  rows for same-position and proxy candidates of the 6 instruments missing by rsID
  proxy_candidates_queried.txt, same_position_queries.txt
  ld_proxy_table_EUR_r2_ge_0.8.tsv           proxy candidates (r2 >= 0.8, 1000 Genomes EUR reference panel, 503 individuals)
  cad_outcome_for_52_instruments.tsv         final outcome table used by both the Python and R analyses (50 rows)
  rc_52_instruments_exposure.tsv             the 52 exposure instruments (Richardson 2022, GCST90092943)
Anyone can re-create inputs/ from the public file with scripts s11_extract_cad.py (Python).

Exposure: Richardson TG, et al. Characterising metabolomic signatures of lipid-modifying therapies through drug target Mendelian randomisation. PLoS Biol 2022;20(2):e3001547. doi:10.1371/journal.pbio.3001547 (GWAS Catalog GCST90092943; remnant cholesterol, UK Biobank NMR, n = 115,082, European ancestry).
Sample overlap: the 2015 CAD meta-analysis predates UK Biobank genotype data release, so it should not contain UK Biobank participants.

Substitutions for instruments absent from the CAD file by rsID (see cad_substitutes_used.tsv):
  rs115478735  proxy rs507666 (r2 0.988)
  rs838146     proxy rs838147 (r2 1.0)
  rs72848251   same chromosome:position and same allele pair, listed as rs115129597 in the CAD file
  rs72843670   same chromosome:position and same allele pair, listed as rs115366102 in the CAD file
  rs574454081, rs532754706: not found, no usable proxy -> dropped
  rs1293261: palindromic with intermediate allele frequency -> dropped at harmonisation (both Python and R)
Result: 49 of 52 instruments used.
