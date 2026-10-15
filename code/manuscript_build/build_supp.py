import sys; sys.path.insert(0,'.')
from lib import *; from data import *
import pandas as pd
FIG='/home/claude/remnant-cholesterol-nafld-mr/results/figures_v2/'; OUT='/home/claude/remnant-cholesterol-nafld-mr/results/manuscript_v2/'
D=Doc(10)
D.para('Supplemental Digital Content',align='left',size=10,italic=True,after=2); D.rule()
D.para('**Mendelian randomization study of remnant cholesterol in relation to non-alcoholic fatty liver disease**',align='left',size=17,after=6)
D.para('Alireza Soleymani T, Mahdi Akbarzadeh, Danial Habibi',align='left',size=11,after=10)
D.para('This document accompanies the main manuscript. It contains Supplemental Digital Content 1 (supplementary methods, run history and the STROBE-MR checklist) and Supplemental Digital Content 2 (supplementary tables S1–S12 and supplementary figures S1–S6). The complete tables, including all variant-level data, are in the Excel workbook Supplementary_Tables.xlsx (Supplemental Digital Content 3), and all code, console logs and intermediate files are in the public repository https://github.com/AliiReezaa/remnant-cholesterol-nafld-mr.',size=10)
D.table(('Table of contents','Overview.'),['Item','Content'],[
 ['SDC 1, Part A','Supplementary methods 1–8: data provenance, instrument selection flow, formulas, harmonisation, MR-PRESSO and extended estimator settings, validation, run history, software'],
 ['SDC 1, Part B','STROBE-MR checklist with location in the manuscript'],
 ['SDC 2, Tables','S1 instruments; S2 harmonised variant data; S3 all MR estimates; S4 heterogeneity and pleiotropy; S5 leave-one-out; S6 proxies; S7 positive control; S8 extended estimators; S9 heterogeneity (I², Egger, NOME); S10 influence diagnostics; S11 radial outliers; S12 cross-dataset concordance'],
 ['SDC 2, Figures','S1 graphical abstract; S2–S3 diagnostics for Fairfield and FinnGen; S4 instrument strength, heterogeneity and influence; S5 radial MR; S6 positive-control diagnostics'],
 ['SDC 3','Supplementary_Tables.xlsx']],widths=[1.3,5.4],size=9)
D.h1('Supplemental Digital Content 1')
D.h2('Part A. Supplementary methods')
D.para('**1. Data sources and provenance.** All data are public summary statistics. Table SM1 lists each source, accession and the details recorded when the data were extracted. The CAD file was downloaded on 15 October 2026 from the GWAS Catalog FTP in its genome build 37 harmonised form (239,499,711 bytes, 9,455,778 variants) and its MD5 checksum (0d296695219d8f4867bee625ac031ff0) matched the checksum published by EMBL-EBI. Case and control counts for the CAD meta-analysis (60,801 and 123,504) were taken from the study README on the GRASP site and not from the paper text, which could not be retrieved, and should be checked against the paper before submission.')
D.table(('Table SM1','Data sources used in the study.'),['Role','Source','Accession / identifier','Sample','Notes'],[
 ['Exposure: remnant cholesterol','Richardson et al, PLoS Biol 2022 (UK Biobank NMR, Nightingale)','GCST90092943 (OpenGWAS)','115,082, European ancestry','Inverse rank-normalised, SD units; BOLT-LMM'],
 ['Outcome 1: NAFLD','Ghodsian et al, Cell Rep Med 2021','GCST90091033','8,434 cases / 770,180 controls','eMERGE, FinnGen freeze 4, UK Biobank, Estonian Biobank; EHR-defined'],
 ['Outcome 2: NAFLD','Fairfield et al, Hepatol Commun 2022','GCST90054782','4,761 cases / 373,227 controls','UK Biobank; consensus diagnostic codes'],
 ['Outcome 3: NAFLD','FinnGen release 5','finn-b-NAFLD (OpenGWAS)','894 cases / 217,898 controls','FinnGen'],
 ['Positive control: CAD','Nikpay et al, Nat Genet 2015 (CARDIoGRAMplusC4D 1000G)','GCST003116','60,801 cases / 123,504 controls','Mixed ancestry, mainly European; build 37 harmonised file'],
 ['LD reference','1000 Genomes Phase 3 European panel','–','503 individuals','Clumping (r² < 0.001, 10,000 kb) and proxies (r² ≥ 0.8, 250 kb)']],widths=[1.2,1.8,1.2,1.2,1.3],size=8)
D.para('**2. Instrument selection flow.** 8,873 variants with P < 5 × 10^^−8^^ for remnant cholesterol (after removal of duplicated identifiers) → 7,593 present in the 1000 Genomes European panel (1,280 absent and not clumpable) → 52 independent variants after clumping (r² < 0.001, 10,000 kb, ranked by P). In the outcome datasets: Ghodsian, 3 variants missing without proxy and 6 variants replaced by proxies (Table S6), 1 palindromic variant removed, 48 used; Fairfield, 0 missing, 1 palindromic variant removed, 51 used; FinnGen, 3 missing without proxy, 3 proxies, 1 palindromic variant removed, 48 used. Positive control: 46 found by rsID, 2 by position, 2 by proxy, 2 absent, 1 palindromic removed, 49 used.')
D.para('**3. Formulas.** Per-variant F statistic: F = (β/SE)^^2^^. Variance explained by a variant, R^^2^^ = 2p(1−p)β^^2^^ with β in SD units; the overall F statistic was F = R^^2^^(N − k − 1) / ((1 − R^^2^^)k), with R^^2^^ the summed variance explained, k the number of variants and N the exposure sample size (159.2 for 48 variants, R^^2^^ = 6.23%, N = 115,082). Wald ratio: β~~Y~~/β~~X~~ with standard error SE~~Y~~/|β~~X~~|. IVW estimate: a regression of β~~Y~~ on β~~X~~ weighted by 1/SE~~Y~~^^2^^ with the intercept fixed at zero and multiplicative random-effects standard errors (residual standard error floored at one). Steiger test: for each variant the variance explained in the exposure was compared with that in the outcome; the latter was approximated on the latent logistic scale as 2p(1−p)β^^2^^ / (2p(1−p)β^^2^^ + π^^2^^/3) without a prevalence assumption. Minimum detectable effect with 80% power: |log OR| = (z~~0.975~~ + z~~0.80~~) × SE~~IVW~~, i.e. 2.80 × SE~~IVW~~. Bias from complete sample overlap relative to the confounded association: approximately 1/F~~overall~~.')
D.para('**4. Harmonisation.** Effects were aligned to the RC-increasing allele. Proxies were selected as the variant with the highest r^^2^^ ≥ 0.8 within 250 kb in the 1000 Genomes European panel and alleles were aligned by LD phase. Palindromic variants with effect allele frequency between 0.42 and 0.58 were excluded; other palindromic variants were aligned by allele frequency. No allele mismatches were found.')
D.para('**5. MR-PRESSO and extended estimator settings.** MR-PRESSO was run with the MRPRESSO R package with 10,000 simulations, outlier and distortion tests, a significance threshold of 0.05 for outlier detection and seed 2026; the raw and outlier-corrected estimates, the global test and the outlier list are in the repository (results/mrpresso_results.csv). Extended estimators were run with the MendelianRandomization package (maximum likelihood, MR-Lasso, contamination mixture, debiased IVW, IVW fixed and random effects) with default settings, with mr.raps (Tukey and Huber losses with overdispersion) and with RadialMR (IVW radial, modified second-order weights; outliers flagged at the Bonferroni-adjusted 5% level of the individual Q statistic). The random seed was 2026. Influence diagnostics used Cook’s distance (threshold 4/n) and studentized residuals from the weighted regression of β~~Y~~ on β~~X~~ through the origin.')
D.para('**6. Validation of the Python implementation.** Before use, the Python code was compared with R output of a previous analysis of these data: IVW and MR-Egger estimates and standard errors, the Egger intercept, Cochran’s Q and the MR-PRESSO global statistic (RSSobs 186.53 in both) and outlier-corrected estimate were reproduced. For the positive control, the Python and R (TwoSampleMR) analyses gave identical IVW and MR-Egger estimates and the same weighted median to the third decimal; the weighted mode point estimate (1.408 in Python, 1.393 in R) and the bootstrap standard errors of the simple mode differ slightly, and the reason for the difference in the weighted mode was not investigated further (Table S7).')
D.para('**7. Run history.** Table SM2 lists every analysis run, including runs that were superseded, so that the reported results can be traced.')
D.table(('Table SM2','Run history of the analyses.'),['Run','Analysis','Software','Outcome / note'],[
 ['1','Primary NAFLD analysis, 52-variant set, three datasets','Python (custom TwoSampleMR-equivalent code)','Results reported in main text'],
 ['2','MR-PRESSO confirmation','R MRPRESSO 1.0, 10,000 simulations, seed 2026','Same outliers as Python; FinnGen corrected P .250 (R) vs .24 (Python); R values reported'],
 ['3','Sensitivity sets: APOE region excluded; earlier 48-variant set','Python','Reported in Table 4'],
 ['4','Positive control (CAD), Python','Python','Reported (Table S7)'],
 ['5','Positive control, R run 1 (first run)','R 4.4.3, TwoSampleMR 0.6.22','Superseded: logging only; kept in repository'],
 ['6','Positive control, R run 2 (logged, unseeded)','R 4.4.3','Superseded: bootstrap standard errors for mode estimators varied between runs; kept in repository'],
 ['7','Positive control, R run 3 (final, seed 2026)','R 4.4.3','Reported in Table 6'],
 ['8','Extended estimators, influence, radial, I² (three datasets)','R 4.4.3; MendelianRandomization 0.10.0, mr.raps 0.2, RadialMR 1.2.1','Reported in Table 5 and Tables S8–S11'],
 ['9','Manuscript figures and concordance statistics','Python 3 (matplotlib, scipy)','Script code/s13_manuscript_figures.py']],widths=[0.4,2.6,2.0,1.7],size=8)
D.para('**8. Software.** Python 3 (numpy, pandas, scipy, matplotlib); R 4.4.3 (2025-02-28) with TwoSampleMR 0.6.22, MRPRESSO 1.0, MendelianRandomization 0.10.0, mr.raps 0.2 and RadialMR 1.2.1. The study was not pre-registered.')
D.h2('Part B. STROBE-MR checklist')
chk=[
['1','Identify MR as the design in the title or abstract','Title; Abstract'],
['2','Background and rationale, including causal plausibility and why MR is suitable','1. Introduction'],
['3','Specific objectives and causal hypotheses','1. Introduction (final paragraph)'],
['4a','Setting, design, underlying population, dates','2.1; Table 1'],
['4b','Participants, eligibility, sample size, power calculations','2.1; 2.7; Table 1'],
['4c','Measurement, quality control and selection of genetic variants','2.2'],
['4d','Methods and diagnostic criteria for exposure and outcomes','2.1 (NMR remnant cholesterol; NAFLD case definitions of the source GWAS)'],
['4e','Ethics approval and consent','Information block; Acknowledgments (public summary data; no new approval)'],
['5','State the three IV assumptions','2.1'],
['6a–c','Quantitative variables, genetic variants and weights, estimator and covariates','2.1–2.3'],
['6d','Handling of missing data','2.2 (proxies, dropped variants)'],
['6e','Multiple testing','2.3 (no adjustment; one relation in three datasets)'],
['7','Methods to assess the IV assumptions','2.2 (F, Steiger); 2.4; 2.8'],
['8','Sensitivity and additional analyses','2.3–2.8'],
['9a','Software, packages, versions, settings','2.9; SDC 1 Part A (8)'],
['9b','Pre-registration','2.1 (not pre-registered)'],
['10a','Numbers at each study stage and reasons for exclusion','2.2; Fig. 1; SDC 1 Part A (2)'],
['10b','Summary statistics for exposure and outcomes','Table 1'],
['10c','Heterogeneity across meta-analysed studies','Not applicable (summary statistics used as published)'],
['10d','Similarity of variant–exposure associations across samples; overlapping individuals','Overlap: 2.7, 3.6. Similarity of exposure associations across samples was not tested separately (all samples European ancestry); stated as a limitation'],
['11a','Variant–exposure and variant–outcome associations','Fig. 3, 4; Tables S1, S2'],
['11b','MR estimates with uncertainty on an interpretable scale','3.1; Table 2; Fig. 2'],
['11c','Absolute risk translation','Not done (null result)'],
['11d','Plots (forest, scatter)','Fig. 2, 4, 5; Figs. S2, S3'],
['12a, b','Assessment of assumptions; I², Q','3.3; Table 3; Tables S4, S9'],
['13a, b','Robustness to assumption violations; other sensitivity analyses','3.4, 3.5; Tables 4, 5'],
['13c','Assessment of causal direction (bidirectional MR)','Not done; Steiger filter used; stated as a limitation'],
['13d','Comparison with non-MR estimates','4. Discussion (published observational studies)'],
['13e','Additional plots (leave-one-out)','Fig. 5C; Figs. S2, S3'],
['14','Key results with reference to objectives','4. Discussion (first paragraph)'],
['15','Limitations: IV assumptions, bias, imprecision','4. Discussion (limitations paragraph)'],
['16a–c','Cautious interpretation, mechanisms, clinical relevance','4. Discussion; 5. Conclusion'],
['17','Generalizability','4. Discussion (seventh limitation)'],
['18','Funding','Information block'],
['19','Data and code sharing','Data and code availability'],
['20','Conflicts of interest','Information block']]
D.table(('Table SM3','STROBE-MR checklist (item wording paraphrased).'),['Item','Recommendation','Location'],chk,widths=[0.55,3.4,2.8],size=8)
D.pagebreak()
D.h1('Supplemental Digital Content 2')
D.h2('Supplementary tables')
# S1
s1=pd.read_excel('/home/claude/remnant-cholesterol-nafld-mr/results/Supplementary_Tables.xlsx',sheet_name='S1 Instruments')
for c_ in ['EAF','beta_RC_SD','SE_RC','F_stat','R2']: s1[c_]=pd.to_numeric(s1[c_],errors='coerce')
nn=lambda v,fmt:('–' if pd.isna(v) else format(v,fmt))
rows=[[r.SNP,r.chr,f"{int(float(r.pos)):,}",r.effect_allele+'/'+r.other_allele,nn(r.EAF,'.3f'),nn(r.beta_RC_SD,'.3f'),nn(r.SE_RC,'.3f'),f"{r.P_exposure:.1e}" if r.P_exposure>0 else '< 1e-300',nn(r.F_stat,'.0f'),nn(100*r.R2,'.2f'),''.join(['G' if r.used_Ghodsian=='yes' else '–','F' if r.used_Fairfield=='yes' else '–','N' if r.used_FinnGen=='yes' else '–'])] for r in s1.itertuples()]
D.table(('Table S1','The 52 remnant cholesterol instruments.'),['SNP','Chr','Position (build 37)','EA/OA','EAF','β (SD)','SE','P','F','R² (%)','Used*'],rows,widths=[0.95,0.35,0.9,0.5,0.45,0.5,0.45,0.65,0.4,0.45,0.5],size=7,note='*G = Ghodsian, F = Fairfield, N = FinnGen; – = not used (missing without proxy or palindromic). EA = effect (RC-raising) allele; OA = other allele; EAF = effect allele frequency. Variant-level outcome data are in Supplementary_Tables.xlsx (S2).')
D.para('**S2. Harmonised variant-level data** (exposure and outcome effect, standard error, P value, Wald ratio, Steiger direction) for the three NAFLD datasets and the positive control are in Supplementary_Tables.xlsx (sheets S2 and S7) and in the repository folder analysis/primary_clumped52.',size=9)
# S3 summary of all methods x analyses
rr=[]
for o in (G,F,N):
    for an in ['all SNPs','Steiger-filtered','MR-PRESSO outliers removed']:
        for m in ['ivw','egger','wm','sm','wmo']:
            r=g(o,m,an); rr.append([NAME[o],an,{'ivw':'IVW','egger':'MR-Egger','wm':'Weighted median','sm':'Simple mode','wmo':'Weighted mode'}[m],r['n'],f"{r['OR']:.3f}",f"{r['lo']:.3f}–{r['hi']:.3f}",pf(r['p'])])
D.table(('Table S3','All MR estimates for the primary 52-variant set (Python), by analysis set and method.'),['Dataset','Analysis set','Method','SNPs','OR','95% CI','P'],rr,widths=[0.8,1.7,1.2,0.5,0.6,1.1,0.8],size=7.5,note='The MR-PRESSO-corrected IVW estimates reported in the main text (R package) are in Table 4. Estimates for the APOE-excluded and earlier 48-variant sets are in Supplementary_Tables.xlsx and the repository.')
# S4
s4=[]
for ds,o in DS.items():
    h=het.loc[ds]; s4.append([ds,S[o]['n_used'],f"{h.Q_IVW:.1f}",f"{h.Q_IVW_p:.1e}",f"{h.Q_Egger:.1f}",f"{h.Q_Egger_p:.1e}",f"{S[o]['egger_int']:.4f} ({S[o]['egger_int_se']:.4f})",pf(S[o]['egger_int_p']),f"{S[o]['steiger_ivw_or']:.3f}",len(S[o]['steiger_removed'])])
D.table(('Table S4','Heterogeneity, pleiotropy and Steiger filtering.'),['Dataset','SNPs','Q IVW','P','Q Egger','P','Egger intercept (SE)','P','Steiger-filtered IVW OR','Removed'],s4,widths=[0.8,0.4,0.55,0.7,0.6,0.7,1.1,0.5,0.8,0.55],size=8)
# S5 LOO
s5=[]
for ds,o in DS.items():
    l=L[o]; l=l[l.SNP_left_out!='All'].copy(); l['OR']=np.exp(l.b); l['dist']=(l.b-g(o)['b']).abs(); top=l.sort_values('dist',ascending=False).head(3)
    s5.append([ds,len(l),f"{l.OR.min():.3f}–{l.OR.max():.3f}",f"{l.pval.min():.3f}",', '.join(f"{r.SNP_left_out} ({r.OR:.3f})" for r in top.itertuples())])
D.table(('Table S5','Leave-one-out summary for the IVW estimate.'),['Dataset','Runs','OR range','Smallest P','Three most influential omitted variants (OR when omitted)'],s5,widths=[0.8,0.5,1.0,0.7,3.7],size=8,note='Every single-omission estimate is in Supplementary_Tables.xlsx (S5).')
# S6 proxies
pxt=pd.read_excel('/home/claude/remnant-cholesterol-nafld-mr/results/Supplementary_Tables.xlsx',sheet_name='S6 Proxies used')
D.table(('Table S6','Proxy variants used in the NAFLD outcome datasets.'),['Outcome','Target variant','Proxy','r²','Distance (bp)'],[[r.outcome.capitalize() if r.outcome!='finngen' else 'FinnGen',r.target,r.proxy,f"{r.r2:.3f}",f"{r.dist_bp:,}"] for r in pxt.itertuples()],widths=[1.2,1.5,1.5,0.8,1.0],size=8,note='Three variants (rs7534572, rs574454081, rs532754706) had no proxy and were dropped from the Ghodsian and FinnGen analyses.')
# S7 positive control full
d=pd.read_csv('/home/claude/remnant-cholesterol-nafld-mr/positive_control/positive_control_python_vs_R.csv') if False else pd.read_excel('/home/claude/remnant-cholesterol-nafld-mr/results/Supplementary_Tables.xlsx',sheet_name='S7 CAD Python vs R')
rows=[[r.method,f"{r.OR_python:.3f} ({r.OR_lo_python:.3f}–{r.OR_hi_python:.3f})",f"{r.pval_python:.1e}",f"{r.OR_R:.3f} ({r.OR_lo_R:.3f}–{r.OR_hi_R:.3f})",f"{r.pval_R:.1e}"] for r in d.itertuples()]
D.table(('Table S7','Positive control (CAD): Python and R estimates side by side (49 variants).'),['Method','Python OR (95% CI)','P','R OR (95% CI)','P'],rows,widths=[1.8,1.7,0.8,1.7,0.8],size=8,note=f"Cochran’s Q {pcS['Q']:.1f} (48 df), P = {pcS['Q_p']:.1e}; Egger intercept {pcS['egger_int']:.4f} (SE {pcS['egger_int_se']:.4f}), P = {pcS['egger_int_p']:.3f}. MR-PRESSO (R): global P < .0001; 6 outliers (rs115478735, rs12151108, rs12740374, rs1500188, rs3777411, rs653178); corrected OR 1.492 (1.348–1.651); distortion test P = .175. Python: 5 outliers, corrected OR {pcS['presso_corr_or']:.3f}. Steiger removed rs653178 (IVW OR {pcS['steiger_ivw_or']:.3f}). Leave-one-out IVW OR 1.51–1.62, all P < 10^^−8^^. Substitutions for missing instruments: rs115478735 → rs507666 (r² 0.988); rs838146 → rs838147 (r² 1.0); rs72848251 → rs115129597 (same position); rs72843670 → rs115366102 (same position); rs574454081 and rs532754706 dropped. R run 1 and run 2 are kept in the repository.")
# S8 extended
rows=[]
for r in ext.itertuples():
    rows.append([r.dataset,r.method,r.nsnp,f"{r.OR:.3f}",f"{r.OR_lo:.3f}–{r.OR_hi:.3f}",('–' if pd.isna(r.pval) else pf(r.pval)),('' if pd.isna(r.note) else str(r.note).replace('radial outliers:','outliers:'))])
D.table(('Table S8','Extended estimators for all datasets.'),['Dataset','Estimator','SNPs','OR','95% CI','P','Note'],rows,widths=[0.7,2.0,0.4,0.55,1.0,0.6,1.5],size=7,note='The contamination mixture model does not provide a P value; its interval is a likelihood-based confidence interval. MR-Lasso and the contamination mixture model also report the number of variants judged valid.')
# S9 I2
rows=[[ds,h.nsnp,f"{h.Q_IVW:.1f}",f"{100*h.I2_IVW:.1f}",f"{h.Q_Egger:.1f}",f"{100*h.I2_GX:.2f}",f"{h.Egger_intercept:.4f} ({h.Egger_intercept_se:.4f})",pf(h.Egger_intercept_p),f"{h.phi_random_effects:.2f}"] for ds,h in het.iterrows()]
D.table(('Table S9','Heterogeneity, I² and no-measurement-error (NOME) statistics.'),['Dataset','SNPs','Q IVW','I² (%)','Q Egger','I²~~GX~~ (%)','Egger intercept (SE)','P','φ'],rows,widths=[0.8,0.45,0.6,0.55,0.65,0.7,1.3,0.5,0.45],size=8,note='φ = multiplicative random-effects residual standard error (values above 1 indicate overdispersion).')
# S10 influence
fl=infl[infl.cooks_d>infl.cook_threshold_4_over_n].sort_values(['dataset','cooks_d'],ascending=[True,False])
D.table(('Table S10','Variants exceeding the Cook’s distance threshold (4/n).'),['Dataset','SNP',"Cook’s distance",'Studentized residual','Leverage','Threshold'],[[r.dataset,r.SNP,f"{r.cooks_d:.3f}",f"{r.student_resid:.2f}",f"{r.hat:.3f}",f"{r.cook_threshold_4_over_n:.3f}"] for r in fl.itertuples()],widths=[0.9,1.2,1.0,1.2,0.8,0.8],size=8,note='Full diagnostics for all variants are in Supplementary_Tables.xlsx.')
# S11 radial
rows=[]
for ds in DS:
    r=pd.read_csv(f'/home/claude/remnant-cholesterol-nafld-mr/extended_analyses/R_out/radial_data_{ds}.csv'); rows.append([ds,int((r.Outliers=='Outlier').sum()),', '.join(r[r.Outliers=='Outlier'].SNP)])
D.table(('Table S11','Radial MR outliers.'),['Dataset','Outliers','Variants'],rows,widths=[0.9,0.7,5.1],size=8)
# S12 concordance
D.table(('Table S12','Cross-dataset concordance of single-variant Wald ratios.'),['Pair','Variants','Spearman ρ','P','Same sign (%)'],[[f"{r.a} vs {r.b}",r.n,f"{r.spearman:.3f}",pf(r.p),f"{100*r.sign_concordance:.1f}"] for r in conc.itertuples()],widths=[2.8,0.8,0.9,0.9,1.0],size=8,note='Datasets share participants, so concordance is expected to exceed that of independent replications.')
D.h2('Supplementary figures')
D.figure(FIG+'graphical_abstract.png','Graphical abstract.','Fifty-two genetic instruments for remnant cholesterol were tested against three NAFLD datasets (null) and against coronary artery disease as a positive control (strong association).',width=6.4,label='Figure S1')
D.figure(FIG+'figS_diagnostics_Fairfield.png','Diagnostics for remnant cholesterol and NAFLD in the Fairfield (UK Biobank) dataset.','Panel A: scatter plot with IVW, MR-Egger, weighted median and weighted mode fits; panel B: funnel plot; panel C: leave-one-out; panel D: single-variant forest plot.',width=5.8,label='Figure S2')
D.figure(FIG+'figS_diagnostics_FinnGen.png','Diagnostics for remnant cholesterol and NAFLD in the FinnGen release 5 dataset.','Panels as in Fig. S2. The wide intervals reflect the small number of cases (894).',width=5.8,label='Figure S3')
D.figure(FIG+'fig9_strength_heterogeneity_influence.png','Instrument strength, heterogeneity and influence.','Panel A: distribution of per-variant F statistics (dashed line: F = 10). Panel B: I² of the IVW estimate and I²~~GX~~ of the MR-Egger framework in each dataset. Panel C: studentized residuals against Cook’s distance (dashed line: 4/n); the points above the line are rs58542926 (TM6SF2) and rs112875651 (TRIB1) in all datasets.',width=6.6,label='Figure S4')
D.figure(FIG+'fig8_radial.png','Radial MR plots.','Each point is a variant plotted by inverse-standard-error weight and weighted ratio estimate; the line is the radial IVW estimate and crosses mark variants flagged as outliers.',width=6.6,label='Figure S5')
D.figure(FIG+'figS_diagnostics_positive_control.png','Diagnostics for the positive control (remnant cholesterol and coronary artery disease).','Panel A: scatter plot (circled points are MR-PRESSO outliers in the Python analysis); panel B: funnel plot; panel C: leave-one-out analysis, in which no single omission removes the association; panel D: single-variant forest plot.',width=5.8,label='Figure S6')
D.reflist('References (Supplemental Digital Content)') if D.order else None
D.footer('Supplemental Digital Content • Soleymani T et al.')
D.save(OUT+'RC_NAFLD_MR_supplementary_v2.docx'); print('ok',D.fig,D.tab)
