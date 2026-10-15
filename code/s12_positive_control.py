"""POSITIVE CONTROL: remnant cholesterol (GCST90092943) -> coronary artery disease (Nikpay 2015). Same pipeline as s4.
usage: s4_analysis.py <instrument_id_file> <outdir>"""
import sys, os; sys.path.insert(0, os.path.expanduser('~/scripts'))
import numpy as np, pandas as pd, json
from scipy import stats
import mr_lib as M
W = os.path.expanduser('~/work'); H = os.path.expanduser('~/mnt/Fatty Liver')
ids = [x.strip() for x in open(sys.argv[1]) if x.strip()]; OUT = sys.argv[2]; os.makedirs(OUT, exist_ok=True)

N_EXP = 115082                                  # Richardson 2022, GCST90092943 (GWAS Catalog)
OUTCOMES = {                                    # positive control: CAD, Nikpay 2015 CARDIoGRAMplusC4D (GCST003116)
 'CAD_Nikpay2015': ('cad_px', 60801, 123504, 'Nikpay 2015 CARDIoGRAMplusC4D 1000G, 60,801 cases / 123,504 controls'),
}
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
def is_pal(a, b): return {a, b} in ({'A', 'T'}, {'C', 'G'})

# ---- exposure rows for the chosen instruments
sig = pd.read_csv(f'{W}/exp_sig_pos.tsv', sep='\t').set_index('variant_id')
exp = sig.loc[[i for i in ids if i in sig.index]].reset_index()
exp = exp.rename(columns={'variant_id': 'SNP', 'effect_allele': 'ea_x', 'other_allele': 'oa_x', 'effect_allele_frequency': 'eaf_x',
                          'beta': 'bx', 'standard_error': 'sex', 'p_value': 'p_x'})
exp['F'] = (exp.bx / exp.sex)**2
exp['r2_x'] = exp.F / (exp.F + N_EXP - 2)

def harmonise(e, o):
    rows = []; log = {'missing_in_outcome': [], 'palindromic_removed': [], 'allele_mismatch': []}
    om = o.drop_duplicates('rsid').set_index('rsid')
    for r in e.itertuples():
        if r.SNP not in om.index: log['missing_in_outcome'].append(r.SNP); continue
        q = om.loc[r.SNP]; ea, oa, qea, qoa = r.ea_x, r.oa_x, str(q.ea), str(q.oa)
        beta, eaf = q.beta, q.eaf
        pal = len(ea) == 1 and len(oa) == 1 and is_pal(ea, oa)
        if (ea, oa) == (qea, qoa): pass
        elif (ea, oa) == (qoa, qea): beta, eaf = -beta, 1 - eaf
        else:
            if len(ea) == 1 and len(oa) == 1 and not pal and (ea, oa) == (COMP.get(qea), COMP.get(qoa)): pass
            elif len(ea) == 1 and len(oa) == 1 and not pal and (ea, oa) == (COMP.get(qoa), COMP.get(qea)): beta, eaf = -beta, 1 - eaf
            else: log['allele_mismatch'].append(r.SNP); continue
        if pal:                                    # TwoSampleMR action 2: drop ambiguous palindromes, orient the rest by frequency
            if 0.42 <= r.eaf_x <= 0.58 or 0.42 <= eaf <= 0.58: log['palindromic_removed'].append(r.SNP); continue
            if (r.eaf_x - 0.5) * (eaf - 0.5) < 0:
                if (ea, oa) == (qea, qoa): beta, eaf = -beta, 1 - eaf   # strand/orientation resolved by frequency
        rows.append(dict(SNP=r.SNP, chr=r.chr, pos=r.pos, ea=ea, oa=oa, eaf_x=r.eaf_x, bx=r.bx, sex=r.sex, p_x=r.p_x, F=r.F,
                         r2_x=r.r2_x, by=beta, sey=q.se, p_y=q.p, eaf_y=eaf))
    return pd.DataFrame(rows), log

def steiger(d, ncase, ncontrol):
    neff = 4 / (1 / ncase + 1 / ncontrol)
    g = 2 * d.eaf_y * (1 - d.eaf_y) * d.by**2
    ve = g / (g + np.pi**2 / 3)                    # variance explained on the latent logistic scale
    r_o = np.sqrt(ve); r_x = np.sqrt(d.r2_x)
    z = (np.arctanh(r_x) - np.arctanh(r_o)) / np.sqrt(1 / (N_EXP - 3) + 1 / (neff - 3))
    return (d.r2_x > ve).values, 2 * stats.norm.sf(np.abs(z)), neff, ve.values

summary = []; allres = []
for lab, (key, nca, nco, desc) in OUTCOMES.items():
    o = pd.read_csv(f'{W}/out_{key}' + os.environ.get('SUF','') + '.tsv', sep='\t')
    d, log = harmonise(exp, o)
    d['ratio'] = d.by / d.bx; d['ratio_se'] = d.sey / d.bx.abs()
    steig_ok, steig_p, neff, ve = steiger(d, nca, nco)
    d['steiger_dir'] = steig_ok; d['steiger_p'] = steig_p; d['r2_y_latent'] = ve
    bx, by, sex, sey = d.bx.values, d.by.values, d.sex.values, d.sey.values
    res = M.run_all(bx, by, sex, sey, nboot=1000, seed=7)
    e = M.egger(bx, by, sey); iv = M.ivw(bx, by, sey)
    pr = M.presso(bx, by, sex, sey, nsim=10000, seed=7)
    keep = d.steiger_dir.values
    st = M.ivw(bx[keep], by[keep], sey[keep]) if keep.sum() > 2 else None
    st_all = M.run_all(bx[keep], by[keep], sex[keep], sey[keep], nboot=500, seed=7)
    loo = M.leave_one_out(bx, by, sey, d.SNP.tolist())
    # strength, power, overlap bias
    R2 = d.r2_x.sum(); k = len(d); Foverall = (N_EXP - k - 1) / k * R2 / (1 - R2)
    b_fix, se_fix = M.ivw_fixed(bx, by, sey); mde_b = (1.959964 + 0.841621) * max(se_fix, iv['se'])
    out = dict(outcome=lab, desc=desc, ncase=nca, ncontrol=nco, neff=neff, n_instruments_in=len(exp), n_used=k,
               missing=len(log['missing_in_outcome']), palindromic_removed=len(log['palindromic_removed']),
               mean_F=float(d.F.mean()), min_F=float(d.F.min()), R2=float(R2), F_overall=float(Foverall),
               ivw_b=iv['b'], ivw_se=iv['se'], ivw_p=iv['pval'], ivw_or=float(np.exp(iv['b'])),
               ivw_or_lo=float(np.exp(iv['b'] - 1.959964 * iv['se'])), ivw_or_hi=float(np.exp(iv['b'] + 1.959964 * iv['se'])),
               Q=iv['Q'], Q_df=iv['Q_df'], Q_p=iv['Q_pval'], egger_int=e['intercept'], egger_int_se=e['intercept_se'], egger_int_p=e['intercept_pval'],
               presso_global_p=pr['global_p'], presso_outliers=[d.SNP.iloc[i] for i in pr['outlier_idx']],
               presso_corr_or=float(np.exp(pr['corrected']['b'])), presso_corr_p=pr['corrected']['pval'], presso_corr_n=pr['corrected']['nsnp'],
               steiger_removed=[d.SNP.iloc[i] for i in np.where(~keep)[0]], steiger_ivw_or=(float(np.exp(st['b'])) if st else None),
               steiger_ivw_p=(st['pval'] if st else None),
               mde_or_80pct_power=float(np.exp(mde_b)), mde_or_80pct_power_protective=float(np.exp(-mde_b)),
               overlap_relative_bias_if_full=float(1 / Foverall), log=log)
    summary.append(out)
    rows = []
    for m, v in res.items():
        rows.append(dict(outcome=lab, analysis='all SNPs', method=m, nsnp=v['nsnp'], b=v['b'], se=v['se'], pval=v['pval'],
                         OR=np.exp(v['b']), OR_lo=np.exp(v['b'] - 1.959964 * v['se']), OR_hi=np.exp(v['b'] + 1.959964 * v['se'])))
    for m, v in st_all.items():
        rows.append(dict(outcome=lab, analysis='Steiger-filtered', method=m, nsnp=v['nsnp'], b=v['b'], se=v['se'], pval=v['pval'],
                         OR=np.exp(v['b']), OR_lo=np.exp(v['b'] - 1.959964 * v['se']), OR_hi=np.exp(v['b'] + 1.959964 * v['se'])))
    ke = pr['keep']
    pres = M.run_all(bx[ke], by[ke], sex[ke], sey[ke], nboot=500, seed=7)
    for m, v in pres.items():
        rows.append(dict(outcome=lab, analysis='MR-PRESSO outliers removed', method=m, nsnp=v['nsnp'], b=v['b'], se=v['se'], pval=v['pval'],
                         OR=np.exp(v['b']), OR_lo=np.exp(v['b'] - 1.959964 * v['se']), OR_hi=np.exp(v['b'] + 1.959964 * v['se'])))
    allres += rows
    d.to_csv(f'{OUT}/harmonised_{lab}.tsv', sep='\t', index=False)
    pd.DataFrame(loo, columns=['SNP_left_out', 'b', 'se', 'pval']).to_csv(f'{OUT}/leave_one_out_{lab}.tsv', sep='\t', index=False)
pd.DataFrame(allres).to_csv(f'{OUT}/mr_results.tsv', sep='\t', index=False)
json.dump(summary, open(f'{OUT}/summary.json', 'w'), indent=1, default=str)
for s in summary:
    print(s['outcome'], '| used', s['n_used'], '/', s['n_instruments_in'], '| IVW OR %.3f (%.3f-%.3f) p=%.3f' % (s['ivw_or'], s['ivw_or_lo'], s['ivw_or_hi'], s['ivw_p']),
          '| Q p=%.2g | Egger int p=%.3f | PRESSO global p=%.4f outliers=%s | MDE OR=%.2f' % (s['Q_p'], s['egger_int_p'], s['presso_global_p'], s['presso_outliers'], s['mde_or_80pct_power']))
