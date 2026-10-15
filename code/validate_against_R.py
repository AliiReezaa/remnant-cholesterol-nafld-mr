"""Validate the Python MR code against the original R output for the same data (dataPresso, GCST90091033)."""
import sys, os; sys.path.insert(0, os.path.expanduser('~/rc_nafld_mr/code'))
import pandas as pd, numpy as np, mr_lib as M
d = pd.read_csv(os.path.expanduser('~/rc_nafld_mr/inputs/original_R_output/dataPresso_ebi-a-GCST90091033.txt'), sep='\t')
bx, by, sex, sey = d['beta.exposure'].values, d['beta.outcome'].values, d['se.exposure'].values, d['se.outcome'].values
print('n SNPs:', len(d))
R = {  # copied from the original R log (res_Presso, GCST90091033, 44 SNPs)
 'MR Egger': (-0.04231471, 0.08015224, 0.6003267), 'Weighted median': (-0.03902608, 0.06785680, 0.5652074),
 'Inverse variance weighted': (-0.01531351, 0.05249555, 0.7705079), 'Simple mode': (-0.04656432, 0.11244334, 0.6808494),
 'Weighted mode': (-0.02270068, 0.06586977, 0.7320527)}
res = M.run_all(bx, by, sex, sey)
print(f"{'method':28s} {'b_py':>9s} {'b_R':>9s} | {'se_py':>8s} {'se_R':>8s} | {'p_py':>7s} {'p_R':>7s}")
for k, v in res.items():
    print(f"{k:28s} {v['b']:9.5f} {R[k][0]:9.5f} | {v['se']:8.5f} {R[k][1]:8.5f} | {v['pval']:7.4f} {R[k][2]:7.4f}")
e = M.egger(bx, by, sey); i = M.ivw(bx, by, sey)
print('Egger intercept  py: %.6f se %.6f p %.4f | R: 0.002344596 se 0.00522188 p 0.6557438' % (e['intercept'], e['intercept_se'], e['intercept_pval']))
print('Q IVW py: %.4f (p %.4f) | R: 52.19335 (p 0.1588615)   Q Egger py: %.4f | R: 51.94402' % (i['Q'], i['Q_pval'], e['Q']))
