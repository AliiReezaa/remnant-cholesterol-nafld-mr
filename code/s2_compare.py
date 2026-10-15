import pandas as pd, numpy as np, os
W=os.path.expanduser('~/rc_nafld_mr/work')
e=pd.read_csv(f'{W}/exp_sig_pos.tsv',sep='\t')
orig=set(open(f'{W}/orig_union.txt').read().split())
print('orig union',len(orig),'present in sig table:',len(orig & set(e.variant_id)))
print(e[e.variant_id.isin(orig)][['variant_id','chr','pos','p_value']].sort_values('p_value').head(3))
pos=e.dropna(subset=['pos']).copy(); pos['pos']=pos.pos.astype(int); pos['chr']=pos.chr.astype(str)
def greedy(df,win):
    leads=[]
    for c,g in df.groupby('chr'):
        g=g.sort_values('p_value'); taken=[]
        for r in g.itertuples():
            if all(abs(r.pos-t)>win for t in taken): taken.append(r.pos); leads.append(r.variant_id)
    return set(leads)
for win in [100_000,250_000,500_000,1_000_000,2_500_000,10_000_000]:
    L=greedy(pos,win)
    print(f'window ±{win//1000:>6} kb: leads={len(L):3d}  orig SNPs that are leads={len(orig&L):2d}/{len(orig)}  leads not in orig={len(L-orig)}')
