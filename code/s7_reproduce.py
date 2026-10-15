import sys, os, itertools; sys.path.insert(0, os.path.expanduser('~/rc_nafld_mr/code'))
import numpy as np, pandas as pd, mr_lib as M
W=os.path.expanduser('~/rc_nafld_mr/work'); R=os.path.expanduser('~/rc_nafld_mr/inputs/original_R_output')
o=pd.read_csv(f'{R}/dataPresso_ebi-a-GCST90091033.txt',sep='\t')
print(o.columns.tolist()[:20], len(o))
h=pd.read_csv(os.path.expanduser('~/rc_nafld_mr/analysis/primary_clumped52/harmonised_GCST90091033.tsv'),sep='\t').set_index('SNP')
o=o.set_index('SNP')
common=[s for s in o.index if s in h.index]; print('orig 44 in my harmonised:',len(common))
c=pd.DataFrame({'bx_o':o.loc[common,'beta.exposure'],'bx_m':h.loc[common,'bx'],'by_o':o.loc[common,'beta.outcome'],'by_m':h.loc[common,'by'],'sey_o':o.loc[common,'se.outcome'],'sey_m':h.loc[common,'sey']})
print('max abs diff bx %.2e by %.2e sey %.2e'%((c.bx_o-c.bx_m).abs().max(),(c.by_o-c.by_m).abs().max(),(c.sey_o-c.sey_m).abs().max()))
print(c[((c.by_o-c.by_m).abs()>1e-4)|((c.bx_o-c.bx_m).abs()>1e-4)].round(4))
print('orig SNPs missing in mine:',[s for s in o.index if s not in h.index])
# try to rebuild the original 47: 44 + 3 of my SNPs not in the 44
extra=[s for s in h.index if s not in o.index]; print('my extra SNPs:',extra)
bx0,by0,se0=o['beta.exposure'].values,o['beta.outcome'].values,o['se.outcome'].values
best=[]
for tri in itertools.combinations(extra,3):
    bx=np.r_[bx0,h.loc[list(tri),'bx'].values]; by=np.r_[by0,h.loc[list(tri),'by'].values]; se=np.r_[se0,h.loc[list(tri),'sey'].values]
    r=M.ivw(bx,by,se); best.append((abs(r['b']+0.006722717)+abs(r['se']-0.09016707),tri,r['b'],r['se'],r['Q']))
best.sort(); 
for x in best[:3]: print(x)

tri=['rs112875651','rs58542926','rs4665972']
bx=np.r_[bx0,h.loc[tri,'bx'].values]; by=np.r_[by0,h.loc[tri,'by'].values]; se=np.r_[se0,h.loc[tri,'sey'].values]
sx=np.r_[o['se.exposure'].values,h.loc[tri,'sex'].values]
pr=M.presso(bx,by,sx,se,nsim=10000,seed=1)
print({k:(v if k not in('keep',) else None) for k,v in pr.items() if k!='corrected'})
print('corrected',pr['corrected'])
