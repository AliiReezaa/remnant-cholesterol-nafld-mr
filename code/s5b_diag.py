import numpy as np, pandas as pd, os
exec(open(os.path.expanduser('~/rc_nafld_mr/code/s5_clump.py')).read().split("keep=[]")[0])
clump=pd.read_csv(W+'/instr_clumped.tsv',sep='\t',dtype={'chr':str}).set_index('rsid')
orig=set(open(W+'/orig_union.txt').read().split())
def r2(a,b): ga,gb=G[a],G[b]; return (ga@gb)**2/((ga@ga)*(gb@gb))
b=bim.set_index('rsid')
for s in ['rs111278137','rs117310449','rs118147862','rs41289512','rs61679753']:
    if s not in b.index: print(s,'NOT IN REFERENCE'); continue
    near=[c for c in clump.index if clump.loc[c,'chr']==b.loc[s,'chr'] and abs(clump.loc[c,'pos']-b.loc[s,'pos'])<=10_000_000]
    print(s,'chr',b.loc[s,'chr'],'pos',b.loc[s,'pos'],'p=%.2g'%b.loc[s,'p'],'| nearby clumped:',[(c,'p=%.2g'%clump.loc[c,'p'],'r2=%.4f'%r2(s,c)) for c in near])
print('--- new 9 (in clumped, not in 48):')
for s in sorted(set(clump.index)-orig): print(s,clump.loc[s,'chr'],clump.loc[s,'pos'],'p=%.2g'%clump.loc[s,'p'])
