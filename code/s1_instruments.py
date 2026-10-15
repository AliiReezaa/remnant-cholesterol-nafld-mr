import pandas as pd, numpy as np, gzip, sys
import os
H=os.path.expanduser('~/rc_nafld_mr/inputs'); W=os.path.expanduser('~/rc_nafld_mr/work')
e=pd.read_csv(f'{W}/exp_sig.tsv',sep='\t')
print('sig rows',len(e),'unique ids',e.variant_id.nunique())
e=e.drop_duplicates('variant_id',keep='first')
print('after dedupe (keep first, as TwoSampleMR does):',len(e),'(original 8873)')
# positions from Ghodsian file (chr,pos,rsID)
g=pd.read_csv(f'{H}/GCST90091033_buildGRCh37.tsv.gz',sep='\t',usecols=['variant_id','chromosome','base_pair_location'])
g=g.drop_duplicates('variant_id')
m=e.merge(g,on='variant_id',how='left')
print('with positions from Ghodsian:',m.chromosome.notna().sum(),'of',len(m))
m.to_csv(f'{W}/exp_sig_pos.tsv',sep='\t',index=False)
