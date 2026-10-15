import pandas as pd, numpy as np, os, sys
H=os.path.expanduser('~/rc_nafld_mr/inputs'); W=os.path.expanduser('~/rc_nafld_mr/work')
which=sys.argv[1]
sig=pd.read_csv(f'{W}/exp_sig_pos.tsv',sep='\t'); ids=set(open(W+'/proxy_ids.txt').read().split())
if which=='ghodsian':
    out=[]
    for ch in pd.read_csv(f'{H}/GCST90091033_buildGRCh37.tsv.gz',sep='\t',chunksize=1_000_000,
                          usecols=['variant_id','chromosome','base_pair_location','effect_allele','other_allele','effect_allele_frequency','beta','standard_error','p_value']):
        out.append(ch[ch.variant_id.isin(ids)])
    o=pd.concat(out); o.columns=['rsid','chr','pos','ea','oa','eaf','beta','se','p']
elif which=='fairfield':
    out=[]
    for ch in pd.read_csv(f'{H}/GCST90054782_buildGRCh37.tsv',sep='\t',chunksize=1_000_000,
                          usecols=['variant_id','chromosome','base_pair_location','effect_allele','other_allele','effect_allele_frequency','odds_ratio','ci_lower','ci_upper','p_value','info']):
        out.append(ch[ch.variant_id.isin(ids)])
    o=pd.concat(out)
    o['beta']=np.log(o.odds_ratio); o['se']=(np.log(o.ci_upper)-np.log(o.ci_lower))/(2*1.959964)
    o=o[['variant_id','chromosome','base_pair_location','effect_allele','other_allele','effect_allele_frequency','beta','se','p_value','info']]
    o.columns=['rsid','chr','pos','ea','oa','eaf','beta','se','p','info']
elif which=='finngen':
    out=[]
    for ch in pd.read_csv(f'{H}/finn-b-NAFLD.vcf.gz',sep='\t',comment='#',header=None,chunksize=1_000_000,
                          usecols=[0,1,2,3,4,8,9],names=['chr','pos','rsid','ref','alt','fmt','val'],dtype=str):
        ch=ch.rename(columns={0:'chr',1:'pos',2:'rsid',3:'ref',4:'alt',8:'fmt',9:'val'}) if 0 in ch.columns else ch
        sub=ch[ch.rsid.isin(ids)]
        out.append(sub)
    o=pd.concat(out)
    f=o.val.str.split(':',expand=True)
    o['beta']=f[0].astype(float); o['se']=f[1].astype(float); o['p']=10**(-f[2].astype(float)); o['eaf']=f[3].astype(float)
    o=o.rename(columns={'alt':'ea','ref':'oa'})[['rsid','chr','pos','ea','oa','eaf','beta','se','p']]
print(which,len(o),'rows matched of',len(ids),'significant exposure SNPs; unique rsids',o.rsid.nunique())
o.to_csv(f'{W}/out_proxy_{which}.tsv',sep='\t',index=False)
