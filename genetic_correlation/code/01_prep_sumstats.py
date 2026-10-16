"""Convert each GWAS Catalog harmonised file to a munge_sumstats-ready table (SNP A1 A2 BETA P FRQ [INFO]).
Z = BETA/SE is computed at full precision and re-expressed as P so that munge_sumstats recovers Z exactly
(some catalog files store P with only 2 significant figures).
A1 = effect allele, A2 = other allele, exactly as in the original (non-hm_) columns of each file.
Fairfield uses the ORIGINAL author file (the harmonised file lacks the INFO column) and provides OR and 95% CI only: BETA = ln(OR), SE = (ln(UL)-ln(LL))/(2*1.959964), as in the MR pipeline (code/s3_extract.py)."""
import pandas as pd, numpy as np, sys, os
from scipy.stats import norm
S='/home/claude/rg_work/sumstats'; O='/home/claude/rg_work/munge_in'; os.makedirs(O,exist_ok=True)
SPEC={
 'RC':      dict(f='RC_GCST90092943.h.tsv.gz',      cols=['variant_id','effect_allele','other_allele','effect_allele_frequency','beta','standard_error']),
 'Ghodsian':dict(f='Ghodsian_GCST90091033.h.tsv.gz',cols=['variant_id','effect_allele','other_allele','effect_allele_frequency','beta','standard_error']),
 'CAD':     dict(f='CAD_GCST003116.h.tsv.gz',       cols=['variant_id','effect_allele','other_allele','effect_allele_frequency','beta','standard_error']),
 'Fairfield':dict(f='Fairfield_GCST90054782_buildGRCh37.tsv',cols=['variant_id','effect_allele','other_allele','effect_allele_frequency','odds_ratio','ci_lower','ci_upper']),
}
def run(name):
    sp=SPEC[name]; parts=[]; n0=0
    for ch in pd.read_csv(f'{S}/{sp["f"]}',sep='\t',usecols=sp['cols']+(['info'] if name=='Fairfield' else []),chunksize=1_000_000,dtype={'variant_id':str,'effect_allele':str,'other_allele':str}):
        n0+=len(ch)
        ch=ch[ch.variant_id.astype(str).str.match(r'^rs\d+$',na=False)]
        if name=='Fairfield':
            ch=ch[(ch.odds_ratio>0)&(ch.ci_lower>0)&(ch.ci_upper>0)]
            ch['beta']=np.log(ch.odds_ratio); ch['standard_error']=(np.log(ch.ci_upper)-np.log(ch.ci_lower))/(2*1.959964)
        ch=ch.dropna(subset=['beta','standard_error','effect_allele','other_allele'])
        ch=ch[ch.standard_error>0]
        z=ch.beta/ch.standard_error
        out=pd.DataFrame({'SNP':ch.variant_id,'A1':ch.effect_allele.str.upper(),'A2':ch.other_allele.str.upper(),
                          'BETA':ch.beta,'P':2*norm.sf(np.abs(z)),'FRQ':ch.effect_allele_frequency})
        if name=='Fairfield': out['INFO']=ch['info']
        parts.append(out)
    d=pd.concat(parts); d=d[(d.P>0)]
    d.to_csv(f'{O}/{name}.tsv.gz',sep='\t',index=False,float_format='%.17g')
    print(f'{name}: rows read {n0:,}; rs-ID rows with valid beta/SE/P kept {len(d):,}; unique SNP {d.SNP.nunique():,}',flush=True)
for n in (sys.argv[1:] or SPEC): run(n)
