import pandas as pd, os
H=os.path.expanduser('~/rc_nafld_mr/inputs'); W=os.path.expanduser('~/rc_nafld_mr/work')
dups=set(open(f'{W}/exp_dup_ids.txt').read().split())
# first-instance rule used by TwoSampleMR::read_exposure_data: keep the first row for each duplicated SNP
seen=set(); keep_first={}
rows=0
for ch in pd.read_csv(f'{H}/GCST90092943_buildGRCh37.tsv.gz',sep='\t',chunksize=1_000_000,dtype={'p_value':str}):
    ch['rowi']=range(rows,rows+len(ch)); rows+=len(ch)
    sub=ch[ch.variant_id.isin(dups)]
    for r in sub.itertuples():
        if r.variant_id not in keep_first: keep_first[r.variant_id]=r.rowi
print('rows',rows)
# significant rows with row index
sig=[]
rows=0
for ch in pd.read_csv(f'{H}/GCST90092943_buildGRCh37.tsv.gz',sep='\t',chunksize=1_000_000,dtype={'p_value':str}):
    ch['rowi']=range(rows,rows+len(ch)); rows+=len(ch)
    ch['p']=pd.to_numeric(ch.p_value,errors='coerce')
    s=ch[ch.p<5e-8]
    sig.append(s)
sig=pd.concat(sig)
def ok(r):
    return (r.variant_id not in keep_first) or keep_first[r.variant_id]==r.rowi
sig=sig[[ok(r) for r in sig.itertuples()]]
print('significant after first-instance dedupe:',len(sig),'(original run: 8873)')
sig=sig.drop(columns=['p_value']).rename(columns={'p':'p_value'})
# positions: Ghodsian + Fairfield + FinnGen VCF
pos={}
g=pd.read_csv(f'{H}/GCST90091033_buildGRCh37.tsv.gz',sep='\t',usecols=['variant_id','chromosome','base_pair_location'])
for r in g.itertuples(index=False): pos.setdefault(r[0],(str(r[1]),int(r[2])))
need=set(sig.variant_id)-set(pos)
print('missing positions after Ghodsian:',len(need))
f=pd.read_csv(f'{H}/GCST90054782_buildGRCh37.tsv',sep='\t',usecols=['variant_id','chromosome','base_pair_location'])
for r in f.itertuples(index=False):
    if r[0] in need: pos.setdefault(r[0],(str(r[1]),int(r[2])))
need=set(sig.variant_id)-set(pos); print('missing after Fairfield:',len(need))
import gzip
with gzip.open(f'{H}/finn-b-NAFLD.vcf.gz','rt') as fh:
    for line in fh:
        if line[0]=='#': continue
        p=line.split('\t',4)
        if p[2] in need: pos.setdefault(p[2],(p[0],int(p[1])))
need=set(sig.variant_id)-set(pos); print('missing after FinnGen VCF:',len(need))
sig['chr']=[pos.get(v,(None,None))[0] for v in sig.variant_id]; sig['pos']=[pos.get(v,(None,None))[1] for v in sig.variant_id]
sig.drop(columns=['rowi']).to_csv(f'{W}/exp_sig_pos.tsv',sep='\t',index=False)
print(sig.head(3))
