import numpy as np, pandas as pd, sys, os
W=os.path.expanduser('~/rc_nafld_mr/work'); L=W+'/ld'
sig=pd.read_csv(W+'/exp_sig_pos.tsv',sep='\t')
cand=set(sig.variant_id)
rows=[]
with open(L+'/EUR.bim') as f:
    for i,line in enumerate(f):
        p=line.split('\t')
        if p[1] in cand: rows.append((i,p[0],p[1],int(p[3]),p[4],p[5].strip()))
bim=pd.DataFrame(rows,columns=['idx','chr','rsid','pos','a1','a2']).drop_duplicates('rsid')
print('candidates',len(cand),'in reference',len(bim))
nind=sum(1 for _ in open(L+'/EUR.fam')); nb=(nind+3)//4
mm=np.memmap(L+'/EUR.bed',dtype=np.uint8,mode='r',offset=3)
lut=np.array([2,np.nan,1,0],dtype=np.float32)   # count of A1 allele (bim a1)
def geno(i):
    b=np.asarray(mm[i*nb:(i+1)*nb]); c=np.stack([(b>>s)&3 for s in (0,2,4,6)],1).reshape(-1)[:nind]
    return lut[c]
G={}
m=sig.set_index('variant_id')
bim['p']=bim.rsid.map(m.p_value)
bim=bim.sort_values('p').reset_index(drop=True)
for r in bim.itertuples(): 
    g=geno(r.idx); mu=np.nanmean(g); g=np.where(np.isnan(g),mu,g); G[r.rsid]=g-mu
keep=[]; removed=set()
for r in bim.itertuples():
    if r.rsid in removed: continue
    if r.p>=5e-8: break
    keep.append(r.rsid); removed.add(r.rsid)
    gi=G[r.rsid]; vi=gi@gi
    near=bim[(bim.chr==r.chr)&(abs(bim.pos-r.pos)<=10_000_000)]
    for q in near.itertuples():
        if q.rsid in removed: continue
        gj=G[q.rsid]; vj=gj@gj
        if vi==0 or vj==0: continue
        r2=(gi@gj)**2/(vi*vj)
        if r2>0.001: removed.add(q.rsid)
out=bim[bim.rsid.isin(keep)][['rsid','chr','pos','p']]
print('clumped index SNPs',len(out))
out.to_csv(W+'/instr_clumped.tsv',sep='\t',index=False)
open(W+'/instr_clumped_ids.txt','w').write('\n'.join(out.rsid)+'\n')
orig=set(open(W+'/orig_union.txt').read().split())
print('overlap with provisional 48:',len(orig&set(out.rsid)),' in clumped not in 48:',len(set(out.rsid)-orig),' in 48 not clumped:',sorted(orig-set(out.rsid)))
