import numpy as np, pandas as pd, os
W=os.path.expanduser('~/rc_nafld_mr/work'); L=W+'/ld'
targets=['rs247616','rs7534572','rs115478735','rs574454081','rs72848251','rs56325564','rs12804228','rs838146','rs532754706']
sig=pd.read_csv(W+'/exp_sig_pos.tsv',sep='\t').set_index('variant_id')
nind=503; nb=(nind+3)//4
mm=np.memmap(L+'/EUR.bed',dtype=np.uint8,mode='r',offset=3)
lut=np.array([2,np.nan,1,0],dtype=np.float32)
def geno(i):
    b=np.asarray(mm[i*nb:(i+1)*nb]); c=np.stack([(b>>s)&3 for s in (0,2,4,6)],1).reshape(-1)[:nind]
    g=lut[c]; mu=np.nanmean(g); return np.where(np.isnan(g),mu,g)-mu
tset=set(targets); tpos={}
with open(L+'/EUR.bim') as f:
    for i,line in enumerate(f):
        p=line.split('\t')
        if p[1] in tset: tpos[p[1]]=(p[0],int(p[3]))
rows=[]
with open(L+'/EUR.bim') as f:
    for i,line in enumerate(f):
        p=line.split('\t')
        for t,(c,ps) in tpos.items():
            if p[0]==c and abs(int(p[3])-ps)<=250000:
                rows.append((i,p[0],p[1],int(p[3]),p[4],p[5].strip())); break
bim=pd.DataFrame(rows,columns=['idx','chr','rsid','pos','a1','a2']).drop_duplicates('rsid'); print('bim window rows',len(bim),'targets found',len(tpos))
res=[]
for t in targets:
    r=bim[bim.rsid==t]
    if r.empty: print(t,'NOT in reference'); continue
    r=r.iloc[0]
    win=bim[(bim.chr==r.chr)&(abs(bim.pos-r.pos)<=250000)&(bim.rsid!=t)&(bim.a1.str.len()==1)&(bim.a2.str.len()==1)&(bim.rsid.str.startswith('rs'))]
    gt=geno(r.idx); vt=gt@gt
    for q in win.itertuples():
        gq=geno(q.idx); vq=gq@gq
        if vq==0: continue
        c=(gt@gq)/np.sqrt(vt*vq)
        if c*c>=0.8: res.append(dict(target=t,proxy=q.rsid,ppos=q.pos,r2=c*c,sign=int(np.sign(c)),t_a1=r.a1,t_a2=r.a2,p_a1=q.a1,p_a2=q.a2,tpos=r.pos))
px=pd.DataFrame(res); px.to_csv(W+'/proxy_table.tsv',sep='\t',index=False)
print(px.groupby('target').size())
open(W+'/proxy_ids.txt','w').write('\n'.join(sorted(set(px.proxy)))+'\n')
print(len(set(px.proxy)),'proxy ids')
