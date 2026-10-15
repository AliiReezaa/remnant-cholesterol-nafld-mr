import pandas as pd, numpy as np, os
W=os.path.expanduser('~/rc_nafld_mr/work')
px=pd.read_csv(W+'/proxy_table.tsv',sep='\t')
sig=pd.read_csv(W+'/exp_sig_pos.tsv',sep='\t').set_index('variant_id')
COMP={'A':'T','T':'A','C':'G','G':'C'}
used=[]
for key in ['ghodsian','finngen']:
    base=pd.read_csv(f'{W}/out_{key}.tsv',sep='\t'); pr=pd.read_csv(f'{W}/out_proxy_{key}.tsv',sep='\t').drop_duplicates('rsid').set_index('rsid')
    new=[]
    for t in px.target.unique():
        if t in set(base.rsid): continue
        e=sig.loc[t]; ea,oa=e.effect_allele,e.other_allele
        cand=px[px.target==t].sort_values('r2',ascending=False)
        for c in cand.itertuples():
            if c.proxy not in pr.index: continue
            q=pr.loc[c.proxy]
            # palindromic proxy with ambiguous frequency is unusable
            if {c.p_a1,c.p_a2} in ({'A','T'},{'C','G'}) and 0.42<=q.eaf<=0.58: continue
            # orient target alleles to the reference panel strand
            if {ea,oa}=={c.t_a1,c.t_a2}: tea=ea
            elif len(ea)==1 and {COMP.get(ea),COMP.get(oa)}=={c.t_a1,c.t_a2}: tea=COMP[ea]
            else: continue
            same = (tea==c.t_a1)
            linked = (c.p_a1 if c.sign>0 else c.p_a2) if same else (c.p_a2 if c.sign>0 else c.p_a1)
            qea,qoa=str(q.ea),str(q.oa)
            if qea==linked: b,f=q.beta,q.eaf
            elif qoa==linked: b,f=-q.beta,1-q.eaf
            elif COMP.get(qea)==linked: b,f=q.beta,q.eaf
            elif COMP.get(qoa)==linked: b,f=-q.beta,1-q.eaf
            else: continue
            row=dict(rsid=t,chr=e.chr,pos=e.pos,ea=ea,oa=oa,eaf=(f if tea==ea else 1-f) if False else f,beta=b,se=q.se,p=q.p)
            if 'info' in base.columns: row['info']=np.nan
            new.append(row); used.append(dict(outcome=key,target=t,proxy=c.proxy,r2=round(c.r2,3),dist_bp=abs(int(c.ppos)-int(c.tpos))))
            break
    pd.concat([base,pd.DataFrame(new)]).to_csv(f'{W}/out_{key}_px.tsv',sep='\t',index=False)
    print(key,'proxies added:',len(new))
u=pd.DataFrame(used); u.to_csv(W+'/proxies_used.tsv',sep='\t',index=False); print(u.to_string(index=False))
