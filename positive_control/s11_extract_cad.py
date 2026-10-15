"""Positive control, step 1: pull the 52 remnant-cholesterol instruments (and r2>=0.8 proxies where missing) out of the
Nikpay 2015 CARDIoGRAMplusC4D CAD GWAS (GCST003116, harmonised build 37). Writes ~/work/out_cad_px.tsv in the same
format as the NAFLD outcome tables (rsid,chr,pos,ea,oa,eaf,beta,se,p)."""
import os, numpy as np, pandas as pd, subprocess, sys
W=os.path.expanduser('~/work'); PC=os.path.expanduser('~/pc')
F=os.path.expanduser('~/mnt/Fatty Liver/positive_control_CAD/cad_nikpay2015_GCST003116_build37.f.tsv.gz')
ids=[x.strip() for x in open(W+'/instr_clumped_ids.txt') if x.strip()]
sig=pd.read_csv(W+'/exp_sig_pos.tsv',sep='\t').set_index('variant_id')
px=pd.read_csv(W+'/proxy_table.tsv',sep='\t')
raw=pd.read_csv(PC+'/cad_instr_raw.tsv',sep='\t')
base=raw[raw.variant_id.isin(ids)]
miss=[i for i in ids if i not in set(base.variant_id)]
print('direct hits',base.variant_id.nunique(),'missing',miss)
# step A: same position (rsID merges) for missing targets
posq={f"{int(sig.loc[i].chr)}:{int(sig.loc[i].pos)}":i for i in miss}
# step B: proxy candidates
cand=set(px[px.target.isin(miss)].proxy)
open(PC+'/need.txt','w').write('\n'.join(sorted(cand))+'\n')
with open(PC+'/posq.txt','w') as f: f.write('\n'.join(posq)+'\n')

# ---- build the outcome table
ex=pd.read_csv(PC+'/cad_extra_raw.tsv',sep='\t')
def conv(r): return dict(ea=r.effect_allele,oa=r.other_allele,eaf=r.effect_allele_frequency,beta=r.beta,se=r.standard_error,p=r.p_value)
COMP={'A':'T','T':'A','C':'G','G':'C'}
rows=[]
for r in base.itertuples():
    rows.append(dict(rsid=r.variant_id,chr=r.chromosome,pos=r.base_pair_location,**conv(r)))
used=[]
exi=ex.set_index('variant_id')
for t in miss:
    e=sig.loc[t]; ea,oa=e.effect_allele,e.other_allele
    # (A) same chr:pos with identical allele pair
    m=ex[(ex.chromosome==int(e.chr))&(ex.base_pair_location==int(e.pos))]
    done=False
    for q in m.itertuples():
        if {q.effect_allele,q.other_allele}=={ea,oa}:
            rows.append(dict(rsid=t,chr=e.chr,pos=e.pos,**conv(q))); used.append(dict(target=t,via='same position, '+q.variant_id,r2=1.0)); done=True; break
    if done: continue
    # (B) best r2 proxy (EUR reference panel), same logic as s6b
    for c in px[px.target==t].sort_values('r2',ascending=False).itertuples():
        if c.proxy not in exi.index: continue
        q=exi.loc[c.proxy]
        if {c.p_a1,c.p_a2} in ({'A','T'},{'C','G'}) and 0.42<=q.effect_allele_frequency<=0.58: continue
        if {ea,oa}=={c.t_a1,c.t_a2}: tea=ea
        elif len(ea)==1 and {COMP.get(ea),COMP.get(oa)}=={c.t_a1,c.t_a2}: tea=COMP[ea]
        else: continue
        same=(tea==c.t_a1)
        linked=(c.p_a1 if c.sign>0 else c.p_a2) if same else (c.p_a2 if c.sign>0 else c.p_a1)
        qea,qoa=q.effect_allele,q.other_allele
        if qea==linked or COMP.get(qea)==linked: b,f=q.beta,q.effect_allele_frequency
        elif qoa==linked or COMP.get(qoa)==linked: b,f=-q.beta,1-q.effect_allele_frequency
        else: continue
        rows.append(dict(rsid=t,chr=e.chr,pos=e.pos,ea=ea,oa=oa,eaf=f,beta=b,se=q.standard_error,p=q.p_value))
        used.append(dict(target=t,via='proxy '+c.proxy,r2=round(c.r2,3))); done=True; break
    if not done: used.append(dict(target=t,via='NOT FOUND',r2=np.nan))
out=pd.DataFrame(rows); out.to_csv(W+'/out_cad_px.tsv',sep='\t',index=False)
pd.DataFrame(used).to_csv(PC+'/cad_substitutes_used.tsv',sep='\t',index=False)
print(pd.DataFrame(used).to_string(index=False)); print('outcome rows',len(out))
