import sys, os; sys.path.insert(0, os.path.expanduser('~/rc_nafld_mr/code'))
import numpy as np, pandas as pd, json, mr_lib as M
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
D=os.path.expanduser('~/rc_nafld_mr/analysis/primary_clumped52'); F=D+'/figures'; os.makedirs(F,exist_ok=True)
LAB={'GCST90091033':'Ghodsian 2021 (8,434 cases)','GCST90054782':'Fairfield 2022 UK Biobank (4,761 cases)','finn-b-NAFLD':'FinnGen R5 (894 cases)'}
HEP={'rs58542926':'TM6SF2','rs4665972':'GCKR region','rs112875651':'TRIB1 region'}
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
summ={s['outcome']:s for s in json.load(open(D+'/summary.json'))}
fig,axs=plt.subplots(1,3,figsize=(13,4.2))
for ax,lab in zip(axs,LAB):
    d=pd.read_csv(f'{D}/harmonised_{lab}.tsv',sep='\t'); bx,by,sex,sey=d.bx.values,d.by.values,d.sex.values,d.sey.values
    sg=np.sign(bx); bxo=np.abs(bx); byo=by*sg
    ax.errorbar(bxo,byo,yerr=1.96*sey,fmt='none',ecolor='#c9ced6',lw=.7,zorder=1)
    hep=d.SNP.isin(HEP); ax.scatter(bxo[~hep],byo[~hep],s=14,c='#2b5c8a',zorder=2)
    ax.scatter(bxo[hep],byo[hep],s=26,c='#c8553d',zorder=3,label='hepatic-fat loci')
    for s,x,y in zip(d.SNP[hep],bxo[hep],byo[hep]): ax.annotate(HEP[s],(x,y),fontsize=7,xytext=(3,3),textcoords='offset points')
    iv=M.ivw(bx,by,sey); e=M.egger(bx,by,sey); xs=np.linspace(0,bxo.max()*1.05,2)
    ax.plot(xs,iv['b']*xs,c='#1b1f24',lw=1.4,label='IVW'); ax.plot(xs,e['intercept']+e['b']*xs,c='#6b8e23',lw=1,ls='--',label='MR-Egger')
    ax.axhline(0,c='#999',lw=.5); ax.set_title(LAB[lab],fontsize=9); ax.set_xlabel('SNP effect on remnant cholesterol (SD)'); ax.set_ylabel('SNP effect on NAFLD (log OR)')
    if lab=='GCST90091033': ax.legend(frameon=False,fontsize=7)
plt.tight_layout(); plt.savefig(F+'/fig_scatter.png',dpi=220); plt.close()
fig,axs=plt.subplots(1,3,figsize=(13,4))
for ax,lab in zip(axs,LAB):
    d=pd.read_csv(f'{D}/harmonised_{lab}.tsv',sep='\t'); iv=M.ivw(d.bx.values,d.by.values,d.sey.values)
    ratio=d.by/d.bx; prec=d.bx.abs()/d.sey
    hep=d.SNP.isin(HEP); ax.scatter(ratio[~hep],prec[~hep],s=14,c='#2b5c8a'); ax.scatter(ratio[hep],prec[hep],s=26,c='#c8553d')
    ax.axvline(iv['b'],c='#1b1f24',lw=1); ax.set_xlabel('Single-SNP Wald ratio (log OR per SD)'); ax.set_ylabel('Precision (|bx|/se)'); ax.set_title(LAB[lab],fontsize=9)
plt.tight_layout(); plt.savefig(F+'/fig_funnel.png',dpi=220); plt.close()
fig,axs=plt.subplots(1,3,figsize=(13,7.2))
for ax,lab in zip(axs,LAB):
    l=pd.read_csv(f'{D}/leave_one_out_{lab}.tsv',sep='\t'); d=pd.read_csv(f'{D}/harmonised_{lab}.tsv',sep='\t'); iv=M.ivw(d.bx.values,d.by.values,d.sey.values)
    l=l.sort_values('b').reset_index(drop=True); y=np.arange(len(l))
    ax.errorbar(l.b,y,xerr=1.96*l.se,fmt='o',ms=3,c='#2b5c8a',ecolor='#9db4cc',lw=.8)
    ax.set_yticks(y); ax.set_yticklabels(l.SNP_left_out,fontsize=5.5); ax.axvline(0,c='#999',lw=.6); ax.axvline(iv['b'],c='#c8553d',lw=.8,ls='--')
    ax.set_xlabel('IVW log OR, leaving SNP out'); ax.set_title(LAB[lab],fontsize=9)
plt.tight_layout(); plt.savefig(F+'/fig_leave_one_out.png',dpi=220); plt.close()
r=pd.read_csv(D+'/mr_results.tsv',sep='\t'); r=r[r.analysis.isin(['all SNPs','MR-PRESSO outliers removed'])]
fig,ax=plt.subplots(figsize=(8,6)); y=0; yt=[]; yl=[]
for lab in LAB:
    for an,col in [('all SNPs','#2b5c8a'),('MR-PRESSO outliers removed','#c8553d')]:
        for m in ['Inverse variance weighted','Weighted median','MR Egger','Weighted mode']:
            q=r[(r.outcome==lab)&(r.analysis==an)&(r.method==m)].iloc[0]
            ax.errorbar(q.OR,y,xerr=[[q.OR-q.OR_lo],[q.OR_hi-q.OR]],fmt='o',ms=4,c=col,lw=1)
            yt.append(y); yl.append(f'{m.replace("Inverse variance weighted","IVW")} ({int(q.nsnp)}){"" if an=="all SNPs" else " *"}'); y+=1
        y+=0.3
    ax.text(ax.get_xlim()[1] if False else 0.02,y-0.3,'',fontsize=1); y+=0.6
ax.axvline(1,c='#999',lw=.7); ax.set_xscale('log'); ax.set_yticks(yt); ax.set_yticklabels(yl,fontsize=6.5); ax.invert_yaxis()
ax.set_xlabel('Odds ratio of NAFLD per 1-SD higher remnant cholesterol (log scale)'); ax.set_title('Blue: all instruments. Orange: after MR-PRESSO outlier removal (*)\nOrder: Ghodsian, Fairfield, FinnGen',fontsize=8)
plt.tight_layout(); plt.savefig(F+'/fig_forest_methods.png',dpi=220); plt.close()
# summary table
rows=[]
for lab,s in summ.items():
    rows.append({'Outcome':LAB[lab],'Cases':s['ncase'],'Controls':s['ncontrol'],'SNPs used':s['n_used'],'Mean F':round(s['mean_F']),'R2 exposure':round(s['R2'],3),
      'IVW OR (95% CI)':'%.2f (%.2f-%.2f)'%(s['ivw_or'],s['ivw_or_lo'],s['ivw_or_hi']),'IVW p':round(s['ivw_p'],3),'Cochran Q p':'%.1e'%s['Q_p'],'Egger intercept p':round(s['egger_int_p'],3),
      'PRESSO global p':'<0.001' if s['presso_global_p']<=0.001 else round(s['presso_global_p'],3),'PRESSO outliers':', '.join(s['presso_outliers']),
      'IVW OR after PRESSO':round(s['presso_corr_or'],2),'p after PRESSO':round(s['presso_corr_p'],3),'Detectable OR (80% power)':round(s['mde_or_80pct_power'],2)})
pd.DataFrame(rows).to_csv(D+'/table_main_summary.csv',index=False); print(pd.DataFrame(rows).T.to_string())
