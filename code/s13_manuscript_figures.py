"""Manuscript figure set (v2). Reads only existing analysis outputs; no new estimation except plotting-level summaries."""
import json, numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from scipy import stats
R='/home/claude/remnant-cholesterol-nafld-mr/'; A=R+'analysis/'; O=R+'results/figures_v2/'
mpl.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','axes.titlesize':10,'savefig.dpi':300,'figure.dpi':100})
INK='#1f2a37'; GREY='#8a94a3'; LIGHT='#eef1f5'
DS={'GCST90091033':('Ghodsian','#c2410c'),'GCST90054782':('Fairfield','#1d4ed8'),'finn-b-NAFLD':('FinnGen','#047857')}
NAMES={'GCST90091033':'Ghodsian 2021','GCST90054782':'Fairfield 2022','finn-b-NAFLD':'FinnGen R5'}
GENE={'rs58542926':'TM6SF2','rs4665972':'GCKR','rs112875651':'TRIB1','rs2326077':'rs2326077','rs429358':'APOE','rs7412':'APOE','rs72654473':'APOE','rs12151108':'rs12151108'}
res=pd.read_csv(A+'primary_clumped52/mr_results.tsv',sep='\t')
S={s['outcome']:s for s in json.load(open(A+'primary_clumped52/summary.json'))}
ext=pd.read_csv(R+'extended_analyses/R_out/extended_estimates.csv')
het=pd.read_csv(R+'extended_analyses/R_out/heterogeneity_I2.csv')
infl=pd.read_csv(R+'extended_analyses/R_out/influence_diagnostics.csv')
H={o:pd.read_csv(A+f'primary_clumped52/harmonised_{o}.tsv',sep='\t') for o in DS}
L={o:pd.read_csv(A+f'primary_clumped52/leave_one_out_{o}.tsv',sep='\t') for o in DS}
def orient(h):
    h=h.copy(); s=np.sign(h.bx); h['bx']=h.bx*s; h['by']=h.by*s; h['ratio']=h.ratio*s; return h
def save(fig,name): fig.savefig(O+name+'.png',bbox_inches='tight',facecolor='white'); plt.close(fig)
def box(ax,x,y,w,h,text,fc,ec=None,fs=9,bold=False,tc=INK):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.01,rounding_size=0.015',fc=fc,ec=ec or fc,lw=1.2))
    ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=fs,color=tc,fontweight='bold' if bold else 'normal',linespacing=1.35)
def arrow(ax,a,b,c=GREY):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=14,lw=1.5,color=c))

# ---------- Fig 1 design ----------
fig,ax=plt.subplots(figsize=(11,6.2)); ax.axis('off'); ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.text(0.5,0.97,'Study design: two-sample Mendelian randomization of remnant cholesterol and NAFLD',ha='center',fontsize=12,fontweight='bold',color=INK)
box(ax,0.02,0.60,0.20,0.26,'Exposure\nRemnant cholesterol\n(Richardson 2022)\nGCST90092943\nN = 115,082','#fde7d9','#c2410c',9)
box(ax,0.02,0.20,0.20,0.30,'Instrument selection\n52 independent variants\nP < 5e-8, r$^2$ < 0.001\n(10 Mb window)\nmean F = 150','#fff3e8','#c2410c',9)
box(ax,0.40,0.62,0.22,0.22,'Outcome 1\nGhodsian 2021\n8,434 cases / 770,180 controls','#fde7d9','#c2410c',8.5)
box(ax,0.40,0.36,0.22,0.22,'Outcome 2\nFairfield 2022\n4,761 cases / 373,227 controls','#dbe7fb','#1d4ed8',8.5)
box(ax,0.40,0.10,0.22,0.22,'Outcome 3\nFinnGen R5\n894 cases / 217,898 controls','#d5f0e6','#047857',8.5)
box(ax,0.74,0.62,0.24,0.24,'Primary estimator\nIVW (multiplicative\nrandom effects)','#eef1f5',GREY,9,True)
box(ax,0.74,0.33,0.24,0.26,'Sensitivity analyses\nMR-Egger, weighted median,\nmodes, MR-PRESSO, radial,\nSteiger, leave-one-out,\nML, Lasso, conmix, dIVW, RAPS','#eef1f5',GREY,8.5)
box(ax,0.74,0.04,0.24,0.24,'Validation\nPositive control:\nRC $\\rightarrow$ coronary artery disease\n(Nikpay 2015)','#e9defa','#6d28d9',8.5)
arrow(ax,(0.22,0.73),(0.40,0.73)); arrow(ax,(0.22,0.40),(0.40,0.47)); arrow(ax,(0.22,0.33),(0.40,0.21))
arrow(ax,(0.62,0.73),(0.74,0.74)); arrow(ax,(0.62,0.47),(0.74,0.46)); arrow(ax,(0.62,0.21),(0.74,0.16))
arrow(ax,(0.12,0.60),(0.12,0.50))
ax.text(0.01,-0.03,'Assumptions: (1) relevance, (2) independence, (3) exclusion restriction. Instruments harmonised on effect allele; palindromic variants with ambiguous frequency removed.',fontsize=7.5,color=GREY)
save(fig,'fig1_study_design')

# ---------- Fig 2 forest (all methods incl extended) ----------
order=['Inverse variance weighted','MR Egger','Weighted median','Weighted mode','Simple mode']
lab={'Inverse variance weighted':'IVW','MR Egger':'MR-Egger','Weighted median':'Weighted median','Weighted mode':'Weighted mode','Simple mode':'Simple mode'}
emap={'Ghodsian':'GCST90091033','Fairfield':'GCST90054782','FinnGen':'finn-b-NAFLD'}
fig,axs=plt.subplots(1,3,figsize=(13,6.2),sharey=True)
for ax,(o,(n,c)) in zip(axs,DS.items()):
    rows=[]
    for m in order:
        r=res[(res.outcome==o)&(res.analysis=='all SNPs')&(res.method==m)].iloc[0]; rows.append((lab[m],r.OR,r.OR_lo,r.OR_hi,r.nsnp,'p'))
    r=res[(res.outcome==o)&(res.analysis=='MR-PRESSO outliers removed')&(res.method=='Inverse variance weighted')].iloc[0]; rows.append(('IVW, MR-PRESSO corrected',r.OR,r.OR_lo,r.OR_hi,r.nsnp,'s'))
    r=res[(res.outcome==o)&(res.analysis=='Steiger-filtered')&(res.method=='Inverse variance weighted')].iloc[0]; rows.append(('IVW, Steiger filtered',r.OR,r.OR_lo,r.OR_hi,r.nsnp,'s'))
    for _,e in ext[ext.dataset==n].iterrows(): rows.append((e.method,e.OR,e.OR_lo,e.OR_hi,e.nsnp,'e'))
    ys=np.arange(len(rows))[::-1]
    for y,(l,orr,lo,hi,k,t) in zip(ys,rows):
        col={'p':c,'s':INK,'e':GREY}[t]
        ax.plot([lo,hi],[y,y],color=col,lw=2); ax.plot(orr,y,'D' if t=='p' else 'o',color=col,ms=6 if t=='p' else 5,mec='white')
        ax.text(1.03,y,f'{orr:.2f} ({lo:.2f}-{hi:.2f})',va='center',fontsize=7.5,color=INK,transform=ax.get_yaxis_transform())
        if o=='GCST90091033': ax.text(-0.0,y,l,ha='right',va='center',fontsize=8,transform=ax.get_yaxis_transform())
    ax.axvline(1,color=INK,lw=.8,ls='--'); ax.set_xscale('log'); ax.set_xlim(.3,2.0); ax.set_xticks([.4,.6,1,1.5]); ax.set_xticklabels(['0.4','0.6','1.0','1.5']); ax.minorticks_off()
    ax.set_title(f'{NAMES[o]}\n{int(S[o]["ncase"]):,} cases',color=c); ax.set_yticks([]); ax.spines['left'].set_visible(False)
    ax.set_xlabel('Odds ratio per 1 SD higher RC (log scale)'); ax.set_xlim(.3,2.0)
axs[0].set_ylim(-0.7,len(rows)-0.3)
fig.suptitle('Remnant cholesterol and NAFLD: estimates across three outcome datasets and 18 estimator variants',fontweight='bold',fontsize=12,y=1.0)
fig.subplots_adjust(left=.14,wspace=.42); save(fig,'fig2_forest_all_methods')

# ---------- Fig 3 scatter ----------
fig,axs=plt.subplots(1,3,figsize=(14,4.8))
for ax,(o,(n,c)) in zip(axs,DS.items()):
    h=orient(H[o]); out=S[o]['presso_outliers']
    ax.errorbar(h.bx,h.by,yerr=1.96*h.sey,fmt='none',ecolor='#d5d9df',lw=.8,zorder=1)
    ax.errorbar(h.bx,h.by,xerr=1.96*h.sex,fmt='none',ecolor='#d5d9df',lw=.8,zorder=1)
    ax.scatter(h.bx,h.by,s=22,color=c,alpha=.85,zorder=3,edgecolor='white',lw=.4)
    ho=h[h.SNP.isin(out)]; ax.scatter(ho.bx,ho.by,s=60,facecolor='none',edgecolor=INK,lw=1.3,zorder=4)
    x=np.linspace(0,h.bx.max()*1.05,10)
    for m,ls,col in [('Inverse variance weighted','-',c),('MR Egger','--',INK),('Weighted median',':',GREY)]:
        r=res[(res.outcome==o)&(res.analysis=='all SNPs')&(res.method==m)].iloc[0]
        icpt=S[o]['egger_int'] if m=='MR Egger' else 0
        ax.plot(x,icpt+r.b*x,ls=ls,color=col,lw=1.8,label=m.replace('Inverse variance weighted','IVW'))
    offs={0:(5,5),1:(5,-9),2:(-38,4),3:(5,-9)}
    for i,(_,r) in enumerate(ho.iterrows()): ax.annotate(GENE.get(r.SNP,r.SNP),(r.bx,r.by),xytext=offs[i%4],textcoords='offset points',fontsize=7.5,color=INK,fontweight='bold')
    ax.axhline(0,color=GREY,lw=.6); ax.set_title(NAMES[o],color=c)
    ax.set_xlabel('SNP effect on remnant cholesterol (SD)'); ax.set_ylabel('SNP effect on NAFLD (log OR)')
    ax.legend(frameon=False,fontsize=7.5,loc='lower right')
fig.suptitle('Variant-level association scatter (circled: MR-PRESSO outliers)',fontweight='bold',fontsize=12,y=1.02); fig.tight_layout(); save(fig,'fig3_scatter_all')

# ---------- Fig 4 four-panel diagnostics, one figure per dataset ----------
def diag(o,name):
    n,c=DS[o]; h=orient(H[o]).copy(); r_ivw=res[(res.outcome==o)&(res.analysis=='all SNPs')&(res.method=='Inverse variance weighted')].iloc[0]
    fig,axs=plt.subplots(2,2,figsize=(11,9))
    ax=axs[0,0]; ax.scatter(h.bx,h.by,s=24,color=c,edgecolor='white',lw=.4); x=np.linspace(0,h.bx.max()*1.05,10)
    for m,ls,col in [('Inverse variance weighted','-',c),('MR Egger','--',INK),('Weighted median',':',GREY),('Weighted mode','-.','#6d28d9')]:
        r=res[(res.outcome==o)&(res.analysis=='all SNPs')&(res.method==m)].iloc[0]; ax.plot(x,(S[o]['egger_int'] if m=='MR Egger' else 0)+r.b*x,ls=ls,color=col,lw=1.8,label=m.replace('Inverse variance weighted','IVW'))
    ax.axhline(0,color=GREY,lw=.6); ax.legend(frameon=False,fontsize=8); ax.set_xlabel('SNP effect on RC'); ax.set_ylabel('SNP effect on NAFLD'); ax.set_title('A  Scatter plot')
    ax=axs[0,1]; prec=1/h.ratio_se; ax.scatter(h.ratio,prec,s=24,color=c,edgecolor='white',lw=.4)
    ax.axvline(r_ivw.b,color=c,lw=1.8); ax.axvline(res[(res.outcome==o)&(res.analysis=='all SNPs')&(res.method=='MR Egger')].iloc[0].b,color=INK,ls='--',lw=1.5)
    yy=np.linspace(prec.min()*.8,prec.max()*1.02,50); ax.fill_betweenx(yy,r_ivw.b-1.96/yy,r_ivw.b+1.96/yy,color=LIGHT,zorder=0)
    ax.set_xlim(-3,3); ax.set_xlabel('Wald ratio (log OR per SD)'); ax.set_ylabel('Precision (1/SE)'); ax.set_title('B  Funnel plot')
    ax=axs[1,0]; l=L[o].copy(); l=l[l.SNP_left_out!='All'].sort_values('b').reset_index(drop=True)
    ax.errorbar(l.b,np.arange(len(l)),xerr=1.96*l.se,fmt='o',ms=3,color=c,ecolor='#cfd5dc',lw=.8); ax.axvline(r_ivw.b,color=INK,ls='--',lw=1)
    ax.axvline(0,color=GREY,lw=.6); ax.set_yticks([]); ax.set_xlabel('IVW log OR with variant omitted'); ax.set_ylabel(f'{len(l)} leave-one-out runs, sorted'); ax.set_title('C  Leave-one-out')
    ax=axs[1,1]; h2=h.sort_values('ratio').reset_index(drop=True)
    ax.errorbar(h2.ratio,np.arange(len(h2)),xerr=1.96*h2.ratio_se,fmt='o',ms=3,color=c,ecolor='#cfd5dc',lw=.8)
    ax.axvline(r_ivw.b,color=INK,ls='--',lw=1); ax.axvline(0,color=GREY,lw=.6); ax.set_xlim(-4,4); ax.set_yticks([]); ax.set_xlabel('Single-variant Wald ratio (log OR)'); ax.set_title('D  Forest plot of single variants')
    fig.suptitle(f'Diagnostics for remnant cholesterol and NAFLD, {NAMES[o]}',fontweight='bold',fontsize=12); fig.tight_layout(); save(fig,name)
diag('GCST90091033','fig4_diagnostics_Ghodsian'); diag('GCST90054782','figS_diagnostics_Fairfield'); diag('finn-b-NAFLD','figS_diagnostics_FinnGen')

# ---------- Fig 5 variant-level lollipop + concordance heatmap ----------
base=orient(H['GCST90091033'])[['SNP','bx']].rename(columns={'bx':'bx0'})
Z=pd.DataFrame({'SNP':sorted(set().union(*[set(H[o].SNP) for o in DS]))}).merge(base,how='left')
for o in DS:
    h=orient(H[o]); Z=Z.merge(h[['SNP','ratio','ratio_se']].rename(columns={'ratio':'r_'+o,'ratio_se':'s_'+o}),how='left')
Z['bx_any']=Z.SNP.map(pd.concat([orient(H[o]).set_index('SNP').bx for o in DS]).groupby(level=0).first())
Z=Z.sort_values('bx_any',ascending=False)
fig,(a1,a2)=plt.subplots(1,2,figsize=(13,9),sharey=True,gridspec_kw={'width_ratios':[1.1,1.0],'wspace':.04})
y=np.arange(len(Z))[::-1]
a1.barh(y,Z.bx_any,color='#c2410c',alpha=.85,height=.7); a1.set_yticks(y); a1.set_yticklabels([{'APOE':'APOE region'}.get(GENE.get(s,s),GENE.get(s,s)) for s in Z.SNP],fontsize=6); a1.set_xlabel('Effect on RC (SD per allele, risk-increasing orientation)'); a1.set_title('A  Instrument strength by variant')
M=np.array([[ (Z['r_'+o]/Z['s_'+o]).iloc[i] if pd.notna(Z['r_'+o].iloc[i]) else np.nan for o in DS] for i in range(len(Z))])
im=a2.pcolormesh(np.arange(4),np.append(y+.5,y[-1]-.5)[::-1],np.ma.masked_invalid(M)[::-1],cmap='RdBu_r',vmin=-3,vmax=3); a2.set_facecolor('#e5e7eb'); a2.set_xticks([.5,1.5,2.5]); a2.set_xticklabels([NAMES[o] for o in DS]); a2.tick_params(left=False); a2.set_title('B  Single-variant Wald-ratio z score by dataset')
cb=fig.colorbar(im,ax=a2,shrink=.5,pad=.02); cb.set_label('z (red = RC raises NAFLD risk)')
fig.suptitle('Instrument architecture and cross-dataset concordance',fontweight='bold',fontsize=12); fig.tight_layout(); save(fig,'fig5_variants_heatmap')
# concordance stats
conc=[]
for a,b in [('GCST90091033','GCST90054782'),('GCST90091033','finn-b-NAFLD'),('GCST90054782','finn-b-NAFLD')]:
    m=Z[['r_'+a,'r_'+b]].dropna(); rho,p=stats.spearmanr(m.iloc[:,0],m.iloc[:,1]); sg=(np.sign(m.iloc[:,0])==np.sign(m.iloc[:,1])).mean()
    conc.append(dict(a=NAMES[a],b=NAMES[b],n=len(m),spearman=rho,p=p,sign_concordance=sg))
pd.DataFrame(conc).to_csv(O+'variant_concordance.csv',index=False)

# ---------- Fig 6 power ----------
fig,axs=plt.subplots(1,2,figsize=(12,4.6))
ax=axs[0]; orr=np.linspace(1.0,2.2,200)
for o,(n,c) in DS.items():
    se=S[o]['ivw_se']; z=np.log(orr)/se; pw=stats.norm.cdf(z-1.96)+stats.norm.cdf(-z-1.96); ax.plot(orr,pw,color=c,lw=2.2,label=f'{NAMES[o]} (MDE {S[o]["mde_or_80pct_power"]:.2f})')
ax.axhline(.8,color=GREY,ls='--',lw=1); ax.set_xlabel('True odds ratio per 1 SD higher RC'); ax.set_ylabel('Power (alpha = 0.05)'); ax.legend(frameon=False); ax.set_title('A  Power curves from observed IVW standard errors')
ax=axs[1]
for o,(n,c) in DS.items():
    r=res[(res.outcome==o)&(res.analysis=='all SNPs')&(res.method=='Inverse variance weighted')].iloc[0]
    ax.errorbar(r.OR,S[o]['ncase'],xerr=[[r.OR-r.OR_lo],[r.OR_hi-r.OR]],fmt='D',color=c,ms=8,capsize=4,lw=2)
    ax.text(r.OR_hi+.04,S[o]['ncase'],NAMES[o],va='center',fontsize=8)
ax.axvline(1,color=INK,ls='--',lw=.8); ax.set_yscale('log'); ax.set_yticks([1000,2000,5000,8000]); ax.set_yticklabels(['1,000','2,000','5,000','8,000']); ax.minorticks_off(); ax.set_xlabel('IVW odds ratio'); ax.set_ylabel('Number of cases (log)'); ax.set_xlim(.4,1.8); ax.set_title('B  Precision rises with case numbers')
fig.tight_layout(); save(fig,'fig6_power')

# ---------- Fig 7 positive control vs NAFLD ----------
pc=pd.read_csv(R+'positive_control/python_out/mr_results.tsv',sep='\t'); print(pc.columns.tolist(), pc.head(8).to_string())
pcs=json.load(open(R+'positive_control/python_out/summary.json'))

pcR=pd.read_csv(R+'positive_control/R_out/mr_results_TwoSampleMR.csv')
pl={'Inverse variance weighted (multiplicative random effects)':'IVW','MR Egger':'MR-Egger','Weighted median':'Weighted median','Simple mode':'Simple mode','Weighted mode':'Weighted mode'}
fig,ax=plt.subplots(figsize=(9.5,6)); y=0; yt=[];yl=[]
groups=[('Positive control: CAD (Nikpay 2015)','#6d28d9',[(pl[r.method],r.OR,r.OR_lo,r.OR_hi) for _,r in pcR.iterrows()])]
for o,(n,c) in DS.items():
    groups.append((f'NAFLD: {NAMES[o]}',c,[(lab[m],*res[(res.outcome==o)&(res.analysis=='all SNPs')&(res.method==m)].iloc[0][['OR','OR_lo','OR_hi']]) for m in order]))
for g,c,rows in groups:
    ax.text(.22,y,g,fontweight='bold',color=c,fontsize=9.5,va='center'); y-=1
    for l,orr,lo,hi in rows:
        ax.plot([lo,hi],[y,y],color=c,lw=2.2); ax.plot(orr,y,'D',color=c,mec='white',ms=6); ax.text(.22,y,'   '+l,fontsize=8,va='center'); ax.text(3.3,y,f'{orr:.2f} ({lo:.2f}-{hi:.2f})',fontsize=8,va='center'); y-=1
    y-=.6
ax.axvline(1,color=INK,ls='--',lw=.8); ax.set_xscale('log'); ax.set_xlim(.2,2.8); ax.set_xticks([.5,1,1.5,2]); ax.set_xticklabels(['0.5','1.0','1.5','2.0']); ax.minorticks_off(); ax.set_yticks([]); ax.spines['left'].set_visible(False)
ax.set_xlabel('Odds ratio per 1 SD higher remnant cholesterol'); ax.set_title('Same 52 instruments: a clear signal for coronary disease, none for NAFLD',fontsize=11)
save(fig,'fig7_positive_control_contrast')

# ---------- Fig 8 radial ----------
fig,axs=plt.subplots(1,3,figsize=(14,4.8))
for ax,(o,(n,c)) in zip(axs,DS.items()):
    rd=pd.read_csv(R+f'extended_analyses/R_out/radial_data_{n}.csv'); r=res[(res.outcome==o)&(res.analysis=='all SNPs')&(res.method=='Inverse variance weighted')].iloc[0]
    out=rd.Outliers=='Outlier'; ax.scatter(rd.Wj[~out],rd.BetaWj[~out]*rd.Wj[~out],s=22,color=c,edgecolor='white',lw=.4)
    ax.scatter(rd.Wj[out],rd.BetaWj[out]*rd.Wj[out],s=55,color=INK,marker='X',zorder=4)
    x=np.linspace(0,rd.Wj.max()*1.05,10); ax.plot(x,r.b*x,color=c,lw=1.8)
    ax.set_xlabel('Inverse-SE weight $W_j$'); ax.set_ylabel('$W_j\\hat\\beta_j$'); ax.set_title(f'{NAMES[o]}: {int(out.sum())} radial outliers',color=c)
fig.suptitle('Radial MR plots (X: outliers flagged at Bonferroni-adjusted Q$_j$)',fontweight='bold',y=1.02); fig.tight_layout(); save(fig,'fig8_radial')

# ---------- Fig 9 heterogeneity & F ----------
fig,axs=plt.subplots(1,3,figsize=(14,4.4))
ax=axs[0]; allF=pd.concat([H[o].assign(o=o) for o in DS]).drop_duplicates('SNP'); ax.hist(np.log10(allF.F),bins=18,color='#c2410c',alpha=.85,edgecolor='white'); ax.axvline(np.log10(10),color=INK,ls='--'); ax.set_xlabel('log10 F statistic'); ax.set_ylabel('Variants'); ax.set_title('A  Instrument strength (F > 10 for all)')
ax=axs[1]; xs=np.arange(3); w=.35
ax.bar(xs-w/2,het.I2_IVW,w,color=[DS[o][1] for o in DS],label='I$^2$ (IVW heterogeneity)'); ax.bar(xs+w/2,het.I2_GX,w,color=[DS[o][1] for o in DS],alpha=.35,hatch='//',label='I$^2_{GX}$ (NOME)')
ax.set_xticks(xs); ax.set_xticklabels(het.dataset); ax.set_ylim(0,1.3); ax.legend(frameon=False,fontsize=8,loc='upper left',ncol=1); ax.set_title('B  Heterogeneity and NOME')
ax=axs[2]
for o,(n,c) in DS.items():
    d=infl[infl.dataset==n]; ax.scatter(d.student_resid,d.cooks_d,s=22,color=c,alpha=.8,label=n,edgecolor='white',lw=.3)
ax.axhline(1/12,color=INK,ls='--',lw=.8); ax.set_yscale('log'); ax.set_xlabel('Studentized residual'); ax.set_ylabel("Cook's distance"); ax.legend(frameon=False,fontsize=8); ax.set_title("C  Influence diagnostics")
fig.tight_layout(); save(fig,'fig9_strength_heterogeneity_influence')

# ---------- Graphical abstract ----------
fig,ax=plt.subplots(figsize=(11,4.6)); ax.axis('off'); ax.set_xlim(0,1); ax.set_ylim(0,1)
box(ax,.01,.25,.2,.5,'52 genetic instruments\nfor remnant cholesterol\nmean F = 150\nR$^2$ = 6.2%','#fde7d9','#c2410c',10)
arrow(ax,(.21,.5),(.30,.5))
box(ax,.30,.25,.22,.5,'3 NAFLD GWAS\n14,089 cases in total\n(overlapping cohorts)','#eef1f5',GREY,10)
arrow(ax,(.52,.5),(.61,.5))
box(ax,.61,.52,.37,.34,'NAFLD: null\nOR 1.05 / 1.00 / 0.85 (all P > 0.4)\nall estimator variants agree','#d5f0e6','#047857',10,True)
box(ax,.61,.12,.37,.34,'Positive control: CAD\nOR 1.58 (1.38-1.81), P = 5e-11\ninstruments are valid','#e9defa','#6d28d9',10,True)
ax.text(.5,.95,'Genetically higher remnant cholesterol is not associated with NAFLD',ha='center',fontsize=13,fontweight='bold',color=INK)
save(fig,'graphical_abstract')
print('ok')
