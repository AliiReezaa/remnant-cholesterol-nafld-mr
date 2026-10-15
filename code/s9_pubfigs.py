import sys, os, json
sys.path.insert(0, os.path.expanduser('~/rc_nafld_mr/code'))
import numpy as np, pandas as pd, mr_lib as M
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.patches import FancyBboxPatch
A=os.path.expanduser('~/rc_nafld_mr/analysis'); D=A+'/primary_clumped52'; F=A+'/figures_pub'; os.makedirs(F,exist_ok=True)
LAB={'GCST90091033':'Ghodsian 2021 (8,434 cases)','GCST90054782':'Fairfield 2022, UK Biobank (4,761 cases)','finn-b-NAFLD':'FinnGen R5 (894 cases)'}
SH={'GCST90091033':'Ghodsian','GCST90054782':'Fairfield','finn-b-NAFLD':'FinnGen'}
BL,OR_,GR='#2b5c8a','#c8553d','#6b6f76'
plt.rcParams.update({'font.size':9,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False})
r=pd.read_csv(D+'/mr_results.tsv',sep='\t')
na=pd.read_csv(A+'/sens_noAPOE/mr_results.tsv',sep='\t'); o48=pd.read_csv(A+'/sens_original48/mr_results.tsv',sep='\t')
def get(df,o,an,me):
    x=df[(df.outcome==o)&(df.analysis==an)&(df.method==me)]
    return x.iloc[0] if len(x) else None
# ---------- Fig 2 forest
rows=[]
for o in LAB:
    rows.append(('H',LAB[o]))
    spec=[(r,'all SNPs','Inverse variance weighted','IVW, all variants (primary)',BL),
          (r,'MR-PRESSO outliers removed','Inverse variance weighted','IVW, MR-PRESSO outliers removed',OR_),
          (r,'Steiger-filtered','Inverse variance weighted','IVW, Steiger-filtered',OR_),
          (na,'all SNPs','Inverse variance weighted','IVW, APOE region excluded',OR_),
          (o48,'all SNPs','Inverse variance weighted','IVW, original 48-variant set',OR_),
          (r,'all SNPs','Weighted median','Weighted median, all variants',GR),
          (r,'all SNPs','MR Egger','MR-Egger, all variants',GR),
          (r,'all SNPs','Weighted mode','Weighted mode, all variants',GR)]
    for df,an,me,lbl,c in spec:
        x=get(df,o,an,me)
        if x is None: x=get(df,o,'all SNPs' if df is not r else an,me)
        if x is None: continue
        rows.append(('R',lbl,x.OR,x.OR_lo,x.OR_hi,x.pval,int(x.nsnp),c,lbl.startswith('IVW, all')))
fig,ax=plt.subplots(figsize=(9.2,8.2)); y=0; yt=[]; yl=[]; ys=[]
for rw in rows:
    if rw[0]=='H':
        y-=0.8; ax.text(0.37,y,rw[1],fontweight='bold',fontsize=9,va='center'); y-=1; continue
    _,lbl,orr,lo,hi,p,n,c,bold=rw
    ax.plot([lo,hi],[y,y],c=c,lw=1.6); ax.plot(orr,y,'s' if bold else 'o',c=c,ms=7 if bold else 5)
    yt.append(y); yl.append(f'{lbl} (n={n})')
    ax.text(3.0,y,f'{orr:.2f} ({lo:.2f} to {hi:.2f})   P={p:.2f}',va='center',fontsize=8)
    y-=1
ax.axvline(1,c='#555',lw=.8); ax.set_xscale('log'); ax.set_xlim(0.35,2.2); ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xticks([0.4,0.6,0.8,1,1.25,1.6,2]); ax.set_xticklabels(['0.4','0.6','0.8','1.0','1.25','1.6','2.0'])
ax.set_yticks(yt); ax.set_yticklabels(yl,fontsize=8); ax.set_ylim(y+0.3,0.2)
ax.set_xlabel('Odds ratio of NAFLD per 1-SD higher genetically predicted remnant cholesterol (log scale)')
ax.spines['left'].set_visible(False); ax.tick_params(axis='y',length=0)
plt.subplots_adjust(left=.32,right=.72,top=.98,bottom=.07); plt.savefig(F+'/fig2_forest_all_analyses.png',dpi=250); plt.close()
# ---------- Fig 3 scatter
HEP={'rs58542926':'TM6SF2','rs4665972':'GCKR region','rs112875651':'TRIB1 region'}
OFF={'rs58542926':(6,-3),'rs4665972':(-8,8),'rs112875651':(6,6)}
fig,axs=plt.subplots(1,3,figsize=(13.5,4.4))
for ax,o in zip(axs,LAB):
    d=pd.read_csv(f'{D}/harmonised_{o}.tsv',sep='\t'); bx,by,sey=d.bx.values,d.by.values,d.sey.values
    sg=np.sign(bx); bxo=np.abs(bx); byo=by*sg
    ax.errorbar(bxo,byo,yerr=1.96*sey,fmt='none',ecolor='#d3d7dd',lw=.7,zorder=1)
    h=d.SNP.isin(HEP).values
    ax.scatter(bxo[~h],byo[~h],s=14,c=BL,zorder=2)
    ax.scatter(bxo[h],byo[h],s=30,c=OR_,zorder=3,label='Labelled variants')
    for s,x,yy in zip(d.SNP[h],bxo[h],byo[h]):
        dx,dy=OFF[s]; ha='right' if dx<0 else 'left'
        ax.annotate(HEP[s],(x,yy),fontsize=7.5,xytext=(dx,dy),textcoords='offset points',ha=ha)
    iv=M.ivw(bx,by,sey); e=M.egger(bx,by,sey); xs=np.linspace(0,bxo.max()*1.05,2)
    ax.plot(xs,iv['b']*xs,c='#1b1f24',lw=1.4,label='IVW'); ax.plot(xs,e['intercept']+e['b']*xs,c='#6b8e23',lw=1.1,ls='--',label='MR-Egger')
    ax.axhline(0,c='#aaa',lw=.5); ax.set_title(LAB[o],fontsize=9.5); ax.set_xlabel('Variant effect on remnant cholesterol (SD)'); ax.set_ylabel('Variant effect on NAFLD (log OR)')
    if o=='GCST90091033': ax.legend(frameon=False,fontsize=7.5,loc='lower right')
plt.tight_layout(); plt.savefig(F+'/fig3_scatter.png',dpi=250); plt.close()
# ---------- Fig 5 what each dataset can exclude
S={s['outcome']:s for s in json.load(open(D+'/summary.json'))}
fig,ax=plt.subplots(figsize=(8,3.4))
for i,o in enumerate(LAB):
    s=S[o]; y=-i
    ax.axvspan(1/s['mde_or_80pct_power'] if False else s['mde_or_80pct_power_protective'],s['mde_or_80pct_power'],ymin=0,ymax=1,alpha=0) 
    ax.plot([s['mde_or_80pct_power_protective'],s['mde_or_80pct_power']],[y-.28,y-.28],c=GR,lw=5,alpha=.35,solid_capstyle='butt')
    ax.plot([s['ivw_or_lo'],s['ivw_or_hi']],[y,y],c=BL,lw=2); ax.plot(s['ivw_or'],y,'s',c=BL,ms=7)
    ax.text(2.2,y-.28,f"80% power detects: {s['mde_or_80pct_power_protective']:.2f} to {s['mde_or_80pct_power']:.2f}",va='center',fontsize=7,color='#555')
    ax.text(2.2,y,f"IVW {s['ivw_or']:.2f} ({s['ivw_or_lo']:.2f} to {s['ivw_or_hi']:.2f})",va='center',fontsize=8)
ax.axvline(1,c='#555',lw=.8); ax.set_xscale('log'); ax.set_xlim(.45,2.1); ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xticks([.5,.7,1,1.4,2]); ax.set_xticklabels(['0.5','0.7','1.0','1.4','2.0'])
ax.set_yticks([0,-1,-2]); ax.set_yticklabels([SH[o]+f" ({S[o]['ncase']:,} cases)" for o in LAB]); ax.spines['left'].set_visible(False); ax.tick_params(axis='y',length=0)
ax.set_xlabel('Odds ratio per 1-SD (log scale). Blue: IVW 95% CI. Grey bar: range of effects detectable with 80% power')
plt.subplots_adjust(left=.2,right=.72,bottom=.22,top=.97); plt.savefig(F+'/fig5_power_precision.png',dpi=250); plt.close()
# ---------- Fig 1 instrument flow
fig,ax=plt.subplots(figsize=(8.6,4.6)); ax.axis('off'); ax.set_xlim(0,10); ax.set_ylim(0,6)
def box(x,y,w,h,t,c='#eef2f7'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.04',fc=c,ec='#5b6b7d',lw=1)); ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=8)
def arr(x1,y1,x2,y2): ax.annotate('',(x2,y2),(x1,y1),arrowprops=dict(arrowstyle='->',color='#444'))
box(1.3,5.0,7.4,.8,'Remnant cholesterol GWAS (UK Biobank NMR, n up to 115,082)\n8,873 variants with P < 5e-8 (after removing duplicate identifiers)')
box(1.3,3.8,7.4,.7,'7,593 variants present in the 1000 Genomes EUR reference panel\n(1,280 absent from the panel)')
box(1.3,2.6,7.4,.7,'LD clumping (r2 < 0.001, 10,000 kb): 52 independent variants\n(proxies r2 >= 0.8 used where a variant was missing in an outcome)')
arr(5,5.0,5,4.55); arr(5,3.8,5,3.35)
for x,t in [(.2,'Ghodsian\n48 variants used'),(3.6,'Fairfield\n51 variants used'),(7.0,'FinnGen\n48 variants used')]:
    box(x,.9,2.8,.8,t); arr(5,2.6,x+1.4,1.75)
ax.text(5,.35,'Sensitivity sets: MR-PRESSO outliers removed, Steiger-filtered, APOE region excluded, original 48-variant set',ha='center',fontsize=8,style='italic')
plt.savefig(F+'/fig1_instrument_selection.png',dpi=250,bbox_inches='tight'); plt.close()
print('done')
