import json, numpy as np, pandas as pd
R='/home/claude/remnant-cholesterol-nafld-mr/'; A=R+'analysis/'
G,F,N='GCST90091033','GCST90054782','finn-b-NAFLD'
NAME={G:'Ghodsian',F:'Fairfield',N:'FinnGen'}
res=pd.read_csv(A+'primary_clumped52/mr_results.tsv',sep='\t')
apoe=pd.read_csv(A+'sens_noAPOE/mr_results.tsv',sep='\t'); o48=pd.read_csv(A+'sens_original48/mr_results.tsv',sep='\t')
S={s['outcome']:s for s in json.load(open(A+'primary_clumped52/summary.json'))}
SA={s['outcome']:s for s in json.load(open(A+'sens_noAPOE/summary.json'))}
S48={s['outcome']:s for s in json.load(open(A+'sens_original48/summary.json'))}
ext=pd.read_csv(R+'extended_analyses/R_out/extended_estimates.csv'); het=pd.read_csv(R+'extended_analyses/R_out/heterogeneity_I2.csv').set_index('dataset')
infl=pd.read_csv(R+'extended_analyses/R_out/influence_diagnostics.csv'); presso=pd.read_csv(R+'results/mrpresso_results.csv').set_index('dataset')
conc=pd.read_csv(R+'results/figures_v2/variant_concordance.csv')
H={o:pd.read_csv(A+f'primary_clumped52/harmonised_{o}.tsv',sep='\t') for o in NAME}
L={o:pd.read_csv(A+f'primary_clumped52/leave_one_out_{o}.tsv',sep='\t') for o in NAME}
pcS=json.load(open(R+'positive_control/python_out/summary.json'))[0]; pcR=pd.read_csv(R+'positive_control/R_out/mr_results_TwoSampleMR.csv')
pcH=pd.read_csv(R+'positive_control/R_out/heterogeneity.csv'); pcE=pd.read_csv(R+'positive_control/R_out/egger_intercept.csv'); pcP=pd.read_csv(R+'positive_control/R_out/mrpresso_main.csv')
M={'ivw':'Inverse variance weighted','egger':'MR Egger','wm':'Weighted median','sm':'Simple mode','wmo':'Weighted mode'}
def g(o,m='ivw',an='all SNPs',df=None):
    d=(res if df is None else df); r=d[(d.outcome==o)&(d.analysis==an)&(d.method==M[m])].iloc[0]; return dict(OR=r.OR,lo=r.OR_lo,hi=r.OR_hi,p=r.pval,n=int(r.nsnp),b=r.b,se=r.se)
def loo(o):
    l=L[o]; l=l[l.SNP_left_out!='All']; OR=np.exp(l.b); return OR.min(),OR.max(),l.pval.min()
def ex(ds,method):
    r=ext[(ext.dataset==ds)&(ext.method==method)].iloc[0]; return dict(OR=r.OR,lo=r.OR_lo,hi=r.OR_hi,p=r.pval,n=int(r.nsnp),note=r.note)
def pr(ds):
    r=presso.loc[ds]; b,sd=r.corrected_beta,r.corrected_sd; return dict(OR=np.exp(b),lo=np.exp(b-1.96*sd),hi=np.exp(b+1.96*sd),p=r.corrected_p,n=int(r.n)-len(r.outliers.split(';')),out=r.outliers)
DS={'Ghodsian':G,'Fairfield':F,'FinnGen':N}
