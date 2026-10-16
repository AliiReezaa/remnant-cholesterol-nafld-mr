"""Parse LDSC logs into tidy CSV tables (h2 and rg). Reads /home/claude/rg_work/ldsc_out/*.log and writes results/*.csv"""
import re, glob, os, sys, pandas as pd
O='/home/claude/rg_work/ldsc_out'; R=os.path.join(os.path.dirname(__file__),'..','results'); os.makedirs(R,exist_ok=True)
def name(p): return os.path.basename(p).replace('.sumstats.gz','')
def parse_rg(path):
    t=open(path).read()
    call=re.search(r'--rg (\S+)',t).group(1).split(',')
    names=[name(c) for c in call]; rows=[]
    blocks=re.split(r'Computing rg for phenotype \d+/\d+',t)[1:]
    for i,b in enumerate(blocks,start=1):
        g=lambda pat: (re.search(pat,b))
        rg=g(r'Genetic Correlation: (\S+) \((\S+)\)'); z=g(r'Z-score: (\S+)'); p=g(r'\nP: (\S+)')
        gi=g(r'Genetic Covariance\n-+\n.*?\nMean z1\*z2: \S+\nIntercept: (\S+) \((\S+)\)')
        gc=g(r'Genetic Covariance\n-+\nTotal \w+ scale gencov: (\S+) \((\S+)\)')
        h1=re.search(r'Heritability of phenotype 1\n-+\nTotal (\w+) scale h2: (\S+) \((\S+)\)',t)
        h2=g(r'Heritability of phenotype \d+/\d+\n-+\nTotal (\w+) scale h2: (\S+) \((\S+)\)')
        i2=g(r'Heritability of phenotype \d+/\d+\n-+\n.*?Intercept: (\S+) \((\S+)\)')
        rows.append(dict(run=os.path.basename(path)[:-4],trait1=names[0],trait2=names[i],
            rg=float(rg.group(1)) if rg and rg.group(1)!='nan' else float('nan'),se=float(rg.group(2)) if rg and rg.group(2)!='nan' else float('nan'),
            z=float(z.group(1)) if z else float('nan'),p=float(p.group(1)) if p else float('nan'),
            h2_trait2=float(h2.group(2)) if h2 else float('nan'),h2_trait2_se=float(h2.group(3)) if h2 else float('nan'),h2_trait2_scale=h2.group(1) if h2 else '',
            h2_int_trait2=float(i2.group(1)) if i2 else float('nan'),gcov_int=float(gi.group(1)) if gi else float('nan'),gcov_int_se=float(gi.group(2)) if gi else float('nan')))
    return rows
def parse_h2(path):
    t=open(path).read(); m=re.search(r'Total (\w+) scale h2: (\S+) \((\S+)\)',t)
    mc=re.search(r'Mean Chi\^2: (\S+)',t); lg=re.search(r'Lambda GC: (\S+)',t); it=re.search(r'Intercept: (\S+) \((\S+)\)',t)
    return dict(trait=os.path.basename(path)[3:-4],scale=m.group(1),h2=float(m.group(2)),h2_se=float(m.group(3)),z_h2=float(m.group(2))/float(m.group(3)),mean_chi2=float(mc.group(1)),lambda_gc=float(lg.group(1)),intercept=float(it.group(1)),intercept_se=float(it.group(2)))
rg=[r for f in sorted(glob.glob(f'{O}/rg_*.log'))+sorted(glob.glob(f'{O}/sens_*.log')) if 'Genetic Correlation:' in open(f).read() for r in parse_rg(f)]
pd.DataFrame(rg).to_csv(f'{R}/ldsc_rg_all_runs.csv',index=False)
pd.DataFrame([parse_h2(f) for f in sorted(glob.glob(f'{O}/h2_*.log'))]).to_csv(f'{R}/ldsc_h2.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
print(pd.DataFrame(rg)[['run','trait1','trait2','rg','se','z','p','gcov_int','h2_trait2','h2_int_trait2']].round(4).to_string(index=False))
print(); print(pd.read_csv(f'{R}/ldsc_h2.csv').round(4).to_string(index=False))
