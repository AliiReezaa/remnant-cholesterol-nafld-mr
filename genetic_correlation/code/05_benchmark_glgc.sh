#!/bin/bash
# External benchmark of the pipeline: GLGC (Willer 2013, GCST002216 TG, GCST002222 LDL) vs CAD (Nikpay 2015) and vs RC / NAFLD.
# Literature expectation (LD Hub / published LDSC): rg(TG,CAD) and rg(LDL,CAD) clearly positive (~0.3 to 0.5).
W=/home/claude/rg_work; P=/home/claude/ldsc_env/bin/python; L=/home/claude/ldsc_env/bin/ldsc.py; MU=/home/claude/ldsc_env/bin/munge_sumstats.py
$P - <<'PY'
import pandas as pd, numpy as np
from scipy.stats import norm
W='/home/claude/rg_work'
for n in ['TG','LDL']:
    d=pd.read_csv(f'{W}/sumstats/GLGC_{n}.txt.gz',sep='\t',usecols=['rsid','A1','A2','beta','se','N','Freq.A1.1000G.EUR'])
    d=d.dropna(); d=d[d.se>0]
    out=pd.DataFrame({'SNP':d.rsid,'A1':d.A1.str.upper(),'A2':d.A2.str.upper(),'BETA':d.beta,'P':2*norm.sf(np.abs(d.beta/d.se)),'N':d.N,'FRQ':d['Freq.A1.1000G.EUR']})
    out=out[out.P>0]; out.to_csv(f'{W}/munge_in/GLGC_{n}.tsv.gz',sep='\t',index=False,float_format='%.17g'); print(n,len(out))
PY
for n in TG LDL; do
  $P $MU --sumstats $W/munge_in/GLGC_$n.tsv.gz --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats BETA,0 --N-col N --frq FRQ --merge-alleles $W/ref/w_hm3_1000G.snplist --chunksize 500000 --out $W/munged/GLGC_$n > /dev/null 2>&1 &
done; wait
LD=$W/ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.; M=$(cat $W/ref/M_hm3_noMHC.txt)
$P $L --rg $W/munged/GLGC_TG.sumstats.gz,$W/munged/CAD.sumstats.gz,$W/munged/GLGC_LDL.sumstats.gz,$W/munged/RC.sumstats.gz,$W/munged/Ghodsian.sumstats.gz,$W/munged/Fairfield.sumstats.gz --ref-ld-chr $LD --w-ld-chr $LD --M $M --out $W/ldsc_out/rg_GLGC_TG_vs_all > /dev/null 2>&1
$P $L --rg $W/munged/GLGC_LDL.sumstats.gz,$W/munged/CAD.sumstats.gz,$W/munged/RC.sumstats.gz,$W/munged/Ghodsian.sumstats.gz,$W/munged/Fairfield.sumstats.gz --ref-ld-chr $LD --w-ld-chr $LD --M $M --out $W/ldsc_out/rg_GLGC_LDL_vs_all > /dev/null 2>&1
