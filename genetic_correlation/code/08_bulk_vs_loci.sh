#!/bin/bash
# Is the near-zero RC vs CAD rg driven by the large-effect loci or by the polygenic bulk?
# Build SNP sets: (a) "bulk" = all HM3 SNPs more than 1 Mb from any RC genome-wide-significant SNP (P<5e-8);
#                 (b) "loci"  = SNPs within 1 Mb of an RC genome-wide-significant SNP.
# Same SNP set applied to every trait, then LDSC rg (one-step, free intercept) run on each set.
W=/home/claude/rg_work; P=/home/claude/ldsc_env/bin/python; L=/home/claude/ldsc_env/bin/ldsc.py; S=$W/munged; O=$W/ldsc_out; mkdir -p $W/split
$P - <<'PY'
import pandas as pd, numpy as np, glob
W='/home/claude/rg_work'
ld=pd.concat([pd.read_csv(f,sep='\t') for f in sorted(glob.glob(f'{W}/ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.*.l2.ldscore.gz'))])[['CHR','SNP','BP']]
rc=pd.read_csv(f'{W}/munged/RC.sumstats.gz',sep='\t').dropna(subset=['Z']).merge(ld,on='SNP')
gws=rc[rc.Z.abs()>5.45]   # |Z|>5.45 ~ P<5e-8
near=np.zeros(len(ld),bool); ld=ld.reset_index(drop=True)
for c,g in gws.groupby('CHR'):
    pos=g.BP.values; m=(ld.CHR==c).values; b=ld.BP.values[m]
    near[m]=(np.abs(b[:,None]-pos[None,:]).min(axis=1)<=1_000_000)
print('HM3 SNPs',len(ld),'within 1Mb of an RC GWS SNP',int(near.sum()),'bulk',int((~near).sum()),'RC GWS SNPs',len(gws))
sets={'bulk':set(ld.SNP[~near]),'loci':set(ld.SNP[near])}
for t in ['RC','CAD','GLGC_LDL','GLGC_TG','Ghodsian','Fairfield']:
    d=pd.read_csv(f'{W}/munged/{t}.sumstats.gz',sep='\t')
    for k,s in sets.items(): d[d.SNP.isin(s)].to_csv(f'{W}/split/{t}_{k}.sumstats.gz',sep='\t',index=False,compression='gzip')
PY
WL=$W/ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.; M=$(cat $W/ref/M_hm3_noMHC.txt)
for k in bulk loci; do
  X="$W/split/RC_$k.sumstats.gz,$W/split/CAD_$k.sumstats.gz,$W/split/GLGC_LDL_$k.sumstats.gz,$W/split/GLGC_TG_$k.sumstats.gz,$W/split/Ghodsian_$k.sumstats.gz,$W/split/Fairfield_$k.sumstats.gz"
  $P $L --rg $X --ref-ld-chr $WL --w-ld-chr $WL --M $M --out $O/rg_RC_vs_all_$k > /dev/null 2>&1 &
  Y="$W/split/GLGC_LDL_$k.sumstats.gz,$W/split/CAD_$k.sumstats.gz,$W/split/GLGC_TG_$k.sumstats.gz"
  $P $L --rg $Y --ref-ld-chr $WL --w-ld-chr $WL --M $M --out $O/rg_LDL_vs_CAD_TG_$k > /dev/null 2>&1 &
done; wait
