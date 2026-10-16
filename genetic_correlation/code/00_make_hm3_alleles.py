"""Build the HapMap3 allele list (SNP, A1, A2) used by munge_sumstats --merge-alleles.
Source: S-LDSC reference files, Zenodo 10.5281/zenodo.7768714 (1000G Phase 3 EUR).
 - regression SNP set  = SNPs in 1000G_Phase3_weights_hm3_no_MHC (HapMap3, MHC excluded)
 - alleles             = A1/A2 from 1000G_Phase3_frq (1000G.EUR.QC.*.frq)
Also writes the SNP count M used for the HM3-only LD score reference."""
import pandas as pd, glob, gzip, os
R='/home/claude/rg_work/ref'
w=pd.concat([pd.read_csv(f,sep='\t',usecols=['CHR','SNP','BP']) for f in sorted(glob.glob(f'{R}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.*.l2.ldscore.gz'))])
fq=pd.concat([pd.read_csv(f,sep=r'\s+') for f in sorted(glob.glob(f'{R}/1000G_Phase3_frq/1000G.EUR.QC.*.frq'))])
print('weights SNPs',len(w),'unique',w.SNP.nunique()); print('frq SNPs',len(fq),'unique',fq.SNP.nunique())
m=w[['SNP']].drop_duplicates().merge(fq[['SNP','A1','A2']].drop_duplicates('SNP'),on='SNP',how='inner')
print('HM3 list with alleles',len(m),'of',w.SNP.nunique())
m.to_csv(f'{R}/w_hm3_1000G.snplist',sep='\t',index=False)
open(f'{R}/M_hm3_noMHC.txt','w').write(str(w.SNP.nunique())+'\n')
