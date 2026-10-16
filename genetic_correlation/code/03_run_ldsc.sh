#!/bin/bash
# LD score regression: SNP heritability (h2) and genetic correlation (rg).
# Primary reference = HapMap3 regression-SNP LD scores (1000G Phase 3 EUR, MHC excluded; Zenodo 7768714),
# used both as --ref-ld and --w-ld, M = 1,187,349 (analogue of the classic eur_w_ld_chr set).
W=/home/claude/rg_work; P=/home/claude/ldsc_env/bin/python; L=/home/claude/ldsc_env/bin/ldsc.py
LD=$W/ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
M=$(cat $W/ref/M_hm3_noMHC.txt); S=$W/munged; O=$W/ldsc_out; mkdir -p $O
common="--ref-ld-chr $LD --w-ld-chr $LD --M $M"
# sample prevalence = cases/(cases+controls); population prevalence for liability-scale h2 (observed-scale h2 is also reported)
# Ghodsian 8434/778614=0.010832 ; Fairfield 4761/377988=0.012595 ; CAD 60801/184305=0.329888
# NAFLD (clinically recorded) population prevalence set equal to sample prevalence (cohort-based case ascertainment); CAD 0.05
$P $L --h2 $S/RC.sumstats.gz        $common --out $O/h2_RC        > /dev/null 2>&1
$P $L --h2 $S/Ghodsian.sumstats.gz  $common --samp-prev 0.010832 --pop-prev 0.010832 --out $O/h2_Ghodsian  > /dev/null 2>&1
$P $L --h2 $S/Fairfield.sumstats.gz $common --samp-prev 0.012595 --pop-prev 0.012595 --out $O/h2_Fairfield > /dev/null 2>&1
$P $L --h2 $S/CAD.sumstats.gz       $common --samp-prev 0.329888 --pop-prev 0.05     --out $O/h2_CAD       > /dev/null 2>&1
# RC against each outcome (RC first)
$P $L --rg $S/RC.sumstats.gz,$S/Ghodsian.sumstats.gz,$S/Fairfield.sumstats.gz,$S/CAD.sumstats.gz $common \
   --samp-prev nan,0.010832,0.012595,0.329888 --pop-prev nan,0.010832,0.012595,0.05 --out $O/rg_RC_vs_all > /dev/null 2>&1
# between the two NAFLD GWAS (concordance check)
$P $L --rg $S/Ghodsian.sumstats.gz,$S/Fairfield.sumstats.gz $common --out $O/rg_Ghodsian_vs_Fairfield > /dev/null 2>&1
