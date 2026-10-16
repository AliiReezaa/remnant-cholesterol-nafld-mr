#!/bin/bash
# Sensitivity of the RC genetic correlations to LDSC settings. RC is always the first trait.
W=/home/claude/rg_work; P=/home/claude/ldsc_env/bin/python; L=/home/claude/ldsc_env/bin/ldsc.py; S=$W/munged; O=$W/ldsc_out; mkdir -p $O
WL=$W/ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.; FL=$W/ref/LDscore/LDscore.; M=$(cat $W/ref/M_hm3_noMHC.txt)
T="$S/RC.sumstats.gz,$S/CAD.sumstats.gz,$S/Ghodsian.sumstats.gz,$S/Fairfield.sumstats.gz"
SP="--samp-prev nan,0.329888,0.010832,0.012595 --pop-prev nan,0.05,0.010832,0.012595"
# A: no chi2 truncation of large-effect loci
$P $L --rg $T --ref-ld-chr $WL --w-ld-chr $WL --M $M --chisq-max 100000 $SP --out $O/sens_A_chisqmax_off > /dev/null 2>&1 &
# B: full 1000G (all-SNP) LD scores as regressor, HM3 weights for regression weighting (M_5_50 read from files)
$P $L --rg $T --ref-ld-chr $FL --w-ld-chr $WL $SP --out $O/sens_B_full1000G_ref > /dev/null 2>&1 &
# C: two-step estimator for the intercepts (cutoff 30)
$P $L --rg $T --ref-ld-chr $WL --w-ld-chr $WL --M $M --two-step 30 $SP --out $O/sens_C_twostep > /dev/null 2>&1 &
# D: genetic-covariance intercept constrained to 0 (valid for RC vs CAD, no participant overlap expected)
$P $L --rg $T --ref-ld-chr $WL --w-ld-chr $WL --M $M --intercept-gencov 0,0,0 $SP --out $O/sens_D_gcov_int0 > /dev/null 2>&1 &
wait
