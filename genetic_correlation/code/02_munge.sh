#!/bin/bash
# munge_sumstats.py (LDSC v2.0.x, bulik/ldsc, PyPI ldsc 2.0.1) for each trait. HM3 alleles from Zenodo 7768714 1000G reference.
W=/home/claude/rg_work; P=/home/claude/ldsc_env/bin/python; M=/home/claude/ldsc_env/bin/munge_sumstats.py
mkdir -p $W/munged
common="--snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats BETA,0 --frq FRQ --merge-alleles $W/ref/w_hm3_1000G.snplist --chunksize 500000"
run(){ name=$1; shift; $P $M --sumstats $W/munge_in/$name.tsv.gz $common "$@" --out $W/munged/$name > /dev/null 2>&1; }
run RC        --N 115082 &
run Ghodsian  --N-cas 8434 --N-con 770180 &
run Fairfield --N-cas 4761 --N-con 373227 --info INFO &
run CAD       --N-cas 60801 --N-con 123504 &
wait
