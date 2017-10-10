#!/bin/sh
#$ -cwd
#$ -S /bin/bash
#$ -l h_vmem=80G,h_rt=24:00:00
#$ -pe threaded 6

module load plink
gcta=/home/software/gcta64/gcta64

### select the subset of SNPs to calculate the GRM on (genotyped, imputed)
grm=$1

plink --bfile stroke.all.consensus.AFR.EUR.QC2.$grm \
    --maf 0.01 \
    --exclude ./stroke.strata.frequencies.filtered.snps.txt \ ## a list of imputation error SNPs to be dropped
    --make-grm-bin ibc3 \
    --out stroke.all.consensus.AFR.EUR.QC2.$grm \
    --allow-no-sex \
    --memory 80000 \
    --threads 6
