#!/bin/sh
#$ -cwd
#$ -S /bin/bash
#$ -l h_rt=03:00:00
#$ -pe threaded 15

module load gwas/plink/1.9
gcta=/home/software/gcta64/gcta64

##### Command-line arguments that point to (1) which GRM to use, and (2) which phenotype to analyze
grm=$1
pheno=$2

date +"%m/%d/%Y %H:%M:%S"

##### Run GCTA
##### Estimate trait prevalence at 0.01 (--prevalence)
##### Point to sex (binary) covariate (.covar file), and PCs (continuous covariates, .qcovar file)
$gcta --reml --grm /hpc/cog_bioinf/ridder/users/sara/stroke/sign/grms/stroke.all.consensus.AFR.EUR.QC2.$grm \
    --prevalence 0.01 \
    --pheno ./pheno/stroke.all.consensus.AFR.EUR.QC2.$pheno.all.cntrl.gcta.pheno \
    --covar ./stroke.all.consensus.AFR.EUR.QC2.covar \
    --qcovar ./stroke.all.consensus.AFR.EUR.QC2.qcovar \
    --out ./stroke.all.consensus.AFR.EUR.QC2.$grm.$pheno.all.cntrl.0.03 \
    --thread-num 15

date +"%m/%d/%Y %H:%M:%S"

