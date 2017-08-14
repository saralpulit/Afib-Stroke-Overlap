#!/bin/sh
#$ -cwd
#$ -S /bin/bash
#$ -pe threaded 3
#$ -l h_rt=06:00:00

bolt=/home/software/BOLT-LMM_v2.2/bolt
SECONDS=0

echo $(date +%d/%m/%Y\ %H:%M:%S)

grm=/home/stroke/sign/genotyping_merged/stroke.all.consensus.AFR.EUR.QC2

#### GRM subsets (i.e., options for the GRM): genotyped, genotyped.pca, imputed_to_hardcall, imputed_to_hardcall.pca

grm_subset=$1
pheno=$2

#### Point to the grm (--bfile), the phenotype file and phenotype column (--phenoFile, --phenoCol), supporting files, and run REML.
#### Include PCs and sex as covariates

$bolt --bfile=$grm.$grm_subset \
    --phenoFile=./pheno/stroke.all.consensus.AFR.EUR.QC2.$pheno.all.cntrl.boltlmm.pheno \
    --phenoCol=$pheno \
    --geneticMapFile /hpc/local/CentOS7/hers_en/software/BOLT-LMM_v2.2/tables/genetic_map_hg19.txt.gz \
    --LDscoresFile /hpc/local/CentOS7/hers_en/software/BOLT-LMM_v2.2/tables/LDSCORE.1000G_EUR.tab.gz \
    --reml \
    --covarFile ./stroke.all.consensus.AFR.EUR.QC2.covariates.txt \
    --covarCol SEX \
    --qCovarCol PC{1:10} \
    --numThreads 3 \
    --maxModelSnps 1200000

duration=$SECONDS

echo "CPU time $pheno: $(($duration / 60)) min $((duration % 60)) sec"
echo $(date +%d/%m/%Y\ %H:%M:%S)
