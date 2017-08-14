#!/bin/sh
#$ -cwd
#$ -S /bin/bash
#$ -pe threaded 6
#$ -l h_rt=10:00:00

bolt=/home/software/BOLT-LMM_v2.2/bolt
SECONDS=0

echo $(date +%d/%m/%Y\ %H:%M:%S)

#### First, merge together (a) the GRM set of SNPs and (b) the set of SNPs from the GRS, +/- 100kb from these SNPs
window=$1
info=$2

grm=/home/stroke/sign/genotyping_merged/stroke.all.consensus.AFR.EUR.QC2.imputed_to_hardcall
grs=/home/stroke/afib/geneticRiskScore/sensitivityAnalysis/hardcalls/geneticRiskScore.stroke.all.consensus.AFR.EUR.QC.grs$window.maf.0.01.info.$info

#### Run BOLT-REML with two variance components
#### Note that the --modelSnps file has TWO columns.
#### Column 1: rsID, Column 2: a 0/1 indicator if the SNP is in the first variance component or the second
pheno=$3

#### This will result in estimated variance for both components
$bolt --bfile=./geneticRiskScore.stroke.all.consensus.AFR.EUR.QC.grs$window.maf.0.01.info$info.grm.2.merged \
    --phenoFile=./stroke.all.consensus.AFR.EUR.QC2.$pheno.all.cntrl.boltlmm.pheno \
    --phenoCol=$pheno \
    --geneticMapFile /home/software/BOLT-LMM_v2.2/tables/genetic_map_hg19.txt.gz \
    --LDscoresFile /home/software/BOLT-LMM_v2.2/tables/LDSCORE.1000G_EUR.tab.gz \
    --reml \
    --modelSnps=./geneticRiskScore.stroke.all.consensus.AFR.EUR.QC.grs$window.maf.0.01.info.$info.grm.2.merged.modelSnps.txt \
    --numThreads 6 \
    --covarFile ./stroke.all.consensus.AFR.EUR.QC2.covariates.txt \
    --covarCol SEX \
    --qCovarCol PC{1:10} \
    --maxModelSnps 2000000

duration=$SECONDS

echo "CPU time $pheno: $(($duration / 60)) min $((duration % 60)) sec"
echo $(date +%d/%m/%Y\ %H:%M:%S)
