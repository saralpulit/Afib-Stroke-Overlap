#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -pe threaded 2
#$ -l h_rt=12:00:00,h_vmem=10G

chr=$1

gtool=/home/software/gtool/gtool
dosages=/home/stroke/SiGNimpute2/chr$chr

module load gwas/plink/1.9

##################
#### BOLT-LMM ####
##################

## Select a GRM to use
## either genotyped or imputed_to_hardcall
grm=$2
pheno=$3

bolt=/home/software/BOLT-LMM_v2.2/bolt
homedir=/home/stroke/UnionIntersect/boltlmm

SECONDS=0

echo $(date +%d/%m/%Y\ %H:%M:%S)

## Run the GWAS in BOLT-LMM
## --bfile points to the 'GRM' SNPs, --impute2* commands point to the imputed dosages
$bolt --bfile $homedir/$grm/stroke.all.consensus.AFR.EUR.QC2.$grm.pca \
    --phenoFile /home/stroke/afib/stroke.all.consensus.AFR.EUR.QC2.phenotypes.subtypes.txt \
    --phenoCol $pheno \
    --geneticMapFile /home/software/BOLT-LMM_v2.2/tables/genetic_map_hg19.txt.gz \
    --LDscoresFile /home/software/BOLT-LMM_v2.2/tables/LDSCORE.1000G_EUR.tab.gz \
    --numThreads 3 \
    --lmm \
    --statsFileImpute2Snps ./stroke.all.consensus.AFR.EUR.QC2.chr$chr.imputed.boltlmm.afib.$pheno.nonGRM.$grm.results.pcs.txt \
    --statsFile  ./stroke.all.consensus.AFR.EUR.QC2.chr$chr.imputed.boltlmm.afib.$pheno.GRM.$grm.results.pcs.txt\
    --impute2FileList /home/stroke/SiGNimpute2/dosage.lists/stroke.all.consensus.AFR.EUR.QC2.chr$chr.imputed.dosage.list \
    --impute2FidIidFile /home/stroke/SiGNimpute2/stroke.all.consensus.AFR.EUR.QC2.imputed.samples \
    --impute2MinMAF 0.01 \
    --maxMissingPerSnp 0.05 \
    --qCovarCol PC{1:10} \
    --covarFile /home/stroke/afib/stroke.all.consensus.AFR.EUR.QC2.covariates.txt \
    --covarCol SEX \

duration=$SECONDS

echo "CPU time $pheno: $(($duration / 60)) min $((duration % 60)) sec"
echo $(date +%d/%m/%Y\ %H:%M:%S)
