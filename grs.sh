#!/bin/sh
#$ -cwd

#### This is a call to PLINK (v1.90b) to calculate the GRS
#### Note that this can be buggy with dosage files!!
#### It is best to create a .gen file that has *exactly* the SNPs you wish to analyze, and then calculate the GRS
#### ./Allchr_clump.r2_0.5.p_0.0001.Sara.SiGN.txt is a file that contains (a) the rsID, (b) the risk allele for that SNP, and 
#### (c) the desired weight for that SNP in the GRS.

module load plink
lst=$1

plink --dosage ./gen/geneticRiskScore.stroke.all.consensus.AFR.EUR.QC.$lst.gen.gz \
    noheader skip0=1 skip1=1 format=3 \
    --fam geneticRiskScore.stroke.all.consensus.AFR.EUR.QC.CCScCEmajor.stroke.vs.cntrl.all.pheno.fam \
    --map ./gen/geneticRiskScore.stroke.all.consensus.AFR.EUR.QC.$lst.map \
    --score ./Allchr_clump.r2_0.5.p_0.0001.Sara.SiGN.txt header \
    --out ./eneticRiskScore.stroke.all.consensus.AFR.EUR.QC.dosage.$lst

