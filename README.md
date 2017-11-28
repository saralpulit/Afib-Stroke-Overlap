# Afib-Stroke-Overlap
Code and calls to software that accompanies manuscript on the overlap of atrial fibrillation genetics and cardioembolic stroke. Additional data relevant to this project but that is not stored here, can be downloaded at Zenodo. This includes:

## Data

### (1) Downloadable data

GWAS summary-level results of atrial fibrillation (broadly defined) vs. all controls: 
   - https://doi.org/10.5281/zenodo.1035871  

GWAS summary-level results of atrial fibrillation (strictly defined) vs. all controls:
   - https://doi.org/10.5281/zenodo.1035873

Merged GWAS summary-level results from the AFGen GWAS, educational attainment GWAS, and stroke subtype and atrial fibrillation GWAS performed in the SiGN data. These data were used to generate the genetic correlation results:
   - https://doi.org/10.5281/zenodo.1067751 

### (2) Data in this repository

#### (1) SiGN.sampleIdentifiers.continental_group.txt.gz

   - Sample identifiers of those individuals included in analyses for the paper. All samples are of either African (AFR) or European (EUR) ancestry, as determined through principal component analysis. This file also includes the acronym identifier for the cohort from which the sampel was drawn (please see the *Supplementary Information* of the manuscript for more detail).
   
#### (2) GeneticRiskScore.SNP.information.txt

   - The SNP identifiers, effect allele for that SNP, and weight given to that SNP (as derived through a previous genome-wide association study) to cacluate the genetic risk scores presented in the paper.
   
#### (3) SuppTable4.partA.SiGN.AFGen.trait.correlations.txt

   - Supplementary Table 4 (part A): correlation between SNP z-scores in the Atrial Fibrillation Genetics (AFGen) Consortium GWAS of atrial fibrillation and SNP z-scores for ischemic stroke subtypes and atrial fibrillation as measured in SiGN.
   - Columns are:
      - the z-score threshold in AFGen used as a cutoff to calculate Pearson's r (correlation)
      - all following columns contain the Pearson's correlation between the AFGen SNPs and the column named in the column (e.g., atrial fibrillation (in SiGN), educational attainment, cardioembolic stroke, etc)
      - more details on the table formats can be found in the supplementary information of the manuscript.
 
#### (4) SuppTable4.partB.SiGN.AFGen.trait.correlations.drop-pitx2-zfhx3.txt

   - Supplementary Table 4 (part B): same as in Part A, but this time excluding 2Mb up- and downstream of the *PITX2* (chromosome 4) and *ZFHX3* (chromosome 16) loci

#### (5) stroke.all.consensus.AFR.EUR.QC2.chr16.zfhx3.snps and stroke.all.consensus.AFR.EUR.QC2.chr4.pitx2.snps

   - A list of SNPs in the PITX2 and ZFHX3 loci excluded from the summary-level GWAS before generating Supplementary Table 4 (part B)
   
## Scripts in this repository:

#### (1) bolt-lmm.gwas.sh

   - A call to BOLT-LMM to run a GWAS in the SiGN stroke data, for various phenotypes
    
#### (2) run.gcta.h2.sh

   - A call to GCTA to calculate heritability
    
#### (3) run.bolt-lmm.h2.sh

   - A call to BOLT-LMM to calculate heritability
    
#### (4) grs.sh

   - A call to PLINK to calculate a genetic risk score (GRS)
    
#### (5) bolt.reml.sh

   - A call to BOLT-LMM to calculate variance explained by the GRS
    
#### (6) h2.obs2liab.R

   - An R script to convert h2 on the observed scale to the liability scale

#### (7) afib.GRS.R

   - An R script for calculating the association between the GRS and various stroke subtypes.


