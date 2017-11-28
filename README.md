# Afib-Stroke-Overlap
Code and calls to software that accompanies manuscript on the overlap of atrial fibrillation genetics and cardioembolic stroke. Additional data relevant to this project but that is not stored here, can be downloaded at Zenodo. This includes:

## Data

### (1) Downloadable data

GWAS summary-results of atrial fibrillation (broadly defined) vs. all controls: https://doi.org/10.5281/zenodo.1035871  
GWAS summary-results of atrial fibrillation (strictly defined) vs. all controls: https://doi.org/10.5281/zenodo.1035873

### (2) Data in this repository

#### (1) SiGN.sampleIdentifiers.continental_group.txt.gz

   - Sample identifiers of those individuals included in analyses for the paper. All samples are of either African (AFR) or European (EUR) ancestry, as determined through principal component analysis. This file also includes the acronym identifier for the cohort from which the sampel was drawn (please see the *Supplementary Information* of the manuscript for more detail).
   
#### (2) GeneticRiskScore.SNP.information.txt

   - The SNP identifiers, effect allele for that SNP, and weight given to that SNP (as derived through a previous genome-wide association study) to cacluate the genetic risk scores presented in the paper.

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


