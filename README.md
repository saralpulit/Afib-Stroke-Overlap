# Afib-Stroke-Overlap
Code and calls to software that accompanies manuscript on the overlap of atrial fibrillation genetics and cardioembolic stroke. Additional data relevant to this project but that is not stored here, can be downloaded at Zenodo. This includes:

GWAS summary-results of atrial fibrillation (broadly defined) vs. all controls: https://doi.org/10.5281/zenodo.1035871  
GWAS summary-results of atrial fibrillation (strictly defined) vs. all controls: https://doi.org/10.5281/zenodo.1035873

## In this repository:

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
