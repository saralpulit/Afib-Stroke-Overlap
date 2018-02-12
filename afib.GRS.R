#### Read in the data
#### The file contains: family ID, individual ID, sample sex (1=male, 2=female), PCs 1 - 10
#### Several columns containing clinical covariates in 0/1 format (0 = absent, 1 = present)
#### A series of columns containing phenotypes (coded as 1/2, 1=control, 2=case)
#### A column containing the GRS
data <- read.table("geneticRiskScore.stroke.all.consensus.AFR.EUR.QC.all.traits.GRS.data.clean.race.ethnicity.FINAL.txt.gz",header=T)

### Separate dataset (data2) for EUR only; note that this was done to check the signal in the European-ancestry samples only
### Requires a column called 'continent' which indicates ancestry group of the sample
data2 <- subset(data,data$Continent==1)

#### Covariate names correspond to column names in the data frame
#### Here is the base formula (phenotype and covariates)
covariates <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + as.factor(SEX)"
#### Here are the clinical covariates
clinical <- "age + as.factor(hypertension) + as.factor(diabetesMelitus) + as.factor(CAD) + as.factor(smokingFormer) + as.factor(smokingNever)"

#### Run the GRS with just PCs and sex
#### Then run the GRS with the clinical covariates included

#### Create a z-transformed GRS to calculate beta per 1 SD increase in the GRS
grs.z <- (data$SCORESUMfiltered.info.0.8.frq.0.01-mean(data$SCORESUMfiltered.info.0.8.frq.0.01))/sd(data$SCORESUMfiltered.info.0.8.frq.0.01)

data <- cbind(data,grs.z)

# calculate stats using standardized GRS? set to 0 to turn off
mygrs.z = 1

#### Set up an empty matrix to store the results
results.all <- matrix(ncol=9,nrow=28,data=NA)

#### Name the columns
#### b0/se0/p0 refer to beta, se and p-value for a model that only includes PCs and sex
#### b1/se1/p1 refer to beta, se and p-value for a model that includes PCs, sex, and clinical covariates
colnames(results.all) <- c("Pheno","Cases","Controls","b0","se0","p0","b1","se1","p1")

#### set up a counter to tick through all the subtypes
counter <- 1

#### This is my (messy) way of ticking through the columns in the order I wanted to see results.
#### This can also simply be done across the relevant columns containing phenotypes
for(i in c(26,27,40,41,52,53,28,29,42,43,50,51,30,31,54,55,34:39,46,47,58,59)) {
	
	pheno <- colnames(data)[i]
	cases <- table(data[,pheno]==2)["TRUE"][[1]]
	controls <- table(data[,pheno]==1)["TRUE"][[1]]
	
	if( mygrs.z == 0 ) {
		
		### The GRS without the clinical covariates
		f0 <- formula(paste("as.factor(", pheno, ") ~ SCORESUMfiltered.info.0.8.frq.0.01 +", covariates, sep=""))

		### The GRS with the clinical covariates included
		f1 <- formula(paste("as.factor(", pheno, ") ~ SCORESUMfiltered.info.0.8.frq.0.01 +", covariates, " + ", clinical, sep=""))  
	

	} else if( mygrs.z == 1 ) {
		
		### The GRS without the clinical covariates
		f0 <- formula(paste("as.factor(", pheno, ") ~ grs.z +", covariates, sep=""))
	
	
		### The GRS with the clinical covariates included
		f1 <- formula(paste("as.factor(", pheno, ") ~ grs.z +", covariates, " + ", clinical, sep=""))  

	}
	
	### Extract the beta/p-value of SCORESUM
	b0 <- coef(summary(glm(f0, data=data, family="binomial")))[2,1]
	se0 <- coef(summary(glm(f0, data=data, family="binomial")))[2,2]
	p0 <- coef(summary(glm(f0, data=data, family="binomial")))[2,4]

	### Extract the beta/p-value of SCORESUM
	b1 <- coef(summary(glm(f1, data=data, family="binomial")))[2,1]
	se1 <- coef(summary(glm(f1, data=data, family="binomial")))[2,2]
	p1 <- coef(summary(glm(f1, data=data, family="binomial")))[2,4]
	
	### Store the results
	tmp <- c(pheno,cases,controls,b0,se0,p0,b1,se1,p1)
	results.all[counter,]  <- tmp
			
	counter <- counter + 1

}

#################

### Repeat the analysis, but only in 'data2' (i.e., Europeans only)

results.eur <- matrix(ncol=9,nrow=26,data=NA)
colnames(results.eur) <- c("Pheno","Cases","Controls","b0","se0","p0","b1","se1","p1")

counter <- 1

for(i in c(26,27,40,41,52,53,28,29,42,43,50,51,34:39,46,47,58,59)) {
		
	pheno <- colnames(data)[i]
	cases <- table(data[,pheno]==2)["TRUE"][[1]]
	controls <- table(data[,pheno]==1)["TRUE"][[1]]
	
	### The GRS without the clinical covariates
	f0 <- formula(paste("as.factor(", pheno, ") ~ SCORESUMfiltered.info.0.8.frq.0.01 +", covariates, sep=""))
	
	### Extract the beta/p-value of SCORESUM
	b0 <- coef(summary(glm(f0, data=data2, family="binomial")))[2,1]
	se0 <- coef(summary(glm(f0, data=data2, family="binomial")))[2,2]
	p0 <- coef(summary(glm(f0, data=data2, family="binomial")))[2,4]
	
	### The GRS with the clinical covariates included
	f1 <- formula(paste("as.factor(", pheno, ") ~ SCORESUMfiltered.info.0.8.frq.0.01 +", covariates, " + ", clinical, sep=""))  
	
	### Extract the beta/p-value of SCORESUM
	b1 <- coef(summary(glm(f1, data=data2, family="binomial")))[2,1]
	se1 <- coef(summary(glm(f1, data=data2, family="binomial")))[2,2]
	p1 <- coef(summary(glm(f1, data=data2, family="binomial")))[2,4]

	### Store the results
	tmp <- c(pheno,cases,controls,b0,se0,p0,b1,se1,p1)
	results.eur[counter,]  <- tmp
	counter <- counter + 1
	
}

#### Write out the results
write.table(results.all,file="geneticRiskScore.stroke.all.consensus.AFR.EUR.QC.all.traits.GRS.filtered.info.0.8.frq.0.01.results.txt",quote=F,row.names=F,col.names=T)
write.table(results.eur,file="geneticRiskScore.stroke.all.consensus.AFR.EUR.QC.all.traits.GRS.filtered.info.0.8.frq.0.01.results.eur.txt",quote=F,row.names=F,col.names=T)
