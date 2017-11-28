data <- read.table("SiGN.AFGen.correlations.gwas.zscores.merged.txt.gz",header=T)

## Data header (data is downloadable; please see README for this GitHub repository)
## SNP AFGen.EUR.Z EduYrs.Z afib.broad.Z allstroke.Z CCScCEmajor.Z CCScCRYPTCE.Z CCScINCUNC.Z ... [19 columns total]

## Calculate the correlation matrix on the full set of data (all phenotypes vs all phenotypes)
## Only calculate correlations across pairs where both z-score values are present
cor.matrix <- cor(data[,2:19],use="complete.obs")

## Subset the correlation matrix down to just the CCS Causative phenotypes + AFGen + Afib (SiGN) + Education Years (as a null)
## Do a bit of reordering so the columns are in the order we want

plot.fig.matrix <- cor.matrix[c(4,5,8,9,10,6,7,3,1,2),c(4,5,8,9,10,6,7,3,1,2)]

################
### Figure 3 ###
################

## Plot out this correlation matrix
## We place these two triangles together in powerpoint and change the coloring of the printed correlation values in Illustrator to generate figure 3
library(corrplot)
mycol <- colorRampPalette(c("dodgerblue4","white","sienna2")) ## color shading

## This plots out a triangle of color-shaded correlations
corrplot(as.matrix(plot.fig.matrix),col=mycol(200),tl.cex=0.5,number.cex=0.5,type="upper")

## This plots out a triangle with the actual correlation (r) values printed
corrplot(as.matrix(plot.fig.matrix),col=mycol(200),tl.cex=0.5,number.cex=0.5,type="lower",method="number")

#####################
### AFGen vs SiGN ###
#####################

## Note that the correlation between AFGen and SiGN Afib is weak (0.07). Let's investigate
## Using the AFGen z-scores as a "reference," let's see how the correlation changes if we filter SNPs based on the z-score in AFGen

result.afgen <- matrix(ncol=17,nrow=21,data=0); colnames(result.afgen) <- colnames(data)[3:19]

for (pheno in 3:19) {
	for(i in seq(0,5,0.25)) {
	
	## Subset the data based on the z-scores in AFGen; keep AFGen and the pheno we are interested in
	tmp.data <- subset(data,data$AFGen.EUR.Z >= i | data$AFGen.EUR.Z <= -i)[,c(2,pheno)]
	
	## for each z-score threshold, calculate the correlation
	result.afgen[i/0.25+1,pheno-2]<- cor(tmp.data[,1],tmp.data[,2],use="complete.obs")
	
	}
}

result.afgen <- cbind(seq(0,5,0.25),result.afgen); colnames(result.afgen)[1] <- "Z.threshold"
write.table(result.afgen,file="SiGN.AFGen.correlations.AFGen_as_baseline.txt",quote=F,row.names=F,col.names=T)

## Recalculate the correlations after dropping +/- 2Mb from the PITX2 and ZFHX3 loci (to make sure those two loci aren't solely driving the correlation)

pitx2.snps <- read.table("stroke.all.consensus.AFR.EUR.QC2.chr4.pitx2.snps",header=F)
zfhx3.snps <- read.table("stroke.all.consensus.AFR.EUR.QC2.chr16.zfhx3.snps",header=F)
drop.snps <- rbind(pitx2.snps,zfhx3.snps); colnames(drop.snps) <- "SNP"

result.afgen.drop <- matrix(ncol=17,nrow=21,data=0); colnames(result.afgen.drop) <- colnames(data)[3:19]

for (pheno in 3:19) {
	for(i in seq(0,5,0.25)) {
	
	## Subset the data based on the z-scores in AFGen; keep AFGen and the pheno we are interested in
	tmp.data <- subset(data,data$AFGen.EUR.Z >= i | data$AFGen.EUR.Z <= -i)[,c(1,2,pheno)]
	
	## Subset tmp.data again, this time dropping the SNPs in PITX2 and ZFHX3 loci
	tmp.data <- tmp.data[! tmp.data$SNP %in% drop.snps$SNP, ]
	
	## for each z-score threshold, calculate the correlation
	result.afgen.drop[i/0.25+1,pheno-2]<- cor(tmp.data[,2],tmp.data[,3],use="complete.obs")
	
	}
}

result.afgen.drop <- cbind(seq(0,5,0.25),result.afgen.drop); colnames(result.afgen.drop)[1] <- "Z.threshold"
write.table(result.afgen.drop,file="SiGN.AFGen.correlations.AFGen_as_baseline.drop-pitx2-zfhx3.txt",quote=F,row.names=F,col.names=T)


####################
### Plot results ###
####################

### Multi-panel figure; correlation (pearson r) values on top, and scatterplots of underlying data on the bottom
par(mfrow=c(2,3),mar=c(4,4,1,1))

result.afgen <- read.table("SiGN.AFGen.correlations.AFGen_as_baseline.txt",header=T)

## alternatively, plot the results without PITX2 and ZFHX3
## result.afgen <- read.table("SiGN.AFGen.correlations.AFGen_as_baseline.drop-pitx2-zfhx3.txt",header=T)

##
## (1) Set up the plot and axes
##

plot(1,col="white",xlim=c(0,5),ylim=c(-0.4,1),las=1,axes=F,xlab="",ylab="")
lines(x=c(0,5),y=c(0,0),lty="dashed",lwd=0.5,col="grey30") ## line at r = 0
axis(side=1,at=seq(0,5,0.5),las=1,labels=T,cex.axis=0.9); mtext(side=1, line=2.5, "Z-score threshold in AFGen", cex=0.6)
axis(side=2,at=seq(-0.4,1,0.2),las=1,labels=T,cex.axis=0.9); mtext(side=2, line=2.5, "Correlation (r) with AFGen", cex=0.6)

## Plot out the "null" (Years of education) first
lines(result.afgen[,1],result.afgen[,2],lwd=5,col="grey90")
points(result.afgen[,1],result.afgen[,2],pch=21,cex=0.5,bg="grey90",lwd=0.5)

## Plot out SiGN afib correlation
lines(result.afgen[,1],result.afgen[,3],lwd=5,col="lightgoldenrod1")
points(result.afgen[,1],result.afgen[,3],pch=21,cex=0.5,bg="lightgoldenrod1",lwd=0.5) ## afib.broad in SiGn
box()

## Add a legend
legend("topleft",legend=c("Atrial fibrillation (SiGN)", "Educational attainment"),lwd=3,col=c("lightgoldenrod1","grey90"),cex=0.7,bty="n")

##
## (2) Next plot -- the primary stroke subtypes
##

## Set up the plot and axes
plot(1,col="white",xlim=c(0,5),ylim=c(-0.4,1),las=1,axes=F,xlab="",ylab="")
lines(x=c(0,5),y=c(0,0),lty="dashed",lwd=0.5,col="grey30") ## line at r = 0
axis(side=1,at=seq(0,5,0.5),las=1,labels=T,cex.axis=0.9); mtext(side=1, line=2.5, "Z-score threshold in AFGen", cex=0.6)
axis(side=2,at=seq(-0.4,1,0.2),las=1,labels=T,cex.axis=0.9); mtext(side=2, line=2.5, "Correlation (r) with AFGen", cex=0.6)

## And add the three primary stroke subtypes
lines(result.afgen[,1],result.afgen[,9],lwd=5,col="thistle2") ## CCScSAO
points(result.afgen[,1],result.afgen[,9],pch=21,cex=0.5,bg="thistle2",lwd=0.5) ## CCScSAO

lines(result.afgen[,1],result.afgen[,8],lwd=5,col="slateblue3") ## CCScLAA
points(result.afgen[,1],result.afgen[,8],pch=21,cex=0.5,bg="slateblue3",lwd=0.5) ## CCScLAA

lines(result.afgen[,1],result.afgen[,5],lwd=5,col="darkseagreen3") ## CCScCE
points(result.afgen[,1],result.afgen[,5],pch=21,cex=0.5,bg="darkseagreen3",lwd=0.5) ## CCScCE
box()

## Add a legend
legend("topleft",legend=c("Cardioembolic stroke", "Large artery atherosclerosis","Small artery occlusion"),lwd=3,col=c("darkseagreen3","slateblue3","thistle2"),cex=0.7,bty="n")


##
## (3) Last plot -- undetermined subtypes
##

## Set up the plot and axes
plot(1,col="white",xlim=c(0,5),ylim=c(-0.4,1),las=1,axes=F,xlab="",ylab="")
lines(x=c(0,5),y=c(0,0),lty="dashed",lwd=0.5,col="grey30") ## line at r = 0
axis(side=1,at=seq(0,5,0.5),las=1,labels=T,cex.axis=0.9); mtext(side=1, line=2.5, "Z-score threshold in AFGen", cex=0.6)
axis(side=2,at=seq(-0.4,1,0.2),las=1,labels=T,cex.axis=0.9); mtext(side=2, line=2.5, "Correlation (r) with AFGen", cex=0.6)

## undetermined subtypes
lines(result.afgen[,1],result.afgen[,7],lwd=5,col="cornflowerblue") ## CCScINCUNC
points(result.afgen[,1],result.afgen[,7],pch=21,cex=0.5,bg="cornflowerblue",lwd=0.5) ## CCScINCUNC

lines(result.afgen[,1],result.afgen[,6],lwd=5,col="royalblue4") ## CCScCRYPTCE
points(result.afgen[,1],result.afgen[,6],pch=21,cex=0.5,bg="royalblue4",lwd=0.5) ## CCScCRYPTCE

lines(result.afgen[,1],result.afgen[,10],lwd=5,col="lightsteelblue1") ## CCScUNDETER
points(result.afgen[,1],result.afgen[,10],pch=21,cex=0.5,bg="lightsteelblue1",lwd=0.5) ## CCScUNDETER
box()

## Add a legend
legend("topleft",legend=c("All undetermined", "Incomplete/unclassified","Cryptogenic/cardioembolic minor"),lwd=3,col=c("lightsteelblue1","cornflowerblue","royalblue4"),cex=0.7,bty="n")

#########################
## Visualize the underlying data for the correlation calculation

## Take a random subset of the null SNPs (otherwise there are way too many points to plot)
random.subset <- subset(data,data$AFGen.EUR.Z < 2 & data$AFGen.EUR.Z > -2)[seq(1,7346541,20),] ## This takes a random 5% of the "null" data (taking 10% total doesn't make a visual difference)
plot.data <- subset(data,data$AFGen.EUR.Z > 2 | data$AFGen.EUR.Z < -2) ## now get the "non-null" data and combine
plot.data <- rbind(plot.data,random.subset)

## Or, take a random subset of the data, also excluding PITX2 and ZFHX3 (uncomment this if you want to drop these two loci)
# plot.data <- plot.data[! plot.data$SNP %in% drop.snps$SNP, ]

## First, plot AFGen vs afib in SiGN
plot(plot.data$AFGen.EUR.Z, plot.data$EduYrs.Z,pch=16,col="grey90",xlab="",ylab="",main="",cex=1.1,las=1,xlim=c(-20,20),ylim=c(-10,10))
points(plot.data$AFGen.EUR.Z,plot.data$afib.broad.Z,pch=16,col="lightgoldenrod1",cex=1.1,las=1)
mtext(side=1,line=2,text="Z-score in AFGen", cex=0.6); mtext(side=2,line=2,text="Z-score in 2nd GWAS", cex=0.6)

## Add a legend
legend("topleft",legend=c("Atrial fibrillation (SiGN)", "Educational attainment"),lwd=3,col=c("lightgoldenrod1","grey90"),cex=0.7,bty="n")

## Next plot the correlations for the primary subtypes
plot(plot.data$AFGen.EUR.Z,plot.data$CCScSAO.Z,pch=16,col="thistle2",cex=1.1,las=1,xlab="",ylab="",main="",xlim=c(-20,20),ylim=c(-10,10))
points(plot.data$AFGen.EUR.Z, plot.data$CCScLAA.Z,pch=16,col="slateblue3",cex=1.1,las=1)
points(plot.data$AFGen.EUR.Z,plot.data$CCScCEmajor.Z,pch=16,col="darkseagreen3",cex=1.1,las=1)
mtext(side=1,line=2,text="Z-score in AFGen", cex=0.6); mtext(side=2,line=2,text="Z-score in 2nd GWAS", cex=0.6)

## Add a legend
legend("topleft",legend=c("Cardioembolic stroke", "Large artery atherosclerosis","Small artery occlusion"),lwd=3,col=c("darkseagreen3","slateblue3","thistle2"),cex=0.7,bty="n")

## Last, plot the undetermined categories
plot(plot.data$AFGen.EUR.Z, plot.data$CCScINCUNC.Z,pch=16,col="cornflowerblue",xlab="",ylab="",main="",cex=1.1,las=1,xlim=c(-20,20),ylim=c(-10,10))
points(plot.data$AFGen.EUR.Z,plot.data$CCScCRYPTCE.Z,pch=16,col="royalblue4",cex=1.1,las=1)
points(plot.data$AFGen.EUR.Z,plot.data$CCScUNDETER.Z,pch=16,col="lightsteelblue1",cex=1.1,las=1)
mtext(side=1,line=2,text="Z-score in AFGen", cex=0.6); mtext(side=2,line=2,text="Z-score in 2nd GWAS", cex=0.6)

## Add a legend
legend("topleft",legend=c("All undetermined", "Incomplete/unclassified","Cryptogenic/cardioembolic minor"),lwd=3,col=c("lightsteelblue1","cornflowerblue","royalblue4"),cex=0.7,bty="n")








