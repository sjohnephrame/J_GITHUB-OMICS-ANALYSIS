rnaseq2 <- read.csv("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\RNASeq2.txt", header=T, sep="\t")
rnaseq2

#Read in file
rnaseq <- read.csv("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\RNASeq2.txt", header=T, sep="\t")
#rnaseq
names(rnaseq)
#Choose only saline TMG
df <- rnaseq[ -c(3:11) ]
df2 <- df[ -c(8:10) ]
#Change col names
colnames(df2) <- c('GENES','TMG1','TMG2','TMG3','SALINE1','SALINE2','SALINE3')
#RELOCATE SALINE TO FIRST 3 COLUMNS
library(dplyr)
TMG2 <- df2 %>%
  relocate(SALINE1, .before = TMG1)%>%
  relocate(SALINE2, .before = TMG2)%>%
  relocate(SALINE3, .before = TMG3)

TMG2  

#TO SEPARATE ENSEMBL ID AND GENE NAME

library(dplyr)
library(tidyr)

TMG<-TMG2 %>% separate(GENES, c('ENSEMBLID', 'GENENAME'))
TMG

#detach(TMG)                      # Detaching first data frame

#Edge R:  differential expression analysis of digital gene expression data
library(edgeR)
library(GEOquery)


###---------------------------
### (i) edgeR 
###---------------------------
### Genes filltering based on counts
#dat <- dat[rowSums(dat)!=0,]     ## Remove if gene counts is zero for all samples
#TMG<-as.data.frame(TMG) #eror

#plotMDS(TMG, main="Multipledimensional scaling plot (MDS)") # error, 

#so creaated dge list

y1<-DGEList(counts=as.matrix(TMG[,3:8]), genes=TMG[,1:2])
y1
dim(y1)
#[1] 55298     6

o <- order(rowSums(y1$counts),decreasing=TRUE) #doeS not seem to be in proper decreasing order...
y1 <- y1[o,]
y1
dim(y1)
#[1] 55298     6

d <- duplicated(y1$genes$GENENAME)
y1 <- y1[!d,]
nrow(y1)
dim(y1)
#[1] 53314     6 after removing duplicates

y1 #53000 from 55000 after removing duplicates

#y1 <- calcNormFactors(y1)
y1<-calcNormFactors(y1, method="TMM")
y1
dim(y1)
#[1] 53314     6

#y1<-calcNormFactors(y1, method="upperquartile")

y1$samples
#y1$genes
#y1$counts

## Data exploration
### MDS (Multidimensional Scaling) plot:

plotMDS(y1, main="Multipledimensional scaling plot (MDS)") #TMG3 is grouping

dim(y1) #before removing count 0
#[1] 53314     6

## Remove if gene counts is zero for all samples

y1 <- y1[rowSums(y1$counts)!=0,]     
y1
dim(y1)
#[1] 25277     6

y1$samples
y1$genes
y1$counts
dim(y1)
#[1] 25277     6

#MDS AFTER REMOVING GENE COUNTS =0
plotMDS(y1, main="Multipledimensional scaling plot (MDS2)")

### Mean Variance plot by gene
###-----------------------------------------------
y1
mean.x <- apply(y1,1,mean)
var.x <- apply(y1,1,var) 
plot(log10(mean.x),log10(var.x),pch=20, xlab="Mean (in log 10 scale)",ylab="Variance (in log 10 scale)",
     xlim=c(0,9),ylim=c(0,9),main="Variance vs Mean")
abline(0,1,col=3,lwd=2)

### Barplot
###-----------------------------------------------
y2<-as.matrix(y1)
y2
#lib.size <- colSums(y1)
lib.size <- colSums(y2)   
y2
barplot(lib.size,xaxt="n",xlab="Study samples",ylab="Library size values",main="Barplot of library size",
        col=c(rep("lightgreen",1),rep("lightcoral",1)),sub="(Red line represents mean)")
abline(h=mean(lib.size),lwd=2,col="red")
legend("bottomright",c("SALINE","TMG"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	

### Boxplot (In log10 scale)
###-----------------------------------------------
boxplot(x=as.list(as.data.frame(log10(y2+1))),xlab="Samples",ylab="Values in log10 scale",main="Boxplot across study samples",
        col=c(rep("lightgreen",1),rep("lightcoral",1)))
legend("topright",c("SALINE","TMG"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	

### Density plots by sample
###-----------------------------------------------
for (i in 1:ncol(y2)){
  if (i==1){
    plot(density(log10(y2[,1])),main="Density plot across study samples",xlab="Subjects",col="gray",ylim=c(0,0.6))}
  else {
    den <- density(log10(y2[,i]))
    lines(den$x,den$y,col="gray")}
}

#with filtering nothing was significant. So also tried without filtering.







#Filter and keep only genes with cpm>=2

## edgeR DE analysis 
##-------------------------------------------------------------------------------------
keep <- rowSums(cpm(as.matrix(y1))>1)>=2     
# At least 2 samples have to have cpm > 1.
dat.filtered <- y1[keep,]
dim(dat.filtered)
#[1] 12711     6

dat.filtered

#rm(keep)
#d <- DGEList(counts=as.matrix(dat.filtered), lib.size=colSums(dat.filtered), group=c(rep("con",3),rep("case",3)))
dim(dat.filtered)
#d <- calcNormFactors(d, method="TMM")        ## Calculates normalization factors
dat.filtered <- estimateDisp(dat.filtered)	  ## Calculates genewise dispersion parameter adjusted using bayesian empirical method(Recommended, edgeR)
dat.filtered$samples$lib.size <- colSums(dat.filtered$counts)
dat.filtered
### PLOTS 

#after filtering and retaining only those rows 
#in the original count data y1 that have at least two counts per million values 
#greater than 1.
###-----------------------------------------------
dat.filtered
dim(dat.filtered)
#[1] 12711     6

#MDS AFTER Filtering and keeping only genes with cpm>=2
plotMDS(dat.filtered, main="Multipledimensional scaling plot (MDS after filtering)")

### Mean Variance plot by gene
###-----------------------------------------------
dat.filtered
dim(dat.filtered)
#[1] 12711     6

mean.x <- apply(dat.filtered,1,mean)
var.x <- apply(dat.filtered,1,var) 
plot(log10(mean.x),log10(var.x),pch=20, xlab="Mean (in log 10 scale)",ylab="Variance (in log 10 scale)",
     xlim=c(0,9),ylim=c(0,9),main="Variance vs Mean")
abline(0,1,col=3,lwd=2)

### Barplot
###-----------------------------------------------
dat.filtered2<-as.matrix(dat.filtered)
dat.filtered2
dim(dat.filtered2)
#[1] 12711     6

lib.size <- colSums(dat.filtered2)   
barplot(lib.size,xaxt="n",xlab="Study samples",ylab="Library size values",main="Barplot of library size",
        col=c(rep("lightgreen",1),rep("lightcoral",1)),sub="(Red line represents mean)")
abline(h=mean(lib.size),lwd=2,col="red")
legend("bottomright",c("SALINE","TMG"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	

### Boxplot (In log10 scale)
###-----------------------------------------------
boxplot(x=as.list(as.data.frame(log10(dat.filtered2+1))),xlab="Samples",ylab="Values in log10 scale",main="Boxplot across study samples",
        col=c(rep("lightgreen",1),rep("lightcoral",1)))
legend("topright",c("SALINE","TMG"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	

### Density plots by sample
###-----------------------------------------------
for (i in 1:ncol(dat.filtered2)){
  if (i==1){
    plot(density(log10(dat.filtered2[,1])),main="Density plot across study samples",xlab="Subjects",col="gray",ylim=c(0,0.6))}
  else {
    den <- density(log10(dat.filtered2[,i]))
    lines(den$x,den$y,col="gray")}
}


## BCV 
plotBCV(dat.filtered)


#The plotBCV function is commonly used in bioinformatics for assessing 
#the Biological Coefficient of Variation (BCV) in RNA-Seq data. 
#BCV is a measure of how much the expression levels of genes vary between biological replicates.

#DIFFERENTIAL EXPRESSION ANALYSIS

dat.filtered
dat.filtered2

#For the example provided

rawdata<-read.csv("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\Tuch_data.csv", header=T, check.names=FALSE, stringsAsFactors=FALSE)
dim(rawdata)
head(rawdata)
y<-DGEList(counts=rawdata[,4:9], genes=rawdata[,1:3])
y
#The study by Tuch et al. [22] was undertaken a few years ago,
# so not all of the RefSeq IDs
#provided by match RefSeq IDs currently in use. 
#We retain only those transcripts with IDs
#in the current NCBI annotation, which is provided by the 
#org.HS.eg.db package:

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)

idfound <- y$genes$RefSeqID %in% mappedRkeys(org.Hs.egREFSEQ)
y <- y[idfound,]
table(idfound)
y

#We add Entrez Gene IDs to the annotation
egREFSEQ <- toTable(org.Hs.egREFSEQ)
head(egREFSEQ)

m <- match(y$genes$RefSeqID, egREFSEQ$accession)
y$genes$EntrezGene <- egREFSEQ$gene_id[m]

## Now use Entrez Gene IDs to find gene symbol
egSYMBOL <- toTable(org.Hs.egSYMBOL)
head(egSYMBOL)
m <- match(y$genes$EntrezGene, egSYMBOL$gene_id)
y$genes$Symbol <- egSYMBOL$symbol[m]
head(y$genes)

#Different RefSeq transcripts for the same gene symbol count #predominantly the same reads.
#So we keep one transcript for each gene symbol. We choose the #transcript with highest
#overall count:
y
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$Symbol)
y <- y[!d,]
nrow(y)
#[1] 10522

#Normally we would also filter lowly expressed genes. For this #data, all transcripts already
#have at least 50 reads for all samples of at least one of the #tissues types.

#Recompute the library sizes:
y$samples$lib.size <- colSums(y$counts)

#Use Entrez Gene IDs as row names:
rownames(y$counts) <- rownames(y$genes) <- y$genes$EntrezGene
y$genes$EntrezGene <- NULL

## normalization
#TMM normalization is applied to this dataset to account for #compositional difference between
#the libraries as opposed to "upperquartile".
y <- calcNormFactors(y)
y2<-calcNormFactors(y, method="upperquartile")

y$samples

## Data exploration
plotMDS(y)

y


y
dat.filtered
dat.filtered2
#dat.filtered2<-as.matrix(dat.filtered)

#change rownames to gene names
rownames(dat.filtered$counts) <- rownames(dat.filtered$genes) <- dat.filtered$genes$GENENAME
dat.filtered
dat.filtered2<-as.matrix(dat.filtered)
dat.filtered2
y

## Design Matrix
Patient <- factor(c(8,8,33,33,51,51))
Tissue <- factor(c("N","T","N","T","N","T"))
data.frame(Sample=colnames(y),Patient,Tissue)
design <- model.matrix(~Patient+Tissue)
rownames(design) <- colnames(y)
design

## Design Matrix
GROUP <- factor(c(1,1,2,2,3,3))
TREATMENT <- factor(c("SALINE","TMG","SALINE","TMG","SALINE","TMG"))
data.frame(Sample=colnames(dat.filtered),GROUP,TREATMENT)
design1 <- model.matrix(~GROUP+TREATMENT)
rownames(design1) <- colnames(dat.filtered)
design1

##Estimating common dispersion followed by tagwise

y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
#Disp = 0.15948 , BCV = 0.3993 

dat.filtered <- estimateGLMCommonDisp(dat.filtered, design1, verbose=TRUE)
#Disp = 0.03273 , BCV = 0.1809 


#The square root of the common dispersion gives 
#the coefficient of variation of biological
#variation. Here the common dispersion is found to be 0.16, 
#so the coefficient of biological
#variation is around 0.4. Then we estimate gene-wise dispersion 
#estimates, allowing a possible trend with averge count size:

y <- estimateGLMTrendedDisp(y, design)
y2 <- estimateGLMTagwiseDisp(y, design)  

dat.filtered <- estimateGLMTrendedDisp(dat.filtered, design1)
dat.filtered1 <- estimateGLMTagwiseDisp(dat.filtered, design1)  


##tagwise method = This function implements the empirical Bayes 
#strategy proposed by McCarthy et al (2012) for estimating the 
#tagwise negative binomial dispersions. The experimental 
#conditions are specified by design matrix allowing for multiple 
#explanatory factors. The empirical Bayes posterior is 
#implemented as a conditional likelihood with tag-specific 
#weights, and the conditional likelihood is computed using Cox-
#Reid approximate conditional likelihood (Cox and Reid, 1987).
#The prior degrees of freedom determines the weight given to the 
#global dispersion trend. The larger the prior degrees of 
#freedom, the more the tagwise dispersions are squeezed towards 
#the global trend.

plot(y$trended.dispersion, y2$tagwise.dispersion)
abline(0,1, col="red")


plot(dat.filtered$trended.dispersion, dat.filtered1$tagwise.dispersion)
abline(0,1, col="red")

#Now proceed to determine differentially expressed genes. Fit #genewise glms. 
# will be fitting model adjusting for patient (paired).
#If dispersion=NULL will be extracted from y, with order of #precedence: genewise dispersion, trended dispersions, common 
#dispersion.

fit <- glmFit(y, design) ## will use trend
#Conduct likelihood ratio tests for tumour vs normal tissue differences and show the top
#genes:
lrt <- glmLRT(fit)

topTags(lrt)


fit1 <- glmFit(dat.filtered, design1) ## will use trend
#Conduct likelihood ratio tests for tumour vs normal tissue differences and show the top
#genes:
lrt1 <- glmLRT(fit1)

topTags(lrt1)


#data for top tags
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]

summary(de<-decideTestsDGE(lrt))


o <- order(lrt1$table$PValue)
cpm(dat.filtered)[o[1:10],]

summary(de<-decideTestsDGE(lrt1))
dim(dat.filtered)

detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")

detags1 <- rownames(dat.filtered)[as.logical(de)]
plotSmear(lrt1, de.tags=detags1)
abline(h=c(-1, 1), col="blue")

#Thus no significan genes with 1 OR 2 WEEKS TMG TREATMENT

###-----------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
#TEST USING EXACT METHOD
#--------------------------------------------------------------------------------


#Gene Set Enrichment analysis using fgsea package in R
# https://bioconductor.org/packages/release/bioc/html/fgsea.html
# http://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
# https://github.com/ctlab/fgsea

####---------------------------------------------
#### fgsea Example 
####---------------------------------------------
rm(list=ls())
 if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
 BiocManager::install("fgsea")
library(data.table)
library(fgsea)
library(ggplot2)
library(edgeR)
library(RColorBrewer)


### Load Data example
load("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\GSEAexample.RData")
class(dat)
head(dat)
colnames(dat)
dim(dat)
# [1] 20382    22
dat <- dat[1:round(nrow(dat)/4),]   ## Select 1/4th of the genes for demonstration
dim(dat)                            ## Using all data run time will be longer
dat

#mydata

rnaseq <- read.csv("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\RNASeq2.txt", header=T, sep="\t")
#rnaseq
names(rnaseq)
#Choose only saline TMG
df <- rnaseq[ -c(3:11) ]
df2 <- df[ -c(8:10) ]
names(df2)
df2<-as.data.frame(df2)
df2
#Change col names
colnames(df2) <- c('GENES','TMG1','TMG2','TMG3','SALINE1','SALINE2','SALINE3')
#RELOCATE SALINE TO FIRST 3 COLUMNS
df2
library(dplyr)
TMG2 <- df2 %>%
  relocate(SALINE1, .before = TMG1)%>%
  relocate(SALINE2, .before = TMG1)%>%
  relocate(SALINE3, .before = TMG1)

TMG2  

#TO SEPARATE ENSEMBL ID AND GENE NAME

library(dplyr)
library(tidyr)

TMG1<-TMG2 %>% separate(GENES, c('ENSEMBLID', 'GENENAME'))
TMG1
TMG<-as.data.frame(TMG1)
TMG<-TMG[,-2]
# Assign the 'ID' column as row names
rownames(TMG) <- TMG$ENSEMBLID

# Remove the 'ID' column from the data frame (if needed)
TMG$ENSEMBLID <- NULL
TMG
dim(TMG)
#[1] 55298     6

dat
dim(dat)
#[1] 5096   22

###---------------------------
### DE analysis using edgeR 
###---------------------------
### Genes filltering based on counts
dat <- dat[rowSums(dat)!=0,]     ## Remove if gene counts is zero for all samples
dim(dat)
#[1] 5096   22

TMG <- TMG[rowSums(TMG)!=0,]     ## Remove if gene counts is zero for all samples
dim(TMG)
#[1] 25705     6

### MDS (Multidimensional Scaling) plot:
plotMDS(dat, main="Multipledimensional scaling plot (MDS)")
plotMDS(TMG, main="Multipledimensional scaling plot (MDS)")


### Mean Variance plot by gene
###-----------------------------------------------
mean.x <- apply(dat,1,mean)
var.x <- apply(dat,1,var) 
plot(log10(mean.x),log10(var.x),pch=20, xlab="Mean (in log 10 scale)",ylab="Variance (in log 10 scale)",
     xlim=c(0,9),ylim=c(0,9),main="Variance vs Mean")
abline(0,1,col=3,lwd=2)


mean.x <- apply(TMG,1,mean)
var.x <- apply(TMG,1,var) 
plot(log10(mean.x),log10(var.x),pch=20, xlab="Mean (in log 10 scale)",ylab="Variance (in log 10 scale)",
     xlim=c(0,9),ylim=c(0,9),main="Variance vs Mean for TMG data")
abline(0,1,col=3,lwd=2)



### Barplot
###-----------------------------------------------
lib.size <- colSums(dat)   
barplot(lib.size,xaxt="n",xlab="Study samples",ylab="Library size values",main="Barplot of library size",
        col=c(rep("lightgreen",11),rep("lightcoral",11)),sub="(Red line represents mean)")
abline(h=mean(lib.size),lwd=2,col="red")
legend("topright",c("Normal","Tumor"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	


lib.size <- colSums(TMG)   
barplot(lib.size,xaxt="n",xlab="Study samples",ylab="Library size values",main="Barplot of library size",
        col=c(rep("lightgreen",3),rep("lightcoral",3)),sub="(Red line represents mean)")
abline(h=mean(lib.size),lwd=2,col="red")
legend("topright",c("SALINE","TMG"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	


### Boxplot (In log10 scale)
###-----------------------------------------------
boxplot(x=as.list(as.data.frame(log10(dat+1))),xlab="Samples",ylab="Values in log10 scale",main="Boxplot across study samples",
        col=c(rep("lightgreen",11),rep("lightcoral",11)))
legend("topright",c("Normal","Tumor"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	


boxplot(x=as.list(as.data.frame(log10(TMG+1))),xlab="Samples",ylab="Values in log10 scale",main="Boxplot across study samples",
        col=c(rep("lightgreen",3),rep("lightcoral",3)))
legend("topright",c("SALINE","TMG"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	


### Density plots by sample
###-----------------------------------------------

for (i in 1:ncol(dat)){
  if (i==1){
    plot(density(log10(dat[,1])),main="Density plot across study samples",xlab="Subjects",col="gray",ylim=c(0,0.8))}
  else {
    den <- density(log10(dat[,i]))
    lines(den$x,den$y,col="gray")}
}



for (i in 1:ncol(TMG)){
  if (i==1){
    plot(density(log10(TMG[,1])),main="Density plot across study samples for TMG data",xlab="Subjects",col="gray",ylim=c(0,0.8))}
  else {
    den <- density(log10(TMG[,i]))
    lines(den$x,den$y,col="gray")}
}


## edgeR DE analysis 
##-------------------------------------------------------------------------------------
keep <- rowSums(cpm(as.matrix(dat))>1)>=2     
# At least 2 samples have to have cpm > 1.
data.filtered <- dat[keep,]
dim(data.filtered)
rm(keep)
d <- DGEList(counts=as.matrix(data.filtered), lib.size=colSums(data.filtered), group=c(rep("Normal",11),rep("Tumor",11)))
dim(d)
d <- calcNormFactors(d, method="TMM")        ## Calculates normalization factors
d <- estimateDisp(d)                   		 ## Calculates genewise dispersion parameter adjusted using bayesian empirical method(Recommended, edgeR)
d

keep1 <- rowSums(cpm(as.matrix(TMG))>1)>=2     
# At least 2 samples have to have cpm > 1.
data.filtered1 <- TMG[keep1,]
dim(data.filtered1)
#[1] 12775     6
rm(keep1)
d1 <- DGEList(counts=as.matrix(data.filtered1), lib.size=colSums(data.filtered1), group=c(rep("SALINE",3),rep("TMG",3)))
dim(d1)
d1
d1 <- calcNormFactors(d1, method="TMM")        ## Calculates normalization factors
d1
d1 <- estimateDisp(d1)                   		 ## Calculates genewise dispersion parameter adjusted using bayesian empirical method(Recommended, edgeR)
d1



## BCV 
plotBCV(d)
plotBCV(d1)

de.test <- exactTest(d,pair=c("Normal","Tumor")) ## First entry under pair is the baseline 
de.test.FDR <- topTags(de.test,n=Inf,adjust.method="BH", sort.by="PValue") 
head(de.test.FDR$table)
summary(de <- decideTestsDGE(de.test, p=0.05, adjust="BH"))  ## Counts up and down regulated genes	

de.test1 <- exactTest(d1,pair=c("SALINE","TMG")) ## First entry under pair is the baseline 
de.test1.FDR <- topTags(de.test1,n=Inf,adjust.method="BH", sort.by="PValue") 
head(de.test1.FDR$table)
summary(de1 <- decideTestsDGE(de.test1, p=0.05, adjust="BH"))  ## Counts up and down regulated genes	
dim(d1)

### MDS (Multidimensional Scaling) plot:
plotMDS(d1, main="Multipledimensional scaling plot (MDS)")


### Mean Variance plot by gene
###-----------------------------------------------

mean.x <- apply(d1,1,mean)
var.x <- apply(d1,1,var) 
plot(log10(mean.x),log10(var.x),pch=20, xlab="Mean (in log 10 scale)",ylab="Variance (in log 10 scale)",
     xlim=c(0,9),ylim=c(0,9),main="Variance vs Mean for TMG data")
abline(0,1,col=3,lwd=2)



### Barplot
###-----------------------------------------------
d2<-as.matrix(d1)
lib.size <- colSums(d2) 
d2
d1
barplot(lib.size,xaxt="n",xlab="Study samples",ylab="Library size values",main="Barplot of library size",
        col=c(rep("lightgreen",3),rep("lightcoral",3)),sub="(Red line represents mean)")
abline(h=mean(lib.size),lwd=2,col="red")
legend("topright",c("SALINE","TMG"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	


### Boxplot (In log10 scale)
###-----------------------------------------------
boxplot(x=as.list(as.data.frame(log10(dat+1))),xlab="Samples",ylab="Values in log10 scale",main="Boxplot across study samples",
        col=c(rep("lightgreen",11),rep("lightcoral",11)))
legend("topright",c("Normal","Tumor"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	


boxplot(x=as.list(as.data.frame(log10(d2+1))),xlab="Samples",ylab="Values in log10 scale",main="Boxplot across study samples",
        col=c(rep("lightgreen",3),rep("lightcoral",3)))
legend("topright",c("SALINE","TMG"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	


### Density plots by sample
###-----------------------------------------------

for (i in 1:ncol(dat)){
  if (i==1){
    plot(density(log10(dat[,1])),main="Density plot across study samples",xlab="Subjects",col="gray",ylim=c(0,0.8))}
  else {
    den <- density(log10(dat[,i]))
    lines(den$x,den$y,col="gray")}
}



for (i in 1:ncol(d2)){
  if (i==1){
    plot(density(log10(d2[,1])),main="Density plot across study samples for TMG data",xlab="Subjects",col="gray",ylim=c(0,0.8))}
  else {
    den <- density(log10(d2[,i]))
    lines(den$x,den$y,col="gray")}
}


# ### plotSmear
# ###---------------------------------------------------------
detags <- rownames(d)[as.logical(de)]
plotSmear(de.test, de.tags=detags,main="Smear Plot",sub="(The horizontal blue lines show 4-fold changes)",cex=1) 
abline(h = c(-2, 2), col = "blue")

detags1 <- rownames(d1)[as.logical(de1)]
plotSmear(de.test1, de.tags=detags1,main="Smear Plot",sub="(The horizontal blue lines show 4-fold changes)",cex=1) 
abline(h = c(-2, 2), col = "blue")


### Volcano plot
###---------------------------------------------------------
d.volcano <-  de.test.FDR$table[,c("logFC","PValue")]
par(mar=c(5,4,4,5))
plot(d.volcano$logFC, -log(d.volcano$PValue,10), main="",pch=20, cex=2,xlab=expression(log[2]~fold~change), ylab=expression(-log[10]~pvalue),
     xlim=c(min(d.volcano$logFC)-1,max(d.volcano$logFC)+1))
title("Volcano plot")


# Log2 fold change and p-value cutoff
lfc <- 2
pval <- 0.05 
# Selecting interesting genes
sigGenes <- (abs(d.volcano$logFC)> lfc & -log(d.volcano$PValue,10) > -log10(pval))   
# Identifying the selected genes
points(d.volcano[sigGenes,]$logFC,-log(d.volcano[sigGenes,]$PValue,10),pch=20,col="red",cex=2)
abline(h=-log10(pval),col="green3",lty=2)
abline(v=c(-lfc,lfc),col="blue",lty=2)
mtext(paste("pval = ",round(pval,2)),side=4,at=-log10(pval),cex=0.8,line=0.5,las=1)
mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)


d.volcano1 <-  de.test1.FDR$table[,c("logFC","PValue")]
par(mar=c(5,4,4,5))
plot(d.volcano1$logFC, -log(d.volcano1$PValue,10), main="",pch=20, cex=2,xlab=expression(log[2]~fold~change), ylab=expression(-log[10]~pvalue),
     xlim=c(min(d.volcano1$logFC)-1,max(d.volcano1$logFC)+1))
title("Volcano plot for TMG data")

# Log2 fold change and p-value cutoff
lfc1 <- 2
pval1 <- 0.05 
# Selecting interesting genes
sigGenes1 <- (abs(d.volcano1$logFC)> lfc1 & -log(d.volcano1$PValue,10) > -log10(pval1))   
# Identifying the selected genes
points(d.volcano1[sigGenes1,]$logFC,-log(d.volcano1[sigGenes1,]$PValue,10),pch=20,col="red",cex=2)
abline(h=-log10(pval1),col="green3",lty=2)
abline(v=c(-lfc1,lfc1),col="blue",lty=2)
mtext(paste("pval = ",round(pval1,2)),side=4,at=-log10(pval1),cex=0.8,line=0.5,las=1)
mtext(c(paste("-",lfc1,"fold"),paste("+",lfc1,"fold")),side=3,at=c(-lfc1,lfc1),cex=1,line=0.2)


d.edgeR.result <- de.test.FDR$table
head(d.edgeR.result)

d1.edgeR.result <- de.test1.FDR$table
head(d1.edgeR.result)



# Map Ensembl gene IDs to symbol. First create a mapping table.
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensg.id <- rownames(d.edgeR.result)
d.gene.id <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values=ensg.id, mart= mart)
dim(d.gene.id)
head(d.gene.id)
d.gene.id <- subset(d.gene.id,hgnc_symbol!="")  # Remove the ENSG ids that do not have hgnc_symbol 
rownames(d.gene.id) <- d.gene.id$ensembl_gene_id
dim(d.gene.id)



library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
ensg.id <- rownames(d1.edgeR.result)
d.gene.id <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"), values=ensg.id, mart=mart)
#attributes=mgi_symbol for mice
dim(d.gene.id)
head(d.gene.id)
d.gene.id <- subset(d.gene.id,mgi_symbol!="")  # Remove the ENSG ids that do not have mgi_symbol 
rownames(d.gene.id) <- d.gene.id$ensembl_gene_id
dim(d.gene.id)
head(d.gene.id)

listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
ensembl
datasets <- listDatasets(ensembl)
head(datasets)
searchDatasets(mart = ensembl, pattern = "mmusculus")


###Combine the datasets together
common.genes <- intersect(rownames(d.edgeR.result),rownames(d.gene.id))
length(common.genes)
length(unique(common.genes))

common.genes <- intersect(rownames(d1.edgeR.result),rownames(d.gene.id))
length(common.genes)
length(unique(common.genes))

# Merge gene id data with the results data
d.tmp1 <- d.gene.id[common.genes,]
d.tmp2 <- d.edgeR.result[common.genes,]   
all(rownames(d.tmp1)==rownames(d.tmp2))

d.tmp1 <- d.gene.id[common.genes,]
d.tmp2 <- d1.edgeR.result[common.genes,]   
all(rownames(d.tmp1)==rownames(d.tmp2))

# TRUE
d.tmp <- cbind(d.tmp1,d.tmp2)
head(d.tmp)
d.tmp <- d.tmp[order(d.tmp$PValue),] ## Sort the data in the order of significance
head(d.tmp)

d.tmp <- cbind(d.tmp1,d.tmp2)
head(d.tmp)
d.tmp <- d.tmp[order(d.tmp$PValue),] ## Sort the data in the order of significance
head(d.tmp)


# Prepare the gene list
gene.list <- d.tmp$logFC             # rank-ordered gene list
names(gene.list) <- d.tmp$hgnc_symbol

gene.list <- d.tmp$logFC             # rank-ordered gene list
names(gene.list) <- d.tmp$mgi_symbol

# Barplot of ranked fold changes
barplot(sort(gene.list, decreasing = T),axisnames=FALSE,main="Plot of ranked gene Fold changes")

barplot(sort(gene.list, decreasing = T),axisnames=FALSE,main="Plot of ranked gene Fold changes")

###GO analysis
library(org.Hs.eg.db)
library(edgeR)
library(GO.db)
d.go <- d.tmp
#d.go.DE <- subset(d.go,FDR<0.05) 
d.go.DE <- subset(d.go,PValue<0.05) 
d.entrez.id <- mapIds(org.Hs.eg.db, keys=d.go.DE$hgnc_symbol,column="ENTREZID",keytype="SYMBOL")
length(d.entrez.id)
head(d.entrez.id)
all(d.go.DE$hgnc_symbol==names(d.entrez.id)) 
go.test <- goana(d.entrez.id,species="Hs")
go.results <- topGO(go.test, sort = "DE", number = Inf)
head(go.results)
sum(go.results$P.DE<10^(-5))


library(org.Mm.eg.db)
library(edgeR)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")

library(GO.db)
d.go <- d.tmp
#d.go.DE <- subset(d.go,FDR<0.05) 
d.go.DE <- subset(d.go,PValue<0.05) 
#attributes=mgi_symbol for mice
d.entrez.id <- mapIds(org.Mm.eg.db, keys=d.go.DE$mgi_symbol,column="ENTREZID",keytype="SYMBOL")
length(d.entrez.id)
head(d.entrez.id)
all(d.go.DE$mgi_symbol==names(d.entrez.id)) 
go.test <- goana(d.entrez.id,species="Mm")
go.results <- topGO(go.test, sort = "DE", number = Inf)
head(go.results)
sum(go.results$P.DE<10^(-5))

###KEGG analysis
d.kegg <- d.tmp
#d.kegg.DE <- subset(d.kegg,FDR<0.05) 
d.kegg.DE <- subset(d.kegg,PValue<0.05)
all(d.kegg.DE$hgnc_symbol==names(d.entrez.id)) 

kegg.test <- kegga(d.entrez.id,species="Hs")
kegg.results <- topKEGG(kegg.test, sort = "DE", number = Inf)
head(kegg.results)
sum(kegg.results$P.DE<10^(-5))


d.kegg <- d.tmp
#d.kegg.DE <- subset(d.kegg,FDR<0.05) 
d.kegg.DE <- subset(d.kegg,PValue<0.05)
all(d.kegg.DE$mgi_symbol==names(d.entrez.id)) 
d.entrez.id <- mapIds(org.Mm.eg.db, keys=d.go.DE$mgi_symbol,column="ENTREZID",keytype="SYMBOL")
length(d.entrez.id)
head(d.entrez.id)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("KEGGREST")
library(KEGGREST)
kegg.test <- kegga(d.entrez.id,species="Mm")
kegg.results <- topKEGG(kegg.test, sort = "DE", number = Inf)
head(kegg.results)
sum(kegg.results$P.DE<10^(-5))



### GSEA analysis
# Load All gene sets file downloaded from Broad Institute 
# The following website contains the gene set collection or the complete Molecular Signatures Database (MSigDB)  
# http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
#  
all.gene.sets <- gmtPathways("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\msigdb.v7.4.symbols.gmt")
class(all.gene.sets)
length(all.gene.sets)
all.gene.sets[1:2]
# Show first a few pathways, and within those, show only the first few genes. 
library(tidyverse)
all.gene.sets %>% head() %>% lapply(head)

### Now run fgsea 
fgseaRes <- fgsea(pathways = all.gene.sets,stats = gene.list,minSize=15,maxSize=500,eps=0)
head(fgseaRes)
head(fgseaRes[order(pval), ])
sum(fgseaRes[, padj < 0.05])
# Make a few Enrichment Plots
#Plot1
plotEnrichment(all.gene.sets[["GAO_LARGE_INTESTINE_ADULT_CI_MESENCHYMAL_CELLS"]],gene.list) + labs(title="GAO_LARGE_INTESTINE_ADULT_CI_MESENCHYMAL_CELLS")
#Plot2
plotEnrichment(all.gene.sets[["CHEN_METABOLIC_SYNDROM_NETWORK"]],gene.list) + labs(title="CHEN_METABOLIC_SYNDROM_NETWORK")
tail(fgseaRes[order(pval), ])
#Plot3
plotEnrichment(all.gene.sets[["ZNF610_TARGET_GENES"]],gene.list) + labs(title="ZNF610_TARGET_GENES")
#Plot4
plotEnrichment(all.gene.sets[["ZSCAN5DP_TARGET_GENES"]],gene.list) + labs(title="ZSCAN5DP_TARGET_GENES")
fgseaRes[order(pval), ][15:20,]
#Plot5
plotEnrichment(all.gene.sets[["GOBP_RESPONSE_TO_BIOTIC_STIMULUS"]],gene.list) + labs(title="GOBP_RESPONSE_TO_BIOTIC_STIMULUS")

# Make a table plot for a bunch of selected pathways:
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(all.gene.sets[topPathways], gene.list, fgseaRes,gseaParam = 0.5)



### GSEA analysis
# Load All gene sets file downloaded from Broad Institute 
# The following website contains the gene set collection or the complete Molecular Signatures Database (MSigDB)  
# http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
#  
#mydata

all.gene.sets <- gmtPathways("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\m2.all.v2023.1.Mm.symbols.gmt")        #m2curated gene sets
class(all.gene.sets)
length(all.gene.sets)
all.gene.sets[1:2]
# Show first a few pathways, and within those, show only the first few genes. 
library(tidyverse)
all.gene.sets %>% head() %>% lapply(head)

### Now run fgsea 
fgseaRes <- fgsea(pathways = all.gene.sets,stats = gene.list,minSize=15,maxSize=500,eps=0)
head(fgseaRes)
head(fgseaRes[order(pval), ])
sum(fgseaRes[, padj < 0.05])
# Make a few Enrichment Plots
#Plot1
plotEnrichment(all.gene.sets[["ICHIBA_GRAFT_VERSUS_HOST_DISEASE_D7_UP"]],gene.list) + labs(title="ICHIBA_GRAFT_VERSUS_HOST_DISEASE_D7_UP")
#Plot2
plotEnrichment(all.gene.sets[["ICHIBA_GRAFT_VERSUS_HOST_DISEASE_35D_UP"]],gene.list) + labs(title="ICHIBA_GRAFT_VERSUS_HOST_DISEASE_35D_UP")
#Plot3
plotEnrichment(all.gene.sets[["ZEMEK_IMMUNE_CHECKPOINT_BLOCKADE_OVARIAN_CANCER_OVERLAP_UP"]],gene.list) + labs(title="ZEMEK_IMMUNE_CHECKPOINT_BLOCKADE_OVARIAN_CANCER_OVERLAP_UP")
#Plot4
plotEnrichment(all.gene.sets[["MARKEY_RB1_ACUTE_LOF_DN"]],gene.list) + labs(title="MARKEY_RB1_ACUTE_LOF_DN")



tail(fgseaRes[order(pval), ])
#Plot5
plotEnrichment(all.gene.sets[["REACTOME_VESICLE_MEDIATED_TRANSPORT"]],gene.list) + labs(title="REACTOME_VESICLE_MEDIATED_TRANSPORT")
#Plot6
plotEnrichment(all.gene.sets[["REACTOME_VLDLR_INTERNALISATION_AND_DEGRADATION"]],gene.list) + labs(title="REACTOME_VLDLR_INTERNALISATION_AND_DEGRADATION")
#Plot7
plotEnrichment(all.gene.sets[["SANSOM_APC_MYC_TARGETS"]],gene.list) + labs(title="SANSOM_APC_MYC_TARGETS")

fgseaRes[order(pval), ][15:20,]
#Plot8
plotEnrichment(all.gene.sets[["PLASARI_TGFB1_TARGETS_10HR_UP"]],gene.list) + labs(title="PLASARI_TGFB1_TARGETS_10HR_UP")
#Plot9
plotEnrichment(all.gene.sets[["WP_ELECTRON_TRANSPORT_CHAIN"]],gene.list) + labs(title="WP_ELECTRON_TRANSPORT_CHAIN")
#Plot10
plotEnrichment(all.gene.sets[["MORI_LARGE_PRE_BII_LYMPHOCYTE_DN"]],gene.list) + labs(title="MORI_LARGE_PRE_BII_LYMPHOCYTE_DN")

# Make a table plot for a bunch of selected pathways:
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(all.gene.sets[topPathways], gene.list, fgseaRes,gseaParam = 0.5)

#Top pathways UP
#Plot11
plotEnrichment(all.gene.sets[["LIM_MAMMARY_STEM_CELL_UP"]],gene.list) + labs(title="LIM_MAMMARY_STEM_CELL_UP")
#Plot12
plotEnrichment(all.gene.sets[["WP_ESC_PLURIPOTENCY_PATHWAYS"]],gene.list) + labs(title="WP_ESC_PLURIPOTENCY_PATHWAYS")
#Plot13
plotEnrichment(all.gene.sets[["PLASARI_TGFB1_TARGETS_10HR_UP"]],gene.list) + labs(title="PLASARI_TGFB1_TARGETS_10HR_UP")

#Top pathways DOWN

#Plot14
plotEnrichment(all.gene.sets[["ICHIBA_GRAFT_VERSUS_HOST_DISEASE_D7_UP"]],gene.list) + labs(title="ICHIBA_GRAFT_VERSUS_HOST_DISEASE_D7_UP")
#Plot15
plotEnrichment(all.gene.sets[["ICHIBA_GRAFT_VERSUS_HOST_DISEASE_35D_UP"]],gene.list) + labs(title="ICHIBA_GRAFT_VERSUS_HOST_DISEASE_35D_UP")
#Plot16
plotEnrichment(all.gene.sets[["ZEMEK_IMMUNE_CHECKPOINT_BLOCKADE_OVARIAN_CANCER_OVERLAP_UP"]],gene.list) + labs(title="ZEMEK_IMMUNE_CHECKPOINT_BLOCKADE_OVARIAN_CANCER_OVERLAP_UP")

#################################################################################################

#Clustering

load("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\DatCluster.RData")   
dim(dat)
#[1]    499 118209
# Note: Clinical data : columns 1: 75
#		Molecular data: column 76 onwards 
dim(dat.expr)  # Top 1000 most associated genes with survival

#mydata

rnaseq <- read.csv("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\RNASeq2.txt", header=T, sep="\t")
#rnaseq
names(rnaseq)
#Choose only saline TMG
df <- rnaseq[ -c(3:11) ]
df2 <- df[ -c(8:10) ]
names(df2)
df2<-as.data.frame(df2)
df2
#Change col names
colnames(df2) <- c('GENES','TMG1','TMG2','TMG3','SALINE1','SALINE2','SALINE3')
#RELOCATE SALINE TO FIRST 3 COLUMNS
df2
library(dplyr)
TMG2 <- df2 %>%
  relocate(SALINE1, .before = TMG1)%>%
  relocate(SALINE2, .before = TMG1)%>%
  relocate(SALINE3, .before = TMG1)

TMG2  

#TO SEPARATE ENSEMBL ID AND GENE NAME

library(dplyr)
library(tidyr)

TMG1<-TMG2 %>% separate(GENES, c('ENSEMBLID', 'GENENAME'))
TMG1
TMG<-as.data.frame(TMG1)
TMG<-TMG[,-2]
# Assign the 'ID' column as row names
rownames(TMG) <- TMG$ENSEMBLID

# Remove the 'ID' column from the data frame (if needed)
TMG$ENSEMBLID <- NULL
TMG
dim(TMG)
#[1] 55298     6

dat
dim(dat)
#[1]    499 118210

TMG <- TMG[rowSums(TMG)!=0,]     ## Remove if gene counts is zero for all samples
dim(TMG)
#[1] 25705     6

### MDS (Multidimensional Scaling) plot:
plotMDS(TMG, main="Multipledimensional scaling plot (MDS)")


keep1 <- rowSums(cpm(as.matrix(TMG))>1)>=2     
# At least 2 samples have to have cpm > 1.
data.filtered1 <- TMG[keep1,]
dim(data.filtered1)
#[1] 12775     6
rm(keep1)
dim(data.filtered1)

#d1 <- DGEList(counts=as.matrix(data.filtered1), lib.size=colSums(data.filtered1), group=c(rep("SALINE",3),rep("TMG",3)))
dim(d1)
d1
d1 <- calcNormFactors(d1, method="TMM")        ## Calculates normalization factors
d1
d1 <- estimateDisp(d1)                   		 ## Calculates genewise dispersion parameter adjusted using bayesian empirical method(Recommended, edgeR)
d1



## BCV 
plotBCV(d)
plotBCV(d1)

de.test <- exactTest(d,pair=c("Normal","Tumor")) ## First entry under pair is the baseline 
de.test.FDR <- topTags(de.test,n=Inf,adjust.method="BH", sort.by="PValue") 
head(de.test.FDR$table)
summary(de <- decideTestsDGE(de.test, p=0.05, adjust="BH"))  ## Counts up and down regulated genes	

de.test1 <- exactTest(d1,pair=c("SALINE","TMG")) ## First entry under pair is the baseline 
de.test1.FDR <- topTags(de.test1,n=Inf,adjust.method="BH", sort.by="PValue") 
head(de.test1.FDR$table)
summary(de1 <- decideTestsDGE(de.test1, p=0.05, adjust="BH"))  ## Counts up and down regulated genes	
dim(d1)


## 1. Hierarchical clustering method
##----------------------------------------------------------------------------------------------
dd <- as.matrix(dat.expr)

dd <- as.matrix(data.filtered1)
dd <- as.matrix(d1)

require(graphics)
d <- dist(dd, method = "euclidean") # distance matrix 
h.clust <- hclust(d, method="complete") 
str(h.clust)

plot(h.clust,labels = F,xlab="",sub="") # display dendogram
cluster.mem <- cutree(h.clust, k=3) # cut tree into 3 clusters
# draw dendogram with red borders around the 3 clusters 
rect.hclust(h.clust, k=3, border="red")
table(cluster.mem)

plot(h.clust,labels = F,xlab="",sub="") # display dendogram
cluster.mem <- cutree(h.clust, k=2) # cut tree into 3 clusters
# draw dendogram with red borders around the 3 clusters 
rect.hclust(h.clust, k=2, border="red")
table(cluster.mem)

plot(h.clust,labels = F,xlab="",sub="") # display dendogram
cluster.mem <- cutree(h.clust, k=5) # cut tree into 2 clusters
# draw dendogram with red borders around the 2 clusters 
rect.hclust(h.clust, k=5, border="red")
table(cluster.mem)

plot(h.clust,labels = F,xlab="",sub="") # display dendogram
cluster.mem <- cutree(h.clust, k=10) # cut tree into 2 clusters
# draw dendogram with red borders around the 2 clusters 
rect.hclust(h.clust, k=10, border="red")
table(cluster.mem)

plot(h.clust,labels = F,xlab="",sub="") # display dendogram
cluster.mem <- cutree(h.clust, k=20) # cut tree into 2 clusters
# draw dendogram with red borders around the 2 clusters 
rect.hclust(h.clust, k=20, border="red")
table(cluster.mem)

plot(h.clust,labels = F,xlab="",sub="") # display dendogram
cluster.mem <- cutree(h.clust, k=50) # cut tree into 2 clusters
# draw dendogram with red borders around the 2 clusters 
rect.hclust(h.clust, k=50, border="red")
table(cluster.mem)

plot(h.clust,labels = F,xlab="",sub="") # display dendogram
cluster.mem <- cutree(h.clust, k=100) # cut tree into 2 clusters
# draw dendogram with red borders around the 2 clusters 
rect.hclust(h.clust, k=100, border="red")
table(cluster.mem)

#Not doing this bcoz not clinical data

### Assess the clusters using CoxPH analysis among the clusters
all(names(cluster.mem)==dat$bcr_patient_barcode)
# [1] TRUE
D <- cbind(cluster.mem,dat[,1:75]) # Combine the cluster.id with clinical part of data
D[1:5,1:5]


all(names(cluster.mem)==data.filtered1$rowname)
# [1] TRUE
D <- cbind(cluster.mem,data.filtered1) # Combine the cluster.id with clinical part of data
D[1:5,1:5]


### Log rank test:
###-------------------------------------------------------------------------------------
library(survival)
lr.test <- survdiff(Surv(time.ttr, recurrence == "YES")~as.factor(cluster.mem), data = D, na.action=na.exclude)
lr.test
pval <- pchisq(lr.test$chisq,1, lower.tail=F)
pval
### Kaplan Meier plot of the two clusters  
for (i in 1:length(unique(sort(cluster.mem))))
{
  ii <- unique(sort(cluster.mem))[i]
  ddd <- D[D$cluster.mem==ii,]
  f <- survfit(Surv(time.ttr, recurrence == "YES")~1, data = ddd)
  if (i == 1){plot(f,col="black",conf.int=FALSE,xlab="Time to Recurrence",ylab="Survival",main="",lwd=2)
    clust.name <- "Cluster 1"
    clr <- i}
  else       {lines(f,col=i,conf.int=F,lwd=2)
    clust.name <- c(clust.name,paste("Cluster",i))
    clr <- c(clr,i)}
  text(1500,1,paste("p-value = ",format(pval, digits=2,scientific = T)))
}
legend("topright", clust.name, text.col = clr, cex=0.9)

##########################################################################################


## 2. k-means clustering method
##-----------------------------------


#rm(list=setdiff(ls(),c("dat.expr","dat")))   
#ls()
#dd <- as.matrix(dat.expr)

#dd <- as.matrix(data.filtered1)

dim(dd)
dd[1:5,1:5]
## A plot of the within groups sum of squares by number of clusters extracted can help determine 
## the appropriate number of clusters. 
wss <- (nrow(dd)-1)*sum(apply(dd,2,var))
for (i in 2:8) wss[i] <- sum(kmeans(dd, centers=i)$withinss)
plot(1:8, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

k.clust <- kmeans(dd,3)
str(k.clust)
cluster.mem <- k.clust$cluster
cluster.mem[1:5]
table(cluster.mem)

k.clust <- kmeans(dd,2)
str(k.clust)
cluster.mem <- k.clust$cluster
cluster.mem[1:5]
table(cluster.mem)

k.clust <- kmeans(dd,5)
str(k.clust)
cluster.mem <- k.clust$cluster
cluster.mem[1:5]
table(cluster.mem)

k.clust <- kmeans(dd,8)
str(k.clust)
cluster.mem <- k.clust$cluster
cluster.mem[1:5]
table(cluster.mem)

k.clust <- kmeans(dd,10)
str(k.clust)
cluster.mem <- k.clust$cluster
cluster.mem[1:5]
table(cluster.mem)

k.clust <- kmeans(dd,12)
str(k.clust)
cluster.mem <- k.clust$cluster
cluster.mem[1:5]
table(cluster.mem)

### Assess the clusters using CoxPH analysis among the clusters
###
all(names(cluster.mem)==dat$bcr_patient_barcode)
# [1] TRUE
D <- cbind(cluster.mem,dat[,1:75]) # Combine the cluster.id with clinical part of data
D[1:5,1:5]

all(names(cluster.mem)==data.filtered1$rowname)
# [1] TRUE
D2 <- cbind(cluster.mem,data.filtered1) # Combine the cluster.id with clinical part of data
D2[1:5,1:5]


### Log rank test:
###-------------------------------------------------------------------------------------
library(survival)
lr.test <- survdiff(Surv(time.ttr, recurrence == "YES")~as.factor(cluster.mem), data = D, na.action=na.exclude)
lr.test
pval <- pchisq(lr.test$chisq,1, lower.tail=F)
pval
for (i in 1:length(unique(sort(cluster.mem))))
{
  ii <- unique(sort(cluster.mem))[i]
  ddd <- D[D$cluster.mem==ii,]
  f <- survfit(Surv(time.ttr, recurrence == "YES")~1, data = ddd)
  if (i == 1){plot(f,col="black",conf.int=FALSE,xlab="Time to Recurrence",ylab="Survival",main="",lwd=2)
    clust.name <- "Cluster 1"
    clr <- i}
  else       {lines(f,col=i,conf.int=F,lwd=2)
    clust.name <- c(clust.name,paste("Cluster",i))
    clr <- c(clr,i)}
  text(1300,1,paste("p-value = ",format(pval, digits=2,scientific = T)))
}
legend("topright", clust.name, text.col = clr, cex=0.9)

dd
heatmap(dd, main = "SALINE Vs TMG", xlab = "Samples", ylab = "DEG")
#https://www.biostars.org/p/374551/








####################################################################################################


#retry using glmfit
























##Estimating common dispersion followed by tagwise

y <- estimateGLMCommonDisp(y, design, verbose=TRUE)

#The square root of the common dispersion gives 
#the coefficient of variation of biological
#variation. Here the common dispersion is found to be 0.16, 
#so the coefficient of biological
#variation is around 0.4. Then we estimate gene-wise dispersion 
#estimates, allowing a possible trend with averge count size:

y <- estimateGLMTrendedDisp(y, design)
y1 <- estimateGLMTagwiseDisp(y, design)  

##tagwise method = This function implements the empirical Bayes 
#strategy proposed by McCarthy et al (2012) for estimating the 
#tagwise negative binomial dispersions. The experimental 
#conditions are specified by design matrix allowing for multiple 
#explanatory factors. The empirical Bayes posterior is 
#implemented as a conditional likelihood with tag-specific 
#weights, and the conditional likelihood is computed using Cox-
#Reid approximate conditional likelihood (Cox and Reid, 1987).
#The prior degrees of freedom determines the weight given to the 
#global dispersion trend. The larger the prior degrees of 
#freedom, the more the tagwise dispersions are squeezed towards 
#the global trend.

plot(y$trended.dispersion, y1$tagwise.dispersion)
abline(0,1, col="red")

#Now proceed to determine differentially expressed genes. Fit #genewise glms. 
# will be fitting model adjusting for patient (paired).
#If dispersion=NULL will be extracted from y, with order of #precedence: genewise dispersion, trended dispersions, common 
#dispersion.

fit <- glmFit(y, design) ## will use trend
#Conduct likelihood ratio tests for tumour vs normal tissue differences and show the top
#genes:
lrt <- glmLRT(fit)

topTags(lrt)

#data for top tags
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]

summary(de<-decideTestsDGE(lrt))

detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")



###-----------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------


  
  
###################################################################################################




#Edge R:  differential expression analysis of digital gene expression data


####---------------------------------------------
#### edgeR Example 
####---------------------------------------------
#dir1<-"C:/Users/pchalise/Desktop/Genomics teaching material/BIOS 855/rpgm_chalise/RNAseqData"
#dir1<-"C:/Users/Sophiya John/Desktop/J_STATISTICS/J_STATISTICAL GENOMICS BIOS 799/edgeR_Example.csv"
dir1<-read.csv("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\edgeR_Example.csv")
attach(dir1)

#Commands to run edge R
#Opening R in the environment

#need to type in anything that has 
#source("http://bioconductor.org/biocLite.R")

#biocLite()
#installing edgeR
biocLite("edgeR") 
#biocLite("GEOquery")
#Using Case Studies from edge R User's Guide (Chapter 4)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("GEOquery")

#Using edge R in the environment
library(edgeR)
library(GEOquery)

#gse10782 <- getGEO('GSE10782',GSEMatrix=TRUE)
gsm1 <- getGEO("GSM272105") #Gene Expression Omnibus (GEO) is a database repository of high throughput gene expression data and hybridization arrays, chips, microarrays.
gsm2 <- getGEO("GSM272106")
gsm3 <- getGEO("GSM272318")
gsm4 <- getGEO("GSM272319")

##Pulling off counts and combining 4 samples data into one matrix
tmp1<-Table(gsm1) #TPM (transcripts per kilobase million) Counts per length of transcript (kb) per million reads mapped. 
names(tmp1)<-c("seq", "sample1.c", "sample1.tpm")
tmp1<-tmp1[,-3]
tmp2<-Table(gsm2)
names(tmp2)<-c("seq", "sample2.c", "sample2.tpm")
tmp2<-tmp2[,-3]
tmp3<-Table(gsm3)
names(tmp3)<-c("seq", "sample3.c", "sample3.tpm")
tmp3<-tmp3[,-3]
tmp4<-Table(gsm4)
names(tmp4)<-c("seq", "sample4.c", "sample4.tpm")
tmp4<-tmp4[,-3]

x<-merge(tmp1, tmp2, by.all="seq")
x<-merge(x, tmp3, by.all="seq")
x<-merge(x, tmp4, by.all="seq")

x
## Data needs to be numeric
x[,2]<-as.numeric(x[,2])
x[,3]<-as.numeric(x[,3])
x[,4]<-as.numeric(x[,4])
x[,5]<-as.numeric(x[,5])

## Put seq as rowname
dimnames(x)[[1]]<-x$seq

#Put the counts and other information into a DGEList object:
y <- DGEList(counts=x[,2:5], group=c("DCLK", "WT", "DCLK", "WT"))
y

colnames(y) <- c("DCLK1","WT1","DCLK2","WT2")
y

#For this dataset there were over 66000 unique tags sequenced, most of which have a very
#small number of counts in total across all libraries. We want to keep tags that are expressed
#in at least one of wild-type or transgenic mice. In either case, the tag should be expressed
#in at least two libraries. We seek tags that achieve one count per million for at least two
#libraries:
keep <- rowSums(cpm(y) > 1) >= 2
table(keep)
#keep
#FALSE  TRUE 
#15952 50051 
y <- y[keep,]
y
#Having filtered, reset the library sizes:
y$samples$lib.size <- colSums(y$counts)

#Normalization: we align the upper-quartiles of the counts-per-million between the
#libraries:

y <- calcNormFactors(y,method="upperquartile")

#Data normalization is a crucial preliminary step in analyzing genomic datasets. 
#The goal  of normalization is to remove global variation to make readings 
#across different experiments comparable.

#plot showing the sample relations based on multidimensional scaling:
plotMDS(y, method="bcv")

#First we estimate the common dispersion to get an idea 
#of the overall degree of inter-library
#variability in the data:
y<- estimateCommonDisp(y, verbose=TRUE)
#Disp = 0.02821 , BCV = 0.1679 

#The biological coefficient of variation (BCV) is the square root of the common dispersion.
#Generally it is important to allow tag-specific dispersion estimates, so we go on to compute
#empirical Bayes moderated tagwise dispersion estimates. The trend is turned off as it
#is not usually required for SAGE data:

#This function implements the empirical Bayes strategy proposed by Robinson and Smyth (2007) for 
#estimating the tagwise negative binomial dispersions. 
#The prior values for the dispersions are determined by a global trend. The individual tagwise 
#dispersions are then squeezed towards this trend. The prior degrees of freedom determines the 
#weight given to the prior. The larger the prior degrees of freedom, the more the tagwise 
#dispersions are squeezed towards the global trend. If the number of libraries is large, the 
#prior becomes less important and the tagwise dispersion are determined more by the individual 
#tagwise data.

#If trend="none", then the prior dispersion is just a constant, the common dispersion. Otherwise, 
#the trend is determined by a moving average (trend="movingave") or loess smoother applied to the 
#tagwise conditional log-likelihood. method="loess" applies a loess curve of degree 0 as 
#implemented in loessByCol.

y<- estimateTagwiseDisp(y, trend="none")
#y1<-estimateTagwiseDisp(y, trend="movingave") ##  usually used for RNA-seq data
#The following plot displays the estimates:
plotBCV(y)

## Plotting mean vs variance
tag.mean<-rowMeans(y$counts)
#tag.mean<-apply(y$counts, 1, mean, na.rm=T)
tag.var<-apply(y$counts, 1, var, na.rm=T)

#remove genes with no variation
ind<-tag.var == 0
tag.mean<-tag.mean[!ind]
tag.var<-tag.var[!ind]

plot(log10(tag.mean), log10(tag.var), xlab = "Gene Means, log10", ylab = "Gene Variance, log10", pch=18)
abline(0,1, col="red", lwd=2)
lines(lowess(log10(tag.mean), log10(tag.var), f =0.01), col="blue", lwd=2) 
lines(smooth.spline(log10(tag.mean), log10(tag.var), df=20), col="purple", lwd=2) 

## Differential Analysis using exact conditional test
et <- exactTest(y, pair=c("WT","DCLK"))

## Top 10 genes (tags)
topTags(et)

#The following table shows the individual counts per million for the top ten tags. edgeR
#chooses tags that both have large fold changes and are consistent between replicates:
detags <- rownames(topTags(et)$table)
cpm(y)[detags, order(y$samples$group)]

#The total number of differentiallly expressed genes at FDR< 0.05:
de <- decideTestsDGE(et, p=0.05)
summary(de)

#A smearplot displays the log-fold changes with the DE genes highlighted:
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h = c(-2, 2), col = "blue")

#The smear plot shows the relationship between the log FC and mean of normalized 
#counts. Grey points represent genes with non-significant changes in expression, 
#whereas red points represent genes that are significantly differentially expressed

## Outputting all results
results1<-topTags(et, n = 50051, adjust.method="BH",sort.by="p", p.value=1)
## Volcano plot of log(FC) vs -log10(pvalue)

plot(results1$table$logFC, -log10(results1$table$PValue), pch=18, xlab="log(FC)", ylab="-log10(p)")
abline(v=-2, col="blue")
abline(v=2, col="blue")
abline(h=10, col="red", lty=2)



####---------------------------------------------
#### glmfit Example 
####---------------------------------------------

#Edge R:  differential expression analysis of digital gene expression data

dir1<-"C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\Tuch_data.csv"
setwd(dir1)

#Commands to run edge R
#Opening R in the environment
#need to type in anything that has 
#source("http://bioconductor.org/biocLite.R")
#installing edgeR
#biocLite("edgeR") 

#Using Case Studies from edge R User's Guide (Chapter 4)
#Using edge R in the environment
library(edgeR)

rawdata<-read.csv("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\Tuch_data.csv", header=T, check.names=FALSE, stringsAsFactors=FALSE)
dim(rawdata)
head(rawdata)
y<-DGEList(counts=rawdata[,4:9], genes=rawdata[,1:3])
y
#The study by Tuch et al. [22] was undertaken a few years ago,
# so not all of the RefSeq IDs
#provided by match RefSeq IDs currently in use. 
#We retain only those transcripts with IDs
#in the current NCBI annotation, which is provided by the 
#org.HS.eg.db package:

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)

idfound <- y$genes$RefSeqID %in% mappedRkeys(org.Hs.egREFSEQ)
y <- y[idfound,]
table(idfound)
y

#We add Entrez Gene IDs to the annotation
egREFSEQ <- toTable(org.Hs.egREFSEQ)
head(egREFSEQ)

m <- match(y$genes$RefSeqID, egREFSEQ$accession)
y$genes$EntrezGene <- egREFSEQ$gene_id[m]

## Now use Entrez Gene IDs to find gene symbol
egSYMBOL <- toTable(org.Hs.egSYMBOL)
head(egSYMBOL)
m <- match(y$genes$EntrezGene, egSYMBOL$gene_id)
y$genes$Symbol <- egSYMBOL$symbol[m]
head(y$genes)

#Different RefSeq transcripts for the same gene symbol count #predominantly the same reads.
#So we keep one transcript for each gene symbol. We choose the #transcript with highest
#overall count:
y
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$Symbol)
y <- y[!d,]
nrow(y)
#[1] 10522

#Normally we would also filter lowly expressed genes. For this #data, all transcripts already
#have at least 50 reads for all samples of at least one of the #tissues types.

#Recompute the library sizes:
y$samples$lib.size <- colSums(y$counts)

#Use Entrez Gene IDs as row names:
rownames(y$counts) <- rownames(y$genes) <- y$genes$EntrezGene
y$genes$EntrezGene <- NULL

## normalization
#TMM normalization is applied to this dataset to account for #compositional difference between
#the libraries as opposed to "upperquartile".
y <- calcNormFactors(y)
y1<-calcNormFactors(y, method="upperquartile")

y$samples

## Data exploration
plotMDS(y)
y


## Design Matrix
Patient <- factor(c(8,8,33,33,51,51))
Tissue <- factor(c("N","T","N","T","N","T"))
data.frame(Sample=colnames(y),Patient,Tissue)
design <- model.matrix(~Patient+Tissue)
rownames(design) <- colnames(y)

##Estimating common dispersion followed by tagwise

y <- estimateGLMCommonDisp(y, design, verbose=TRUE)

#The square root of the common dispersion gives 
#the coefficient of variation of biological
#variation. Here the common dispersion is found to be 0.16, 
#so the coefficient of biological
#variation is around 0.4. Then we estimate gene-wise dispersion 
#estimates, allowing a possible trend with averge count size:

y <- estimateGLMTrendedDisp(y, design)
y1 <- estimateGLMTagwiseDisp(y, design)  

##tagwise method = This function implements the empirical Bayes 
#strategy proposed by McCarthy et al (2012) for estimating the 
#tagwise negative binomial dispersions. The experimental 
#conditions are specified by design matrix allowing for multiple 
#explanatory factors. The empirical Bayes posterior is 
#implemented as a conditional likelihood with tag-specific 
#weights, and the conditional likelihood is computed using Cox-
#Reid approximate conditional likelihood (Cox and Reid, 1987).
#The prior degrees of freedom determines the weight given to the 
#global dispersion trend. The larger the prior degrees of 
#freedom, the more the tagwise dispersions are squeezed towards 
#the global trend.

plot(y$trended.dispersion, y1$tagwise.dispersion)
abline(0,1, col="red")

#Now proceed to determine differentially expressed genes. Fit #genewise glms. 
# will be fitting model adjusting for patient (paired).
#If dispersion=NULL will be extracted from y, with order of #precedence: genewise dispersion, trended dispersions, common 
#dispersion.

fit <- glmFit(y, design) ## will use trend
#Conduct likelihood ratio tests for tumour vs normal tissue differences and show the top
#genes:
lrt <- glmLRT(fit)

topTags(lrt)

#data for top tags
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]

summary(de<-decideTestsDGE(lrt))

detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")



###-----------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------


### Example DE analysis using two methods
### (i) edgeR   --- Negative binomial distribution assumption
###(ii) metagenomeSeq  --- Zero inflated gaussian distribution
###-------------------------------------------------------------------------------------------------
library(edgeR)
library(RColorBrewer)

### Load Data
dat <- read.csv("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\edgeR_Example.csv",header=T,row.names=1)
class(dat)
head(dat)
colnames(dat)
# [1] "con1"  "con2"  "con3"  "case1" "case2" "case3"
head(dat)
dim(dat)

###---------------------------
### (i) edgeR 
###---------------------------
### Genes filltering based on counts
dat <- dat[rowSums(dat)!=0,]     ## Remove if gene counts is zero for all samples

### MDS (Multidimensional Scaling) plot:
plotMDS(dat, main="Multipledimensional scaling plot (MDS)")

### Mean Variance plot by gene
###-----------------------------------------------
mean.x <- apply(dat,1,mean)
var.x <- apply(dat,1,var) 
plot(log10(mean.x),log10(var.x),pch=20, xlab="Mean (in log 10 scale)",ylab="Variance (in log 10 scale)",
     xlim=c(0,9),ylim=c(0,9),main="Variance vs Mean")
abline(0,1,col=3,lwd=2)

### Barplot
###-----------------------------------------------
lib.size <- colSums(dat)   
barplot(lib.size,xaxt="n",xlab="Study samples",ylab="Library size values",main="Barplot of library size",
        col=c(rep("lightgreen",3),rep("lightcoral",3)),sub="(Red line represents mean)")
abline(h=mean(lib.size),lwd=2,col="red")
legend("topright",c("control","case"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	

### Boxplot (In log10 scale)
###-----------------------------------------------
boxplot(x=as.list(as.data.frame(log10(dat+1))),xlab="Samples",ylab="Values in log10 scale",main="Boxplot across study samples",
        col=c(rep("lightgreen",3),rep("lightcoral",3)))
legend("topright",c("control","case"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	

### Density plots by sample
###-----------------------------------------------
for (i in 1:ncol(dat)){
  if (i==1){
    plot(density(log10(dat[,1])),main="Density plot across study samples",xlab="Subjects",col="gray",ylim=c(0,0.6))}
  else {
    den <- density(log10(dat[,i]))
    lines(den$x,den$y,col="gray")}
}

## edgeR DE analysis 
##-------------------------------------------------------------------------------------
keep <- rowSums(cpm(as.matrix(dat))>1)>=2     
# At least 2 samples have to have cpm > 1.
dat.filtered <- dat[keep,]
dim(dat.filtered)
rm(keep)
d <- DGEList(counts=as.matrix(dat.filtered), lib.size=colSums(dat.filtered), group=c(rep("con",3),rep("case",3)))
dim(d)
d <- calcNormFactors(d, method="TMM")        ## Calculates normalization factors
d <- estimateDisp(d)                   		 ## Calculates genewise dispersion parameter adjusted using bayesian empirical method(Recommended, edgeR)

## BCV 
plotBCV(d)

de.test <- exactTest(d,pair=c("con","case")) ## First entry under pair is the baseline 
de.test.FDR <- topTags(de.test,n=Inf,adjust.method="BH", sort.by="PValue") 
head(de.test.FDR$table)
summary(de <- decideTestsDGE(de.test, p=0.01, adjust="BH"))  ## Counts up and down regulated genes	

# ### plotSmear
# ###---------------------------------------------------------
detags <- rownames(d)[as.logical(de)]
plotSmear(de.test, de.tags=detags,main="Smear Plot",sub="(The horizontal blue lines show 4-fold changes)",cex=1) 
abline(h = c(-2, 2), col = "blue")

### Volcano plot
###---------------------------------------------------------
d.volcano <-  de.test.FDR$table[,c("logFC","PValue")]
par(mar=c(5,4,4,5))
plot(d.volcano$logFC, -log(d.volcano$PValue,10), main="",pch=20, cex=2,xlab=expression(log[2]~fold~change), ylab=expression(-log[10]~pvalue),
     xlim=c(min(d.volcano$logFC)-1,max(d.volcano$logFC)+1))
title("Volcano plot")

# Log2 fold change and p-value cutoff
lfc <- 2
pval <- 0.01 
# Selecting interesting genes
sigGenes <- (abs(d.volcano$logFC)> lfc & -log(d.volcano$PValue,10) > -log10(pval))   
# Identifying the selected genes
points(d.volcano[sigGenes,]$logFC,-log(d.volcano[sigGenes,]$PValue,10),pch=20,col="red",cex=2)
abline(h=-log10(pval),col="green3",lty=2)
abline(v=c(-lfc,lfc),col="blue",lty=2)
mtext(paste("pval = ",round(pval,2)),side=4,at=-log10(pval),cex=0.8,line=0.5,las=1)
mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)

d.edgeR.result <- de.test.FDR$table


###---------------------------
library(metagenomeSeq)
###---------------------------
dat <- dat[rowSums(dat)!=0,]     ## Remove if gene counts is zero for all samples

## metagenomeSeq DE analysis
##-------------------------------------------------------------------------------------
d <- newMRexperiment(dat)   ### Creates MR experiment object 
p <- cumNormStatFast(d)     # percentile for cumumative sum (css) normalization factor, default is used here
d <- cumNorm(d, p = p)      # Cumulative sum scaling normalization
group.Indicator <- c("gr1","gr1","gr1","gr2","gr2","gr2")
mod <- model.matrix(~1 + group.Indicator)
DE.fit <- fitFeatureModel(d, mod)    	# Zero-inflated log normal model 
#DE.fit <- fitZig(d, mod)    			# Zero-inflated Gaussian mixture model
DEresult.metagenomeSeq <- MRtable(DE.fit,number=nrow(d),adjustMethod = "BH")   
DEresult.metagenomeSeq <- DEresult.metagenomeSeq[order(DEresult.metagenomeSeq$adjPvalues),]
head(DEresult.metagenomeSeq)

### Volcano plot
###---------------------------------------------------------
d.volcano <-  DEresult.metagenomeSeq[,c("logFC","pvalues")]
d.volcano <-  subset(d.volcano,!is.na(logFC))
d.volcano$pvalues[d.volcano$pvalues==0] <- .Machine$double.eps   # replace 0 by small number to show the points in plot
par(mar=c(5,4,4,5))
plot(d.volcano$logFC, -log(d.volcano$pvalues,10), main="",pch=20, cex=2,xlab=expression(log[2]~fold~change), ylab=expression(-log[10]~pvalue),
     #		xlim=c(-6,6))
     xlim=c(min(d.volcano$logFC)-1,max(d.volcano$logFC)+1))
title("Volcano plot")
# Log2 fold change and p-value cutoff
lfc <- 2
pval <- 0.01 
# Selecting interesting genes
sigGenes <- (abs(d.volcano$logFC)> lfc & -log(d.volcano$pvalues,10) > -log10(pval))   
# Indicating the selected genes
points(d.volcano[sigGenes,]$logFC,-log(d.volcano[sigGenes,]$pvalues,10),pch=20,col="red",cex=2)
abline(h=-log10(pval),col="green3",lty=2)
abline(v=c(-lfc,lfc),col="blue",lty=2)
mtext(paste("pval = ",round(pval,2)),side=4,at=-log10(pval),cex=0.8,line=0.5,las=1)
mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)

d.meta.result <- DEresult.metagenomeSeq


### Common DE genes identified by both methods
### Venn diagram
library(limma)
d1 <- d.edgeR.result
d2 <- d.meta.result
d1$edgeR <- ifelse(d1$FDR<0.05,1,0)
d2$metagenomeSeq <- ifelse(d2$adjPvalues<0.05,1,0)
d2 <- d2[rownames(d1),]          # Order as edgeR results 
all(rownames(d1)==rownames(d2))
dd <- cbind(d1$edgeR,d2$metagenomeSeq)  ## Note: dd has been ordered as edgeR results
dd[is.na(dd)] <- 0
rownames(dd) <- rownames(d1)
colnames(dd) <- c("edgeR","metagenomeSeq")
vennDiagram(vennCounts(dd),circle.col=c("blue","red"),counts.col=NULL,main="")
