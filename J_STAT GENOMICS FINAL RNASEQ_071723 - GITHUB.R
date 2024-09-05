
#--------------------------------------------------------------------------------
#TEST USING EXACT METHOD
#--------------------------------------------------------------------------------


#Gene Set Enrichment analysis using fgsea package in R
# https://bioconductor.org/packages/release/bioc/html/fgsea.html
# http://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
# https://github.com/ctlab/fgsea

####---------------------------------------------
#### fgsea  
####---------------------------------------------
rm(list=ls())

#Read in file
rnaseq <- read.csv("C:\\Users\\sophi\\OneDrive\\Desktop\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\RNASeq2.txt", header=T, sep="\t")

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("fgsea")
library(data.table)
library(fgsea)
library(ggplot2)
library(edgeR)
library(RColorBrewer)


#mydata

#rnaseq <- read.csv("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\RNASeq2.txt", header=T, sep="\t")
head(rnaseq)
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


###---------------------------
### DE analysis using edgeR - This was used for final data
###---------------------------

TMG <- TMG[rowSums(TMG)!=0,]     ## Remove if gene counts is zero for all samples
dim(TMG)
#[1] 25705     6

#Plots before normalization:

### MDS (Multidimensional Scaling) plot:
plotMDS(TMG, main="Multipledimensional scaling plot (MDS)")


### Mean Variance plot by gene
###-----------------------------------------------

mean.x <- apply(TMG,1,mean)
var.x <- apply(TMG,1,var) 
plot(log10(mean.x),log10(var.x),pch=20, xlab="Mean (in log 10 scale)",ylab="Variance (in log 10 scale)",
     xlim=c(0,9),ylim=c(0,9),main="Variance vs Mean for TMG data")
abline(0,1,col=3,lwd=2)



### Barplot
###-----------------------------------------------

lib.size <- colSums(TMG)   
barplot(lib.size,xaxt="n",xlab="Study samples",ylab="Library size values",main="Barplot of library size",
        col=c(rep("lightgreen",3),rep("lightcoral",3)),sub="(Red line represents mean)")
abline(h=mean(lib.size),lwd=2,col="red")
legend("topright",c("SALINE","TMG"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	


### Boxplot (In log10 scale)
###-----------------------------------------------

boxplot(x=as.list(as.data.frame(log10(TMG+1))),xlab="Samples",ylab="Values in log10 scale",main="Boxplot across study samples",
        col=c(rep("lightgreen",3),rep("lightcoral",3)))
legend("topright",c("SALINE","TMG"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	


### Density plots by sample
###-----------------------------------------------


for (i in 1:ncol(TMG)){
  if (i==1){
    plot(density(log10(TMG[,1])),main="Density plot across study samples for TMG data",xlab="Subjects",col="gray",ylim=c(0,0.8))}
  else {
    den <- density(log10(TMG[,i]))
    lines(den$x,den$y,col="gray")}
}


## edgeR DE analysis 
##-------------------------------------------------------------------------------------

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
plotBCV(d1)

de.test1 <- exactTest(d1,pair=c("SALINE","TMG")) ## First entry under pair is the baseline 
de.test1.FDR <- topTags(de.test1,n=Inf,adjust.method="BH", sort.by="PValue") 
head(de.test1.FDR$table)
summary(de1 <- decideTestsDGE(de.test1, p=0.05, adjust="BH"))  ## Counts up and down regulated genes	
dim(d1)

#Plots after TMM normalization

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

boxplot(x=as.list(as.data.frame(log10(d2+1))),xlab="Samples",ylab="Values in log10 scale",main="Boxplot across study samples",
        col=c(rep("lightgreen",3),rep("lightcoral",3)))
legend("topright",c("SALINE","TMG"),lty=c(1,1),lwd=c(4,4),col=c("lightgreen","red"),bty="o",cex=1,text.font=2)	


### Density plots by sample
###-----------------------------------------------

for (i in 1:ncol(d2)){
  if (i==1){
    plot(density(log10(d2[,1])),main="Density plot across study samples for TMG data",xlab="Subjects",col="gray",ylim=c(0,0.8))}
  else {
    den <- density(log10(d2[,i]))
    lines(den$x,den$y,col="gray")}
}


# ### plotSmear
# ###---------------------------------------------------------

detags1 <- rownames(d1)[as.logical(de1)]
plotSmear(de.test1, de.tags=detags1,main="Smear Plot",sub="(The horizontal blue lines show 4-fold changes)",cex=1) 
abline(h = c(-2, 2), col = "blue")


### Volcano plot
###---------------------------------------------------------
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


d1.edgeR.result <- de.test1.FDR$table
head(d1.edgeR.result)



# Map Ensembl gene IDs to symbol. First create a mapping table.
#library(biomaRt)
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#ensg.id <- rownames(d.edgeR.result)
#d.gene.id <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values=ensg.id, mart= mart)
#dim(d.gene.id)
#head(d.gene.id)
#d.gene.id <- subset(d.gene.id,hgnc_symbol!="")  # Remove the ENSG ids that do not have hgnc_symbol 
#rownames(d.gene.id) <- d.gene.id$ensembl_gene_id
#dim(d.gene.id)

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

common.genes <- intersect(rownames(d1.edgeR.result),rownames(d.gene.id))
length(common.genes)
length(unique(common.genes)) #12734

# Merge gene id data with the results data

d.tmp1 <- d.gene.id[common.genes,]
d.tmp2 <- d1.edgeR.result[common.genes,]   
all(rownames(d.tmp1)==rownames(d.tmp2))

# TRUE

d.tmp <- cbind(d.tmp1,d.tmp2)
head(d.tmp)
d.tmp <- d.tmp[order(d.tmp$PValue),] ## Sort the data in the order of significance
head(d.tmp)


# Prepare the gene list

gene.list <- d.tmp$logFC             # rank-ordered gene list
names(gene.list) <- d.tmp$mgi_symbol

# Barplot of ranked fold changes

barplot(sort(gene.list, decreasing = T),axisnames=FALSE,main="Plot of ranked gene Fold changes")

###GO analysis
library(org.Hs.eg.db)
library(edgeR)
library(GO.db)
#d.go <- d.tmp
#d.go.DE <- subset(d.go,FDR<0.05) 
#d.go.DE <- subset(d.go,PValue<0.05) 
#d.entrez.id <- mapIds(org.Hs.eg.db, keys=d.go.DE$hgnc_symbol,column="ENTREZID",keytype="SYMBOL")
#length(d.entrez.id)
#head(d.entrez.id)
#all(d.go.DE$hgnc_symbol==names(d.entrez.id)) 
#go.test <- goana(d.entrez.id,species="Hs")
#go.results <- topGO(go.test, sort = "DE", number = Inf)
#head(go.results)
#sum(go.results$P.DE<10^(-5))


library(org.Mm.eg.db)
library(edgeR)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("GO.db")

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
sum(go.results$P.DE<10^(-5)) #196

###KEGG analysis
d.kegg <- d.tmp
#d.kegg.DE <- subset(d.kegg,FDR<0.05) 
d.kegg.DE <- subset(d.kegg,PValue<0.05)
#all(d.kegg.DE$hgnc_symbol==names(d.entrez.id)) 

#kegg.test <- kegga(d.entrez.id,species="Hs")
#kegg.results <- topKEGG(kegg.test, sort = "DE", number = Inf)
#head(kegg.results)
#sum(kegg.results$P.DE<10^(-5))


d.kegg <- d.tmp
#d.kegg.DE <- subset(d.kegg,FDR<0.05) 
d.kegg.DE <- subset(d.kegg,PValue<0.05)
all(d.kegg.DE$mgi_symbol==names(d.entrez.id)) 
d.entrez.id <- mapIds(org.Mm.eg.db, keys=d.go.DE$mgi_symbol,column="ENTREZID",keytype="SYMBOL")
length(d.entrez.id)
head(d.entrez.id)


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("KEGGREST")
library(KEGGREST)
kegg.test <- kegga(d.entrez.id,species="Mm")
kegg.results <- topKEGG(kegg.test, sort = "DE", number = Inf)
head(kegg.results)
sum(kegg.results$P.DE<10^(-5)) #1

### GSEA analysis
# Load All gene sets file downloaded from Broad Institute 
# The following website contains the gene set collection or the complete Molecular Signatures Database (MSigDB)  
# http://software.broadinstitute.org/gsea/downloads.jsp#msigdb

#Human  
#all.gene.sets <- gmtPathways("C:\\Users\\sophi\\OneDrive\\Desktop\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\J_CLASS DOWNLOADS ALL\\msigdb.v7.4.symbols.gmt")
#class(all.gene.sets)
#length(all.gene.sets)
#all.gene.sets[1:2]
# Show first a few pathways, and within those, show only the first few genes. 
#library(tidyverse)
#all.gene.sets %>% head() %>% lapply(head)

### Now run fgsea 
#fgseaRes <- fgsea(pathways = all.gene.sets,stats = gene.list,minSize=15,maxSize=500,eps=0)
#head(fgseaRes)
#head(fgseaRes[order(pval), ])
#sum(fgseaRes[, padj < 0.05])
# Make a few Enrichment Plots
#Plot1
#plotEnrichment(all.gene.sets[["GAO_LARGE_INTESTINE_ADULT_CI_MESENCHYMAL_CELLS"]],gene.list) + labs(title="GAO_LARGE_INTESTINE_ADULT_CI_MESENCHYMAL_CELLS")
#Plot2
#plotEnrichment(all.gene.sets[["CHEN_METABOLIC_SYNDROM_NETWORK"]],gene.list) + labs(title="CHEN_METABOLIC_SYNDROM_NETWORK")
#tail(fgseaRes[order(pval), ])
#Plot3
#plotEnrichment(all.gene.sets[["ZNF610_TARGET_GENES"]],gene.list) + labs(title="ZNF610_TARGET_GENES")
#Plot4
#plotEnrichment(all.gene.sets[["ZSCAN5DP_TARGET_GENES"]],gene.list) + labs(title="ZSCAN5DP_TARGET_GENES")
#fgseaRes[order(pval), ][15:20,]
#Plot5
#plotEnrichment(all.gene.sets[["GOBP_RESPONSE_TO_BIOTIC_STIMULUS"]],gene.list) + labs(title="GOBP_RESPONSE_TO_BIOTIC_STIMULUS")

# Make a table plot for a bunch of selected pathways:
#topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
#topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
#topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#plotGseaTable(all.gene.sets[topPathways], gene.list, fgseaRes,gseaParam = 0.5)



### GSEA analysis
# Load All gene sets file downloaded from Broad Institute 
# The following website contains the gene set collection or the complete Molecular Signatures Database (MSigDB)  
# http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
#  
#mydata - Mouse

all.gene.sets <- gmtPathways("C:\\Users\\sophi\\OneDrive\\Desktop\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\J_CLASS DOWNLOADS ALL\\m2.all.v2023.1.Mm.symbols.gmt")        #m2curated gene sets
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
sum(fgseaRes[, padj < 0.05]) #64
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

#mydata

#rnaseq <- read.csv("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\RNASeq2.txt", header=T, sep="\t")
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
plotBCV(d1)

de.test1 <- exactTest(d1,pair=c("SALINE","TMG")) ## First entry under pair is the baseline 
de.test1.FDR <- topTags(de.test1,n=Inf,adjust.method="BH", sort.by="PValue") 
head(de.test1.FDR$table)
summary(de1 <- decideTestsDGE(de.test1, p=0.05, adjust="BH"))  ## Counts up and down regulated genes	
dim(d1) # 12775     6


## 1. Hierarchical clustering method
##----------------------------------------------------------------------------------------------

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

######################################################################################################

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

######################################################################################################
#Not doing this bcoz not clinical data
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

######################################################################################################

dd
heatmap(dd, main = "SALINE Vs TMG", xlab = "Samples", ylab = "DEG")
#https://www.biostars.org/p/374551/


####################################################################################################


