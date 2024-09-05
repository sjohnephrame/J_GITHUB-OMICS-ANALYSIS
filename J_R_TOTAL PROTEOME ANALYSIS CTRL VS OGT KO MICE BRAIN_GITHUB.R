#CTRL VS OGT KO MICE BRAINS #WT-OGT_KO-MouseBrain-Mecano-10OutOf50-Fr1-24-Combined

# Clear environment
rm(list = ls())

# Attach librarues
library(readxl)
library(tidyverse)
library(limma)
library(EnvStats)   # used to get geometric means
library(missForest) # Imputation
library(biomaRt) # library for mapping between annotations
library(qvalue)

# Load helper functions
path.to.functions = "C:\\Users\\sophi\\OneDrive\\Desktop\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_TOTAL PROTEOME_R\\J_PROTEOMICS R CODE\\functions2.R"
source(path.to.functions)
extract_string = function(x, k, pos) unlist(lapply(strsplit(x, k), function(y) y[pos]))

#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
#code explanation:

#The extract_string function takes a vector of strings (x), splits each
#string using a specified delimiter (k), and extracts the element at a specified position (pos) 
#from each resulting vector of substrings. The final result is a vector containing the extracted elements.

#The overall result is a list containing the elements extracted from each string at the specified position. 

#function(y) y[pos]: This is an anonymous function (a function without a name) that
#takes a vector y and extracts the element at position pos from that vector.

#lapply - This applies the anonymous function to each element of the list obtained
#from strsplit(x, k). It returns a list where each element is the extracted element
#at position pos from the corresponding vector of substrings.

#strsplit(x, k): This function splits each string in vector x using the specified 
#delimiter k. It returns a list where each element is a vector of substrings obtained by splitting the corresponding string.

#unlist(...): This converts the list of extracted elements back into a single vector.


#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------

#=========================#
# Proteomics Preprocessing ----
#=========================#

dat1k = read_xlsx("C:\\Users\\sophi\\OneDrive\\Desktop\\J_Dr. Slawson projects_ 2023\\J_DR DIAS OGT KO MICE\\WT-OGT_KO-MouseBrain-Mecano-10OutOf50-Fr1-24-Combined (1) - HIGH AND MED FDR.xlsx")
colnames(dat1k)
data2k<-dat1k

colnames(dat1k)
head(dat1k)
#dat1<-dat1[,-(20:23)]
colnames(dat1k)
names(dat1k)[names(dat1k) == 'Gene Symbol'] <- 'Symbol'
colnames(dat1k)
head(dat1k)
data2newk<-dat1k
library(dplyr)
dat1k <- dat1k %>%
  relocate(Symbol, .after = Accession)
colnames(dat1k)
head(dat1k)
data2new1k<-dat1k


dim(dat1k)
#[1] 8476   11 before removing duplicates 

new_dat1<-dat1k
head(new_dat1)
new_dat1<-new_dat1[,-(1)] #Remove FDR
head(new_dat1)
new_dat1<-new_dat1[,-(3)] #Remove Coverage
head(new_dat1)
new_dat1<-new_dat1[,-(9)] #Remove modfn
head(new_dat1) 

dim(new_dat1) #[1] 8476    8#without accession


#add_symbols - adds symbols to refseq id and removes duplicates
head(new_dat1)
new_dat2<-new_dat1
new_dat3<-new_dat1
new_dat3<-new_dat3[,-(2)] #Remove modfn
head(new_dat3)
new_dat2 = add_symbols(new_dat2) # This returns unique gene symbols by numbering any duplicates. Be careful!
new_dat3 = add_symbols(new_dat3) # This returns unique gene symbols by numbering any duplicates. Be careful!
# Duplicates are made unique by appending '.2', '.3', etc. 
# sum(grepl("\\.", dat$Symbol))
head(new_dat1)
head(new_dat2)
head(new_dat3)
new_dat4<-new_dat3
#remove duplicates
d <- duplicated(new_dat4$Symbol)
head(d)
new_dat4 <- new_dat4[!d,]
nrow(new_dat4)
dim(new_dat4)
#[1] 7868   8 after removing duplicates

#but we need like matrix with no symbols

head(new_dat4)
new_dat4<-as.data.frame(new_dat4)
head(new_dat4)
rownames(new_dat4) <- new_dat4$Accession #convert column value to row names
head(new_dat4)
new_dat4.1<-new_dat4
head(new_dat4.1)
new_dat4<-new_dat4[,-(1)] #Remove Accession
head(new_dat4)

# This might look weird, but there will be NaN for some values that had no data
# and this will change them to NA, which is easier to handle.
new_dat4[is.na(new_dat4)] <- NA
dim(new_dat4) #[1] 7868    7 #without accession
new_dat4.1[is.na(new_dat4.1)] <- NA
dim(new_dat4.1) #[1] 7868   8 #with accession

# Assess missingness of proteins

head(new_dat4)
dim(new_dat4) #[1] 7868 7 #without accession, with symbol
new_dat5<-new_dat4[,-1] #without accession, without symbol
head(new_dat4) #with symbols
head(new_dat5) #without symbols
dim(new_dat4) #with symbols [1] 7868   7
dim(new_dat5) #without symbols [1] 7868   6

head(new_dat4) #with symbols
table(apply(new_dat4, 1, function(x) sum(is.na(x)))) #with symbols
# 0    6 
#7171  697 

head(new_dat5) #without symbols
table(apply(new_dat5, 1, function(x) sum(is.na(x)))) #with symbols
#0    6 
#7171  697 

#-------------------
#Code explanation:
#-------------------

#apply(new_dat1, 1, function(x) sum(is.na(x))):

#new_dat1: Your data frame.
#1: Indicates that the function should be applied to rows. 2 for cols
#function(x) sum(is.na(x)): A function that calculates the number of NA values in each row.

#----------------------------------------------------------------------------------------------


#thus, there is no missing value with symbol
#so ideally good to map symbols here before proceeding further

#For now let us omit NA with symbols

# Reduce to proteins with no missing data.

dim(new_dat4) #7868 rows 7 cols in new_dat4 #with symbols
head(new_dat4)
dat_comp1 = na.omit(new_dat4) #with symbols and no missing values
dim(dat_comp1) #with symbols and no missing values
#[1] 7171    7 #19 becauseof Symbol
head(new_dat4)
dim(new_dat4) #7868 rows 7 cols in new_dat1 #with symbols, without omit
dim(dat_comp1) #[1] 7171   7 #7 because of symbol, with omit
mapping1<-dat_comp1 
head(mapping1)
dim(mapping1)  #[1] 7171 7 #7 because of symbol, with omit
#mapping1<-mapping1[,-(2:19)]

#without symbols
head(new_dat5)
dim(new_dat5)
dat_comp2 = na.omit(new_dat5)
dim(dat_comp2)
#[1] 7171 6 #without symbols

#It is same with/without symbols

# log2 transform - to reduce skewness of a measurement variable 
#to make data more symmetrical, which helps it meet the assumptions of statistical models

head(dat_comp2)
dim(dat_comp2) #[1] 7171    6 after removing symbol, without NA


head(new_dat5)
dim(new_dat5) #[1] 7868    6 after remving symbols # with NAs

head(dat_comp2)
dim(dat_comp2) #[1] 7171    6 after removing symbol # without NAs
dat_comp_mat1 = as.matrix(log2(dat_comp2)) # without NAs #log2 transformation
head(dat_comp_mat1)
dim(dat_comp_mat1) #[1] 7171   6 after removing symbol # without NAs #log2 transformation

head(new_dat5)
dim(new_dat5)# with NAs
dat_mat1 = as.matrix(log2(new_dat5)) # with NAs
head(dat_mat1)
dim(dat_mat1) #[1] 7868    6 with NAs_log
dim(dat_comp_mat1) #[1] 7171   6  without NAs_log
head(dat_comp_mat1)
head(dat_mat1)
dim(dat_mat1) #[1] 7868    6 with NAs_log
dim(dat_comp_mat1) #[1] 7171   6  without NAs_log

#density plot - plot is applied to log10 transformed data in code 

head(dat_comp2)#not log transformed, not normalized, omitted NA
dim(dat_comp2)#[1] 7171 6 without NAs_ not log, not normalized
for (i in 1:ncol(dat_comp2)){
  if (i==1){
    plot(density(log10(dat_comp2[,1])),main="Density plot across study samples",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density(log10(dat_comp2[,i]))
    lines(den$x,den$y,col="blue")}
}


# Quantile normalize - to remove any technical variation before differential expression analysis

head(dat_comp_mat1)
dim(dat_comp_mat1) #[1] 7171   6  without NAs_log
dat_comp_mat_quant1 = limma::normalizeQuantiles(dat_comp_mat1) # without NAs, log
head(dat_comp_mat_quant1)
dim(dat_comp_mat_quant1) #7171    6

head(dat_mat1)
dim(dat_mat1) #7868    6
dat_mat_quant1 = limma::normalizeQuantiles(dat_mat1) # with NAs, log
head(dat_mat_quant1)
dim(dat_mat_quant1)# 7868    6 with NAs, log, normalize #9206   18

head(dat_comp_mat_quant1) # without NAs, log, normalize 
dim(dat_comp_mat_quant1) #7171    6 # without NAs, log, normalize #8334   18

head(dat_comp2)
dim(dat_comp2) # 7171    6 without NAs, not log,
dat_comp1_norm = limma::normalizeQuantiles(dat_comp2) # without NAs, not log, normalized
dim(dat_comp1_norm) #7171   6 # without NAs, not log, normalized
head(dat_comp1_norm)

#https://www.youtube.com/watch?v=ecjN6Xpv6SE&t=279s


## Batch Effect Correction ##
# Here is where one should run diagnostics for batch effects (i.e., systematic differences in measurements due to technical factors rather than biological signal), if the experimental design introduces batch effects

#for erk data, do a box plot
# boxplot: intensities of all 16 channels after data preprocessing and normalihttp://127.0.0.1:40129/graphics/17b44cfd-ded5-4327-9e49-56e1f2f75d97.pngzation
head(dat_comp_mat_quant1)
dim(dat_comp_mat_quant1) # without NAs, log, normalize #7171 6
par(mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
boxplot(dat_comp_mat_quant1, main="Boxplot normalized Intensities")


#=================================#
# Differential Expression Analysis ----
#     using limma-trend - # This assumes uniform and normal distribution of data (seen in density plot)
#=================================#

#When the library sizes are quite variable between samples, then the voom approach 
#is theoretically more powerful than limma-trend. 

#Library size could mean one of two things: the total number of reads that were sequenced 
#in the run or the total number of mapped reads

#Read = sequenced fragment of cDNA (obtained from RNA) here ptn

dat1 = dat_comp_mat_quant1
head(dat1)
dim(dat1) #[1] 7171   6 #removed NA, without Symbol, norm, log

# Change column names to not contain "-", since this will conflict with naming syntax in contrast matrix

colnames(dat1)
head(dat1)

# Getting a data.frame with 'Accession' and 'Symbol' columns

head(mapping1)
head(dat1)
dim(dat1) #[1] 7171 6 Without symbol
dim(mapping1) #[1] 7171 7 With Symbol
head(mapping1)
#both dat1 and mapping1 match
#so, try to separate symbol and accession for mapping1

mapping1.1<-mapping1
head(mapping1.1)
dim(mapping1.1) #[1] 7171   7
head(mapping1) #[1] 7171   7
dim(mapping1)
#make dim(mapping1) 8334 2
mapping2<-as.data.frame(cbind(rownames(mapping1),mapping1$Symbol))
dim(mapping2) #7171 2
head(mapping2)
colnames(mapping2) <- c('Accession','Symbol')
colnames(mapping2)
head(mapping2)
dim(mapping2) #[1] 7171    2

#----------------------------------------------------------------------------------------
#Code to add entrez gene id - need for AMEND
#----------------------------------------------------------------------------------------

#library(org.Hs.eg.db)
#dat4<-dat1
#Accession<-rownames(dat4)
#dim(dat4)

#idfound <-  Accession %in% mappedRkeys(org.Hs.egUNIPROT)
#y <- dat4[idfound,]
#table(idfound)
#y

#We add Entrez Gene IDs to the annotation
#egUNIPROT <- toTable(org.Hs.egUNIPROT)
#head(egREFSEQ)

#m <- match(y$genes$RefSeqID, egREFSEQ$accession)
#y$genes$EntrezGene <- egREFSEQ$gene_id[m]

## Now use Entrez Gene IDs to find gene symbol
#egSYMBOL <- toTable(org.Hs.egSYMBOL)
#head(egSYMBOL)
#m <- match(y$genes$EntrezGene, egSYMBOL$gene_id)
#y$genes$Symbol <- egSYMBOL$symbol[m]
#head(y$genes)

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

# Design matrix
head(dat1)
colnames(dat1) <- c('CTRL_1','CTRL_2','CTRL_3','OGTKO_1','OGTKO_2','OGTKO_3')
head(dat1)
txa = extract_string(colnames(dat1), "_", 1)
txa
Xa = model.matrix(~0 + txa)
Xa
colnames(Xa) = sort(unique(txa))
Xa #design matrix

# Contrast matrix
colnames(Xa)
#contrast.matrix <- makeContrasts(group2-group1, group3-group2, group3-group1, levels=design)
contrastsa = c("CTRL-OGTKO")
contrastsa
Xa
cm_data = makeContrasts(contrasts=contrastsa, levels=colnames(Xa))
cm_data #contrast matrix  #Note: Check with Dr. Thompson if it is correct!!!!!!

#Check out this recent paper by the authors which is precisely focused on explaining design matrices.

#Basically, limma (and other software) are focused on using linear models to analyse your data, 
#because these tools let you represent/model your study, taking into account any variables which are
#of interest to you.

#For example, you may want to compare expression of gene X between young and old samples. 
#A t-test will simply compare the means of the groups. However, sex may influence 
#the expression of gene X. A linear model lets you take into account sex when comparing
#the expression between young and old groups.

#Having said this, design matrices are the way that you indicate the study 
#variables to model (e.g., age and sex), and contrast matrices are the way that you 
#indicate which variables you want to test (e.g. age) between which groups.
#https://www.biostars.org/p/9554210/

#Linear models describe a continuous response variable as a function of one or more predictor variables. 
#eg linear regression


# Fitting models. The trend=TRUE argument indicates that the mean-variance trend will be accounted for in eBayes (limma-trend method)

head(mapping2)
head(dat1)
Xa
cm_data
dea = top_table_results(mapping2, dat1, Xa, cm_data, annotate = TRUE, trend = TRUE)
head(dea)
head(dea$tt) #8 columns #log FC 
gene.symbola = extract_string(dea$tt$Symbol, "\\.", 1)
head(gene.symbola)
head(dea$tt$Symbol)
# both match
tmpa = aggregate(dea$tt[,-c(7,8)], by = list(Symbol = gene.symbola), FUN = median) #remove last 2 columns
head(dea$tt)
head(tmpa)
tmpa$Accession = dea$tt$Accession[match(tmpa$Symbol, gene.symbola)] # Get the first listed Accession ID for each gene symbol
head(tmpa)
head(dea$tt)
dea$tt = tmpa
head(dea$tt)
head(tmpa)

##################################################################################################################
######################################################################################################################
#logFCprot - if you need an excel sheet of the differentially expressed values - For heatmap
head(dea$tt)
logFCprotKO<-dea$tt
head(logFCprotKO)
logFCprotKO<-as.data.frame(logFCprotKO)
head(logFCprotKO)
type(logFCprotKO)
library("writexl") 
#install.packages("writexl")
#write_xlsx(logFCprotKO,"C:\\Users\\sophi\\OneDrive\\Desktop\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_TOTAL PROTEOME_R\\logFCprotKO.xlsx")
dim(logFCprotKO)

# Fitting models - Simplified code - The trend=TRUE argument indicates that the mean-variance trend will be accounted for in eBayes (limma-trend method)
#de = top_table_results(mapping, dat, X, cm_dat, annotate = TRUE, trend = TRUE)
#gene.symbol = extract_string(de$tt$Symbol, "\\.", 1)
#tmp = aggregate(de$tt[,-c(11,12)], by = list(Symbol = gene.symbol), FUN = median)
#tmp$Accession = de$tt$Accession[match(tmp$Symbol, gene.symbol)] # Get the first listed Accession ID for each gene symbol
#de$tt = tmp

#head(de$tt)

##################################################################################################################
######################################################################################################################

#MDS - Multi-dimensional scaling plot - similar to PCA
head(dea$tt) #values after comparisons with contrast matrix
head(dat1) #normailzed and filtered values - use this for mds plot
dim(dat1)
labels<-colnames(dat1)
labels
library(ggplot2)
col <- c("blue","blue","blue","orange","orange","orange")
col
#https://support.bioconductor.org/p/100859/
plotMDS(dat1, main="Multipledimensional scaling plot (MDS)",col=col)


#https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html

######################################################################################

# PCA - Principal Component Analysis

labels<-c("CTRL","CTRL","CTRL","OGTKO","OGTKO","OGTKO")
labels
get_pca
get_pca1 <- function(dat, labels, legend_title = 'Treatment') {
  require(ggplot2)
  
  pca_res <- prcomp(t(dat), center = T, scale. = F)
  
  U <- data.frame(pca_res$x)
  
  p <- ggplot(data = U, aes(x = PC1, y = PC2, color = labels )) +
    geom_point(size = 3, alpha = 0.5) + 
    theme_bw() +
    labs(color = legend_title, title="Principal Component Analysis (PCA)" )+
    scale_color_manual(values = c("blue","orange"))+
    theme(plot.title = element_text(hjust=0.5))
  
  return(list(p = p , summ = summary(pca_res)))
} # end get_pca

head(dat1)  
PCA<-get_pca1(dat1,labels)
PCA$p
#https://stackoverflow.com/questions/33640492/change-point-colors-and-color-of-frame-ellipse-around-points
#https://www.geeksforgeeks.org/how-to-change-position-of-ggplot-title-in-r/

##########################################################################################

#do a box plot
# boxplot: intensities of all 16 channels after data preprocessing and normalihttp://127.0.0.1:40129/graphics/17b44cfd-ded5-4327-9e49-56e1f2f75d97.pngzation

col <- c("blue","blue","blue","orange","orange","orange")
col
fill <- c("blue","orange")
fill
#par(mar = c(3, 3, 3, 10), mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2, xpd = TRUE)
boxplot(dat_comp_mat_quant1, main="Boxplot normalized Intensities",col=col)
#legend("right", title="Treatment",legend=c("GFP-0","GFP-10","OGTKD-605-0","OGTKD-605-10","OGTKD-606-0","OGTKD-606-10"), fill=fill, inset = c(-0.155,0))

#inset depends on margin of plot window

##################################################################################################################################

# Volcano plots for p-value and adjusted p-values

colnames(dea$tt) #logFC and pvalues
head(dea$tt)

#CTRL-OGTKO

rx <- c(-1, 1)*max(abs(dea$tt$logFC))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.CTRL.OGTKO), -log10(dea$tt$adj.P.Val.CTRL.OGTKO))))

#https://biostatsquid.com/volcano-plot/

#par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
#par(las=1, xaxs="i", yaxs="i")


#adjusted p value
plot(dea$tt$logFC, -log10(dea$tt$adj.P.Val.CTRL.OGTKO), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Volcano plot of adjusted p-values")

# only adj p-value cutoff
lfc <- 0
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
# Selecting interesting genes
sigGenes1 <- ((dea$tt$logFC)> lfc & -log(dea$tt$adj.P.Val.CTRL.OGTKO,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$logFC)< (-lfc) & -log(dea$tt$adj.P.Val.CTRL.OGTKO,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$logFC,-log(dea$tt[sigGenes1,]$adj.P.Val.CTRL.OGTKO,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$logFC,-log(dea$tt[sigGenes2,]$adj.P.Val.CTRL.OGTKO,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
#abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)

#############################################################################################################################################

#Prep for GSEA

dea
dea2<-dea
dea$tt
head(dea2)
dea2$tt
head(dea2$tt)

#CTRL VS OGT KO
updownOGTTENGFPTEN<-cbind(dea2$tt$logFC,dea2$tt$adj.P.Val.CTRL.OGTKO,dea2$tt$Symbol)
head(updownOGTTENGFPTEN)
head(dea2$tt)
updownOGTTENGFPTEN<-as.data.frame(updownOGTTENGFPTEN)
colnames(updownOGTTENGFPTEN)<-c('logFC','adj.P.Val.CTRL.OGTKO','Symbol')
head(updownOGTTENGFPTEN)
rownames(updownOGTTENGFPTEN)<-updownOGTTENGFPTEN$Symbol
head(updownOGTTENGFPTEN)
updownOGTTENGFPTEN<-updownOGTTENGFPTEN[,-3]
head(updownOGTTENGFPTEN)
updownOGTTENGFPTEN<-as.data.frame(updownOGTTENGFPTEN)
head(updownOGTTENGFPTEN)
#convert e to nmeric so can sort
updownOGTTENGFPTEN$adj.P.Val.CTRL.OGTKO <- as.numeric(formatC(updownOGTTENGFPTEN$adj.P.Val.CTRL.OGTKO, format = "e", digits = 2))
head(updownOGTTENGFPTEN)
sortupdownOGTTENGFPTEN<-as.data.frame(updownOGTTENGFPTEN[order(updownOGTTENGFPTEN$adj.P.Val.CTRL.OGTKO,decreasing = TRUE),]) 
head(sortupdownOGTTENGFPTEN)


###################################################################################################################
###################################################################################################################

#NOT NEEDED FOR GSEA BUT CAN DO FOR KEGG AND GO
#remove genes not significant p>0.05 
remsortupdownOGTTENGFPTEN <- sortupdownOGTTENGFPTEN[!(sortupdownOGTTENGFPTEN$adj.P.Val.CTRL.OGTKO>0.05), ]
head(remsortupdownOGTTENGFPTEN)
remsortupdownOGTTENGFPTEN<-as.data.frame(remsortupdownOGTTENGFPTEN)
head(remsortupdownOGTTENGFPTEN)
dim(remsortupdownOGTTENGFPTEN)
#463 proteins < p=0.05 and FDR Filtration 

# Filter log FC >1.5 and <-1.5 separately
head(remsortupdownOGTTENGFPTEN)

logFCless<-remsortupdownOGTTENGFPTEN
head(logFCless)
logFCless<-as.data.frame(logFCless)
head(logFCless)
type(logFCless)
#logFCless<-as.numeric(logFCless$logFC.OGTTEN.GFPTEN)
#head(logFCless)
logFCless <- transform(logFCless, logFC = as.numeric(logFC))
head(logFCless)
type(logFCless)
#logFCless<-subset(logFCless,(((logFCless$logFC.OGTTEN.GFPTEN)<(-1.5))&((logFCless$logFC.OGTTEN.GFPTEN)>(1.5))))
logFCless<-subset(logFCless,(((logFCless$logFC)<(-1.5))))
dim(logFCless)
#4,2

logFCgreat<-remsortupdownOGTTENGFPTEN
head(logFCgreat)
logFCgreat<-as.data.frame(logFCgreat)
head(logFCgreat)
type(logFCgreat)
#logFCgreat<-as.numeric(logFCgreat$logFC.OGTTEN.GFPTEN)
#head(logFCgreat)
logFCgreat <- transform(logFCgreat, logFC = as.numeric(logFC))
head(logFCgreat)
type(logFCgreat)
#logFCgreat<-subset(logFCgreat,(((logFCgreat$logFC.OGTTEN.GFPTEN)<(-1.5))&((logFCgreat$logFC.OGTTEN.GFPTEN)>(1.5))))
logFCgreat<-subset(logFCgreat,(((logFCgreat$logFC)>(1.5))))


head(logFCgreat)
dim(logFCgreat) #0 2

#>0.5 and <-0.5 combined proteins

logFCgreat<-remsortupdownOGTTENGFPTEN
head(logFCgreat)
logFCgreat<-as.data.frame(logFCgreat)
head(logFCgreat)
type(logFCgreat)
#logFCgreat<-as.numeric(logFCgreat$logFC.OGTTEN.GFPTEN)
#head(logFCgreat)
logFCgreat <- transform(logFCgreat, logFC = as.numeric(logFC))
head(logFCgreat)
type(logFCgreat)
#logFCgreat<-subset(logFCgreat,(((logFCgreat$logFC.OGTTEN.GFPTEN)<(-1.5))&((logFCgreat$logFC.OGTTEN.GFPTEN)>(1.5))))
#logFCgreat<-subset(logFCgreat,(((logFCgreat$logFC.OGTTEN.GFPTEN)>(1.5))))
logFCgreat<-as.data.frame(logFCgreat)
head(logFCgreat)
type(logFCgreat)
logFCgreat<-subset(logFCgreat,logFCgreat$logFC>0.5 | logFCgreat$logFC<(-0.5))
head(logFCgreat)
#logFCgreat <- logFCgreat[logFCgreat$logFC.OGTTEN.GFPTEN > 1.5 | logFCgreat$logFC.OGTTEN.GFPTEN < -1.5, ]
head(logFCgreat)
dim(logFCgreat) #78  2 combined >0.5 and <-0.5

###################################################################################################################
###################################################################################################################

#Prepare protlist for waterfall plot - rank-ordered protein list

#GSEA 606 PROT LIST
head(sortupdownOGTTENGFPTEN)
prot.list <- sortupdownOGTTENGFPTEN$logFC             # rank-ordered protein list
head(prot.list)
prot.list<-as.numeric(prot.list)
head(prot.list)
head(sortupdownOGTTENGFPTEN)
names(prot.list) <- rownames(sortupdownOGTTENGFPTEN)
head(prot.list)

# Waterfall plot - rank-ordered protein list

# Barplot of ranked fold changes (waterfall plot)
#barplot(sort(gene.list, decreasing = T),axisnames=FALSE,main="Plot of ranked gene Fold changes")

head(prot.list)
col <- ifelse(prot.list >0, "orange", "blue")
barplot(sort(prot.list, decreasing = T),axisnames=FALSE,main="Plot of ranked Fold changes of total proteome - OGT - 606 10 min vs. GFP - 10 min",  col=col, border=col, ylim=c(-4,4))

###############################################################################################################################################

###GO analysis

#goana: Gene Ontology or KEGG Pathway Analysis
#In limma: Linear Models for Microarray Data

#topGO package provides tools for testing GO terms while accounting for the topology of the GO graph. Different test statistics and different methods for eliminating local similarities and dependencies between GO terms can be implemented and applied.

#In R, topKEGG is a function that extracts the most significant KEGG pathways from kegga output. 
#kegga function to gather pathway enrichment for my dataset.

library(org.Mm.eg.db)
library(edgeR)
library(GO.db)

head(remsortupdownOGTTENGFPTEN)
d.go1 <- remsortupdownOGTTENGFPTEN
head(d.go1) #subsetted already for p<0.05
#d.go.DE <- subset(d.go,PValue<0.05) 
d.go.DE1<-d.go1
head(d.go.DE1)
d.entrez.id1 <- mapIds(org.Mm.eg.db, keys=rownames(d.go.DE1),column="ENTREZID",keytype="SYMBOL")
length(d.entrez.id1)
head(d.entrez.id1)
all(rownames(d.go.DE1)==names(d.entrez.id1)) 
go.test1 <- goana(d.entrez.id1,species="Mm")
go.results1 <- topGO(go.test1, sort = "DE", number = Inf)
head(go.results1)
sum(go.results1$P.DE<10^(-5)) #633
sum(go.results1$P.DE<0.05) #3488

head(remsortupdownOGTTENGFPTEN)
head(logFCgreat)
dim(logFCgreat) #78  2
d.go2 <- logFCgreat
head(d.go2) #subsetted already for p<0.05, FDR, log FC >1.5 and <-1.5
#d.go.DE <- subset(d.go,PValue<0.05) 
d.go.DE2<-d.go2
head(d.go.DE2)
dim(d.go.DE2)
d.entrez.id2 <- mapIds(org.Mm.eg.db, keys=rownames(d.go.DE2),column="ENTREZID",keytype="SYMBOL")
length(d.entrez.id2)
head(d.entrez.id2)
all(rownames(d.go.DE2)==names(d.entrez.id2)) 
go.test2 <- goana(d.entrez.id2,species="Mm")
go.results2 <- topGO(go.test2, sort = "DE", number = Inf)
head(go.results2)
sum(go.results2$P.DE<10^(-5)) #500
sum(go.results2$P.DE<0.05) #2006

###KEGG analysis

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("KEGGREST")
library(KEGGREST)
#kegg.test <- kegga(d.entrez.id,species="Mm")
#kegg.results <- topKEGG(kegg.test, sort = "DE", number = Inf)
#head(kegg.results)
#sum(kegg.results$P.DE<10^(-5))

d.kegg1 <- remsortupdownOGTTENGFPTEN
#d.kegg.DE <- subset(d.kegg,FDR<0.05) 
#d.kegg.DE <- subset(d.kegg,PValue<0.05)
d.kegg.DE1<-d.kegg1
head(d.kegg.DE1)
dim(d.kegg.DE1)
head(d.entrez.id1)
all(rownames(d.kegg.DE1)==names(d.entrez.id1)) 

kegg.test1 <- kegga(d.entrez.id1,species="Mm")
kegg.results1 <- topKEGG(kegg.test1, sort = "DE", number = Inf)
head(kegg.results1)
sum(kegg.results1$P.DE<10^(-5))#10
sum(kegg.results1$P.DE<0.05) #74

d.kegg2 <- logFCgreat
#d.kegg.DE <- subset(d.kegg,FDR<0.05) 
#d.kegg.DE <- subset(d.kegg,PValue<0.05)
d.kegg.DE2<-d.kegg2
head(d.kegg.DE2)
head(d.entrez.id2)
all(rownames(d.kegg.DE2)==names(d.entrez.id2)) 

kegg.test2 <- kegga(d.entrez.id2,species="Mm")
kegg.results2 <- topKEGG(kegg.test2, sort = "DE", number = Inf)
head(kegg.results2)
sum(kegg.results2$P.DE<10^(-5)) #4
sum(kegg.results2$P.DE<0.05) #39



######################################################################################################################
######################################################################################################################

### J_GSEA analysis for OGT 606 10 MIN VS GFP 10 MIN

# Load All gene sets file downloaded from Broad Institute 
# The following website contains the gene set collection or the complete Molecular Signatures Database (MSigDB)  
# http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
#  
library(fgsea)

#all.gene.sets <- gmtPathways("C:\\Users\\sophi\\OneDrive\\Desktop\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\msigdb.v2023.2.Hs.symbols.gmt") #got lot of warnings
#all.gene.sets <- gmtPathways("C:\\Users\\sophi\\OneDrive\\Desktop\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\msigdb.v7.4.symbols.gmt")
#all.gene.sets <- gmtPathways("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\m2.all.v2023.1.Mm.symbols.gmt")        #m2curated gene sets
all.gene.sets <- gmtPathways("C:\\Users\\sophi\\OneDrive\\Desktop\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\J_CLASS DOWNLOADS ALL\\m2.all.v2023.1.Mm.symbols.gmt")        #m2curated gene sets

class(all.gene.sets)
length(all.gene.sets)
all.gene.sets[1:2]
# Show first a few pathways, and within those, show only the first few genes. 
library(tidyverse)
all.gene.sets %>% head() %>% lapply(head)

head(prot.list)
#head(d.entrez.id1)
#d.entrez.id2<-
### Now run fgsea 
fgseaRes <- fgsea(pathways = all.gene.sets, stats = prot.list, minSize=15, maxSize=500,eps=0)
warnings()
head(fgseaRes)
head(fgseaRes[order(pval), ])
sum(fgseaRes[, padj < 0.05])#258
fgseaRes1<-fgseaRes[order(padj), ]
head(fgseaRes1)


# Make a table plot for a bunch of selected pathways:
topPathwaysUp <- fgseaRes[NES > 0][head(order(padj), n=10), pathway]
topPathwaysDown <- fgseaRes[NES < 0][head(order(padj), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
p<-plotGseaTable(all.gene.sets[topPathways], prot.list, fgseaRes,gseaParam = 0.5, pathwayLabelStyle=list(size=7),
                 valueStyle = list(size=7),
                 axisLabelStyle = list(size=4))
p

#ggsave(p, width=10, height=8, file="table.png")

plotGseaTable(all.gene.sets[topPathwaysUp], prot.list, fgseaRes,gseaParam = 0.5, pathwayLabelStyle=list(size=7),
              valueStyle = list(size=7),
              axisLabelStyle = list(size=6))

plotGseaTable(all.gene.sets[rev(topPathwaysDown)], prot.list, fgseaRes,gseaParam = 0.5, pathwayLabelStyle=list(size=7),
              valueStyle = list(size=7),
              axisLabelStyle = list(size=6))

head(topPathways)

library(ggplot2)
topPathways
# Make a few Enrichment Plots

#Top 4 up

#Plot1
plotEnrichment(all.gene.sets[["REACTOME_CHOLESTEROL_BIOSYNTHESIS"]],prot.list) + labs(title="REACTOME_CHOLESTEROL_BIOSYNTHESIS")
#Plot2
plotEnrichment(all.gene.sets[["WP_ELECTRON_TRANSPORT_CHAIN"]],prot.list) + labs(title="WP_ELECTRON_TRANSPORT_CHAIN")
#Plot3
#plotEnrichment(all.gene.sets[["KEGG_OXIDATIVE_PHOSPHORYLATION"]],prot.list) + labs(title="KEGG_OXIDATIVE_PHOSPHORYLATION")
#Plot4
#plotEnrichment(all.gene.sets[["WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA"]],prot.list) + labs(title="WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA")

#Top 4 down

#Plot1
plotEnrichment(all.gene.sets[["MARKEY_RB1_ACUTE_LOF_DN"]],prot.list) + labs(title="MARKEY_RB1_ACUTE_LOF_DN")
#gseaplot(all.gene.sets[["REACTOME_EUKARYOTIC_TRANSLATION_INITIATION"]],prot.list) + labs(title="REACTOME_EUKARYOTIC_TRANSLATION_INITIATION")

#Plot2
plotEnrichment(all.gene.sets[["CHEN_METABOLIC_SYNDROM_NETWORK"]],prot.list) + labs(title="CHEN_METABOLIC_SYNDROM_NETWORK")
#Plot3
plotEnrichment(all.gene.sets[["LEE_AGING_CEREBELLUM_UP"]],prot.list) + labs(title="LEE_AGING_CEREBELLUM_UP")
#Plot4
plotEnrichment(all.gene.sets[["LIU_OVARIAN_CANCER_TUMORS_AND_XENOGRAFTS_XDGS_DN"]],prot.list) + labs(title="LIU_OVARIAN_CANCER_TUMORS_AND_XENOGRAFTS_XDGS_DN")
#Plot5
plotEnrichment(all.gene.sets[["ICHIBA_GRAFT_VERSUS_HOST_DISEASE_35D_UP"]],prot.list) + labs(title="ICHIBA_GRAFT_VERSUS_HOST_DISEASE_35D_UP")

#----------------------------------------------------------------------------------------------------

#Pretty gsea visualizations    

#library(org.Hs.eg.db)
#install.packages("msigdbr")
library(msigdbr)
library(clusterProfiler)
head(prot.list)

head(prot.list)
type(prot.list)


#GSEA 606 PROT LIST


library(enrichplot)

head(fgseaRes)

topPathwaysUp1 <- fgseaRes[ES > 0][head(order(padj), n=10)]
head(topPathwaysUp1)
topPathwaysUp1<-topPathwaysUp1[,-2]
topPathwaysUp1
topPathwaysUp1<-topPathwaysUp1[,-3]
topPathwaysUp1
topPathwaysUp1<-topPathwaysUp1[,-3]
topPathwaysUp1
topPathwaysUp1<-topPathwaysUp1[,-(4:5)]
topPathwaysUp1
topPathwaysUp1<-as.data.frame(topPathwaysUp1)
topPathwaysUp1

topPathwaysDown1 <- fgseaRes[ES < 0][head(order(padj), n=10)]
topPathwaysDown1<-topPathwaysDown1[,-2]
topPathwaysDown1
topPathwaysDown1<-topPathwaysDown1[,-3]
topPathwaysDown1
topPathwaysDown1<-topPathwaysDown1[,-3]
topPathwaysDown1
topPathwaysDown1<-topPathwaysDown1[,-(4:5)]
topPathwaysDown1
topPathwaysDown1<-as.data.frame(topPathwaysDown1)
topPathwaysDown1

topPathways1 <- rbind(topPathwaysUp1, rev(topPathwaysDown1)) #merge rows
topPathways1
topPathways1<-as.data.frame(topPathways1)
topPathways1

#https://stephenturner.github.io/deseq-to-fgsea/#using_the_fgsea_package
ggplot(topPathways1, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05&NES>0)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()+
  theme(text = element_text(face = "bold"))

#Change legend

