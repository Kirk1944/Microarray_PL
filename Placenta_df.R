library(affy)
library(pvclust)
library(vsn)
library(gplots)
library(limma)
library(Biobase)
library(sva)
library(affyPLM)
library(arrayQualityMetrics)
library(ALLMLL)
library(dplyr)
library(hgu133plus2hsentrezg.db)
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezgcdf)

hgu133plus2hsentrezg()


setwd("D:/R/RStudio/Working directory/SCI_PROJECT")
getwd()

#Anotation
gns <- select(hgu133plus2hsentrezg.db,  probes, c("SYMBOL","ENTREZID","GENENAME"))
gnlst <- tapply(1:nrow(gns), gns$PROBEID, function(x) gns[x,])
gnlst <- lapply(gnlst, function(x) apply(x, 2, function(y) paste(unique(y), collapse = " | ")))
gns <- do.call("rbind", gnlst)
g <- as.data.frame(gns)
g[which(grepl("NA",g$ENTREZID)),]$ENTREZID <- NA
fullList <- fullList[complete.cases(fullList),]
probes <- row.names(fit2$coefficients)


#Loading samples
pd_9984 = read.AnnotatedDataFrame(filename="pdata_9984.txt")
pd_12767 = read.AnnotatedDataFrame(filename="pdata_12767.txt")
pd_18809 = read.AnnotatedDataFrame(filename="pdata_18809.txt")

celnames_9984 <- paste(pd_9984@data$geo_accession,".CEL",sep = '')
celnames_12767 <- paste(pd_12767@data$geo_accession,".CEL",sep = '')
celnames_18809 <- paste(pd_18809@data$geo_accession,".CEL",sep = '')

eset_9984 = ReadAffy(phenoData=pd_9984, sampleNames=pd_9984$geo_accession,filenames = celnames_9984 )
eset_12767 = ReadAffy(phenoData=pd_12767, sampleNames=pd_12767$geo_accession,filenames = celnames_12767)
eset_18809 = ReadAffy(phenoData=pd_18809, sampleNames=pd_18809$geo_accession,filenames = celnames_18809)

#cdfNAME from Brainarray
eset_12767@cdfName<-"hgu133plus2hsentrezg"
eset_9984@cdfName<-"hgu133plus2hsentrezg"
eset_18809@cdfName<-"hgu133plus2hsentrezg"

#RMA
n_affyData9984 <- rma(eset_9984) 
n_affyData12767 <- rma(eset_12767) 
n_affyData18809 <- rma(eset_18809)

#Combat batch effect removing
allGenes <- combine(n_affyData9984, n_affyData12767, n_affyData18809)
batch <- allGenes@phenoData@data$GSE
modcombat = model.matrix(~as.factor(Trimester), data = allGenes)
combat_edata = ComBat(dat=exprs(allGenes), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

#Annotation & diff expression analysis
types = factor(c(rep("First",4),rep("Second",4),rep("Term",4),rep("First",8),rep("Preterm",5),rep("Term",5)))
types
design = model.matrix(~ 0+types)
colnames(design) = levels(types)
fit = lmFit(combat_edata,design)
fit$genes <- g
contrast.matrix = makeContrasts(First-Second, levels=design)
fit2 = contrasts.fit(fit,contrast.matrix)
fit2 = eBayes(fit2)
fullList = topTable(fit2, number=nrow(fit2), sort.by="logFC")

df_genes <- subset(x = fullList, abs(logFC)>1.5)

diff_reg_F_S_1_5 <- df_genes
diff_reg_F_T_1_5 <- df_genes
diff_reg_F_PT_1_5 <- df_genes
diff_reg_S_T_1_5 <- df_genes
diff_reg_S_PT_1_5 <- df_genes
diff_reg_T_PT_1_5 <- df_genes


#Output 
write((diff_reg_F_S_1_5$SYMBOL), "diff_reg_F_S_LogFC_1_5.txt", sep="\n")
write((diff_reg_F_T_1_5$SYMBOL), "diff_reg_F_T_LogFC_1_5.txt", sep="\n")
write((diff_reg_F_PT_1_5$SYMBOL), "diff_reg_F_PT_LogFC_1_5.txt", sep="\n")
write((diff_reg_S_T_1_5$SYMBOL), "diff_reg_S_T_LogFC_1_5.txt", sep="\n")
write((diff_reg_S_PT_1_5$SYMBOL), "diff_reg_S_PT_LogFC_1_5.txt", sep="\n")
write((diff_reg_T_PT_1_5$SYMBOL), "diff_reg_T_PT_LogFC_1_5.txt", sep="\n")


#Junk
df_genes$memlev <- abs(df_genes$logFC/max(df_genes$logFC))
write.table(df_genes[,c(2,11)], file = "allGenes_F_S_all.txt", row.names = FALSE, col.names = FALSE,sep = ",",quote = FALSE)
