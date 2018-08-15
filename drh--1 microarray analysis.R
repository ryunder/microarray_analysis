source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("GEOquery")
biocLite("celegans.db")
biocLite("affy")
biocLite("limma")
biocLite("gcrma")
biocLite("org.Ce.eg.db")
biocLite("annotate")

#Data: GSE41056 - sarkies et al 2013 microarray data

library(limma)
library(affy)
library(celegans.db)
library(Biobase)
library(GEOquery)
library(gcrma)
library(org.Ce.eg.db)
library(annotate)

#Import data from GEO
gset <- getGEOSuppFiles("GSE41056")

#Unpack downloaded data
untar("GSE41056/GSE41056_RAW.tar", exdir = "data")
cels <- list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep="/"), gunzip)
cels <- list.files("data/",pattern = "CEL")

#data

#Load and organize data
ab <- ReadAffy(filenames = cels, celfile.path = "data/")
eset <- rma(ab)
head(exprs(eset))
samplenames <- substring(colnames(eset), 12, nchar(colnames(eset))-14)
colnames(eset) <- samplenames

#MDS plot
plotMDS(eset)

#Naming eset - Probes to gene names
ID <- featureNames(eset)
gene_ids <- getSYMBOL(ID, "celegans")
fData(eset) <- data.frame(Symbol=gene_ids)
head(fData(eset))
dim(eset)
  
#Organize experiment factors
treatment <- as.factor(rep(c("un","inf","un","inf","un","inf") ,c(4,4,3,3,3,3)))
strain <- as.factor(rep(c("N2","JU1580","RDE1"), c(8,6,6)))               
TS <- paste(treatment, strain, sep = ".")
TS <- factor(TS, levels = c("un.N2","inf.N2","un.JU1580","inf.JU1580","un.RDE1","inf.RDE1"))

#Design and organize design matrix
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
design                                    

#Design contrast matrix
contrast_matrix <- makeContrasts(
  infN2.vs.unN2 = inf.N2 - un.N2,
  infJU1580.vs.unJU1580 = inf.JU1580 - un.JU1580,
  Diff = (inf.JU1580 - un.JU1580) - (inf.N2 - un.N2),
  levels = design)

#linear modeling and ebayes fit
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
efit <- eBayes(fit2, trend = T)
results <- decideTests(efit)
summary(results)

#Exploring data
vennDiagram(results[,1:2],circle.col = c("blue","red"), include = "up", 
            cex=0.67, show.include = T)

p.value=0.05
un.N2.vs.inf.N2 <- topTable(efit, coef = 1, n=Inf, p.value = p.value)
un.JU1580.vs.inf.JU1580 <- topTable(efit, coef = 2, n=Inf, p.value = p.value)
diff.JU1580.vs.N2 <- topTable(efit, coef=3, n=Inf, p.value = p.value)
dim(un.N2.vs.inf.N2)
dim(un.JU1580.vs.inf.JU1580)
dim(diff.JU1580.vs.N2)

head(un.N2.vs.inf.N2)
head(un.JU1580.vs.inf.JU1580)
head(diff.JU1580.vs.N2)

#To introduce a logFC cutoff, use `treat` instead of `eBayes`
tfit <- treat(fit2, lfc = log2(1.5), trend=T)
tresults <- decideTests(tfit)
summary(tresults)

vennDiagram(tresults[,1:2],circle.col = c("blue","red"), include = "up", 
            cex=0.67, show.include = T)

p.value=0.05
un.N2.vs.inf.N2.treat <- topTreat(tfit, coef = 1, n=Inf, p.value = p.value)
un.JU1580.vs.inf.JU1580.treat <- topTreat(tfit, coef = 2, n=Inf, p.value = p.value)
diff.JU1580.vs.N2.treat <- topTreat(tfit, coef=3, n=Inf, p.value = p.value)
dim(un.N2.vs.inf.N2.treat)
dim(un.JU1580.vs.inf.JU1580.treat)
dim(diff.JU1580.vs.N2.treat)

head(un.N2.vs.inf.N2.treat)
head(un.JU1580.vs.inf.JU1580.treat)
head(diff.JU1580.vs.N2.treat)

sessionInfo()

