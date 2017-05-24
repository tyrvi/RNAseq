library(pheatmap)
library(reshape)
library(gplots)
library(ops)
library(calibrate)
library(biomaRt)
library(sva)
library(ggplot2)
library(org.Hs.eg.db) # for transferring gene identifiers
library(data.table) # for collapsing transcript RPKMs

gtex <- read.delim("GTEx-3-Tissues-N.txt", sep="\t")
gtexcolumns <- colnames(gtex)
gtex.fpkms <- gtex[,gtexcolumns[2:length(gtexcolumns)]]
colnames(gtex.fpkms)[1] <- "Gene"

tcga <- read.delim("TCGA-3-Tissues-N.txt", sep="\t")
tcgacolumns <- colnames(tcga)
tcga.genes <- tcga[,"Gene"]
tcga[,"Gene"] <- gsub("\\|.*", "", tcga.genes)
tcga.fpkms <- tcga

fpkms <- merge(gtex.fpkms, tcga.fpkms, by="Gene")

temp <- fpkms[!duplicated(fpkms[,1]),]
gene <- temp[,1]
f <- temp[,2:ncol(temp)]
rownames(f) <- gene


data <- read.delim("combined-DataMatrix-QuantNormed.txt", sep = "\t")
datacol <- colnames(data)
# remove breast samples
data[,datacol[grep(".*breast", colnames(data))]] <- list(NULL)
data.genes <- data[,"Gene"]

qn <- data[,2:length(data)]
rownames(qn) <- data.genes

f.nozero <- f[-which(rowMeans(f[,])<=0.01),]
qn.nozero <- qn

f.gtex.bladder <- length(grep("GTEX.*bladder", colnames(f.nozero)))
f.gtex.prostate <- length(grep("GTEX.*prostate", colnames(f.nozero)))
f.gtex.thyroid <- length(grep("GTEX.*thyroid", colnames(f.nozero)))
f.tcga.bladder <- length(grep("TCGA.*bladder", colnames(f.nozero)))
f.tcga.prostate <- length(grep("TCGA.*prostate", colnames(f.nozero)))
f.tcga.thyroid <- length(grep("TCGA.*thyroid", colnames(f.nozero)))
f.gtex.total <- f.gtex.bladder + f.gtex.prostate + f.gtex.thyroid
f.tcga.total <- f.tcga.bladder + f.tcga.prostate + f.tcga.thyroid

qn.gtex.bladder <- length(grep("GTEX.*bladder", colnames(qn.nozero)))
qn.gtex.prostate <- length(grep("GTEX.*prostate", colnames(qn.nozero)))
qn.gtex.thyroid <- length(grep("GTEX.*thyroid", colnames(qn.nozero)))
qn.tcga.bladder <- length(grep("TCGA.*bladder", colnames(qn.nozero)))
qn.tcga.prostate <- length(grep("TCGA.*prostate", colnames(qn.nozero)))
qn.tcga.thyroid <- length(grep("TCGA.*thyroid", colnames(qn.nozero)))
qn.gtex.total <- qn.gtex.bladder + qn.gtex.prostate + qn.gtex.thyroid
qn.tcga.total <- qn.tcga.bladder + qn.tcga.prostate + qn.tcga.thyroid

plot.pca.published <- function(df,x,y,z,l,counts){
  
  p <- prcomp(t(df))
  
  v <- 100*(p$sdev)^2 / sum(p$sdev^2)
  v.x <- v[x]
  v.y <- v[y]
  
  colors <- c(rep("indianred",counts[3]), rep("dodgerblue",counts[4]), rep("forestgreen",counts[5]),
              rep("indianred",counts[6]), rep("dodgerblue",counts[7]), rep("forestgreen",counts[8]))
  
  shapes <- c(rep(16,counts[1]),rep(18,counts[2]))
  
  plot(p$x[,x],p$x[,y],pch=shapes,cex=1.5,col=colors,xlab=paste(paste("PC",x),round(v.x),"% of variance"),ylab=paste(paste("PC",y),round(v.y),"% of variance"),main=paste(z," FPKM \n n=",l))
  
}

# vectors of counts of samples
f.counts <- c(f.gtex.total, f.tcga.total, f.gtex.bladder, f.gtex.prostate, f.gtex.thyroid, f.tcga.bladder, f.tcga.prostate, f.tcga.thyroid)
qn.counts <- c(qn.gtex.total, qn.tcga.total, qn.gtex.bladder, qn.gtex.prostate, qn.gtex.thyroid, qn.tcga.bladder, qn.tcga.prostate, qn.tcga.thyroid)


# fig 1
pheatmap(cor(f.nozero), method="spearman", main="fig. 1b", show_colnames=FALSE, show_rownames=FALSE)

plot.pca.published(f.nozero, 1, 2, "fig. 1c", length(f.nozero[,1]), f.counts)
legend("bottomleft",legend=c("Bladder","Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(16,18),ncol=2, bty="n")

# log transform
pseudo <- 1
f.log <- log2(f.nozero + pseudo)

# fig 2
pheatmap(cor(f.log), method="spearman", main="fig. 2a log2", show_rownames=FALSE, show_colnames=FALSE)

plot.pca.published(f.log, 1, 2, "fig. 2b log2", length(f.nozero[,1]), f.counts)
legend("bottomleft",legend=c("Bladder","Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(16,18),ncol=2, bty="n")

plot.pca.published(f.log, 2, 3, "fig. 2c log2", length(f.nozero[,1]), f.counts)
legend("bottomleft",legend=c("Bladder","Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(16,18),ncol=2, bty="n")

meta <- data.frame(study=c(rep("GTEx",f.counts[1]), rep("TCGA",f.counts[2])), tissue=c(rep("Bladder",f.counts[3]), rep("Prostate",f.counts[4]), rep("Thyroid",f.counts[5]), rep("Bladder",f.counts[6]), rep("Prostate",f.counts[7]), rep("Thyroid",f.counts[8])))
batch <- meta$study
design <- model.matrix(~1, data=meta)
combat <- ComBat(dat=f.log, batch=batch, mod=design, par.prior=TRUE)

# fig 3
pheatmap(cor(combat), method="spearman", main="fig. 3a COMBAT", show_rownames=FALSE, show_colnames=FALSE)

plot.pca.published(combat, 1, 2, "fig. 3b COMBAT", length(f.nozero[,1]), f.counts)
legend("bottomleft",legend=c("Bladder","Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(16,18),ncol=2, bty="n")

# Quantile Normalized
# fig 1
pheatmap(cor(qn.nozero), method="spearman", main="fig. 1b Quantile Normalized", show_colnames=FALSE, show_rownames=FALSE)

plot.pca.published(qn.nozero, 1, 2, "fig. 1c Quantile Normalized", length(qn.nozero[,1]), qn.counts)
legend("bottomleft",legend=c("Bladder","Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(16,18),ncol=2, bty="n")

# log transform
pseudo <- 1
qn.log <- log2(qn.nozero + pseudo)

# fig 2
pheatmap(cor(qn.log), method="spearman", main="fig. 2a log2 Quantile Normalized", show_rownames=FALSE, show_colnames=FALSE)

plot.pca.published(qn.log, 1, 2, "fig. 2b log2 Quantile Normalized", length(qn.nozero[,1]), qn.counts)
legend("bottomleft",legend=c("Bladder","Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(16,18),ncol=2, bty="n")

plot.pca.published(qn.log, 2, 3, "fig. 2c log2 Quantile Normalized", length(qn.nozero[,1]), qn.counts)
legend("bottomleft",legend=c("Bladder","Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(16,18),ncol=2, bty="n")

# remove rows where the row variance is 0 otherwise combat throws "Error in while (change > conv) { : missing value where TRUE/FALSE needed"
# see https://stackoverflow.com/questions/21532998/error-when-using-combat
#qn.combat <- qn.log
qn.combat <- qn.log[-which(rowVars(qn.log[,])==0),]

meta <- data.frame(study=c(rep("GTEx",qn.counts[1]), rep("TCGA",qn.counts[2])), tissue=c(rep("Bladder",qn.counts[3]), rep("Prostate",qn.counts[4]), rep("Thyroid",qn.counts[5]), rep("Bladder",qn.counts[6]), rep("Prostate",qn.counts[7]), rep("Thyroid",qn.counts[8])))
batch <- meta$study
design <- model.matrix(~1, data=meta)
combat <- ComBat(dat=qn.combat, batch=batch, mod=design, par.prior=TRUE)

# fig 3
pheatmap(cor(combat), method="spearman", main="fig. 3a COMBAT Quantile Normalized", show_rownames=FALSE, show_colnames=FALSE)

plot.pca.published(combat, 1, 2, "fig. 3b COMBAT Quantile Normalized", length(qn.combat[,1]), qn.counts)
legend("bottomleft",legend=c("Bladder","Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(16,18),ncol=2, bty="n")
