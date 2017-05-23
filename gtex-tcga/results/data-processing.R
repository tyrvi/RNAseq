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

# Processes and write the GTEx data set to gtexfpkms.txt file
gtex <- read.delim("GTEx-3-Tissues-N.txt", sep="\t")
gtexcolumns <- colnames(gtex)

# random selection of sample is commented out in order to replicate graphs
random.gtex.bladder <- "GTEX.SNOS.0526.SM.4DM54.bladder." #gtexcolumns[sample(grep("bladder", gtexcolumns), size=1)]
random.gtex.prostate <- "GTEX.V1D1.1926.SM.4JBGX.prostate." #gtexcolumns[sample(grep("prostate", gtexcolumns), size=1)]
random.gtex.thyroid <- "GTEX.SIU7.1126.SM.2XCDW.thyroid." #gtexcolumns[sample(grep("thyroid", gtexcolumns), size=1)]

gtex.genes <- gtex[,"Description"]
gtex.fpkms <- data.frame(GENE=gtex.genes, GTEx_bladder=gtex[,c(random.gtex.bladder)], GTEx_prostate=gtex[,c(random.gtex.prostate)], GTEx_thyroid=gtex[,c(random.gtex.thyroid)])

# write to gtex file
#write.table(gtex.fpkms, file="gtexfpkms.txt", quote=FALSE, sep="\t")

tcga <- read.delim("TCGA-3-Tissues-N.txt", sep="\t")
tcgacolumns <- colnames(tcga)

# random selection of sample is commented out in order to replicate graphs
random.tcga.bladder <- "TCGA.K4.A54R.11A.11R.A26T.07.bladder." #tcgacolumns[sample(grep("bladder", tcgacolumns), size=1)]
random.tcga.prostate <- "TCGA.EJ.7794.11A.01R.2118.07.prostate." #tcgacolumns[sample(grep("prostate", tcgacolumns), size=1)]
random.tcga.thyroid <- "TCGA.KS.A41I.11A.11R.A23N.07.thyroid." #tcgacolumns[sample(grep("thyroid", tcgacolumns), size=1)]

tcga.genes <- tcga[,"Gene"]
tcga.genes.nobar <- gsub("\\|.*", "", tcga.genes)
tcga.fpkms <- data.frame(GENE=tcga.genes.nobar, TCGA_bladder=tcga[,c(random.tcga.bladder)], TCGA_prostate=tcga[,c(random.tcga.prostate)], TCGA_thyroid=tcga[,c(random.tcga.thyroid)])

#write.table(tcga.fpkms, file="tcgafpkms.txt", sep="\t")

# merge data based on gene
fpkms <- merge(gtex.fpkms, tcga.fpkms, by="GENE")

# remove duplicate genes
temp <- fpkms[!duplicated(fpkms[,1]),]
gene <- temp[,1]
f <- temp[,2:ncol(temp)]
rownames(f) <- gene

#write.table(f, file="combined-gtex-tcga-bygene.txt", quote=FALSE, sep="\t")

# remove samples where the FPKM is less than or equal to 0.01 in all samples
f.nozero <- f[-which(rowMeans(f[,])<=0.01),]

plot.pca.published <- function(df,x,y,z,l){
  
  p <- prcomp(t(df))
  
  v <- 100*(p$sdev)^2 / sum(p$sdev^2)
  v.x <- v[x]
  v.y <- v[y]
  
  colors <- c("indianred", "dodgerblue", "forestgreen",
              "indianred", "dodgerblue", "forestgreen")
  
  shapes <- c(rep(15,3),rep(16,3))
  
  plot(p$x[,x],p$x[,y],pch=shapes,cex=1.5,col=colors,xlab=paste(paste("PC",x),round(v.x),"% of variance"),ylab=paste(paste("PC",y),round(v.y),"% of variance"),main=paste(z," FPKM \n n=",l))
  
}

# figure 1b
pheatmap(cor(f.nozero), method="spearman", main="fig. 1a")

# figure 1c
plot.pca.published(f.nozero, 1, 2, "fig. 1b", length(f.nozero[,1]))
legend("bottomleft",legend=c("Bladder","Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(15,16,17,8),ncol=2, bty="n")

# figure 2 log transformation
pseudo <- 1
f.log <- log2(f.nozero + pseudo)
# figure 2a
pheatmap(cor(f.log), method="spearman", main="fig. 2a log2")
# figure 2b
plot.pca.published(f.log, 1, 2, "fig. 2b log2", length(f.nozero[,1]))
legend("bottomleft",legend=c("Bladder","Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(15,16,17,8),ncol=2, bty="n")
# figure 2c
plot.pca.published(f.log, 2, 3, "fig. 2c log2", length(f.nozero[,1]))
legend("bottomleft",legend=c("Bladder","Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(15,16,17,8),ncol=2, bty="n")

# figure 3
meta <- data.frame(study=c(rep("GTEx", 3), rep("TCGA", 3)), tissue=c("Bladder", "Prostate", "Thyroid", "Bladder", "Prostate", "Thyroid"))
batch <- meta$study
design <- model.matrix(~1, data=meta)
combat <- ComBat(dat=f.log, batch=batch, mod=design, par.prior=TRUE)
# figure 3a
pheatmap(cor(combat), method="spearman", main="fig. 3a COMBAT")
# figure 3b
plot.pca.published(combat, 1, 2, "fig. 3b COMBAT", length(f.nozero[,1]))
legend("bottomleft",legend=c("Bladder","Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(15,16,17,8),ncol=2, bty="n")
