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

data <- read.delim("combined-DataMatrix-QuantNormed.txt", sep = "\t")
datacol <- colnames(data)

# sample 1
# "GTEX.TMMY.1526.SM.4DXST.bladder."       "GTEX.T6MN.0726.SM.32PML.breast."        "GTEX.OIZH.2026.SM.3NB1M.prostate."      "GTEX.QXCU.0326.SM.2TC63.thyroid."      
# "TCGA.GD.A3OQ.11A.21R.A220.07.bladder."  "TCGA.BH.A0DL.11A.13R.A115.07.breast."   "TCGA.EJ.7789.11A.01R.2118.07.prostate." "TCGA.EL.A3GZ.11A.11R.A20F.07.thyroid."
# sample 2
# "GTEX.S4UY.0926.SM.4AD6O.bladder."       "GTEX.T5JC.2126.SM.32PMO.breast."        "GTEX.S33H.1826.SM.4AD65.prostate."      "GTEX.WHPG.0226.SM.3NMB9.thyroid."      
# "TCGA.BT.A20N.11A.11R.A14Y.07.bladder."  "TCGA.BH.A18R.11A.42R.A12D.07.breast."   "TCGA.EJ.7314.11A.01R.2118.07.prostate." "TCGA.BJ.A2N8.11A.11R.A18C.07.thyroid." 
random.gtex.bladder <- datacol[sample(grep("GTEX.*bladder", colnames(data)), size=1)]
random.gtex.breast <- datacol[sample(grep("GTEX.*breast", colnames(data)), size=1)]
random.gtex.prostate <- datacol[sample(grep("GTEX.*prostate", colnames(data)), size=1)]
random.gtex.thyroid <- datacol[sample(grep("GTEX.*thyroid", colnames(data)), size=1)]

random.tcga.bladder <- datacol[sample(grep("TCGA.*bladder", colnames(data)), size=1)]
random.tcga.breast <- datacol[sample(grep("TCGA.*breast", colnames(data)), size=1)]
random.tcga.prostate <- datacol[sample(grep("TCGA.*prostate", colnames(data)), size=1)]
random.tcga.thyroid <- datacol[sample(grep("TCGA.*thyroid", colnames(data)), size=1)]

#write.table(tcga.fpkms, file="tcgafpkms.txt", sep="\t")

gene <- data[,"Gene"]
gtex.bladder <- data[,random.gtex.bladder]
gtex.breast <- data[,random.gtex.breast]
gtex.prostate <- data[,random.gtex.prostate]
gtex.thyroid <- data[,random.gtex.thyroid]
tcga.bladder <- data[,random.tcga.bladder]
tcga.breast <- data[,random.tcga.breast]
tcga.prostate <- data[,random.tcga.prostate]
tcga.thyroid <- data[,random.tcga.thyroid]

f <- data.frame(GTEx_bladder=gtex.bladder, GTEx_breast=gtex.breast, GTEx_prostate=gtex.prostate, GTEx_thyroid=gtex.thyroid, TCGA_bladder=tcga.bladder, TCGA_breast=tcga.breast, TCGA_prostate=tcga.prostate, TCGA_thyroid=tcga.thyroid)
rownames(f) <- gene

#write.table(qn, file="", quote=FALSE, sep="\t")

# remove samples where the FPKM is less than or equal to 0.01 in all samples
f.nozero <- f[-which(rowMeans(f[,])<=0.01),]

plot.pca.published <- function(df,x,y,z,l){
  
  p <- prcomp(t(df))
  
  v <- 100*(p$sdev)^2 / sum(p$sdev^2)
  v.x <- v[x]
  v.y <- v[y]
  
  colors <- c("indianred", "dodgerblue", "forestgreen", "gold",
              "indianred", "dodgerblue", "forestgreen", "gold")
  
  shapes <- c(rep(15,4),rep(16,4))
  
  plot(p$x[,x],p$x[,y],pch=shapes,cex=1.5,col=colors,xlab=paste(paste("PC",x),round(v.x),"% of variance"),ylab=paste(paste("PC",y),round(v.y),"% of variance"),main=paste(z," FPKM \n n=",l))
  
}

# figure 1b
pheatmap(cor(f.nozero), method="spearman", main="fig. 1b Quantile Normalized")

# figure 1c
plot.pca.published(f.nozero, 1, 2, "fig. 1c Quantile Normalized", length(f.nozero[,1]))
legend("bottomleft",legend=c("Bladder","Breast", "Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen", "gold"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(15,16,17,8),ncol=2, bty="n")

# figure 2 log transformation
pseudo <- 1
f.log <- log2(f.nozero + pseudo)
# figure 2a
pheatmap(cor(f.log), method="spearman", main="fig. 2a Quantile Normalized log2")
# figure 2b
plot.pca.published(f.log, 1, 2, "fig. 2b Quantile Normalized log2", length(f.nozero[,1]))
legend("bottomleft",legend=c("Bladder","Breast", "Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen", "gold"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(15,16,17,8),ncol=2, bty="n")
# figure 2c
plot.pca.published(f.log, 2, 3, "fig. 2c Quantile Normalized log2", length(f.nozero[,1]))
legend("bottomleft",legend=c("Bladder","Breast", "Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen", "gold"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(15,16,17,8),ncol=2, bty="n")

# figure 3
meta <- data.frame(study=c(rep("GTEx", 4), rep("TCGA", 4)), tissue=c("Bladder", "Breast", "Prostate", "Thyroid", "Bladder", "Breast", "Prostate", "Thyroid"))
batch <- meta$study
design <- model.matrix(~1, data=meta)
combat <- ComBat(dat=f.log, batch=batch, mod=design, par.prior=TRUE)
# figure 3a
pheatmap(cor(combat), method="spearman", main="fig. 3a Quantile Normalized COMBAT")
# figure 3b
plot.pca.published(combat, 1, 2, "fig. 3b Quantile Normalized COMBAT", length(f.nozero[,1]))
legend("bottomleft",legend=c("Bladder","Breast", "Prostate","Thyroid"),col=c("indianred", "dodgerblue", "forestgreen", "gold"),cex=1.5,pch=20, bty="n")
legend("top",legend=c("GTEx","TCGA"),col="black",pch=c(15,16,17,8),ncol=2, bty="n")
