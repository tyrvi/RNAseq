# Processes and write the GTEx data set to gtexfpkms.txt file
gtex <- read.delim("GTEx-3-Tissues-N.txt", sep="\t")
gtexcolumns <- colnames(gtex)
# pick a random sample from samples of bladder, prostate and thyroid
random.gtex.bladder <- gtexcolumns[sample(grep("bladder", gtexcolumns), size=1)]
random.gtex.prostate <- gtexcolumns[sample(grep("prostate", gtexcolumns), size=1)]
random.gtex.thyroid <- gtexcolumns[sample(grep("bladder", gtexcolumns), size=1)]

gtex.genes <- gtex[,"Description"]
gtex.fpkms <- data.frame(GENE=gtex.genes, GTEx_bladder=gtex[,c(random.gtex.bladder)], GTEx_prostate=gtex[,c(random.gtex.prostate)], GTEx_thyroid=gtex[,c(random.gtex.thyroid)])

# write to gtex file
#write.table(gtex.fpkms, file="gtexfpkms.txt", quote=FALSE, sep="\t")

tcga <- read.delim("TCGA-3-Tissues-N.txt", sep="\t")
tcgacolumns <- colnames(tcga)

random.tcga.bladder <- tcgacolumns[sample(grep("bladder", tcgacolumns), size=1)]
random.tcga.prostate <- tcgacolumns[sample(grep("prostate", tcgacolumns), size=1)]
random.tcga.thyroid <- tcgacolumns[sample(grep("thyroid", tcgacolumns), size=1)]

tcga.genes <- tcga[,"Gene"]
tcga.genes.nobar <- gsub("\\|.*", "", tcga.genes)
tcga.fpkms <- data.frame(GENE=tcga.genes.nobar, TCGA_bladder=tcga[,c(random.tcga.bladder)], TCGA_prostate=tcga[,c(random.tcga.prostate)], TCGA_thyroid=tcga[,c(random.tcga.thyroid)])

#write.table(tcga.fpkms, file="tcgafpkms.txt", sep="\t")

fpkms <- merge(gtex.fpkms, tcga.fpkms, by="GENE")

write.table(fpkms, file="combined-gtex-tcga-bygene.txt", quote=FALSE, sep="\t")