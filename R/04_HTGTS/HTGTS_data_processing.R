library(dplyr)
library(grid)
library(fields)
library(stringr)
library(mapplots)
library(reshape2)

# CRISPR cut sites (based on hg19)
## SHANK2 -- chr11:70658290
## RARA -- chr17:38491451

# Bin count for modeling -- 250-Kbp bin count
args <- commandArgs(trailingOnly = TRUE)

hg19_coord <- read.csv("hg19_coord.txt", header=T, as.is=T, sep="\t")
binsize <- 250000
bindf <- as.data.frame(matrix(NA, ncol = 5, nrow = 0))
colnames(bindf) <- c("chr", "start", "end", "ct_control", "ct_e2")
for (i in 1:(nrow(hg19_coord)-1)){
  bindf[(nrow(bindf)+1):(nrow(bindf))+floor(hg19_coord$length[i]/binsize),] <- NA
  bindf$chr[is.na(bindf$chr)] <- hg19_coord$chr[i]
  bindf$start[is.na(bindf$start)] <- seq(1, floor(hg19_coord$length[i]/binsize)*binsize + 1, by=binsize)
  bindf$end[is.na(bindf$end)] <- bindf$start[is.na(bindf$end)] + binsize -1
}
condf <- read.csv(paste0("bedfiles/", args[1], "_h", args[2], "_withoutE2.hg19.bed"), header=F, as.is=T, sep="\t")
expdf <- read.csv(paste0("bedfiles/", args[1], "_h", args[2], "_withE2.hg19.bed"), header=F, as.is=T, sep="\t")
for (i in 1:nrow(bindf)){
  bindf$ct_control[i] <- length(unique(condf$V2[condf$V1 == paste0("chr", bindf$chr[i]) & condf$V2 >= bindf$start[i] & condf$V2 <= bindf$end[i]]))
  bindf$ct_e2[i] <- length(unique(expdf$V2[expdf$V1 == paste0("chr", bindf$chr[i]) & expdf$V2 >= bindf$start[i] & expdf$V2 <= bindf$end[i]]))
}
write.table(bindf, paste0("Counts.", args[1], ".", args[2], ".e2.", binsize/1000, "kbp.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

# Collecting data files 
## The output of this script "HTGTS.count.lograt.average.allcells.bait.marked.txt". This is available in the Data folder.
tempdf <- read.csv("htgts_intermediates/Counts.MCF7.SHANK2.e2.nonoverlapping.250kbp.txt", header=T, as.is=T, sep="\t")
colnames(tempdf)[4:5] <- c("mcf7.shank2.ctrl", "mcf7.shank2.e2")
df <- tempdf

tempdf <- read.csv("htgts_intermediates/Counts.T47D.SHANK2.e2.nonoverlapping.250kbp.txt", header=T, as.is=T, sep="\t")
colnames(tempdf)[4:5] <- c("t47d.shank2.ctrl", "t47d.shank2.e2")
df$t47d.shank2.ctrl <- tempdf$t47d.shank2.ctrl
df$t47d.shank2.e2 <- tempdf$t47d.shank2.e2

tempdf <- read.csv("htgts_intermediates/Counts.MCF7.RARA.e2.nonoverlapping.250kbp.txt", header=T, as.is=T, sep="\t")
colnames(tempdf)[4:5] <- c("mcf7.rara.ctrl", "mcf7.rara.e2")
df$mcf7.rara.ctrl <- tempdf$mcf7.rara.ctrl
df$mcf7.rara.e2 <- tempdf$mcf7.rara.e2

tempdf <- read.csv("htgts_intermediates/Counts.T47D.RARA.e2.nonoverlapping.250kbp.txt", header=T, as.is=T, sep="\t")
colnames(tempdf)[4:5] <- c("t47d.rara.ctrl", "t47d.rara.e2")
df$t47d.rara.ctrl <- tempdf$t47d.rara.ctrl
df$t47d.rara.e2 <- tempdf$t47d.rara.e2

df$inclusion <- "yes"
df$inclusion[df$chr=="11" & df$start <= 70658290 + 1000000 & df$end >= 70658290 - 1000000] <- "shank2_bait"
df$inclusion[df$chr=="17" & df$start <= 38491451 + 1000000 & df$end >= 38491451 - 1000000] <- "rara_bait"

df$msr <- log2((df$mcf7.shank2.e2+1)/(df$mcf7.shank2.ctrl+1))
df$tsr <- log2((df$t47d.shank2.e2+1)/(df$t47d.shank2.ctrl+1))
df$mrr <- log2((df$mcf7.rara.e2+1)/(df$mcf7.rara.ctrl+1))
df$trr <- log2((df$t47d.rara.e2+1)/(df$t47d.rara.ctrl+1))
df$msr[df$inclusion == "shank2_bait"] <- NA
df$tsr[df$inclusion == "shank2_bait"] <- NA
df$mrr[df$inclusion == "rara_bait"] <- NA
df$trr[df$inclusion == "rara_bait"] <- NA
df$average <- NA
for (i in 1:nrow(df)){
  df$average[i] <- mean(c(df$msr[i], df$tsr[i], df$mrr[i], df$trr[i]), na.rm=T)
}
write.table(df, paste0("HTGTS.count.lograt.average.allcells.bait.marked.txt"), row.names = F, col.names = T, quote = F, sep = "\t")


## Gene-based counting and ordering, then GSEA -- per individual experiment
args <- commandArgs(trailingOnly = TRUE)

marginlength <- 5000
bindf <- read.csv("genedf.GRCh37.all.genes.txt", header=T, as.is=T, sep="\t")
bindf$ct_control <- NA
bindf$ct_e2 <- NA
condf <- read.csv(paste0("bedfiles/", args[1], "_h", args[2], "_withoutE2.hg19.bed"), header=F, as.is=T, sep="\t")
expdf <- read.csv(paste0("bedfiles/", args[1], "_h", args[2], "_withE2.hg19.bed"), header=F, as.is=T, sep="\t")
for (i in 1:nrow(bindf)){
  bindf$ct_control[i] <- length(unique(condf$V2[condf$V1 == paste0("chr", bindf$chromosome[i]) & condf$V2 >= bindf$start[i] - marginlength & condf$V2 <= bindf$end[i] + marginlength]))
  bindf$ct_e2[i] <- length(unique(expdf$V2[expdf$V1 == paste0("chr", bindf$chromosome[i]) & expdf$V2 >= bindf$start[i] - marginlength & expdf$V2 <= bindf$end[i] + marginlength]))
}
write.table(bindf, paste0("Counts.", args[1], ".", args[2], ".per.gene.e2.ratio.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

geneofi <- "RARA"
geneofi <- "SHANK2"
cellofi <- "T47D"
cellofi <- "MCF7"

tempdf <- read.csv(paste0("../28_HTGTS/Counts.", cellofi, ".", geneofi, ".per.gene.e2.ratio.txt"), header=T, as.is=T, sep="\t")
tempdf$ctrl <- tempdf$ct_control+1
tempdf$e2 <- tempdf$ct_e2+1
tempdf$ratio <- tempdf$e2/tempdf$ctrl
tempdf$lograt <- log2(tempdf$ratio)
tempdf$inclusion <- "yes"
if (geneofi == "SHANK2"){
  tempdf$inclusion[tempdf$chromosome=="11" & tempdf$start <= 70658290 + 1000000 & tempdf$end >= 70658290 - 1000000] <- "bait"
} else if (geneofi == "RARA"){
  tempdf$inclusion[tempdf$chromosome=="17" & tempdf$start <= 38491451 + 1000000 & tempdf$end >= 38491451 - 1000000] <- "bait"
}
for (i in which(tempdf$chromosome %in% c("13", "14", "15", "21", "22"))){
  if (tempdf$start[i] < hg19_coord$centro_end[hg19_coord$chr == tempdf$chromosome[i]]){
    tempdf$inclusion[i] <- "acrocentric"
  }
}
for (i in 1:nrow(tempdf)){
  if (str_sub(tempdf$gene[i], -4) == "-AS1"){
    tempdf$inclusion[i] <- "as1"
  }
}

hist(tempdf$lograt[tempdf$inclusion=="yes"])
tempdf <- tempdf[order(tempdf$lograt, decreasing = T),]
write.table(tempdf[tempdf$inclusion %in% c("yes", "bait") ,c("gene", "lograt", "inclusion")], paste0("htgts_intermediates/lograt.", cellofi, ".", geneofi, ".gene.v2.merge"), row.names = F, col.names = F, quote = F, sep = "\t")

## Pooling the gene-based analysis outcome for all genes -- The output of this code is "HTGTS.lograt.average.allcells.twotarget.gene.v2.corrected.allinfo.txt", which is available in the Data folder.
tempdf <- read.csv("htgts_intermediates/lograt.MCF7.SHANK2.gene.v2.merge", header=F, as.is=T, sep="\t")
tempdf$V2[tempdf$V3 == "bait"] <- NA
colnames(tempdf) <- c("gene", "mcf7.shank2", "inclusion")
df <- tempdf

tempdf <- read.csv("htgts_intermediates/lograt.MCF7.RARA.gene.v2.merge", header=F, as.is=T, sep="\t")
tempdf$V2[tempdf$V3 == "bait"] <- NA
colnames(tempdf) <- c("gene", "mcf7.rara", "inclusion")
table(tempdf$gene == df$gene)
df$mcf7.rara <- tempdf$mcf7.rara

tempdf <- read.csv("htgts_intermediates/lograt.T47D.SHANK2.gene.v2.merge", header=F, as.is=T, sep="\t")
tempdf$V2[tempdf$V3 == "bait"] <- NA
colnames(tempdf) <- c("gene", "t47d.shank2", "inclusion")
table(tempdf$gene == df$gene)
df$t47d.shank2 <- tempdf$t47d.shank2

tempdf <- read.csv("htgts_intermediates/lograt.T47D.RARA.gene.v2.merge", header=F, as.is=T, sep="\t")
tempdf$V2[tempdf$V3 == "bait"] <- NA
colnames(tempdf) <- c("gene", "t47d.rara", "inclusion")
table(tempdf$gene == df$gene)
df$t47d.rara <- tempdf$t47d.rara
df <- df[,colnames(df) != "inclusion"]

df$average <- NA
for (i in 1:nrow(df)){
  df$average[i] <- mean(c(df$mcf7.shank2[i], df$mcf7.rara[i], df$t47d.shank2[i], df$t47d.rara[i]), na.rm=T)
}
df <- df[order(df$average, decreasing = T),]
write.table(df, paste0("HTGTS.lograt.average.allcells.twotarget.gene.v2.corrected.allinfo.txt"), row.names = F, col.names = T, quote = F, sep = "\t")


