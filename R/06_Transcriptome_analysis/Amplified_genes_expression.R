library(dplyr)
library(grid)
library(fields)
library(stringr)
library(plotrix)
library(mapplots)
library(viridis)
library(reshape2)
library(vioplot)

# Expression level of the genes in the amplicon (Supplementary Note)
## Sanger study expression profile
geneofi <- c("ERBB2", "CCND1", "FADD", "PAK1", "RSF1", "USP32", "PPM1D", "ZNF703", "NSD3", "FGFR1", "TRPS1", "MYC", "ZNF217", "QRSL1", "C6orf203")
ebasis <- read.csv("Supplementary Table 7.Transcriptomic.342.txt", header=T, as.is=T, sep="\t")
rnalist <- colnames(ebasis)[6:ncol(ebasis)]
rnalist <- gsub("PR", "PD", gsub(".RNA", "", rnalist))
for (i in 1:length(rnalist)){
  rnalist[i] <- strsplit(strsplit(strsplit(rnalist[i], "a", fixed=T)[[1]][1], "b", fixed=T)[[1]][1], "c", fixed=T)[[1]][1]
}
rnalist <- rnalist[rnalist %in% purdf$study_id[purdf$availability=="yes"]]

df <- NULL
for (i in geneofi){
  if (i %in% c("ERBB2", "CCND1", "RSF1", "USP32", "PPM1D", "ZNF703", "FGFR1", "MYC", "ZNF217")){
    print(i)
    tempdf <- as.data.frame(c(masterdf$study_id[masterdf$study_id %in% drvp$study_id[drvp$gene == i & drvp$driver == "AMP"] & masterdf$study_id %in% rnalist], masterdf$study_id[!(masterdf$study_id %in% drvp$study_id[drvp$gene == i & drvp$driver == "AMP"]) & masterdf$study_id %in% rnalist]))
    colnames(tempdf) <- "study_id"
    tempdf$amplification <- c(rep("yes", length(masterdf$study_id[masterdf$study_id %in% drvp$study_id[drvp$gene == i & drvp$driver == "AMP"] & masterdf$study_id %in% rnalist])), rep("no", length(masterdf$study_id[!(masterdf$study_id %in% drvp$study_id[drvp$gene == i & drvp$driver == "AMP"]) & masterdf$study_id %in% rnalist])))
    tempdf$expression <- NA
    for (j in 1:nrow(tempdf)){
      tempdf$expression[j] <- ebasis[ebasis$Name==i,tempdf$study_id[j]]
    }
    tempdf$gene <- i
    averageval <- mean(tempdf$expression)
    sdval <- sd(tempdf$expression)
    tempdf$zscore <- NA
    for (j in 1:nrow(tempdf)){
      tempdf$zscore[j] <- (tempdf$expression[j] - averageval)/sdval
    }
    if (i == "ERBB2"){
      df <- tempdf
    } else {
      df <- rbind(df, tempdf)
    }
  } else {
    print(i)
    if (i == "FADD"){
      ebidf <- read.csv(paste0("Gistic2/expression/Additional.AMP.11q13.3_6.txt"), header=T, as.is=T, sep="\t")
    } else if (i == "PAK1"){
      ebidf <- read.csv(paste0("Gistic2/expression/Additional.AMP.11q14.1_5.txt"), header=T, as.is=T, sep="\t")
    } else if (i == "USP32"){
      ebidf <- read.csv(paste0("Gistic2/expression/Additional.AMP.11q14.1_5.txt"), header=T, as.is=T, sep="\t")
    } else if (i == "NSD3"){
      i <- "WHSC1L1"
      ebidf <- read.csv(paste0("Gistic2/expression/Additional.AMP.8p11.21_35.txt"), header=T, as.is=T, sep="\t")
    } else if (i == "QRSL1" | i == "C6orf203"){
      ebidf <- read.csv(paste0("Gistic2/expression/Additional.AMP.6q21_15.txt"), header=T, as.is=T, sep="\t")
    } else if (i == "TRPS1"){
      ebidf <- read.csv(paste0("Gistic2/expression/Additional.AMP.8q23.3_10.txt"), header=T, as.is=T, sep="\t")
    }
    tempdf <- as.data.frame(c(masterdf$study_id[masterdf$study_id %in% ebidf$study_id[ebidf$gene == i & ebidf$driver == "AMP"] & masterdf$study_id %in% rnalist], masterdf$study_id[!(masterdf$study_id %in% ebidf$study_id[ebidf$gene == i & ebidf$driver == "AMP"]) & masterdf$study_id %in% rnalist]))
    colnames(tempdf) <- "study_id"
    tempdf$amplification <- c(rep("yes", length(masterdf$study_id[masterdf$study_id %in% ebidf$study_id[ebidf$gene == i & ebidf$driver == "AMP"] & masterdf$study_id %in% rnalist])), rep("no", length(masterdf$study_id[!(masterdf$study_id %in% ebidf$study_id[ebidf$gene == i & ebidf$driver == "AMP"]) & masterdf$study_id %in% rnalist])))
    tempdf$expression <- NA
    for (j in 1:nrow(tempdf)){
      tempdf$expression[j] <- ebasis[ebasis$Name==i,tempdf$study_id[j]]
    }
    tempdf$gene <- i
    averageval <- mean(tempdf$expression)
    sdval <- sd(tempdf$expression)
    tempdf$zscore <- NA
    for (j in 1:nrow(tempdf)){
      tempdf$zscore[j] <- (tempdf$expression[j] - averageval)/sdval
    }
    df <- rbind(df, tempdf)
  }
}
df$gene[df$gene == "WHSC1L1"] <- "NSD3"
df$gene <- factor(df$gene, levels = c("ERBB2", "CCND1", "FADD", "PAK1", "RSF1", "USP32", "PPM1D", "ZNF703", "NSD3", "FGFR1", "TRPS1", "MYC", "ZNF217", "QRSL1", "C6orf203"))

pdf("Boxplot.multiple.genes.basis.pdf", height=4, width=8)
boxplot(df[,"expression"] ~ df[,"gene"], boxfill = NA, border = NA, frame=F, las=2, xlab="", ylab="mRNA expression in BASIS", ylim=c(-8,10)) #invisible boxes - only axes and plot area
boxplot(df[df$amplification=="no", "expression"] ~ df[df$amplification=="no","gene"], xaxt = "n", yaxt = "n", frame=F, outline=F, add = TRUE, boxfill="#003469", 
        boxwex=0.25, at = 1:length(unique(df$gene)) - 0.15, ylab="", ylim=c(-3,13)) #shift these left by -0.15
boxplot(df[df$amplification=="yes", "expression"] ~ df[df$amplification=="yes", "gene"], xaxt = "n", yaxt = "n", frame=F, outline=F, add = TRUE, boxfill="#A51C30", 
        boxwex=0.25, at = 1:length(unique(df$gene)) + 0.15, ylab="", ylim=c(-3,13)) #shift to the right by +0.15
dev.off()

table(df$gene, df$amplification)
prop.table(table(df$gene, df$amplification),1)
for (i in geneofi){
  print(i)
  print(t.test(df$expression[df$gene == i] ~ df$amplification[df$gene == i]))
}
for (i in geneofi){
  print(i)
  print(wilcox.test(df$expression[df$gene == i] ~ df$amplification[df$gene == i]))
}

# CRISPR screen analysis (Supplementary Note)
library(readxl)
ct47d <- read_excel("../30_References/Fei_PNAS_2019_CRISPR/pnas.1908155116.sd01.xlsx")
ct47d <- as.data.frame(ct47d)
head(ct47d[order(ct47d$`T47D_Day21_vs_T47D_Day0|pos|rank`, decreasing = F),], 10)

plot(ct47d$`T47D_Day21_vs_T47D_Day0|neg|rank`, log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`)*(-1), pch=20, col=rgb(0,0,0,.1))

points(ct47d$`T47D_Day21_vs_T47D_Day0|neg|rank`[ct47d$id=="CCND1"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="CCND1"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|neg|rank`[ct47d$id=="ERBB2"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="ERBB2"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|neg|rank`[ct47d$id=="MYC"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="MYC"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|neg|rank`[ct47d$id=="FADD"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="FADD"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|neg|rank`[ct47d$id=="RSF1"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="RSF1"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|neg|rank`[ct47d$id=="PAK1"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="PAK1"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|neg|rank`[ct47d$id=="USP32"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="USP32"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|neg|rank`[ct47d$id=="PPM1D"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="PPM1D"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|neg|rank`[ct47d$id=="ZNF703"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="ZNF703"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|neg|rank`[ct47d$id=="ZNF217"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="ZNF217"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|neg|rank`[ct47d$id=="QRSL1"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="QRSL1"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|neg|rank`[ct47d$id=="WHSC1L1"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="WHSC1L1"])*(-1), pch=19, col=rgb(1,0,0,.5))

plot(ct47d$`T47D_Day21_vs_T47D_Day0|pos|rank`, log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`)*(-1), pch=20, col=rgb(0,0,0,.1))
points(ct47d$`T47D_Day21_vs_T47D_Day0|pos|rank`[ct47d$id=="CCND1"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="CCND1"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|pos|rank`[ct47d$id=="ERBB2"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="ERBB2"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|pos|rank`[ct47d$id=="MYC"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="MYC"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|pos|rank`[ct47d$id=="FADD"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="FADD"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|pos|rank`[ct47d$id=="RSF1"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="RSF1"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|pos|rank`[ct47d$id=="PAK1"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="PAK1"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|pos|rank`[ct47d$id=="USP32"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="USP32"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|pos|rank`[ct47d$id=="PPM1D"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="PPM1D"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|pos|rank`[ct47d$id=="ZNF703"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="ZNF703"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|pos|rank`[ct47d$id=="ZNF217"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="ZNF217"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|pos|rank`[ct47d$id=="QRSL1"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="QRSL1"])*(-1), pch=19, col=rgb(1,0,0,.5))
points(ct47d$`T47D_Day21_vs_T47D_Day0|pos|rank`[ct47d$id=="WHSC1L1"], log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="WHSC1L1"])*(-1), pch=19, col=rgb(1,0,0,.5))

plot(log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`)*(-1), log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`)*(-1), pch=20, col=rgb(0,0,0,.1))
points(log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="CCND1"])*(-1), log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="CCND1"])*(-1), pch=19, col=rgb(1,0,0,.7))
points(log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="ERBB2"])*(-1), log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="ERBB2"])*(-1), pch=19, col=rgb(1,0,0,.7))
points(log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="MYC"])*(-1), log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="MYC"])*(-1), pch=19, col=rgb(1,0,0,.7))
points(log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="TRPS1"])*(-1), log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="TRPS1"])*(-1), pch=19, col=rgb(1,0,0,.7))
points(log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="FADD"])*(-1), log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="FADD"])*(-1), pch=19, col=rgb(1,0,0,.7))
points(log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="RSF1"])*(-1), log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="RSF1"])*(-1), pch=19, col=rgb(1,0,0,.7))
points(log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="PAK1"])*(-1), log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="PAK1"])*(-1), pch=19, col=rgb(1,0,0,.7))
points(log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="USP32"])*(-1), log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="USP32"])*(-1), pch=19, col=rgb(1,0,0,.7))
points(log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="PPM1D"])*(-1), log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="PPM1D"])*(-1), pch=19, col=rgb(1,0,0,.7))
points(log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="ZNF703"])*(-1), log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="ZNF703"])*(-1), pch=19, col=rgb(1,0,0,.7))
points(log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="ZNF217"])*(-1), log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="ZNF217"])*(-1), pch=19, col=rgb(1,0,0,.7))
points(log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="QRSL1"])*(-1), log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="QRSL1"])*(-1), pch=19, col=rgb(1,0,0,.7))
points(log10(ct47d$`T47D_Day21_vs_T47D_Day0|pos|score`[ct47d$id=="WHSC1L1"])*(-1), log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|score`[ct47d$id=="WHSC1L1"])*(-1), pch=19, col=rgb(1,0,0,.7))

plot(ct47d$`T47D_Day21_vs_T47D_Day0|neg|lfc`, log10(ct47d$`T47D_Day21_vs_T47D_Day0|neg|p-value`)*(-1), pch=20, col=rgb(0,0,0,.1))

## DepMap analysis
library(data.table)
depmapgenelist <- unlist(read.csv("../../References/depmap_public_22Q1/CRISPR_gene_effect.csv", header=F, as.is=T, nrows=1), use.names = F)
depmapgenelist[str_sub(depmapgenelist, 1, 5) == "TRPS1"]
depmapgenelist[str_sub(depmapgenelist, 1, 5) == "PAK1 "]
depmapgenelist[str_sub(depmapgenelist, 1, 7) == "MTRES1 "]
df <- cbind(as.data.frame(fread("../../References/depmap_public_22Q1/CRISPR_gene_effect.csv", header=T, stringsAsFactors=F, select = "DepMap_ID")), as.data.frame(fread("../../References/depmap_public_22Q1/CRISPR_gene_effect.csv", header=T, stringsAsFactors=F, select = depmapgenelist[str_sub(depmapgenelist, 1, 5) == "TRPS1"])))

for (i in c("ERBB2", "CCND1", "FADD", "PAK1", "RSF1", "USP32", "PPM1D", "ZNF703", "NSD3", "FGFR1", "TRPS1", "MYC", "ZNF217", "QRSL1", "MTRES1")){
  if (i == "ERBB2"){
    df <- cbind(as.data.frame(fread("../../References/depmap_public_22Q1/CRISPR_gene_effect.csv", header=T, stringsAsFactors=F, select = "DepMap_ID")), as.data.frame(fread("../../References/depmap_public_22Q1/CRISPR_gene_effect.csv", header=T, stringsAsFactors=F, select = depmapgenelist[str_sub(depmapgenelist, 1, nchar(i)+1) == paste0(i, " ")])))
  } else {
    df <- cbind(df, as.data.frame(fread("../../References/depmap_public_22Q1/CRISPR_gene_effect.csv", header=T, stringsAsFactors=F, select = depmapgenelist[str_sub(depmapgenelist, 1, nchar(i)+1) == paste0(i, " ")])))
  }
}
colnames(df)[2:ncol(df)] <- c("ERBB2", "CCND1", "FADD", "PAK1", "RSF1", "USP32", "PPM1D", "ZNF703", "NSD3", "FGFR1", "TRPS1", "MYC", "ZNF217", "QRSL1", "MTRES1")

dms <- read.csv("../../References/depmap_public_22Q1/sample_info.csv", header=T, as.is=T)
dms <- dms[dms$lineage == "breast",]
for (i in c("ERBB2", "CCND1", "FADD", "PAK1", "RSF1", "USP32", "PPM1D", "ZNF703", "NSD3", "FGFR1", "TRPS1", "MYC", "ZNF217", "QRSL1", "MTRES1")){
  dms[,ncol(dms)+1] <- NA
  colnames(dms)[ncol(dms)] <- i
  for (j in 1:nrow(dms)){
    if (nrow(df[df$DepMap_ID == dms$DepMap_ID[j],]) != 0){
      dms[j,ncol(dms)] <- df[,i][df$DepMap_ID==dms$DepMap_ID[j]]
    }
  }
}

dms <- dms[!is.na(dms$ERBB2),]
boxplot(dms$ERBB2 ~ dms$lineage_sub_subtype,las=2)
boxplot(dms$CCND1 ~ dms$lineage_sub_subtype,las=2)
boxplot(dms$FADD ~ dms$lineage_sub_subtype,las=2)
boxplot(dms$PAK1 ~ dms$lineage_sub_subtype,las=2)
boxplot(dms$RSF1 ~ dms$lineage_sub_subtype,las=2)
boxplot(dms$USP32 ~ dms$lineage_sub_subtype,las=2)
boxplot(dms$ZNF703 ~ dms$lineage_sub_subtype,las=2)
boxplot(dms$ZNF217 ~ dms$lineage_sub_subtype,las=2)

## DepMap analysis for all genes in the GISTIC amplicon list
library(data.table)
infolist <- depmapgenelist[depmapgenelist != "DepMap_ID"]
for (i in 1:length(infolist)){
  infolist[i] <- strsplit(infolist[i], " (", fixed=T)[[1]][1]
}
tempdf <- read.csv("Gistic2/genes.in.the.amplicons.gistic.tsv", header=T, as.is=T, sep="\t")
for (w in 1:20){
  if (tempdf$genes[w] == "ZNF703"){
    geneofi <- c("ZNF703", "NSD3", "FGFR1")
  } else if (tempdf$genes[w] == "CCND1"){
    geneofi <- c("CCND1", "FGF3")
  } else if (tempdf$genes[w] == "RPS6KB1;RNFT1;TUBD1;VMP1;MIR21"){
    geneofi <- c("RPS6KB1", "RNFT1", "TUBD1", "VMP1", "MIR21", "PPM1D", "USP32", "APPBP2")
  } else {
    geneofi <- strsplit(tempdf$genes[w], ";", fixed=T)[[1]]
  }
  for (i in geneofi){
    if (i == "ERBB2"){
      df <- cbind(as.data.frame(fread("../../References/depmap_public_22Q1/CRISPR_gene_effect.csv", header=T, stringsAsFactors=F, select = "DepMap_ID")), as.data.frame(fread("../../References/depmap_public_22Q1/CRISPR_gene_effect.csv", header=T, stringsAsFactors=F, select = depmapgenelist[str_sub(depmapgenelist, 1, nchar(i)+1) == paste0(i, " ")])))
      print(paste0(w, "_", tempdf$cytoband[w], "_", i))
    } else {
      if (i == "C6orf203"){
        i <- "MTRES1"
      } else if (i == "WHSC1L1"){
        i <- "NSD3"
      }
      if (i %in% infolist){
        df <- cbind(df, as.data.frame(fread("../../References/depmap_public_22Q1/CRISPR_gene_effect.csv", header=T, stringsAsFactors=F, select = depmapgenelist[str_sub(depmapgenelist, 1, nchar(i)+1) == paste0(i, " ")])))
        print(paste0(w, "_", tempdf$cytoband[w], "_", i))
      }
    }
  }
}
gisticgenelist <- colnames(df)[2:length(colnames(df))]
for (i in 1:length(gisticgenelist)){
  gisticgenelist[i] <- strsplit(gisticgenelist[i], " (", fixed=T)[[1]][1]
}
colnames(df)[2:ncol(df)] <- gisticgenelist

"MTRES1" %in% gisticgenelist
"NSD3" %in% gisticgenelist
"WHSC1L1" %in% gisticgenelist
"C6orf203" %in% gisticgenelist

dms <- read.csv("../../References/depmap_public_22Q1/sample_info.csv", header=T, as.is=T)
dms <- dms[dms$lineage == "breast",]
#table(dms$lineage_molecular_subtype[dms$lineage=="breast"])
#table(dms$lineage_sub_subtype[dms$lineage=="breast"])

for (i in gisticgenelist){
  dms[,ncol(dms)+1] <- NA
  colnames(dms)[ncol(dms)] <- i
  for (j in 1:nrow(dms)){
    if (nrow(df[df$DepMap_ID == dms$DepMap_ID[j],]) != 0){
      dms[j,ncol(dms)] <- df[,i][df$DepMap_ID==dms$DepMap_ID[j]]
    }
  }
}
dms <- dms[!is.na(dms$ERBB2),]
write.table(dms, "depmap.genes.amplicon.effect.tsv", row.names = F, col.names = T, quote = F, sep = "\t")

for (i in gisticgenelist){
  pdf(paste0("Gistic2/plot_depmap/", i, ".depmap.pdf"), height=5, width=4)
  par(mar=c(8.1, 4.1, 4.1, 2.1))
  boxplot(dms[,i] ~ dms$lineage_sub_subtype, las=2, frame=F, outline=F, xlab="", ylab="Gene effect")
  stripchart(dms[,i] ~ dms$lineage_sub_subtype, pch=20, col=rgb(0,0,0,.5), vertical=T, method="jitter", add=T)
  dev.off()
}

fivecolors <- c("#DC143C", "#FFB5C5", "#FF1493", "#FF9912", "#C1C1C1")
geneofi <- c("ERBB2", "TRPS1", "FADD", "QRSL1")
for (i in geneofi){
  pdf(paste0("FIP_revision/gene_effect_score_depmap/", i, ".depmap.persubsubtype.pdf"), height=2.7, width=6)
  par(mar=c(5.1, 9.1, 4.1, 2.1))
  boxplot(dms[,i] ~ dms$lineage_sub_subtype, ylim=c(-1.5,0.5), las=1, frame=F, outline=F, ylab="", xlab="Gene effect", horizontal=T, col=rev(fivecolors))
  stripchart(dms[,i] ~ dms$lineage_sub_subtype, pch=20, col=rgb(0,0,0,.5), vertical=F, method="jitter", add=T)
  dev.off()
}

library(vioplot)
fivecolors <- c("#DC143C", "#FFB5C5", "#FF1493", "#FF9912", "#C1C1C1")
geneofi <- c("ERBB2", "CCND1", "MYC", "TRPS1", "FADD", "QRSL1", "MTRES1", "BEND3", "USP32", "PPM1D")
for (i in geneofi){
  pdf(paste0("FIP_revision/gene_effect_score_depmap/Vioplot.", i, ".depmap.persubsubtype.reproduce.pdf"), height=2.8, width=4.5)
  par(mar=c(5.1, 9.1, 4.1, 2.1))
  vioplot(dms[,i] ~ dms$lineage_sub_subtype, ylim=c(-2.5,0.5), las=1, frame=, ylab="", xlab="Gene effect", horizontal=T, col=rev(fivecolors))
  dev.off()
}

