library(dplyr)
library(grid)
library(fields)
library(stringr)
library(plotrix)
library(mapplots)
library(viridis)
library(reshape2)
library(data.table)
library(vioplot)
library(beeswarm)

# ER activity analysis -- Using 268 transcriptomes in the Sanger study --> 263 samples available in our genome analysis
ebasis <- read.csv("Supplementary Table 7.Transcriptomic.342.txt", header=T, as.is=T, sep="\t") ## We used supplementary table 7 of the reference 15 (Nik-Zainal et al. 2016 Nature).
rnalist <- colnames(ebasis)[6:ncol(ebasis)]
rnalist <- gsub("PR", "PD", gsub(".RNA", "", rnalist))
for (i in 1:length(rnalist)){
  rnalist[i] <- strsplit(strsplit(strsplit(rnalist[i], "a", fixed=T)[[1]][1], "b", fixed=T)[[1]][1], "c", fixed=T)[[1]][1]
}
length(rnalist)
colnames(ebasis)[6:ncol(ebasis)] <- rnalist

table(rnalist %in% purdf$study_id[purdf$availability=="yes"])
rnalist <- rnalist[rnalist %in% purdf$study_id[purdf$availability=="yes"]]

## Quick validation -- ER, PR, HER2 status in the expression profile
rdf <- as.data.frame(matrix(NA, ncol = 8, nrow = length(rnalist)))
colnames(rdf) <- c("study_id", "is.wgd", "er_plot", "pr_plot", "her2_plot", "pam50_plot", "hr_status", "tba_num")
rdf$study_id <- rnalist
for (i in 1:nrow(rdf)){
  rdf$is.wgd[i] <- purdf$is.wgd[purdf$availability == "yes" & purdf$study_id == rdf$study_id[i]]
  rdf$er_plot[i] <- hrdf$er_plot[hrdf$study_id == rdf$study_id[i]]
  rdf$pr_plot[i] <- hrdf$pr_plot[hrdf$study_id == rdf$study_id[i]]
  rdf$her2_plot[i] <- hrdf$her2_plot[hrdf$study_id == rdf$study_id[i]]
  rdf$pam50_plot[i] <- hrdf$pam50_plot[hrdf$study_id == rdf$study_id[i]]
  rdf$hr_status[i] <- hrdf$hr_status[hrdf$study_id == rdf$study_id[i]]
  rdf$tba_num[i] <- oncdf$tba[oncdf$study_id == rdf$study_id[i]]
}

eg <- read.csv("/Users/jjklee/gsea_home/output/apr21/corrected.htgts.twocells.twogenes.GseaPreranked.1650575918912/HALLMARK_ESTROGEN_RESPONSE_EARLY.tsv", header=T, as.is=T, sep="\t")
eg <- eg$SYMBOL
table(eg %in% ebasis$Name)
eg <- eg[eg %in% ebasis$Name]

lg <- read.csv("/Users/jjklee/gsea_home/output/apr21/corrected.htgts.twocells.twogenes.GseaPreranked.1650575918912/HALLMARK_ESTROGEN_RESPONSE_LATE.tsv", header=T, as.is=T, sep="\t")
lg <- lg$SYMBOL
table(lg %in% ebasis$Name)
lg <- lg[lg %in% ebasis$Name]
lg <- lg[!(lg %in% eg)]

eg <- c(eg, lg)

for (i in 1:length(eg)){
  print(paste0(eg[i], " ", i))
  rdf[,ncol(rdf)+1] <- NA
  colnames(rdf)[ncol(rdf)] <- eg[i]
  for (j in 1:nrow(rdf)){
    if (nrow(ebasis[ebasis$Name == eg[i],]) == 1){
      rdf[j,ncol(rdf)] <- ebasis[ebasis$Name == eg[i], which(colnames(ebasis) == rdf$study_id[j])]
    } else {
      rdf[j,ncol(rdf)] <- ebasis[ebasis$Name == eg[i] & ebasis$UNIQID == min(ebasis$UNIQID[ebasis$Name == eg[i]]), which(colnames(ebasis) == rdf$study_id[j])]
    }
  }
}
rdf$ESR1 <- NA
for (j in 1:nrow(rdf)){
  if (nrow(ebasis[ebasis$Name == "ESR1",]) == 1){
    rdf$ESR1[j] <- ebasis[ebasis$Name == "ESR1", which(colnames(ebasis) == rdf$study_id[j])]
  } else {
    rdf$ESR1[j] <- ebasis[ebasis$Name == "ESR1" & ebasis$UNIQID == min(ebasis$UNIQID[ebasis$Name == "ESR1"]), which(colnames(ebasis) == rdf$study_id[j])]
  }
}
write.table(rdf, "rna.expression.basis.263.samples.txt", row.names = F, col.names = T, quote = F, sep = "\t")
rdf <- read.csv("rna.expression.basis.263.samples.txt", header=T, as.is=T, sep="\t")

## Genes associated with differentiation status of breast cancer
rdf$NEUROD1 <- NA
rdf$SYP <- NA
rdf$CHGA <- NA
rdf$NCAM1 <- NA
for (i in c("NEUROD1", "SYP", "CHGA", "NCAM1")){
  for (j in 1:nrow(rdf)){
    if (nrow(ebasis[ebasis$Name == i,]) == 1){
      rdf[j,which(colnames(rdf) == i)] <- ebasis[ebasis$Name == i, which(colnames(ebasis) == rdf$study_id[j])]
    } else {
      rdf[j,which(colnames(rdf) == i)] <- ebasis[ebasis$Name == i & ebasis$UNIQID == min(ebasis$UNIQID[ebasis$Name == i]), which(colnames(ebasis) == rdf$study_id[j])]
    }
  }
}

rdf$ENO2 <- NA
rdf$REST <- NA
rdf$CALCA <- NA
rdf$TP53 <- NA
for (i in c("ENO2", "REST", "CALCA", "TP53")){
  for (j in 1:nrow(rdf)){
    if (nrow(ebasis[ebasis$Name == i,]) == 1){
      rdf[j,which(colnames(rdf) == i)] <- ebasis[ebasis$Name == i, which(colnames(ebasis) == rdf$study_id[j])]
    } else {
      rdf[j,which(colnames(rdf) == i)] <- ebasis[ebasis$Name == i & ebasis$UNIQID == min(ebasis$UNIQID[ebasis$Name == i]), which(colnames(ebasis) == rdf$study_id[j])]
    }
  }
}

rdf$EPCAM <- NA
rdf$GATA3 <- NA
rdf$ALDH2 <- NA
rdf$ZEB1 <- NA
for (i in c("EPCAM", "GATA3", "ALDH2", "ZEB1")){
  for (j in 1:nrow(rdf)){
    if (nrow(ebasis[ebasis$Name == i,]) == 1){
      rdf[j,which(colnames(rdf) == i)] <- ebasis[ebasis$Name == i, which(colnames(ebasis) == rdf$study_id[j])]
    } else {
      rdf[j,which(colnames(rdf) == i)] <- ebasis[ebasis$Name == i & ebasis$UNIQID == min(ebasis$UNIQID[ebasis$Name == i]), which(colnames(ebasis) == rdf$study_id[j])]
    }
  }
}

rdf$ITGA6 <- NA
rdf$KRT14 <- NA
rdf$KRT8 <- NA
rdf$VIM <- NA
for (i in c("ITGA6", "KRT14", "KRT8", "VIM")){
  for (j in 1:nrow(rdf)){
    if (nrow(ebasis[ebasis$Name == i,]) == 1){
      rdf[j,which(colnames(rdf) == i)] <- ebasis[ebasis$Name == i, which(colnames(ebasis) == rdf$study_id[j])]
    } else {
      rdf[j,which(colnames(rdf) == i)] <- ebasis[ebasis$Name == i & ebasis$UNIQID == min(ebasis$UNIQID[ebasis$Name == i]), which(colnames(ebasis) == rdf$study_id[j])]
    }
  }
}

rdf$ALDH1A1 <- NA
rdf$CLDN1 <- NA
rdf$KRT18 <- NA
rdf$KRT19 <- NA
rdf$KIT <- NA
rdf$TP63 <- NA
rdf$CDH3 <- NA
rdf$EGFR <- NA
rdf$KRT6 <- NA
rdf$KRT5 <- NA
for (i in c("ALDH1A1", "CLDN1", "KRT19", "KIT", "TP63", "CDH3", "EGFR", "KRT5")){
  for (j in 1:nrow(rdf)){
    if (nrow(ebasis[ebasis$Name == i,]) == 1){
      rdf[j,which(colnames(rdf) == i)] <- ebasis[ebasis$Name == i, which(colnames(ebasis) == rdf$study_id[j])]
    } else {
      rdf[j,which(colnames(rdf) == i)] <- ebasis[ebasis$Name == i & ebasis$UNIQID == min(ebasis$UNIQID[ebasis$Name == i]), which(colnames(ebasis) == rdf$study_id[j])]
    }
  }
}

for (i in 9:ncol(rdf)){
  pdf(paste0("ER_activity/plots/boxplot.", colnames(rdf)[i], ".erstatus.pdf"), height=5, width=4)
  boxplot(rdf[,i] ~ rdf$er_plot, las=1, ylab = paste0("mRNA expression of ", colnames(rdf)[i]), xlab="", main=t.test(rdf[rdf$er_plot == 1,i], rdf[rdf$er_plot == 2,i])$p.value, frame=F, ylim=c(-8,+15), outline=F, col=c(rgb(219/255, 112/255, 147/255, .5), rgb(154/255, 50/255, 205/255, .5), rgb(0,0,0,.3)))
  stripchart(rdf[,i] ~ rdf$er_plot, pch=19, col=rgb(0,0,0,.2), vertical=T, add=T, method="jitter")
  dev.off()
}

## Defining ER transcriptome activity
tempdf <- as.data.frame(matrix(NA, ncol = 4, nrow = length(eg)))
colnames(tempdf) <- c("gene", "mean.ERp", "mean.ERn", "t.test.p")
tempdf$gene <- eg
for (i in 1:nrow(tempdf)){
  tempdf$mean.ERp[i] <- mean(rdf[rdf$er_plot == 1,which(colnames(rdf) == tempdf$gene[i])], na.rm = T)
  tempdf$mean.ERn[i] <- mean(rdf[rdf$er_plot == 2,which(colnames(rdf) == tempdf$gene[i])], na.rm = T)
  tempdf$t.test.p[i] <- t.test(rdf[rdf$er_plot == 1,which(colnames(rdf) == tempdf$gene[i])], rdf[rdf$er_plot == 2,which(colnames(rdf) == tempdf$gene[i])])$p.value
}
tempdf$significance <- "no"
tempdf$significance[tempdf$t.test.p < 0.05] <- "yes"
tempdf$marker <- "no"
tempdf$marker[tempdf$significance == "yes" & tempdf$mean.ERp > tempdf$mean.ERn] <- "yes"

for (i in which(tempdf$marker == "yes")){
  pdf(paste0("ER_activity/tbaplots/scatter.", tempdf$gene[i], ".tbastatus.pdf"), height=5, width=8)
  par(mfrow=c(1,2))
  boxplot(rdf[,which(colnames(rdf) == tempdf$gene[i])] ~ rdf$er_plot, las=1, ylab = paste0("mRNA expression of ", tempdf$gene[i]), xlab="", main=t.test(rdf[rdf$er_plot == 1,which(colnames(rdf) == tempdf$gene[i])], rdf[rdf$er_plot == 2,which(colnames(rdf) == tempdf$gene[i])])$p.value, frame=F, ylim=c(-8,+15), outline=F, col=c(rgb(219/255, 112/255, 147/255, .5), rgb(154/255, 50/255, 205/255, .5), rgb(0,0,0,.3)))
  stripchart(rdf[,which(colnames(rdf) == tempdf$gene[i])] ~ rdf$er_plot, pch=19, col=rgb(0,0,0,.2), vertical=T, add=T, method="jitter")
  boxplot(rdf[,which(colnames(rdf) == tempdf$gene[i])] ~ rdf$tba_num, las=1, ylab = paste0("mRNA expression of ", tempdf$gene[i]), xlab="", frame=F, ylim=c(-8,+15), outline=F)
  stripchart(rdf[,which(colnames(rdf) == tempdf$gene[i])] ~ rdf$tba_num, pch=19, col=rgb(0,0,0,.2), vertical=T, add=T, method="jitter")
  dev.off()
}

rdf$score10 <- 0
rdf$score30 <- 0
rdf$score50 <- 0
rdf$score70 <- 0
for (i in 9:which(colnames(rdf) == "UNC13B")){
  rdf$score10[rdf[,i] >= as.vector(quantile(rdf[,i], 0.9, na.rm=T)) & !is.na(rdf[,i])] <- rdf$score10[rdf[,i] >= as.vector(quantile(rdf[,i], 0.9, na.rm=T)) & !is.na(rdf[,i])]+1
  rdf$score30[rdf[,i] >= as.vector(quantile(rdf[,i], 0.7, na.rm=T)) & !is.na(rdf[,i])] <- rdf$score30[rdf[,i] >= as.vector(quantile(rdf[,i], 0.7, na.rm=T)) & !is.na(rdf[,i])]+1
  rdf$score50[rdf[,i] >= as.vector(quantile(rdf[,i], 0.5, na.rm=T)) & !is.na(rdf[,i])] <- rdf$score50[rdf[,i] >= as.vector(quantile(rdf[,i], 0.5, na.rm=T)) & !is.na(rdf[,i])]+1
  rdf$score70[rdf[,i] >= as.vector(quantile(rdf[,i], 0.3, na.rm=T)) & !is.na(rdf[,i])] <- rdf$score50[rdf[,i] >= as.vector(quantile(rdf[,i], 0.3, na.rm=T)) & !is.na(rdf[,i])]+1
}
boxplot(rdf$score30 ~ rdf$er_plot, las=1, frame=F, ylab="ER activity score", outline=F, col=c(rgb(219/255, 112/255, 147/255, .5), rgb(154/255, 50/255, 205/255, .5), rgb(0,0,0,.3)), ylim=c(0,300))
stripchart(rdf$score30 ~ rdf$er_plot, pch=19, col=rgb(0,0,0,.2), vertical=T, add=T, method="jitter")
boxplot(rdf$score50 ~ rdf$er_plot, las=1, frame=F, ylab="ER activity score", outline=F, col=c(rgb(219/255, 112/255, 147/255, .5), rgb(154/255, 50/255, 205/255, .5), rgb(0,0,0,.3)), ylim=c(0,300))
stripchart(rdf$score50 ~ rdf$er_plot, pch=19, col=rgb(0,0,0,.2), vertical=T, add=T, method="jitter")
boxplot(rdf$score70 ~ rdf$er_plot, las=1, frame=F, ylab="ER activity score", outline=F, col=c(rgb(219/255, 112/255, 147/255, .5), rgb(154/255, 50/255, 205/255, .5), rgb(0,0,0,.3)), ylim=c(0,300))
stripchart(rdf$score70 ~ rdf$er_plot, pch=19, col=rgb(0,0,0,.2), vertical=T, add=T, method="jitter")

pdf("FIP_revision/boxplot.ER.activity.score.pdf", height=5, width=2)
boxplot(rdf$score50 ~ rdf$er_plot, las=1, frame=F, ylab="ER activity score", outline=F, col=c(rgb(219/255, 112/255, 147/255, .5), rgb(154/255, 50/255, 205/255, .5), rgb(0,0,0,.3)), ylim=c(50,250))
stripchart(rdf$score50 ~ rdf$er_plot, pch=20, col=rgb(0,0,0,.2), vertical=T, add=T, method="jitter")
dev.off()

t.test(rdf$score50[rdf$er_plot == 1], rdf$score50[rdf$er_plot == 2])
table(rdf$er_plot)

pdf("FIP_revision/boxplot.ESR1.pdf", height=5, width=3)
boxplot(rdf$ESR1 ~ rdf$er_plot, las=1, frame=F, ylab="ESR1 mRNA expression", outline=F, col=c(rgb(219/255, 112/255, 147/255, .5), rgb(154/255, 50/255, 205/255, .5), rgb(0,0,0,.3)), ylim=c(-2,12))
stripchart(rdf$ESR1 ~ rdf$er_plot, pch=20, col=rgb(0,0,0,.2), vertical=T, add=T, method="jitter")
dev.off()

pdf("FIP_revision/boxplot.PGR.pdf", height=5, width=3)
boxplot(rdf$PGR ~ rdf$er_plot, las=1, frame=F, ylab="PGR mRNA expression", outline=F, col=c(rgb(219/255, 112/255, 147/255, .5), rgb(154/255, 50/255, 205/255, .5), rgb(0,0,0,.3)), ylim=c(-7,8))
stripchart(rdf$PGR ~ rdf$er_plot, pch=20, col=rgb(0,0,0,.2), vertical=T, add=T, method="jitter")
dev.off()

rdf$tba_plot <- rdf$tba_num
rdf$tba_plot[rdf$tba_plot >= 4] <- 4

pdf("FIP_revision/boxplot.ER.activity.tbamp.linear.regression.pdf", height=5, width=2.5)
boxplot(rdf$score50[rdf$er_plot == 1] ~ rdf$tba_plot[rdf$er_plot == 1], col=RColorBrewer::brewer.pal(9, "YlGnBu"), ylim=c(50,250), las=1, frame=F, outline=F, ylab="ER activity score", xlab="Number of TB amplification events per tumor")
stripchart(rdf$score50[rdf$er_plot == 1] ~ rdf$tba_plot[rdf$er_plot == 1], pch=20, col=rgb(0,0,0,.2), vertical=T, add=T, method="jitter")
abline(lm(rdf$score50[rdf$er_plot == 1] ~ rdf$tba_num[rdf$er_plot == 1]), lwd=2, lty=5, col=rgb(0,0,1,.5))
dev.off()

summary(lm(rdf$score50[rdf$er_plot == 1] ~ rdf$tba_num[rdf$er_plot == 1]))

## Key ER target genes
boxplot(rdf$GREB1[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1])
summary(lm(rdf$GREB1[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1]))

boxplot(rdf$TFF1[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1])
summary(lm(rdf$TFF1[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1])) # Used in the text

boxplot(rdf$CCND1[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1])
summary(lm(rdf$CCND1[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1]))

boxplot(rdf$PGR[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1])
summary(lm(rdf$PGR[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1])) # Used in the text

## Aldehyde dehydrogenases -- associated with stem cell features
summary(lm(rdf$ALDH1A1[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1])) # SIG
summary(lm(rdf$ALDH2[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1])) # SIG

boxplot(rdf$ALDH1A1[rdf$er_plot==1] ~ rdf$tba_plot[rdf$er_plot==1], col=RColorBrewer::brewer.pal(9, "YlGnBu"), ylim=c(-2,8), las=1, frame=F, outline=F, ylab="ASCL1 mRNA expression", xlab="Number of TB amplification events per tumor")
stripchart(rdf$ALDH1A1[rdf$er_plot == 1] ~ rdf$tba_plot[rdf$er_plot == 1], pch=20, col=rgb(0,0,0,.2), vertical=T, add=T, method="jitter")
abline(lm(rdf$ALDH1A1[rdf$er_plot == 1] ~ rdf$tba_num[rdf$er_plot == 1]), lwd=2, lty=5, col=rgb(0,0,1,.5))

boxplot(rdf$ALDH2[rdf$er_plot==1] ~ rdf$tba_plot[rdf$er_plot==1], col=RColorBrewer::brewer.pal(9, "YlGnBu"), ylim=c(-2,8), las=1, frame=F, outline=F, ylab="ASCL1 mRNA expression", xlab="Number of TB amplification events per tumor")
stripchart(rdf$ALDH2[rdf$er_plot == 1] ~ rdf$tba_plot[rdf$er_plot == 1], pch=20, col=rgb(0,0,0,.2), vertical=T, add=T, method="jitter")
abline(lm(rdf$ALDH2[rdf$er_plot == 1] ~ rdf$tba_num[rdf$er_plot == 1]), lwd=2, lty=5, col=rgb(0,0,1,.5))

## Basal marker  
summary(lm(rdf$TP63[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1])) # SIG
summary(lm(rdf$KRT5[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1])) # SIG -- CK5 (basal cytokeratin)
summary(lm(rdf$KRT14[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1])) # SIG -- CK14 (basal cytokeratin)

## Relationship between the extent of TB amplification and ER target gene expression in consideration of multiple hypotheses testing
p.adjust(summary(lm(rdf$PGR[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1]))$coefficients[2,4], method = "fdr", length(eg))
p.adjust(summary(lm(rdf$TFF1[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1]))$coefficients[2,4], method = "fdr", length(eg))
p.adjust(summary(lm(rdf$FOS[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1]))$coefficients[2,4], method = "fdr", length(eg))
p.adjust(summary(lm(rdf$KRT5[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1]))$coefficients[2,4], method = "fdr", length(eg))
p.adjust(summary(lm(rdf$KRT14[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1]))$coefficients[2,4], method = "fdr", length(eg))
p.adjust(summary(lm(rdf$TP63[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1]))$coefficients[2,4], method = "fdr", length(eg))
p.adjust(summary(lm(rdf$ASCL1[rdf$er_plot==1] ~ rdf$tba_num[rdf$er_plot==1]))$coefficients[2,4], method = "fdr", length(eg))

df <- data.frame(matrix(vector(), 0, 7, dimnames=list(c(), c("gene", "all_num", "er_num", "estimate", "rsquared", "pval", "fdr"))), stringsAsFactors=F)
genelist <- c("ESR1", eg)
for (i in genelist){
  df[nrow(df)+1,] <- NA
  df$gene[nrow(df)] <- i
  df$all_num[nrow(df)] <- nrow(rdf) - sum(is.na(rdf[,i]))
  df$er_num[nrow(df)] <- nrow(rdf[rdf$er_plot==1,]) - sum(is.na(rdf[rdf$er_plot==1,i]))
}
df <- df[df$all_num > 100,]
for (i in 1:nrow(df)){
  df$estimate[i] <- summary(lm(rdf[rdf$er_plot==1,df$gene[i]] ~ rdf$tba_num[rdf$er_plot==1]))$coefficients[2,1]
  df$rsquared[i] <- summary(lm(rdf[rdf$er_plot==1,df$gene[i]] ~ rdf$tba_num[rdf$er_plot==1]))$r.squared
  df$pval[i] <- summary(lm(rdf[rdf$er_plot==1,df$gene[i]] ~ rdf$tba_num[rdf$er_plot==1]))$coefficients[2,4]
  df$fdr[i] <- p.adjust(summary(lm(rdf[rdf$er_plot==1,df$gene[i]] ~ rdf$tba_num[rdf$er_plot==1]))$coefficients[2,4], method = "fdr", nrow(df))
}

df[df$fdr < 0.1,]
genelist <- c("ESR1", df$gene[df$fdr < 0.1])
for (i in genelist){
  pdf(paste0("FIP_revision/boxplot.expression.", i, ".tbamp.fdr0.1.pdf"), height=4.5, width=3)
  boxplot(rdf[rdf$er_plot==1, i] ~ rdf$tba_plot[rdf$er_plot==1], col=RColorBrewer::brewer.pal(9, "YlGnBu"), ylim=c(-10,15), las=1, frame=F, outline=F, ylab=paste0(i, " mRNA expression"), xlab="Number of TB amplification events per tumor")
  stripchart(rdf[rdf$er_plot == 1, i] ~ rdf$tba_plot[rdf$er_plot == 1], pch=20, col=rgb(0,0,0,.2), vertical=T, add=T, method="jitter")
  abline(lm(rdf[rdf$er_plot == 1, i] ~ rdf$tba_num[rdf$er_plot == 1]), lwd=2, lty=2, col=rgb(0,0,1,.5))
  dev.off()
}

