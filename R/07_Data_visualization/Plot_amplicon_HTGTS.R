library(dplyr)
library(grid)
library(fields)
library(stringr)
library(plotrix)
library(mapplots)
library(viridis)
library(reshape2)

## Amplicon boundaries in ERBB2 neighborhood (Extended Data Fig. 8)
genedf[genedf$gene == "ERBB2",]
her2list <- unique(amphmf$study_id[amphmf$chromosome == "17" & amphmf$start <= 37856333 & amphmf$end >= 37884915])
hp <- hrdf$study_id[hrdf$study_id %in% her2list & hrdf$er_plot == "positive"]
hn <- hrdf$study_id[hrdf$study_id %in% her2list & hrdf$er_plot == "negative"]
hu <- hrdf$study_id[hrdf$study_id %in% her2list & hrdf$er_plot == "unknown"]

pdf("segments.ERBB2.amplified.regions.zoom.in.pdf", height=4.5, width=8)
posdf <- amphmf[amphmf$chromosome == "17" & amphmf$study_id %in% hp,]
posdf <- posdf[posdf$end >= ((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)-5000000 & posdf$start <= ((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)+5000000,]
posdf <- posdf[order(posdf$end, decreasing = F),]
orderdf <- posdf[posdf$chromosome == "17" & posdf$start <= 37856333 & posdf$end >= 37884915,]
tumorlist <- unique(orderdf$study_id)

plot(c(), xlim=c(((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)-5000000, ((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)+5000000), ylim = c(0,130), frame=F, las=1, xlab = "Positions in chromosome 17 (Mbp)", xaxt='n', ylab=paste0("ERBB2-amplified breast cancers (n=", length(her2list), ")"))
rect(genedf$start[genedf$gene=="ERBB2"], 0, genedf$end[genedf$gene=="ERBB2"], 130, col=rgb(165/255, 28/255, 48/255, 0.3), border = F)
rect(genedf$start[genedf$gene=="RARA"], 0, genedf$end[genedf$gene=="RARA"], 130, col=rgb(165/255, 28/255, 48/255, 0.3), border = F)
rect(genedf$start[genedf$gene=="CCR7"], 0, genedf$end[genedf$gene=="CCR7"], 130, col=rgb(165/255, 28/255, 48/255, 0.3), border = F)
abline(v=((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)-5000000)
abline(v=((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)+5000000)

axis(1, at=c(34000000, 36000000, 38000000, 40000000, 42000000), labels = c(34, 36, 38, 40, 42))
for (i in 1:length(tumorlist)){
  df <- posdf[posdf$study_id == tumorlist[i],]
  for (j in 1:nrow(df)){
    segments(df$start[j], 124-i, df$end[j], 124-i, col=rgb(39/255, 64/255, 139/255, 0.7), lwd=0.8, lend='butt')
  }
}
for (i in 1:nrow(orderdf)){
  if (orderdf$mcj_e[i] != "no"){
    if (orderdf$svclass_e[i] == "TRA"){
      points(orderdf$end[i], 124-i, pch=20, col=rgb(0,0,0,.3))
    } else if (orderdf$svclass_e[i] == "h2hINV" & orderdf$svlen_e[i] <= 5000){
      points(orderdf$end[i], 124-i, pch=17, cex=0.7, col=rgb(0,0,0,.3))
    }
  }
}

negdf <- amphmf[amphmf$chromosome == "17" & amphmf$study_id %in% hn,]
negdf <- negdf[negdf$end >= ((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)-5000000 & negdf$start <= ((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)+5000000,]
negdf <- negdf[order(negdf$end, decreasing = F),]
orderdf <- negdf[negdf$chromosome == "17" & negdf$start <= 37856333 & negdf$end >= 37884915,]
tumorlist <- unique(orderdf$study_id)

for (i in 1:length(tumorlist)){
  df <- negdf[negdf$study_id == tumorlist[i],]
  for (j in 1:nrow(df)){
    segments(df$start[j], 124-length(hp)-i-1, df$end[j], 124-length(hp)-i-1, col=rgb(139/255, 10/255, 80/255, 0.7), lwd=0.8, lend='butt')
  }
}
for (i in 1:nrow(orderdf)){
  if (orderdf$mcj_e[i] != "no"){
    if (orderdf$svclass_e[i] == "TRA"){
      points(orderdf$end[i], 124-length(hp)-i-1, pch=20, col=rgb(0,0,0,.3))
    } else if (orderdf$svclass_e[i] == "h2hINV" & orderdf$svlen_e[i] <= 5000){
      points(orderdf$end[i], 124-length(hp)-i-1, pch=17, cex=0.7, col=rgb(0,0,0,.3))
    }
  }
}

unkdf <- amphmf[amphmf$chromosome == "17" & amphmf$study_id %in% hu,]
unkdf <- unkdf[unkdf$end >= ((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)-5000000 & unkdf$start <= ((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)+5000000,]
unkdf <- unkdf[order(unkdf$end, decreasing = F),]
orderdf <- unkdf[unkdf$chromosome == "17" & unkdf$start <= 37856333 & unkdf$end >= 37884915,]
tumorlist <- unique(orderdf$study_id)

for (i in 1:length(tumorlist)){
  df <- unkdf[unkdf$study_id == tumorlist[i],]
  for (j in 1:nrow(df)){
    segments(df$start[j], 124-length(hp)-length(hn)-i-2, df$end[j], 124-length(hp)-length(hn)-i-2, col=rgb(0, 0, 0, 0.3), lwd=0.8, lend='butt')
  }
}
for (i in 1:nrow(orderdf)){
  if (orderdf$mcj_e[i] != "no"){
    if (orderdf$svclass_e[i] == "TRA"){
      points(orderdf$end[i], 124-length(hp)-length(hn)-i-2, pch=20, col=rgb(0,0,0,.3))
    } else if (orderdf$svclass_e[i] == "h2hINV" & orderdf$svlen_e[i] <= 5000){
      points(orderdf$end[i], 124-length(hp)-length(hn)-i-2, pch=17, cex=0.7, col=rgb(0,0,0,.3))
    }
  }
}

dev.off()


## Plotting the HTGTS breakpoints in the ERBB2 neighborhood (Extended Data Fig. 8)
library(beeswarm)
df <- read.csv("../28_HTGTS/MCF7_hSHANK2_withoutE2.hg19.bed", header=F, as.is=T, sep="\t")
df <- df[df$V1 == "chr17" & df$V2 >= ((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)-5000000 & df$V2 <= ((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)+5000000,]
df$V7 <- "control"
tempdf <- df

df <- read.csv("../28_HTGTS/T47D_hSHANK2_withoutE2.hg19.bed", header=F, as.is=T, sep="\t")
df <- df[df$V1 == "chr17" & df$V2 >= ((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)-5000000 & df$V2 <= ((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)+5000000,]
df$V7 <- "control"
tempdf <- rbind(tempdf, df)

df <- read.csv("../28_HTGTS/MCF7_hSHANK2_withE2.hg19.bed", header=F, as.is=T, sep="\t")
df <- df[df$V1 == "chr17" & df$V2 >= ((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)-5000000 & df$V2 <= ((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)+5000000,]
df$V7 <- "e2"
tempdf <- rbind(tempdf, df)

df <- read.csv("../28_HTGTS/T47D_hSHANK2_withE2.hg19.bed", header=F, as.is=T, sep="\t")
df <- df[df$V1 == "chr17" & df$V2 >= ((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)-5000000 & df$V2 <= ((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)+5000000,]
df$V7 <- "e2"
df <- rbind(tempdf, df)

pdf("beeswarm.htgts.junctions.erbb2.mcf7.t47d.pdf", height=8, width=3)
beeswarm(df$V2 ~ df$V7, method='hex', pwpch=ifelse(df$V6 == "+", 25, 24), pwbg=ifelse(df$V7 == "control", "#FFC1C1", "#CD3278"), pwcol=ifelse(df$V7 == "control", "#FFC1C1", "#CD3278"), las=3, pwcex=ifelse(df$V6 == "+", 0.5, 0.5), ylim=c(((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)-5000000,((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)+5000000), cex=0.4, xlab="", ylab="", yaxt='n')
abline(h=((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)-5000000)
abline(h=((genedf$start[genedf$gene=="ERBB2"]+genedf$end[genedf$gene=="ERBB2"])/2)+5000000)
axis(2, at=c(34000000, 36000000, 38000000, 40000000, 42000000), labels = c(34, 36, 38, 40, 42))
dev.off()

