library(dplyr)
library(grid)
library(fields)
library(stringr)
library(mapplots)
library(reshape2)
library(data.table)
library(MutationTimeR)

# Copy-number synchronicity analysis
tmdf <- purdf[purdf$availability == "yes" & purdf$timer == "yes" & purdf$veclonal == "yes",] # purdf from 01_Creating_metadata
tmdf$nt.wgd <- NA
tmdf$nt.total <- NA
tmdf$time.wgd <- NA
tmdf$n.wgd <- NA
tmdf$n.all <- NA
tmdf$chr.wgd <- NA
tmdf$chr.all <- NA
tmdf$sd.wgd <- NA
tmdf$avg.ci <- NA
tmdf$sd.all <- NA

for (i in 1:nrow(tmdf)){
  df <- read.csv(paste0("../27_EBI/final_data_022222/timeR/wgd_synchronicity/", tmdf$study_id[i], ".wgd.synchronicity.txt"), header=T, as.is=T, sep = "\t")
  tmdf$nt.wgd[i] <- df$x[1]
  tmdf$nt.total[i] <- df$x[2]
  tmdf$time.wgd[i] <- df$x[3]
  tmdf$n.wgd[i] <- df$x[4]
  tmdf$n.all[i] <- df$x[5]
  tmdf$chr.wgd[i] <- df$x[6]
  tmdf$chr.all[i] <- df$x[7]
  tmdf$sd.wgd[i] <- df$x[8]
  tmdf$avg.ci[i] <- df$x[9]
  tmdf$sd.all[i] <- df$x[10]
}

## Synchronous vs asynchronous chromosomal copy gains -- we used the criteria used in the PCAWG analysis (reference: https://gerstung-lab.github.io/PCAWG-11/#9_synchronous_gains)
tmdf$info <- "uninformative"
tmdf$info[tmdf$avg.ci <= 0.5 & tmdf$chr.all > 2] <- "informative"
tmdf$timingclass <- "uninformative"
for (i in 1:nrow(tmdf)){
  if (tmdf$info[i] == "informative" & tmdf$nt.wgd[i]/tmdf$nt.total[i] > 0.75){
    tmdf$timingclass[i] <- "sync"
  } else if (tmdf$info[i] == "informative" & tmdf$nt.wgd[i]/tmdf$nt.total[i] <= 0.75){
    tmdf$timingclass[i] <- "async"
  }
}

tmdf$er_plot <- "unknown"
tmdf$er_plot <- factor(tmdf$er_plot, levels <- c("positive", "negative", "unknown"))
for (i in 1:nrow(tmdf)){
  tmdf$er_plot[i] <- demdf$er_plot[demdf$study_id==tmdf$study_id[i]]
}

plot(tmdf$ploidy[tmdf$availability=="yes"] ~ tmdf$haploidProportion[tmdf$availability=="yes"], col=rev(twocolors)[factor(tmdf$wholeGenomeDuplication[tmdf$availability=="yes"])], pch=20, las=1, xlab="Proportion of genome with LOH", ylab="Ploidy estimate", ylim=c(1,7.5), frame=F, main=paste0("Breast cancers (n=", nrow(tmdf[tmdf$availability=="yes",]), ")"), cex=0.4)
abline(a=2.9, b=-2) # PCAWG WGD classifier based on their function -- .classWgd <- function(ploidy, hom) 2.9 -2*hom <= ploidy
abline(a=2.9, b=-1.8)
abline(a=2.9, b=-1.6)
abline(a=2.85, b=-1.6)
tmdf$is.wgd <- ifelse(2.85 - 1.6*tmdf$haploidProportion <= tmdf$ploidy, "WGD", "ND") # Our classifier for 780 breast cancer cohort.
tmdf$study_id[tmdf$is.wgd == "ND" & tmdf$wholeGenomeDuplication == "true"] # Borderline cases.

## Summary chronicity
timingClass <- factor(paste(tmdf$is.wgd, tmdf$timingclass, sep=" "))
pdf("FIP_revision/pie.synchronicity.chrom.gains.expanded.pdf", height=4.5, width=4)

colTime <- c("#E3CF57","#FF9912","#FFE4C4","#CAE1FF","#104E8B","#5CACEE")
names(colTime) <- levels(timingClass)[c(4,5,6,3,2,1)]
c <- c(RColorBrewer::brewer.pal(9, "Pastel1"),"#DDDDDD")
t <- table(timingClass)[names(colTime)]
pie(t, init.angle=90, labels=paste0(names(t), ",\nn=", t), col=colTime)
par(new=TRUE)
symbols(x=0,y=0,circles=0.4, inches=FALSE, add=TRUE, bg="white")
dev.off()

## Timing of WGD
library(grDevices)
prgn <- RColorBrewer::brewer.pal(11,"PRGn")
set1 <- RColorBrewer::brewer.pal(9,"Set1")
colorRampPalette(set1[c(4,9,3)])(20)

tmdf$time.group <- floor(tmdf$time.wgd/0.05)
tmdf$er_plot <- "unknown"
tmdf$er_plot <- factor(tmdf$er_plot, levels <- c("positive", "negative", "unknown"))
for (i in 1:nrow(tmdf)){
  tmdf$er_plot[i] <- demdf$er_plot[demdf$study_id==tmdf$study_id[i]]
}

args <- c("WGD", "sync", "positive")
args <- c("WGD", "sync", "negative")
args <- c("ND", "sync", "positive")
args <- c("ND", "sync", "negative")

df <- tmdf[tmdf$is.wgd==args[1] & tmdf$timingclass==args[2] & tmdf$er_plot==args[3],]
tempdf <- as.data.frame(c(0:19))
tempdf$count <- 0
colnames(tempdf)[1] <- "time.group"
for (i in 1:nrow(tempdf)){
  tempdf$count[i] <- nrow(df[df$time.group == tempdf$time.group[i],])
}
pdf(paste0("FIP_revision/barplot.chrom.gain.timing.", args[1], ".", args[2], ".", args[3], ".expanded.pdf"), height=3, width=2.7)
barplot(tempdf$count, col=rev(colorRampPalette(set1[c(4,9,3)])(20)), border=NA, space=0, las=1, ylab="Number of patients", xlab="Relative time (fraction of mutations)", ylim=c(0,12), main=paste0("Synchronous chrom gains in ", args[1], " tumors, ER ", args[3]))
axis(side=1, at=seq(0,20,5), labels=c(0,0.25,0.5,0.75,1))
dev.off()





