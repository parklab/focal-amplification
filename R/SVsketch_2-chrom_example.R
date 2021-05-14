library(dplyr)
library(grid)
library(fields)
library(stringr)
options(scipen=999)

sumdf <- read.csv("Summaryinfo.table.1.19.txt", header=T, as.is=T, sep="\t")
cgcv90_coord <- read.csv("cgcv90_coord.txt", header=T, as.is=T, sep="\t")
hg19_coord <- read.csv("hg19_coord.txt", header=T, as.is=T, sep="\t")

i = "DO1281"

emphasis_col = rgb(138/255, 54/255, 15/255, 1)
tumorid <- sumdf$tumor[sumdf$icgc_donor_id==i]
caseid <- sumdf$submitted_specimen_id[sumdf$icgc_donor_id==i]

chri1 = "17"
chri2 = "8"
orientation <- c("+", "+")
geneofi = c("ERBB2", "CCND1", "VMP1-mir21", "MDM2", "MYC", "CCND2", "CCNE1", "RUNX1", "CDK4", "CDK6", "TP53", "PTEN", "RB1", "GNAS", "AR", "FGFR1", "EGFR", "PAK1")

inipos1 <- 0
inipos2 <- 0
endpos1 <- hg19_coord$length[hg19_coord$chr == chri1]
endpos2 <- hg19_coord$length[hg19_coord$chr == chri2]
#endpos1 <- 90000000
#endpos2 <- 40000000
delta <- ceiling(endpos1/20000000) * 20000000

# Combination Plot
pdf(paste0(i, "_", chri1, "_", chri2, "_", sumdf$histology_abbreviation[sumdf$icgc_donor_id == i], orientation[1], orientation[2], ".svsketch.ordered.10kb.pdf"), width=10, height=5)
par(mar=c(5.1,4.1,4.1,4.1)) # default

abscn <- read.csv("DO1281.depth.ratio.absCN.gc.correct.v2.txt", header=F, as.is=T, sep="\t")

#if (bcdf$project_code[bcdf$icgc_donor_id == i] == "BRCA-FR"){
#  abscn <- read.csv(paste0("~/Documents/SV_Enhancers/24_BRCA-FR_consensus/AbsCN_GC_Corrected/", i, ".depth.ratio.absCN.gc.correct.v2.txt"), header=F, as.is=T, sep="\t")
#} else {
#  abscn <- read.csv(paste0("~/Documents/SV_Enhancers/25_BRCA-PCAWG_abscn_final/", i, ".depth.ratio.absCN.gc.correct.v2.txt"), header=F, as.is=T, sep="\t")
#}

abscn$V2 <- abscn$V2 - 4999
abscn <- abscn[,c(1,2,3,4,6)]
colnames(abscn)[5] <- "V5"

abscn <- abscn[abscn$V1 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"),]
abscn <- abscn[!is.na(abscn$V5),]
abssub1 <- subset(abscn, abscn$V1 == chri1 & abscn$V2 >= inipos1 & abscn$V2 <= endpos1)
abssub2 <- subset(abscn, abscn$V1 == chri2 & abscn$V2 >= inipos2 & abscn$V2 <= endpos2)

if (orientation[1] == "-"){
  tempdf <- abssub1
  maxvalue <- max(tempdf$V2)
  tempdf$V6 <- maxvalue + 1 - tempdf$V2
  tempdf <- tempdf[order(tempdf$V6, decreasing = F),]
  tempdf <- tempdf[,c(1,6,3,4,5)]
  colnames(tempdf)[2] <- "V2"
  abssub1 <- tempdf
  value1 <- maxvalue
}
if (orientation[2] == "-"){
  tempdf <- abssub2
  maxvalue <- max(tempdf$V2)
  tempdf$V6 <- maxvalue + 1 - tempdf$V2
  tempdf <- tempdf[order(tempdf$V6, decreasing = F),]
  tempdf <- tempdf[,c(1,6,3,4,5)]
  colnames(tempdf)[2] <- "V2"
  abssub2 <- tempdf
  value2 <- maxvalue
}

absadj <- abssub2
absadj$V2 <- absadj$V2 + delta
absadj <- rbind(abssub1, absadj)

cadf <- read.csv("BreastCancer278.Ampdf.v2.full.txt", header=T, as.is=T, sep="\t")
#if (sumdf$histology_abbreviation[sumdf$icgc_donor_id == i] == "Breast-LobularCA"){
#  cadf <- read.csv(paste0("ampdf_PCAWG_for_figures/ampdf.Breast-AdenoCA.txt"), header=T, as.is=T, sep="\t")
#} else {
#  cadf <- read.csv(paste0("ampdf_PCAWG_for_figures/ampdf.", sumdf$histology_abbreviation[sumdf$icgc_donor_id == i], ".txt"), header=T, as.is=T, sep="\t")
#}

cadf <- cadf[cadf$icgc_donor_id==i & (cadf$chr == chri1 | cadf$chr == chri2),]
cadf$s[cadf$chr == chri2] <- cadf$s[cadf$chr == chri2] + delta
cadf$e[cadf$chr == chri2] <- cadf$e[cadf$chr == chri2] + delta

k = max(absadj$V5, na.rm=T) #k = integer #in case that the y axis should be specified
plot(0, type='n', xlim=c(inipos1, max(absadj$V2)), ylim=c(0,1.1*k), xlab=paste0("Position on chromosome ", chri1, " and ", chri2), ylab="", frame=F, xaxt = "n", yaxt = "n", las=1)
axis(1, at = seq(0, delta + floor(endpos2/20000000)*20000000, by= 20000000), labels = c(seq(0, (ceiling(endpos1/20000000)-1)*20, by = 20), seq(0, floor(endpos2/20000000)*20, by = 20)))
for (j in 1:nrow(cadf)){
  polygon(x = c(cadf$s[j], cadf$s[j], cadf$e[j], cadf$e[j]), y = c(0, max(absadj$V2), max(absadj$V2), 0), col = rgb(255/255, 153/255, 18/255, .3), border = NA)
}

## Plotting Genes of Interest
for (j in 1:length(geneofi)){
  if (geneofi[j] == "VMP1-mir21"){
    if (chri1 == "17"){
      segments(57784863, 1.15*k, 57918698, 1.15*k, lty=1, lwd=3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
    } else if (chri2 == "17"){
      segments(57784863 + delta, 1.15*k, 57918698 + delta, 1.15*k, lty=1, lwd=3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
    }
  } else if (geneofi[j] == "PAK1"){
    if (chri1 == "11"){
      segments(77033060, 1.15*k, 77185018, 1.15*k, lty=1, lwd=3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
    } else if (chri2 == "11"){
      segments(77033060 + delta, 1.15*k, 77185018 + delta, 1.15*k, lty=1, lwd=3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
    }
  } else if (cgcv90_coord$chr[cgcv90_coord$Gene.Symbol == geneofi[j]] == chri1){
    segments(cgcv90_coord$start[cgcv90_coord$Gene.Symbol == geneofi[j]], 1.15*k, cgcv90_coord$end[cgcv90_coord$Gene.Symbol == geneofi[j]], 1.15*k, lty=1, lwd=3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
  } else if (cgcv90_coord$chr[cgcv90_coord$Gene.Symbol == geneofi[j]] == chri2){
    segments(cgcv90_coord$start[cgcv90_coord$Gene.Symbol == geneofi[j]] + delta, 1.15*k, cgcv90_coord$end[cgcv90_coord$Gene.Symbol == geneofi[j]] + delta, 1.15*k, lty=1, lwd=3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
  }
}

## Importing SV Information and Plotting
#isv <- svall[svall$icgc_donor_id == i,]
#isv <- read.csv(paste0("BRCA206_AnnotatedSV/Old/", tumorid, ".pcawg_consensus_1.6.161116.somatic.sv.bedpe.annotated.ctxv3.txt"), header=T, as.is=T, sep="\t")
#isv <- read.csv(paste0("BRCA206_AnnotatedSV/ctx_jabba_intermediate2/", i, ".pcawg_consensus_1.6.161116.somatic.sv.bedpe.annotated.jabba.ctxv3.txt"), header=T, as.is=T, sep="\t")
#isv <- read.csv(paste0("Amp_boundary_fromO2_renewal/", sumdf$histology_abbreviation[sumdf$icgc_donor_id == i], "/SV_amp_boundaries/", i, ".amp_boundaries.txt"), header=T, as.is=T, sep="\t")

paired.read.support.threshold.v4<-function(file=file){
  t=read.table(file=file,header=T,sep="\t",stringsAsFactor=F);
  y1=sort(t$pe_support)
  x=(1:dim(t)[1])/(dim(t)[1])
  y=y1/(max(y1));
  #find the value where its targent value begins to continually larger than 1. 
  is=which(diff(y)/diff(x)>=1);
  i.h=is[which.max(unlist(lapply(is,function(i){
    a=y[i]-x[i];
    y2=x+a;
    sum(y-y2)})))[1]]
  thr=round(y1[i.h-1])
  #if there is no such value, return top 20% value.
  if(length(thr)==0){i.h=round(length(y1)*0.8);thr=round(y1[i.h])} 
  return(thr)
}

isv <- read.csv(paste0("DO1281.amp_boundaries.txt"), header=T, as.is=T, sep="\t")
cutoff <- paired.read.support.threshold.v4(paste0("DO1281.amp_boundaries.txt"))

#if (bcdf$project_code[bcdf$icgc_donor_id == i] == "BRCA-FR"){
#  isv <- read.csv(paste0("Breast-HER2/Result_1mb/SV_amp_boundaries/", i, ".amp_boundaries.txt"), header=T, as.is=T, sep="\t")
#  cutoff <- paired.read.support.threshold.v4(paste0("Breast-HER2/Result_1mb/SV_amp_boundaries/", i, ".amp_boundaries.txt"))
#} else {
#  isv <- read.csv(paste0("Breast/Result_1mb/SV_amp_boundaries/", i, ".amp_boundaries.txt"), header=T, as.is=T, sep="\t")
#  cutoff <- paired.read.support.threshold.v4(paste0("Breast/Result_1mb/SV_amp_boundaries/", i, ".amp_boundaries.txt"))
#}

#cutoff <- paired.read.support.threshold.v4(paste0("svboundary_for_figures/", i, ".amp_boundaries.txt"))

#isv <- read.csv(paste0("svboundary_for_figures/", i, ".amp_boundaries.txt"), header=T, as.is=T, sep="\t")
isv <- isv[isv$chrom1 %in% c(chri1, chri2) | isv$chrom2 %in% c(chri1, chri2),]

if (orientation[1] == "-"){
  isv$start1[isv$chrom1 == chri1] <- value1 - isv$start1[isv$chrom1 == chri1]
  isv$end1[isv$chrom1 == chri1] <- value1 - isv$end1[isv$chrom1 == chri1]
  isv$start2[isv$chrom2 == chri1] <- value1 - isv$start2[isv$chrom2 == chri1]
  isv$end2[isv$chrom2 == chri1] <- value1 - isv$end2[isv$chrom2 == chri1]
  isv$strand1[isv$chrom1 == chri1 & isv$strand1 == "+"] <- "negative"
  isv$strand1[isv$chrom1 == chri1 & isv$strand1 == "-"] <- "positive"
  isv$strand2[isv$chrom2 == chri1 & isv$strand1 == "+"] <- "negative"
  isv$strand2[isv$chrom2 == chri1 & isv$strand1 == "-"] <- "positive"
  isv$strand1[isv$chrom1 == chri1 & isv$strand1 == "positive"] <- "+"
  isv$strand1[isv$chrom1 == chri1 & isv$strand1 == "negative"] <- "-"
  isv$strand2[isv$chrom2 == chri1 & isv$strand1 == "positive"] <- "+"
  isv$strand2[isv$chrom2 == chri1 & isv$strand1 == "negative"] <- "-"
}
if (orientation[2] == "-"){
  isv$start1[isv$chrom1 == chri2] <- value2 - isv$start1[isv$chrom1 == chri2]
  isv$end1[isv$chrom1 == chri2] <- value2 - isv$end1[isv$chrom1 == chri2]
  isv$start2[isv$chrom2 == chri2] <- value2 - isv$start2[isv$chrom2 == chri2]
  isv$end2[isv$chrom2 == chri2] <- value2 - isv$end2[isv$chrom2 == chri2]
  isv$strand1[isv$chrom1 == chri2 & isv$strand1 == "+"] <- "negative"
  isv$strand1[isv$chrom1 == chri2 & isv$strand1 == "-"] <- "positive"
  isv$strand2[isv$chrom2 == chri2 & isv$strand1 == "+"] <- "negative"
  isv$strand2[isv$chrom2 == chri2 & isv$strand1 == "-"] <- "positive"
  isv$strand1[isv$chrom1 == chri2 & isv$strand1 == "positive"] <- "+"
  isv$strand1[isv$chrom1 == chri2 & isv$strand1 == "negative"] <- "-"
  isv$strand2[isv$chrom2 == chri2 & isv$strand1 == "positive"] <- "+"
  isv$strand2[isv$chrom2 == chri2 & isv$strand1 == "negative"] <- "-"
}

isv$start1[isv$chrom1 == chri2] <- isv$start1[isv$chrom1 == chri2] + delta
isv$end1[isv$chrom1 == chri2] <- isv$end1[isv$chrom1 == chri2] + delta
isv$start2[isv$chrom2 == chri2] <- isv$start2[isv$chrom2 == chri2] + delta
isv$end2[isv$chrom2 == chri2] <- isv$end2[isv$chrom2 == chri2] + delta
#isv <- isv[,c(1:19)]

isvintra <- subset(isv, isv$chrom1 == isv$chrom2)
isvintra$group[isvintra$strand1 == "+" & isvintra$strand2 == "-"] = 1 ## DELETION
isvintra$group[isvintra$strand1 == "-" & isvintra$strand2 == "+"] = 2 ## DUPLICATION
isvintra$group[isvintra$strand1 == "+" & isvintra$strand2 == "+"] = 3 ## HEAD-TO-HEAD INVERSION
isvintra$group[isvintra$strand1 == "-" & isvintra$strand2 == "-"] = 4 ## TAIL-TO-TAIL INVERSION
isvinter <- subset(isv, isv$chrom1 != isv$chrom2)
isvinterA <- subset(isvinter, isvinter$chrom1 %in% c(chri1, chri2) & isvinter$chrom2 %in% c(chri1, chri2))
isvinterB <- subset(isvinter, !(isvinter$chrom1 %in% c(chri1, chri2) & isvinter$chrom2 %in% c(chri1, chri2)))

par(new = T)
theta=seq(0,pi, len=100)
isvintra$rad=abs(isvintra$start2-isvintra$start1)/2

deldf <- isvintra[isvintra$group==1,]
dupdf <- isvintra[isvintra$group==2,]
hhidf <- isvintra[isvintra$group==3,]
ttidf <- isvintra[isvintra$group==4,]

w = max(isv$pe_support)
plot(0, type='n', xlim=c(inipos1, max(absadj$V2)), ylim=c(0,1.1*w), frame=F, axes = F, xlab = "", ylab = "")
axis(4, las=1)
mtext("Read support", side=4, line=3)

## DELETIONS
if (length(isvintra$sv_id[isvintra$group==1])!=0){
  for (j in seq(1,length(unique(deldf$sv_id)),1)){
    if (deldf$amp_boundary[j] == "no"){
      svcol <- rgb(0/255,0/255,255/255,.3)
    } else {
      svcol <- emphasis_col
    }
    x = deldf$rad[j]*cos(theta)+((deldf$start1[j]+deldf$start2[j])/2)
    y = w*0.1*sin(theta)+deldf$pe_support[j]
    lines(x,y,col=svcol, xpd=TRUE)
    segments(deldf$start1[j], 0, deldf$start1[j], deldf$pe_support[j], lty=1, xpd=TRUE, col=svcol)
    segments(deldf$start2[j], 0, deldf$start2[j], deldf$pe_support[j], lty=1, xpd=TRUE, col=svcol)
  }
}

## DUPLICATIONS
if (length(isvintra$sv_id[isvintra$group==2])!=0){
  for (j in seq(1,length(unique(dupdf$sv_id)),1)){
    if (dupdf$amp_boundary[j] == "no"){
      svcol <- rgb(0/255,128/255,128/255,.3)
    } else {
      svcol <- emphasis_col
    }
    x = dupdf$rad[j]*cos(theta)+((dupdf$start1[j]+dupdf$start2[j])/2)
    y = -w*0.1*sin(theta)+dupdf$pe_support[j]
    lines(x,y,col=svcol, xpd=TRUE)
    segments(dupdf$start1[j], 0, dupdf$start1[j], dupdf$pe_support[j], lty=1, xpd=TRUE, col=svcol)
    segments(dupdf$start2[j], 0, dupdf$start2[j], dupdf$pe_support[j], lty=1, xpd=TRUE, col=svcol)
  }
}

## HEAD-TO-HEAD INVERSIONS
if (length(isvintra$sv_id[isvintra$group==3])!=0){
  for (j in seq(1,length(unique(hhidf$sv_id)),1)){
    if (hhidf$amp_boundary[j] == "no"){
      svcol <- rgb(220/255,20/255,60/255,.3)
    } else {
      svcol <- emphasis_col
    }
    x = hhidf$rad[j]*cos(theta)+((hhidf$start1[j]+hhidf$start2[j])/2)
    y = w*0.1*sin(theta)+hhidf$pe_support[j]
    lines(x,y,col=svcol, xpd=TRUE)
    segments(hhidf$start1[j], 0, hhidf$start1[j], hhidf$pe_support[j], lty=1, xpd=TRUE, col=svcol)
    segments(hhidf$start2[j], 0, hhidf$start2[j], hhidf$pe_support[j], lty=1, xpd=TRUE, col=svcol)
  }
}

## TAIL-TO-TAIL INVERSIONS
if (length(isvintra$sv_id[isvintra$group==4])!=0){
  for (j in seq(1,length(unique(ttidf$sv_id)),1)){
    if (ttidf$amp_boundary[j] == "no"){
      svcol <- rgb(128/255,128/255,0/255,.3)
    } else {
      svcol <- emphasis_col
    }
    x = ttidf$rad[j]*cos(theta)+((ttidf$start1[j]+ttidf$start2[j])/2)
    y = -w*0.1*sin(theta)+ttidf$pe_support[j]
    lines(x,y,col=svcol, xpd=TRUE)
    segments(ttidf$start1[j], 0, ttidf$start1[j], ttidf$pe_support[j], lty=1, xpd=TRUE, col=svcol)
    segments(ttidf$start2[j], 0, ttidf$start2[j], ttidf$pe_support[j], lty=1, xpd=TRUE, col=svcol)
  }
}

## INTERCHROMOSOMAL TRANSLOCATIONS
isvinterA$rad=abs(isvinterA$start2-isvinterA$start1)/2
if (nrow(isvinterA) != 0){
  for (j in seq(1,nrow(isvinterA),1)){
    if (isvinterA$amp_boundary[j] == "boundary" & isvinterA$pe_support[j] > cutoff){
      svcol <- emphasis_col
    } else {
      svcol <- rgb(85/255,26/255,139/255,.3)
    }
    x = isvinterA$rad[j]*cos(theta)+((isvinterA$start1[j]+isvinterA$start2[j])/2)
    y = w*0.1*sin(theta)+isvinterA$pe_support[j]
    lines(x,y,col=svcol, xpd=TRUE)
    segments(isvinterA$start1[j], 0, isvinterA$start1[j], isvinterA$pe_support[j], lty=1, xpd=TRUE, col=svcol)
    segments(isvinterA$start2[j], 0, isvinterA$start2[j], isvinterA$pe_support[j], lty=1, xpd=TRUE, col=svcol)
  }
}

## INTERCHROMOSOMAL TRANSLOCATIONS OUT OF THE BRIDGE
if (nrow(isvinterB) != 0){
  for (j in seq(1,nrow(isvinterB),1)){
    if (isvinterB$amp_boundary[j] == "boundary" & isvinterB$pe_support[j] > cutoff){
      svcol <- emphasis_col
    } else {
      svcol <- rgb(85/255,26/255,139/255,.3)
    }
    if (isvinterB$chrom1[j] %in% c(chri1, chri2)){
      arrows(isvinterB$start1[j], 0, isvinterB$start1[j], isvinterB$pe_support[j], lty=1, xpd=TRUE, col=svcol, code=2, length = 0.1)
    } else {
      arrows(isvinterB$start2[j], 0, isvinterB$start2[j], isvinterB$pe_support[j], lty=1, xpd=TRUE, col=svcol, code=2, length = 0.1)
    }
  }
}

par(new = T)
k = max(absadj$V5, na.rm=T) #k = integer #in case that the y axis should be specified
plot(0, type='n', xlim=c(inipos1, max(absadj$V2)), ylim=c(0,1.1*k), xlab=paste0("Position on chromosome ", chri1, " and ", chri2), ylab="Absolute CN", frame=F, xaxt = "n", las=1)
points(absadj$V2, absadj$V5, pch=20, col=rgb(0,0,0,.5))

dev.off()
