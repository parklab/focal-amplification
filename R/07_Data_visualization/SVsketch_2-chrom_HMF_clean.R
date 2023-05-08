library(dplyr)
library(grid)
library(fields)
library(stringr)
options(scipen=999)

# Plotting copy number and structural variations for the two chromosomes undergoing translocation-bridge amplification
## Amplified SVs
paired.read.support.threshold.hmf<-function(file=file){
  t=read.table(file=file,header=T,sep="\t",stringsAsFactor=F);
  y1=sort(t$vf)
  x=(1:dim(t)[1])/(dim(t)[1])
  y=y1/(max(y1));
  is=which(diff(y)/diff(x)>=2);
  i.h=is[which.max(unlist(lapply(is,function(i){
    a=y[i]-x[i];
    y2=x+a;
    sum(y-y2)})))[1]]
  thr=round(y1[i.h-1])
  if(length(thr)==0){i.h=round(length(y1)*0.8);thr=round(y1[i.h])} 
  return(thr)
}

## Essential files
cgcv90_coord <- read.csv("cancer_gene_census_v90_coordinates.csv", header=T, as.is=T)
hg19_coord <- read.csv("hg19chr.coordinates.csv", header=T, as.is=T)
amphmf <- read.csv("List.amplicons.summary.allcohort.str.fbi.mech.031022.txt", header=T, as.is=T, sep="\t")

## Gene of interest
geneofi <- c("ERBB2", "CCND1", "FGF3", "RSF1", "PAK1", "ZNF703", "FGFR1", "MYC", "TRPS1", "PPM1D", "USP32", "ZNF217", "MDM2", "MCL1", "PTEN", "CDKN2A", "RB1", "TP53", "BRCA1", "BRCA2", "PALB2", "MAP3K1", "MAP2K4", "PIK3CA", "PIK3R1", "AKT1", "GATA3", "RUNX1", "CBFB", "SF3B1", "CDH1", "KMT2C", "KMT2D", "ARID1A", "SETD2", "QRSL1", "IGF1R")

## Colors
emphasis_col = rgb(138/255, 54/255, 15/255, 1)

## Output directory
svsketchoutputdir <- "" ## Output directory

## Function -- SVsketch
svsketch.hmf <- function(i, chri1, chri2, ori1, ori2, pwt, pht){
  inipos1 <- 0
  inipos2 <- 0
  endpos1 <- hg19_coord$length[hg19_coord$chr == chri1]
  endpos2 <- hg19_coord$length[hg19_coord$chr == chri2]
  delta <- ceiling(endpos1/20000000) * 20000000
  
  cnvdf <- read.csv(paste0("../27_EBI/final_data_022222/acgr/", i, ".purple.cnv.somatic.acgr.tsv"), header=T, as.is=T, sep="\t")
  cnv1 <- subset(cnvdf, cnvdf$chromosome == chri1 & cnvdf$start >= inipos1 & cnvdf$end <= endpos1)
  cnv1 <- cnv1[,c("chromosome", "start", "end", "majorAlleleCopyNumber", "minorAlleleCopyNumber", "loh", "acgr")]
  cnv2 <- subset(cnvdf, cnvdf$chromosome == chri2 & cnvdf$start >= inipos2 & cnvdf$end <= endpos2)
  cnv2 <- cnv2[,c("chromosome", "start", "end", "majorAlleleCopyNumber", "minorAlleleCopyNumber", "loh", "acgr")]
  
  if (ori1 == "-"){
    tempdf <- cnv1
    maxvalue <- max(tempdf$end)
    tempdf$rev_start <- maxvalue + 1 - tempdf$end
    tempdf$rev_end <- maxvalue + 1 - tempdf$start
    tempdf$acgr[tempdf$acgr == "amplicon_start"] <- "tail"
    tempdf$acgr[tempdf$acgr == "amplicon_end"] <- "head"
    tempdf$acgr[tempdf$acgr == "tail"] <- "amplicon_end"
    tempdf$acgr[tempdf$acgr == "head"] <- "amplicon_start"
    tempdf <- tempdf[order(tempdf$rev_start, decreasing = F),]
    tempdf <- tempdf[,c("chromosome", "rev_start", "rev_end", "majorAlleleCopyNumber", "minorAlleleCopyNumber", "loh", "acgr")]
    colnames(tempdf)[2:3] <- c("start", "end")
    cnv1 <- tempdf
    value1 <- maxvalue
  }
  if (ori2 == "-"){
    tempdf <- cnv2
    maxvalue <- max(tempdf$end)
    tempdf$rev_start <- maxvalue + 1 - tempdf$end
    tempdf$rev_end <- maxvalue + 1 - tempdf$start
    tempdf$acgr[tempdf$acgr == "amplicon_start"] <- "tail"
    tempdf$acgr[tempdf$acgr == "amplicon_end"] <- "head"
    tempdf$acgr[tempdf$acgr == "tail"] <- "amplicon_end"
    tempdf$acgr[tempdf$acgr == "head"] <- "amplicon_start"
    tempdf <- tempdf[order(tempdf$rev_start, decreasing = F),]
    tempdf <- tempdf[,c("chromosome", "rev_start", "rev_end", "majorAlleleCopyNumber", "minorAlleleCopyNumber", "loh", "acgr")]
    colnames(tempdf)[2:3] <- c("start", "end")
    cnv2 <- tempdf
    value2 <- maxvalue
  }
  
  cnv2$start <- cnv2$start + delta
  cnv2$end <- cnv2$end + delta
  cnvadj <- rbind(cnv1, cnv2)
  
  ### Importing SV Information and Plotting
  isv <- read.csv(paste0("../27_EBI/final_data_022222/rearrangement_bedpe/",  i, ".purple.sv.vcf.gz.bedpe"), header=T, as.is=T, sep="\t")
  isv$boundary <- "no"
  if (nrow(amphmf[amphmf$study_id == i,]) != 0){
    boundary_list <- unique(c(amphmf$svname_s[amphmf$study_id == i & amphmf$mcj_s != "no"], amphmf$svname_e[amphmf$study_id == i & amphmf$mcj_e != "no"]))
    isv$boundary[isv$name %in% boundary_list] <- "yes"
  }
  isv <- isv[isv$chrom1 %in% c(chri1, chri2) | isv$chrom2 %in% c(chri1, chri2),]
  isv <- isv[!is.na(isv$vf),]
  if (nrow(isv) == 0){
    return(NULL)
  }
  
  if (ori1 == "-"){
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
  if (ori2 == "-"){
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
  
  isvintra <- subset(isv, isv$chrom1 == isv$chrom2)
  isvintra$group[isvintra$strand1 == "+" & isvintra$strand2 == "-"] = 1 ## DELETION
  isvintra$group[isvintra$strand1 == "-" & isvintra$strand2 == "+"] = 2 ## DUPLICATION
  isvintra$group[isvintra$strand1 == "+" & isvintra$strand2 == "+"] = 3 ## HEAD-TO-HEAD INVERSION
  isvintra$group[isvintra$strand1 == "-" & isvintra$strand2 == "-"] = 4 ## TAIL-TO-TAIL INVERSION
  isvinter <- subset(isv, isv$chrom1 != isv$chrom2)
  isvinterA <- subset(isvinter, isvinter$chrom1 %in% c(chri1, chri2) & isvinter$chrom2 %in% c(chri1, chri2))
  isvinterB <- subset(isvinter, !(isvinter$chrom1 %in% c(chri1, chri2) & isvinter$chrom2 %in% c(chri1, chri2)))
  
  
  ### Combination Plot
  pdf(paste0(svsketchoutputdir, i, "_", chri1, "_", chri2, "_", ori1, ori2, ".svsketch.pdf"), width=pwt, height=pht) # All preliminary figures were width 10 height 5
  par(mar=c(5.1,4.1,4.1,4.1)) # default
  
  theta=seq(0,pi, len=100)
  isvintra$rad=abs(isvintra$start2-isvintra$start1)/2
  
  deldf <- isvintra[isvintra$group==1,]
  dupdf <- isvintra[isvintra$group==2,]
  hhidf <- isvintra[isvintra$group==3,]
  ttidf <- isvintra[isvintra$group==4,]
  
  w = max(isv$vf)
  plot(0, type='n', xlim=c(inipos1, max(cnvadj$end)), ylim=c(0,1.1*w), frame=F, axes = F, xlab = "", ylab = "")
  axis(4, las=1)
  mtext("Read support", side=4, line=3)
  
  ### DELETIONS
  if (length(isvintra$name[isvintra$group==1])!=0){
    for (j in seq(1,length(unique(deldf$name)),1)){
      if (deldf$boundary[j] == "yes"){
        svcol <- emphasis_col
      } else {
        svcol <- rgb(0/255,0/255,255/255,.3)
      }
      x = deldf$rad[j]*cos(theta)+((deldf$start1[j]+deldf$start2[j])/2)
      y = w*0.1*sin(theta)+deldf$vf[j]
      lines(x,y,col=svcol, xpd=TRUE, lwd=2/3)
      segments(deldf$start1[j], 0, deldf$start1[j], deldf$vf[j], lty=1, xpd=TRUE, col=svcol, lwd=2/3)
      segments(deldf$start2[j], 0, deldf$start2[j], deldf$vf[j], lty=1, xpd=TRUE, col=svcol, lwd=2/3)
    }
  }
  
  ### DUPLICATIONS
  if (length(isvintra$name[isvintra$group==2])!=0){
    for (j in seq(1,length(unique(dupdf$name)),1)){
      if (dupdf$boundary[j] == "yes"){
        svcol <- emphasis_col
      } else {
        svcol <- rgb(0/255,128/255,128/255,.3)
      }
      x = dupdf$rad[j]*cos(theta)+((dupdf$start1[j]+dupdf$start2[j])/2)
      y = -w*0.1*sin(theta)+dupdf$vf[j]
      lines(x,y,col=svcol, xpd=TRUE, lwd=2/3)
      segments(dupdf$start1[j], 0, dupdf$start1[j], dupdf$vf[j], lty=1, xpd=TRUE, col=svcol, lwd=2/3)
      segments(dupdf$start2[j], 0, dupdf$start2[j], dupdf$vf[j], lty=1, xpd=TRUE, col=svcol, lwd=2/3)
    }
  }
  
  ### HEAD-TO-HEAD INVERSIONS
  if (length(isvintra$name[isvintra$group==3])!=0){
    for (j in seq(1,length(unique(hhidf$name)),1)){
      if (hhidf$boundary[j] == "yes"){
        svcol <- emphasis_col
      } else {
        svcol <- rgb(220/255,20/255,60/255,.3)
      }
      x = hhidf$rad[j]*cos(theta)+((hhidf$start1[j]+hhidf$start2[j])/2)
      y = w*0.1*sin(theta)+hhidf$vf[j]
      lines(x,y,col=svcol, xpd=TRUE, lwd=2/3)
      segments(hhidf$start1[j], 0, hhidf$start1[j], hhidf$vf[j], lty=1, xpd=TRUE, col=svcol, lwd=2/3)
      segments(hhidf$start2[j], 0, hhidf$start2[j], hhidf$vf[j], lty=1, xpd=TRUE, col=svcol, lwd=2/3)
    }
  }
  
  ### TAIL-TO-TAIL INVERSIONS
  if (length(isvintra$name[isvintra$group==4])!=0){
    for (j in seq(1,length(unique(ttidf$name)),1)){
      if (ttidf$boundary[j] == "yes"){
        svcol <- emphasis_col
      } else {
        svcol <- rgb(128/255,128/255,0/255,.3)
      }
      x = ttidf$rad[j]*cos(theta)+((ttidf$start1[j]+ttidf$start2[j])/2)
      y = -w*0.1*sin(theta)+ttidf$vf[j]
      lines(x,y,col=svcol, xpd=TRUE, lwd=2/3)
      segments(ttidf$start1[j], 0, ttidf$start1[j], ttidf$vf[j], lty=1, xpd=TRUE, col=svcol, lwd=2/3)
      segments(ttidf$start2[j], 0, ttidf$start2[j], ttidf$vf[j], lty=1, xpd=TRUE, col=svcol, lwd=2/3)
    }
  }
  
  ### INTERCHROMOSOMAL TRANSLOCATIONS
  isvinterA$rad=abs(isvinterA$start2-isvinterA$start1)/2
  if (nrow(isvinterA) != 0){
    for (j in seq(1,nrow(isvinterA),1)){
      if (isvinterA$boundary[j] == "yes"){
        svcol <- emphasis_col
      } else {
        svcol <- rgb(85/255,26/255,139/255,.3)
      }
      x = isvinterA$rad[j]*cos(theta)+((isvinterA$start1[j]+isvinterA$start2[j])/2)
      y = w*0.1*sin(theta)+isvinterA$vf[j]
      lines(x,y,col=svcol, xpd=TRUE, lwd=2/3)
      segments(isvinterA$start1[j], 0, isvinterA$start1[j], isvinterA$vf[j], lty=1, xpd=TRUE, col=svcol, lwd=2/3)
      segments(isvinterA$start2[j], 0, isvinterA$start2[j], isvinterA$vf[j], lty=1, xpd=TRUE, col=svcol, lwd=2/3)
    }
  }
  
  ### INTERCHROMOSOMAL TRANSLOCATIONS OUT OF THE BRIDGE
  if (nrow(isvinterB) != 0){
    for (j in seq(1,nrow(isvinterB),1)){
      if (isvinterB$boundary[j] == "yes"){
        svcol <- emphasis_col
      } else {
        svcol <- rgb(85/255,26/255,139/255,.3)
      }
      if (isvinterB$chrom1[j] %in% c(chri1, chri2)){
        arrows(isvinterB$start1[j], 0, isvinterB$start1[j], isvinterB$vf[j], lty=1, xpd=TRUE, col=svcol, code=2, length = 0.1, lwd=2/3)
      } else {
        arrows(isvinterB$start2[j], 0, isvinterB$start2[j], isvinterB$vf[j], lty=1, xpd=TRUE, col=svcol, code=2, length = 0.1, lwd=2/3)
      }
    }
  }
  
  par(new = T)
  
  k <- max(cnvadj$majorAlleleCopyNumber)
  plot(0,type='n', xlim=c(inipos1, max(cnvadj$end)), ylim=c(0,1.1*k), frame=F, xaxt="n", las=1, xlab=paste0("Position on chromosome ", chri1, " and ", chri2), ylab="Allele-specific copy number")
  axis(1, at = seq(0, delta + floor(endpos2/20000000)*20000000, by= 20000000), labels = c(seq(0, (ceiling(endpos1/20000000)-1)*20, by = 20), seq(0, floor(endpos2/20000000)*20, by = 20)))
  segments(cnvadj$start[cnvadj$start != cnvadj$end], cnvadj$majorAlleleCopyNumber[cnvadj$start != cnvadj$end], cnvadj$end[cnvadj$start != cnvadj$end], cnvadj$majorAlleleCopyNumber[cnvadj$start != cnvadj$end], col="red", lwd=4/3)
  segments(cnvadj$start[cnvadj$start != cnvadj$end], cnvadj$minorAlleleCopyNumber[cnvadj$start != cnvadj$end]-0.2, cnvadj$end[cnvadj$start != cnvadj$end], cnvadj$minorAlleleCopyNumber[cnvadj$start != cnvadj$end]-0.2, col="blue", lwd=4/3)
  if (nrow(cnvadj[cnvadj$acgr == "amplicon_start",]) != 0){
    for (j in 1:nrow(cnvadj[cnvadj$acgr == "amplicon_start",])){
      polygon(x = c(cnvadj$start[which(cnvadj$acgr == "amplicon_start")][j], cnvadj$start[which(cnvadj$acgr == "amplicon_start")][j], cnvadj$end[which(cnvadj$acgr == "amplicon_end")][j], cnvadj$end[which(cnvadj$acgr == "amplicon_end")][j]), y=c(0, 1.1*k, 1.1*k, 0), col=rgb(255/255, 153/255, 18/255, .3), border = NA)
    }
  }
  if (nrow(cnvadj[cnvadj$acgr == "amplicon_itself",]) != 0){
    for (j in 1:nrow(cnvadj[cnvadj$acgr == "amplicon_itself",])){
      polygon(x = c(cnvadj$start[which(cnvadj$acgr == "amplicon_itself")][j], cnvadj$start[which(cnvadj$acgr == "amplicon_itself")][j], cnvadj$end[which(cnvadj$acgr == "amplicon_itself")][j], cnvadj$end[which(cnvadj$acgr == "amplicon_itself")][j]), y=c(0, 1.1*k, 1.1*k, 0), col=rgb(255/255, 153/255, 18/255, .3), border = NA)
    }
  }
  
  ### Plotting Genes of Interest
  for (j in 1:length(geneofi)){
    if (geneofi[j] == "VMP1-mir21"){
      if (chri1 == "17" & ori1 == "+"){
        segments(57784863, 1.15*k, 57918698, 1.15*k, lty=1, lwd=4/3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
      } else if (chri1 == "17" & ori1 == "-"){
        segments(hg19_coord$length[hg19_coord$chr == chri1] - 57918698, 1.15*k, hg19_coord$length[hg19_coord$chr == chri1] - 57784863, 1.15*k, lty=1, lwd=4/3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
      } else if (chri2 == "17" & ori2 == "+"){
        segments(57784863 + delta, 1.15*k, 57918698 + delta, 1.15*k, lty=1, lwd=4/3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
      } else if (chri2 == "17" & ori2 == "-"){
        segments(hg19_coord$length[hg19_coord$chr == chri2] + delta - 57918698, 1.15*k, hg19_coord$length[hg19_coord$chr == chri2] + delta - 57784863, 1.15*k, lty=1, lwd=4/3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
      }
    } else if (geneofi[j] == "PAK1"){
      if (chri1 == "11" & ori1 == "+"){
        segments(77033060, 1.15*k, 77185018, 1.15*k, lty=1, lwd=4/3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
      } else if (chri1 == "11" & ori1 == "-"){
        segments(hg19_coord$length[hg19_coord$chr == chri1] - 77185018, 1.15*k, hg19_coord$length[hg19_coord$chr == chri1] - 77033060, 1.15*k, lty=1, lwd=4/3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
      } else if (chri2 == "11" & ori2 == "+"){
        segments(77033060 + delta, 1.15*k, 77185018 + delta, 1.15*k, lty=1, lwd=4/3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
      } else if (chri2 == "11" & ori2 == "-"){
        segments(hg19_coord$length[hg19_coord$chr == chri2] + delta - 77185018, 1.15*k, hg19_coord$length[hg19_coord$chr == chri2] + delta - 77033060, 1.15*k, lty=1, lwd=4/3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
      }
    } else if (genedf$chromosome[genedf$gene == geneofi[j]] == chri1 & ori1 == "+"){
      segments(genedf$start[genedf$gene == geneofi[j]], 1.15*k, genedf$end[genedf$gene == geneofi[j]], 1.15*k, lty=1, lwd=4/3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
    } else if (genedf$chromosome[genedf$gene == geneofi[j]] == chri1 & ori1 == "-"){
      segments(hg19_coord$length[hg19_coord$chr == chri1] - genedf$end[genedf$gene == geneofi[j]], 1.15*k, hg19_coord$length[hg19_coord$chr == chri1] - genedf$start[genedf$gene == geneofi[j]], 1.15*k, lty=1, lwd=4/3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
    } else if (genedf$chromosome[genedf$gene == geneofi[j]] == chri2 & ori2 == "+"){
      segments(genedf$start[genedf$gene == geneofi[j]] + delta, 1.15*k, genedf$end[genedf$gene == geneofi[j]] + delta, 1.15*k, lty=1, lwd=4/3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
    } else if (genedf$chromosome[genedf$gene == geneofi[j]] == chri2 & ori2 == "-"){
      segments(hg19_coord$length[hg19_coord$chr == chri2] + delta - genedf$end[genedf$gene == geneofi[j]], 1.15*k, hg19_coord$length[hg19_coord$chr == chri2] + delta - genedf$start[genedf$gene == geneofi[j]], 1.15*k, lty=1, lwd=4/3, xpd=TRUE, col=rgb(0/255, 0/255, 255/255, 1))
    }
  }
  
  ### Plotting Centromere
  if (ori1 == "+"){
    points(hg19_coord$centro_mid[hg19_coord$chr == chri1], (-0.039)*k, pch=19, col=rgb(0,0,0,0.5))
  } else if (ori1 == "-"){
    points(hg19_coord$length[hg19_coord$chr == chri1] - hg19_coord$centro_mid[hg19_coord$chr == chri1], (-0.039)*k, pch=19, col=rgb(0,0,0,0.5))
  }
  if (ori2 == "+"){
    points(hg19_coord$centro_mid[hg19_coord$chr == chri2] + delta, (-0.039)*k, pch=19, col=rgb(0,0,0,0.5))
  } else if (ori2 == "-"){
    points(hg19_coord$length[hg19_coord$chr == chri2] + delta - hg19_coord$centro_mid[hg19_coord$chr == chri2], (-0.039)*k, pch=19, col=rgb(0,0,0,0.5))
  }
  
  dev.off()
}

## LOH plot
lohplot.hmf <- function(i, chri1, chri2, ori1, ori2, pwt){
  inipos1 <- 0
  inipos2 <- 0
  endpos1 <- hg19_coord$length[hg19_coord$chr == chri1]
  endpos2 <- hg19_coord$length[hg19_coord$chr == chri2]
  delta <- ceiling(endpos1/20000000) * 20000000
  
  cnvdf <- read.csv(paste0("../27_EBI/final_data_022222/acgr/", i, ".purple.cnv.somatic.acgr.tsv"), header=T, as.is=T, sep="\t")
  cnv1 <- subset(cnvdf, cnvdf$chromosome == chri1 & cnvdf$start >= inipos1 & cnvdf$end <= endpos1)
  cnv1 <- cnv1[,c("chromosome", "start", "end", "majorAlleleCopyNumber", "minorAlleleCopyNumber", "loh", "acgr")]
  cnv2 <- subset(cnvdf, cnvdf$chromosome == chri2 & cnvdf$start >= inipos2 & cnvdf$end <= endpos2)
  cnv2 <- cnv2[,c("chromosome", "start", "end", "majorAlleleCopyNumber", "minorAlleleCopyNumber", "loh", "acgr")]
  
  if (ori1 == "-"){
    tempdf <- cnv1
    maxvalue <- max(tempdf$end)
    tempdf$rev_start <- maxvalue + 1 - tempdf$end
    tempdf$rev_end <- maxvalue + 1 - tempdf$start
    tempdf$acgr[tempdf$acgr == "amplicon_start"] <- "tail"
    tempdf$acgr[tempdf$acgr == "amplicon_end"] <- "head"
    tempdf$acgr[tempdf$acgr == "tail"] <- "amplicon_end"
    tempdf$acgr[tempdf$acgr == "head"] <- "amplicon_start"
    tempdf <- tempdf[order(tempdf$rev_start, decreasing = F),]
    tempdf <- tempdf[,c("chromosome", "rev_start", "rev_end", "majorAlleleCopyNumber", "minorAlleleCopyNumber", "loh", "acgr")]
    colnames(tempdf)[2:3] <- c("start", "end")
    cnv1 <- tempdf
    value1 <- maxvalue
  }
  if (ori2 == "-"){
    tempdf <- cnv2
    maxvalue <- max(tempdf$end)
    tempdf$rev_start <- maxvalue + 1 - tempdf$end
    tempdf$rev_end <- maxvalue + 1 - tempdf$start
    tempdf$acgr[tempdf$acgr == "amplicon_start"] <- "tail"
    tempdf$acgr[tempdf$acgr == "amplicon_end"] <- "head"
    tempdf$acgr[tempdf$acgr == "tail"] <- "amplicon_end"
    tempdf$acgr[tempdf$acgr == "head"] <- "amplicon_start"
    tempdf <- tempdf[order(tempdf$rev_start, decreasing = F),]
    tempdf <- tempdf[,c("chromosome", "rev_start", "rev_end", "majorAlleleCopyNumber", "minorAlleleCopyNumber", "loh", "acgr")]
    colnames(tempdf)[2:3] <- c("start", "end")
    cnv2 <- tempdf
    value2 <- maxvalue
  }
  
  cnv2$start <- cnv2$start + delta
  cnv2$end <- cnv2$end + delta
  cnvadj <- rbind(cnv1, cnv2)
  
  pdf(paste0(svsketchoutputdir, i, "_", chri1, "_", chri2, "_", ori1, ori2, ".lohplot.pdf"), width=pwt, height=2.5) # All preliminary figures were width 10 height 5
  par(mar=c(5.1,4.1,4.1,4.1)) # default
  
  plot(NA, type='h', xlim=c(inipos1, max(cnvadj$end)), ylim=c(0,1.5), xlab="", ylab="", frame=F, las=1, xaxt='n')
  axis(1, at = seq(0, delta + floor(endpos2/20000000)*20000000, by= 20000000), labels = c(seq(0, (ceiling(endpos1/20000000)-1)*20, by = 20), seq(0, floor(endpos2/20000000)*20, by = 20)))
  for (k in 1:nrow(cnvadj)){
    if (cnvadj$minorAlleleCopyNumber[k] < 0.5){
      segments(cnvadj$start[k], 1, cnvadj$end[k], 1, col=rgb(122/255, 55/255, 139/255, 1), lwd=1)
    }
  }
  
  dev.off()
}

## MutCN plot
vcnplot.hmf <- function(i, chri1, chri2, ori1, ori2, pwt, pht){
  inipos1 <- 0
  inipos2 <- 0
  endpos1 <- hg19_coord$length[hg19_coord$chr == chri1]
  endpos2 <- hg19_coord$length[hg19_coord$chr == chri2]
  delta <- ceiling(endpos1/20000000) * 20000000
  
  cnvdf <- read.csv(paste0("../27_EBI/final_data_022222/acgr/", i, ".purple.cnv.somatic.acgr.tsv"), header=T, as.is=T, sep="\t")
  cnv1 <- subset(cnvdf, cnvdf$chromosome == chri1 & cnvdf$start >= inipos1 & cnvdf$end <= endpos1)
  cnv1 <- cnv1[,c("chromosome", "start", "end", "majorAlleleCopyNumber", "minorAlleleCopyNumber", "loh", "acgr")]
  cnv2 <- subset(cnvdf, cnvdf$chromosome == chri2 & cnvdf$start >= inipos2 & cnvdf$end <= endpos2)
  cnv2 <- cnv2[,c("chromosome", "start", "end", "majorAlleleCopyNumber", "minorAlleleCopyNumber", "loh", "acgr")]
  
  if (ori1 == "-"){
    tempdf <- cnv1
    maxvalue <- max(tempdf$end)
    tempdf$rev_start <- maxvalue + 1 - tempdf$end
    tempdf$rev_end <- maxvalue + 1 - tempdf$start
    tempdf$acgr[tempdf$acgr == "amplicon_start"] <- "tail"
    tempdf$acgr[tempdf$acgr == "amplicon_end"] <- "head"
    tempdf$acgr[tempdf$acgr == "tail"] <- "amplicon_end"
    tempdf$acgr[tempdf$acgr == "head"] <- "amplicon_start"
    tempdf <- tempdf[order(tempdf$rev_start, decreasing = F),]
    tempdf <- tempdf[,c("chromosome", "rev_start", "rev_end", "majorAlleleCopyNumber", "minorAlleleCopyNumber", "loh", "acgr")]
    colnames(tempdf)[2:3] <- c("start", "end")
    cnv1 <- tempdf
    value1 <- maxvalue
  }
  if (ori2 == "-"){
    tempdf <- cnv2
    maxvalue <- max(tempdf$end)
    tempdf$rev_start <- maxvalue + 1 - tempdf$end
    tempdf$rev_end <- maxvalue + 1 - tempdf$start
    tempdf$acgr[tempdf$acgr == "amplicon_start"] <- "tail"
    tempdf$acgr[tempdf$acgr == "amplicon_end"] <- "head"
    tempdf$acgr[tempdf$acgr == "tail"] <- "amplicon_end"
    tempdf$acgr[tempdf$acgr == "head"] <- "amplicon_start"
    tempdf <- tempdf[order(tempdf$rev_start, decreasing = F),]
    tempdf <- tempdf[,c("chromosome", "rev_start", "rev_end", "majorAlleleCopyNumber", "minorAlleleCopyNumber", "loh", "acgr")]
    colnames(tempdf)[2:3] <- c("start", "end")
    cnv2 <- tempdf
    value2 <- maxvalue
  }
  
  cnv2$start <- cnv2$start + delta
  cnv2$end <- cnv2$end + delta
  cnvadj <- rbind(cnv1, cnv2)
  
  vecdf <- read.csv(paste0("../27_EBI/final_data_022222/timing_veclonal/", i, ".sage30.snv.timeR.VEClonal.txt"), header=T, as.is=T, sep="\t")
  vec1 <- subset(vecdf, vecdf$chrom == chri1 & vecdf$pos >= inipos1 & vecdf$pos <= endpos1)
  vec1 <- vec1[,c("chrom", "pos", "vcn", "CLS", "kt", "pVEClonal", "pOtherClonal", "pSClonal")]
  vec2 <- subset(vecdf, vecdf$chrom == chri2 & vecdf$pos >= inipos2 & vecdf$pos <= endpos2)
  vec2 <- vec2[,c("chrom", "pos", "vcn", "CLS", "kt", "pVEClonal", "pOtherClonal", "pSClonal")]
  
  if (ori1 == "-"){
    tempdf <- vec1
    maxvalue <- hg19_coord$length[hg19_coord$chr == chri1]
    tempdf$rev_pos <- maxvalue + 1 - tempdf$pos
    tempdf <- tempdf[order(tempdf$rev_pos, decreasing = F),]
    tempdf <- tempdf[,c("chrom", "rev_pos", "vcn", "CLS", "kt", "pVEClonal", "pOtherClonal", "pSClonal")]
    colnames(tempdf)[2] <- "pos"
    vec1 <- tempdf
    value1 <- maxvalue
  }
  if (ori2 == "-"){
    tempdf <- vec2
    maxvalue <- hg19_coord$length[hg19_coord$chr == chri2]
    tempdf$rev_pos <- maxvalue + 1 - tempdf$pos
    tempdf <- tempdf[order(tempdf$rev_pos, decreasing = F),]
    tempdf <- tempdf[,c("chrom", "rev_pos", "vcn", "CLS", "kt", "pVEClonal", "pOtherClonal", "pSClonal")]
    colnames(tempdf)[2] <- "pos"
    vec2 <- tempdf
    value2 <- maxvalue
  }
  vec2$pos <- vec2$pos + delta
  vecadj <- rbind(vec1, vec2)
  vecadj$type <- 0
  vecadj$type[vecadj$CLS == "clonal [early]"] <- 2
  vecadj$type[vecadj$CLS == "clonal [late]"] <- 3
  vecadj$type[vecadj$CLS == "clonal [NA]"] <- 4
  vecadj$type[vecadj$CLS == "subclonal"] <- 5
  vecadj$type[vecadj$CLS == "clonal [early]" & vecadj$pVEClonal > vecadj$pOtherClonal & vecadj$pVEClonal > vecadj$pSClonal] <- 1
  vecadj$type[vecadj$kt != "not_analyzed"] <- 6
  
  cnvadj <- cnvadj[cnvadj$majorAlleleCopyNumber < ceiling(max(vecadj$vcn)) & cnvadj$minorAlleleCopyNumber < ceiling(max(vecadj$vcn)),]
  
  pdf(paste0(svsketchoutputdir, i, "_", chri1, "_", chri2, "_", ori1, ori2, ".vcnplot.pdf"), width=pwt, height=pht)
  par(mar=c(5.1,4.1,4.1,4.1)) # default
  plot(vecadj$pos, vecadj$vcn, xlim=c(inipos1, (delta + hg19_coord$length[hg19_coord$chr == chri2])), ylim=c(0,ceiling(max(vecadj$vcn))), xlab="", ylab="", frame=F, las=1, xaxt='n', col=clonalcolor[factor(vecadj$type, levels = c(1,2,3,4,5,6))], pch=19, cex=0.5)
  #plot(vecadj$pos, vecadj$vcn, xlim=c(inipos1, (delta + hg19_coord$length[hg19_coord$chr == chri2])), ylim=c(0,ceiling(max(cnvadj$majorAlleleCopyNumber)+max(cnvadj$minorAlleleCopyNumber))), xlab="", ylab="", frame=F, las=1, xaxt='n', col=clonalcolor[factor(vecadj$type, levels = c(1,2,3,4,5,6))], pch=19, cex=0.5)
  axis(1, at = seq(0, delta + floor(endpos2/20000000)*20000000, by= 20000000), labels = c(seq(0, (ceiling(endpos1/20000000)-1)*20, by = 20), seq(0, floor(endpos2/20000000)*20, by = 20)))
  segments(cnvadj$start[cnvadj$start != cnvadj$end], cnvadj$majorAlleleCopyNumber[cnvadj$start != cnvadj$end], cnvadj$end[cnvadj$start != cnvadj$end], cnvadj$majorAlleleCopyNumber[cnvadj$start != cnvadj$end], col=rgb(1,0,0,0.3), lwd=4/3)
  segments(cnvadj$start[cnvadj$start != cnvadj$end], cnvadj$minorAlleleCopyNumber[cnvadj$start != cnvadj$end]-0.2, cnvadj$end[cnvadj$start != cnvadj$end], cnvadj$minorAlleleCopyNumber[cnvadj$start != cnvadj$end]-0.2, col=rgb(0,0,1,0.3), lwd=4/3)
  dev.off()
}


## Example
svsketch.hmf("TCGA-A8-A08S", "11", "17", "+", "-", 4.7, 3.5)
lohplot.hmf("TCGA-A8-A08S", "11", "17", "+", "-", 4.7)


