library(dplyr)
library(grid)
library(fields)
library(stringr)
library(plotrix)
library(mapplots)
library(viridis)
library(reshape2)
options(scipen=999)

sourcedir <- "" # Location of the source directory, where the copy number segment files are located. Raw data is available in the "copy_number_segment" folder at 'http://compbio.med.harvard.edu/TBAmplification/'. 
outputdir <- "" # Location of the output directory.

maxval <- NA
for (w in which(masterdf$availability=="yes")){
  if (masterdf$file_nov26[w] != "no"){
    cnvdf <- read.csv(paste0(sourcedir, masterdf$study_id[w], ".purple.cnv.somatic.tsv"), header=T, as.is=T, sep="\t")
    isv <- read.csv(paste0(gsub("copy_number_segment", "rearrangement_bedpe", sourcedir), masterdf$study_id[w], ".purple.sv.vcf.gz.bedpe"), header=T, as.is=T, sep="\t")
    # Print ID
    print(paste0(masterdf$study_id[w], " ; ", w, " ; ", nrow(cnvdf)))
    
    # Annotation of arms
    cnvdf$arm <- "no"
    arm_idx <- "no"
    for (i in 1:nrow(cnvdf)){
      if (cnvdf$segmentStartSupport[i] == "TELOMERE"){
        arm_idx <- "p"
      } else if (cnvdf$segmentStartSupport[i] == "CENTROMERE"){
        arm_idx <- "q"
      }
      cnvdf$arm[i] <- paste0(cnvdf$chromosome[i], arm_idx)
    }
    cnvdf$seglen <- cnvdf$end - cnvdf$start
    
    # Annotation of CN baseline for each arm
    cnvdf$cnbase <- 2
    for (i in unique(cnvdf$arm)){
      infodf1 <- NULL
      for (k in sort(unique(round(cnvdf$copyNumber[cnvdf$arm == i])))){
        lenvalue <- sum(cnvdf$seglen[cnvdf$arm == i & round(cnvdf$copyNumber) == k], na.rm = T)
        infoline1 <- cbind(i, k, lenvalue)
        infodf1 <- rbind(infodf1, infoline1)
      }
      infodf1 <- as.data.frame(infodf1)
      colnames(infodf1) <- c("arm", "total_cn", "lenvalue")
      infodf1$arm <- as.character(infodf1$arm)
      infodf1$total_cn <- as.numeric(as.character(infodf1$total_cn))
      infodf1$lenvalue <- as.numeric(as.character(infodf1$lenvalue))
      if (min(infodf1$total_cn[infodf1$lenvalue == max(infodf1$lenvalue)]) > 2){
        cnvdf$cnbase[cnvdf$arm == i] <- min(infodf1$total_cn[infodf1$lenvalue == max(infodf1$lenvalue)])
      }
    }
    
    # Annotation of amplified and unamplified regions
    cnvdf$loh <- "no"
    cnvdf$loh[round(cnvdf$minorAlleleCopyNumber) == 0] <- "loh"
    cnvdf$cnstate <- "no"
    for (i in 1:nrow(cnvdf)){
      if (cnvdf$copyNumber[i] > 3*cnvdf$cnbase[i] & cnvdf$seglen[i] != 0){ #  | cnvdf$copyNumber[i] > cnvdf$cnbase[i]+6 -- Size filter (cutting out <1000 sized amplicons)
        cnvdf$cnstate[i] <- "amp"
      } else if (round(cnvdf$copyNumber[i]) < cnvdf$cnbase[i]){
        cnvdf$cnstate[i] <- "below_baseline"
      } else if (round(cnvdf$copyNumber[i]) == cnvdf$cnbase[i]){
        cnvdf$cnstate[i] <- "baseline"
      } else if (round(cnvdf$copyNumber[i]) == cnvdf$cnbase[i]+1){
        cnvdf$cnstate[i] <- "baseline+1"
      }
    }
    for (i in which(cnvdf$seglen == 0 & cnvdf$copyNumber > 3*cnvdf$cnbase)){
      if (i != 1 & i != nrow(cnvdf) & cnvdf$chromosome[i-1] == cnvdf$chromosome[i] & cnvdf$cnstate[i-1] == "amp" & cnvdf$chromosome[i] == cnvdf$chromosome[i+1] & cnvdf$cnstate[i+1] == "amp"){
        cnvdf$cnstate[i] <- "amp"
      }
    }
    
    # Annotation of amplified contiguous genomic regions
    cnvdf$cgr <- "no"
    cnvdf$cgr[cnvdf$cnstate == "amp"] <- "amp"
    cnvdf$cgr[cnvdf$cnstate == "filled_gap"] <- "amp"
    cnvdf$cgr[cnvdf$cnstate == "extended_edge_f"] <- "amp"
    cnvdf$cgr[cnvdf$cnstate == "extended_edge_r"] <- "amp"
    cnvdf$acgr <- cnvdf$cgr
    for (i in 1:nrow(cnvdf)){
      if (i == 1){
        if (cnvdf$cgr[i] == "amp" & cnvdf$cgr[i+1] == "amp"){
          cnvdf$acgr[i] <- "amplicon_start"
        } else if (cnvdf$cgr[i] == "amp" & cnvdf$cgr[i+1] == "no"){
          cnvdf$acgr[i] <- "amplicon_itself"
        }
      } else if (i == nrow(cnvdf)){
        if (cnvdf$cgr[i] == "amp" & cnvdf$cgr[i-1] == "amp"){
          cnvdf$acgr[i] <- "amplicon_end"
        } else if (cnvdf$cgr[i] == "amp" & cnvdf$cgr[i-1] == "no"){
          cnvdf$acgr[i] <- "amplicon_itself"
        }
      } else {
        if (cnvdf$cgr[i] == "amp" & cnvdf$cgr[i-1] == "no" & cnvdf$cgr[i+1] == "no"){
          cnvdf$acgr[i] <- "amplicon_itself"
        } else if (cnvdf$cgr[i] == "amp" & cnvdf$segmentStartSupport[i] != "TELOMERE" & cnvdf$cgr[i-1] == "no" & cnvdf$cgr[i+1] == "amp" | cnvdf$cgr[i] == "amp" & cnvdf$segmentStartSupport[i] == "TELOMERE" & cnvdf$cgr[i+1] == "amp"){
          cnvdf$acgr[i] <- "amplicon_start"
        } else if (cnvdf$cgr[i] == "amp" & cnvdf$segmentEndSupport[i] != "TELOMERE" & cnvdf$cgr[i-1] == "amp" & cnvdf$cgr[i+1] == "no" | cnvdf$cgr[i] == "amp" & cnvdf$segmentEndSupport[i] == "TELOMERE" & cnvdf$cgr[i+1] == "no" ){
          cnvdf$acgr[i] <- "amplicon_end"
        }
      }
    }
    
    # Size filter (cutting out <1000 sized amplicons)
    imp_num <- NA
    for (i in 1:nrow(cnvdf)){
      if (cnvdf$acgr[i] == "amplicon_itself" & cnvdf$seglen[i] < 1000){
        #print(cnvdf$seglen[i])
        cnvdf$acgr[i] <- "no"
      }
      if (cnvdf$acgr[i] == "amplicon_start"){
        imp_num <- i
      } else if (cnvdf$acgr[i] == "amplicon_end"){
        if (cnvdf$chromosome[i] == cnvdf$chromosome[imp_num] & (cnvdf$end[i] - cnvdf$start[imp_num]) < 1000){
          #print(cnvdf$end[i] - cnvdf$start[imp_num])
          cnvdf$acgr[i:imp_num] <- "no"
        }
        imp_num <- NA
      }
    }
    cnvdf <- cnvdf[,colnames(cnvdf) != "cgr"]
    
    # Consideration of Fold-back Inversions at the border
    for (i in 2:(nrow(cnvdf)-1)){
      if (cnvdf$acgr[i] == "amplicon_start" & cnvdf$segmentStartSupport[i] == "INV" & nrow(isv[isv$chrom1 == cnvdf$chromosome[i] & isv$start2 <= cnvdf$start[i] & isv$end2 >= cnvdf$start[i],]) != 0){
        maxval <- max(isv$vf[isv$chrom1 == cnvdf$chromosome[i] & isv$start2 <= cnvdf$start[i] & isv$end2 >= cnvdf$start[i]], na.omit=T)
        if (!is.na(maxval)){
          if (isv$svclass[isv$chrom1 == cnvdf$chromosome[i] & isv$start2 <= cnvdf$start[i] & isv$end2 >= cnvdf$start[i] & isv$vf == maxval] == "t2tINV"){
            k <- which(isv$chrom1 == cnvdf$chromosome[i] & isv$start2 <= cnvdf$start[i] & isv$end2 >= cnvdf$start[i] & isv$vf == maxval)
            if (i != 1 & i != nrow(cnvdf) & cnvdf$chromosome[i-1] == cnvdf$chromosome[i] & cnvdf$start[i-1] >= isv$start1[k] & cnvdf$start[i-1] <= isv$end1[k] & isv$svlen[k] <=5000 & cnvdf$copyNumber[i-1] > 2*cnvdf$cnbase[i]){
              cnvdf$acgr[i] <- "amp"
              cnvdf$acgr[i-1] <- "amplicon_start"
            }
          }
        }
        rm(maxval)
      } else if (cnvdf$acgr[i] == "amplicon_end" & cnvdf$segmentEndSupport[i] == "INV" & nrow(isv[isv$chrom1 == cnvdf$chromosome[i] & isv$start1 <= cnvdf$end[i] & isv$end1 >= cnvdf$end[i],]) != 0){
        maxval <- max(isv$vf[isv$chrom1 == cnvdf$chromosome[i] & isv$start1 <= cnvdf$end[i] & isv$end1 >= cnvdf$end[i]], na.omit=T)
        if (!is.na(maxval)){
          if (isv$svclass[isv$chrom1 == cnvdf$chromosome[i] & isv$start1 <= cnvdf$end[i] & isv$end1 >= cnvdf$end[i] & isv$vf == maxval] == "h2hINV"){
            k <- which(isv$chrom1 == cnvdf$chromosome[i] & isv$start1 <= cnvdf$end[i] & isv$end1 >= cnvdf$end[i] & isv$vf == maxval)
            if (i != 1 & i != nrow(cnvdf) & cnvdf$chromosome[i] == cnvdf$chromosome[i+1] & cnvdf$end[i+1] >= isv$start2[k] & cnvdf$end[i+1] <= isv$end2[k] & isv$svlen[k] <=5000 & cnvdf$copyNumber[i+1] > 2*cnvdf$cnbase[i]){
              cnvdf$acgr[i] <- "amp"
              cnvdf$acgr[i+1] <- "amplicon_end"
            }
          }
        }
        rm(maxval)
      }
    }
    if (sum(cnvdf$cnbase >= 6) != 0){
      print(paste0("high CN baseline (>=6): ", paste(unique(cnvdf$chromosome[cnvdf$cnbase >= 6]), collapse = ";")))
    }
    write.table(cnvdf, paste0(outputdir, masterdf$study_id[w], ".purple.cnv.somatic.acgr.tsv"), row.names=F, col.names=T, quote=F, sep="\t")
  }
}