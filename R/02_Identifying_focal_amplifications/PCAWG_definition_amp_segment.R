library(dplyr)
library(grid)
library(fields)
library(stringr)
library(plotrix)
library(mapplots)
library(viridis)
library(reshape2)
options(scipen=999)

sourcedir <- "" # Location of the source directory, where the copy number segment files are located. Raw data is available in the "CNA_absCN_GC.zip" file at the Data folder. 
outputdir <- "" # Location of the output directory.

sumdf <- read.csv("Summaryinfo.table.1.19.txt", sep="\t", as.is=T, header=T) # "Summaryinfo.table.1.19.txt" available in the Data folder.
sumdf$tempnum <- as.numeric(as.character(gsub("DO", "", sumdf$icgc_donor_id)))
sumdf <- sumdf[order(sumdf$tempnum, decreasing = F),]
sumdf <- sumdf[order(sumdf$project_code, decreasing = F),]
sumdf <- sumdf[,colnames(sumdf) != "tempnum"]

for (w in which(!is.na(sumdf$SV.events))){
  svdf <- read.csv(paste0("/Users/jjklee/Documents/ICGC_PCAWG_calls/SV_calls_1.6/BEDPE/", sumdf$tumor[w], ".pcawg_consensus_1.6.161022.somatic.sv.bedpe"), header=T, as.is=T, sep="\t")
  if (nrow(svdf) != 0){
    svdf$svlen <- NA
    svdf$svlen[svdf$svclass != "TRA"] <- abs((svdf$start2[svdf$svclass != "TRA"]+svdf$end2[svdf$svclass != "TRA"])/2 - (svdf$start1[svdf$svclass != "TRA"]+svdf$end1[svdf$svclass != "TRA"])/2)
  }
  cnvdf <- read.csv(paste0(sourcedir, sumdf$icgc_donor_id[w], ".somatic.cna.abscn.txt"), head=T, as.is=T, sep="\t")
  cnvdf$copyNumber <- cnvdf$depth
  
  # Preprocessing of CNV files (segmentation in centromeres)
  for (i in 1:nrow(hg19_coord)){
    if (nrow(cnvdf[cnvdf$chromosome == hg19_coord$chr[i] & cnvdf$start < hg19_coord$centro_mid[i] & cnvdf$end > hg19_coord$centro_mid[i],]) == 1){
      k <- which(cnvdf$chromosome == hg19_coord$chr[i] & cnvdf$start < hg19_coord$centro_mid[i] & cnvdf$end > hg19_coord$centro_mid[i])
      cnvdf[nrow(cnvdf)+1,] <- cnvdf[k,]
      cnvdf$end[k] <- hg19_coord$centro_mid[i]-1
      cnvdf$start[nrow(cnvdf)] <- hg19_coord$centro_mid[i]
    } 
  }
  cnvdf$chr_num <- cnvdf$chromosome
  cnvdf$chr_num[cnvdf$chr_num=="X"] <- "23"
  cnvdf$chr_num[cnvdf$chr_num=="Y"] <- "24"
  cnvdf$chr_num[cnvdf$chr_num=="MT"] <- "25"
  cnvdf$chr_num <- as.numeric(as.character(cnvdf$chr_num))
  cnvdf <- cnvdf[order(cnvdf$start, decreasing = F),]
  cnvdf <- cnvdf[order(cnvdf$chr_num, decreasing = F),]
  cnvdf <- cnvdf[,colnames(cnvdf) != "chr_num"]
  cnvdf$segmentStartSupport <- "NONE"
  cnvdf$segmentEndSupport <- "NONE"
  for (i in 2:nrow(cnvdf)){
    if (nrow(svdf[svdf$chrom1 == cnvdf$chromosome[i] & svdf$start1 == cnvdf$start[i] | svdf$chrom2 == cnvdf$chromosome[i] & svdf$start2 == cnvdf$start[i],]) != 0){
      k <- which(svdf$chrom1 == cnvdf$chromosome[i] & svdf$start1 == cnvdf$start[i] | svdf$chrom2 == cnvdf$chromosome[i] & svdf$start2 == cnvdf$start[i])
      if (length(k) >1){
        k <- k[svdf$pe_support[k] == max(svdf$pe_support[k])][1]
      }
      cnvdf$segmentStartSupport[i] <- svdf$svclass[k]
      if (cnvdf$chromosome[i] == cnvdf$chromosome[i-1]){
        cnvdf$segmentEndSupport[i-1] <- svdf$svclass[k]
      }
    }
  }
  cnvdf$segmentStartSupport[1] <- "TELOMERE"
  cnvdf$segmentEndSupport[nrow(cnvdf)] <- "TELOMERE"
  for (i in 2:nrow(cnvdf)){
    if (cnvdf$chromosome[i] != cnvdf$chromosome[i-1]){
      cnvdf$segmentStartSupport[i] <- "TELOMERE"
      cnvdf$segmentEndSupport[i-1] <- "TELOMERE"
    }
  }
  for (i in 1:nrow(hg19_coord)){
    k <- which(cnvdf$chromosome == hg19_coord$chr[i] & cnvdf$start == hg19_coord$centro_mid[i])
    cnvdf$segmentStartSupport[k] <- "CENTROMERE"
    cnvdf$segmentEndSupport[k-1] <- "CENTROMERE"
  }
  cnvdf$segmentStartSupport[cnvdf$segmentStartSupport=="TRA"] <- "BND"
  cnvdf$segmentEndSupport[cnvdf$segmentEndSupport=="TRA"] <- "BND"
  cnvdf$segmentStartSupport[cnvdf$segmentStartSupport=="h2hINV"] <- "INV"
  cnvdf$segmentEndSupport[cnvdf$segmentEndSupport=="h2hINV"] <- "INV"
  cnvdf$segmentStartSupport[cnvdf$segmentStartSupport=="t2tINV"] <- "INV"
  cnvdf$segmentEndSupport[cnvdf$segmentEndSupport=="t2tINV"] <- "INV"
  
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
    if (length(infodf1) != 0){ # In PCAWG, some acrocentric arms were copyNumber == 0 per calculation
      infodf1 <- as.data.frame(infodf1)
      colnames(infodf1) <- c("arm", "total_cn", "lenvalue")
      infodf1$arm <- as.character(infodf1$arm)
      infodf1$total_cn <- as.numeric(as.character(infodf1$total_cn))
      infodf1$lenvalue <- as.numeric(as.character(infodf1$lenvalue))
      if (min(infodf1$total_cn[infodf1$lenvalue == max(infodf1$lenvalue)]) > 2){
        cnvdf$cnbase[cnvdf$arm == i] <- min(infodf1$total_cn[infodf1$lenvalue == max(infodf1$lenvalue)])
      }
    }
  }
  
  # Annotation of amplified and unamplified regions
  cnvdf$loh <- "no"
  cnvdf$loh[round(cnvdf$minor_cn) == 0] <- "loh"
  cnvdf$cnstate <- "no"
  for (i in 1:nrow(cnvdf)){
    if (!is.na(cnvdf$copyNumber[i])){ # in PCAWG, some of the segments are copyNumber == NA hence needs this additional step
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
  }
  for (i in which(cnvdf$seglen == 0 & cnvdf$copyNumber > 3*cnvdf$cnbase & !is.na(cnvdf$copyNumber))){
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
  if (nrow(svdf) != 0){
    for (i in 2:(nrow(cnvdf)-1)){
      if (cnvdf$acgr[i] == "amplicon_start" & cnvdf$segmentStartSupport[i] == "INV" & nrow(svdf[svdf$chrom1 == cnvdf$chromosome[i] & svdf$start2 <= cnvdf$start[i] & svdf$end2 >= cnvdf$start[i],]) != 0){
        jsv <- svdf[svdf$chrom1 == cnvdf$chromosome[i] & svdf$start2 <= cnvdf$start[i] & svdf$end2 >= cnvdf$start[i] & svdf$svclass == "t2tINV",]
        if (nrow(jsv) == 1){
          k <- which(svdf$chrom1 == cnvdf$chromosome[i] & svdf$start2 <= cnvdf$start[i] & svdf$end2 >= cnvdf$start[i] & svdf$svclass == "t2tINV")
          if (i != 1 & i != nrow(cnvdf) & cnvdf$chromosome[i-1] == cnvdf$chromosome[i] & cnvdf$start[i-1] >= svdf$start1[k] & cnvdf$start[i-1] <= svdf$end1[k] & svdf$svlen[k] <=5000 & cnvdf$copyNumber[i-1] > 2*cnvdf$cnbase[i]){
            cnvdf$acgr[i] <- "amp"
            cnvdf$acgr[i-1] <- "amplicon_start"
          }
        } else if (nrow(jsv) >1){
          k <- which(svdf$chrom1 == cnvdf$chromosome[i] & svdf$start2 <= cnvdf$start[i] & svdf$end2 >= cnvdf$start[i] & svdf$svclass == "t2tINV" & svdf$pe_support == max(jsv$pe_support))
          if (i != 1 & i != nrow(cnvdf) & cnvdf$chromosome[i-1] == cnvdf$chromosome[i] & cnvdf$start[i-1] >= svdf$start1[k] & cnvdf$start[i-1] <= svdf$end1[k] & svdf$svlen[k] <=5000 & cnvdf$copyNumber[i-1] > 2*cnvdf$cnbase[i]){
            cnvdf$acgr[i] <- "amp"
            cnvdf$acgr[i-1] <- "amplicon_start"
          }
        }
      } else if (cnvdf$acgr[i] == "amplicon_end" & cnvdf$segmentEndSupport[i] == "INV" & nrow(svdf[svdf$chrom1 == cnvdf$chromosome[i] & svdf$start1 <= cnvdf$end[i] & svdf$end1 >= cnvdf$end[i],]) != 0){
        jsv <- svdf[svdf$chrom1 == cnvdf$chromosome[i] & svdf$start1 <= cnvdf$end[i] & svdf$end1 >= cnvdf$end[i] & svdf$svclass == "h2hINV",]
        if (nrow(jsv) == 1){
          k <- which(svdf$chrom1 == cnvdf$chromosome[i] & svdf$start1 <= cnvdf$end[i] & svdf$end1 >= cnvdf$end[i] & svdf$svclass == "h2hINV")
          if (i != 1 & i != nrow(cnvdf) & cnvdf$chromosome[i] == cnvdf$chromosome[i+1] & cnvdf$end[i+1] >= svdf$start2[k] & cnvdf$end[i+1] <= svdf$end2[k] & svdf$svlen[k] <=5000 & cnvdf$copyNumber[i+1] > 2*cnvdf$cnbase[i]){
            cnvdf$acgr[i] <- "amp"
            cnvdf$acgr[i+1] <- "amplicon_end"
          }
        } else if (nrow(jsv) >1){
          k <- which(svdf$chrom1 == cnvdf$chromosome[i] & svdf$start1 <= cnvdf$end[i] & svdf$end1 >= cnvdf$end[i] & svdf$svclass == "h2hINV" & svdf$pe_support == max(jsv$pe_support))
          if (i != 1 & i != nrow(cnvdf) & cnvdf$chromosome[i] == cnvdf$chromosome[i+1] & cnvdf$end[i+1] >= svdf$start2[k] & cnvdf$end[i+1] <= svdf$end2[k] & svdf$svlen[k] <=5000 & cnvdf$copyNumber[i+1] > 2*cnvdf$cnbase[i]){
            cnvdf$acgr[i] <- "amp"
            cnvdf$acgr[i+1] <- "amplicon_end"
          }
        }
      }
    }
  }
  
  exprange <- 1000000
  cnvdf$expanded_cgr <- cnvdf$acgr
  for (i in 2:(nrow(cnvdf)-1)){
    if (cnvdf$acgr[i] == "amplicon_start"){
      curr_chr <- cnvdf$chromosome[i]
      curr_s <- cnvdf$start[i]
      if (curr_s > exprange){
        if (curr_s - exprange < cnvdf$start[which(cnvdf$chromosome == curr_chr)[1]]){ # exceptional case of no report of telomeric segment
          anchor_val <- which(cnvdf$chromosome == curr_chr)[1]
        } else {
          anchor_val <- which(cnvdf$chromosome == curr_chr & cnvdf$end >= (curr_s-exprange) & cnvdf$start <= (curr_s-exprange))
        }
      } else {
        anchor_val <- which(cnvdf$chromosome == curr_chr)[1]
      }
      for (j in (i-1):anchor_val){
        if (cnvdf$acgr[j] != "no"){
          imp_num1 <- i
          break
        } else if (cnvdf$depth[j] <= cnvdf$cnbase[j] & !is.na(cnvdf$depth[j]) & !is.na(cnvdf$cnbase[j])){
          imp_num1 <- j+1
          break
        } else {
          imp_num1 <- i
        }
      }
      if (imp_num1 != i){
        cnvdf$expanded_cgr[imp_num1] <- "amplicon_start"
        cnvdf$expanded_cgr[i] <- "no"
        rm(imp_num1)
      }
    } else if (cnvdf$acgr[i] == "amplicon_end"){
      curr_chr <- cnvdf$chromosome[i]
      curr_e <- cnvdf$end[i]
      if ((max(cnvdf$end[cnvdf$chromosome == curr_chr])-curr_e) > exprange){
        if (curr_e + exprange > cnvdf$end[which(cnvdf$chromosome == curr_chr)[length(which(cnvdf$chromosome == curr_chr))]]){ # exceptional case of no report of telomeric segment
          anchor_val <- which(cnvdf$chromosome == curr_chr)[length(which(cnvdf$chromosome == curr_chr))]
        } else {
          anchor_val <- which(cnvdf$chromosome == curr_chr & cnvdf$end >= (curr_e+exprange) & cnvdf$start <= (curr_e+exprange))
        }
      } else {
        anchor_val <- which(cnvdf$chromosome == curr_chr)[length(which(cnvdf$chromosome == curr_chr))]
      }
      for (j in (i+1):anchor_val){
        if (cnvdf$acgr[j] != "no"){
          imp_num2 <- i
          break
        } else if (cnvdf$depth[j] <= cnvdf$cnbase[j] & !is.na(cnvdf$depth[j]) & !is.na(cnvdf$cnbase[j])){
          imp_num2 <- j-1
          break
        } else {
          imp_num2 <- i
        }
      }
      if (imp_num2 != i){
        cnvdf$expanded_cgr[imp_num2] <- "amplicon_end"
        cnvdf$expanded_cgr[i] <- "no"
        rm(imp_num2)
      }
    } else if (cnvdf$acgr[i] == "amplicon_itself"){
      curr_chr <- cnvdf$chromosome[i]
      curr_s <- cnvdf$start[i]
      if (curr_s > exprange){
        if (curr_s - exprange < cnvdf$start[which(cnvdf$chromosome == curr_chr)[1]]){ # exceptional case of no report of telomeric segment
          anchor_val <- which(cnvdf$chromosome == curr_chr)[1]
        } else {
          anchor_val <- which(cnvdf$chromosome == curr_chr & cnvdf$end >= (curr_s-exprange) & cnvdf$start <= (curr_s-exprange))
        }
      } else {
        anchor_val <- which(cnvdf$chromosome == curr_chr)[1]
      }
      for (j in (i-1):anchor_val){
        if (cnvdf$acgr[j] != "no"){
          imp_num1 <- i
          break
        } else if (cnvdf$depth[j] <= cnvdf$cnbase[j] & !is.na(cnvdf$depth[j]) & !is.na(cnvdf$cnbase[j])){
          imp_num1 <- j+1
          break
        } else {
          imp_num1 <- i
        }
      }
      curr_e <- cnvdf$end[i]
      if ((max(cnvdf$end[cnvdf$chromosome == curr_chr])-curr_e) > exprange){
        if (curr_e + exprange > cnvdf$end[which(cnvdf$chromosome == curr_chr)[length(which(cnvdf$chromosome == curr_chr))]]){ # exceptional case of no report of telomeric segment
          anchor_val <- which(cnvdf$chromosome == curr_chr)[length(which(cnvdf$chromosome == curr_chr))]
        } else {
          anchor_val <- which(cnvdf$chromosome == curr_chr & cnvdf$end >= (curr_e+exprange) & cnvdf$start <= (curr_e+exprange))
        }
      } else {
        anchor_val <- which(cnvdf$chromosome == curr_chr)[length(which(cnvdf$chromosome == curr_chr))]
      }
      for (j in (i+1):anchor_val){
        if (cnvdf$acgr[j] != "no"){
          imp_num2 <- i
          break
        } else if (cnvdf$depth[j] <= cnvdf$cnbase[j] & !is.na(cnvdf$depth[j]) & !is.na(cnvdf$cnbase[j])){
          imp_num2 <- j-1
          break
        } else {
          imp_num2 <- i
        }
      }
      if (imp_num1 != i){
        cnvdf$expanded_cgr[imp_num1] <- "amplicon_start"
        cnvdf$expanded_cgr[i] <- "amplicon_end"
        rm(imp_num1)
      }
      if (imp_num2 != i){
        cnvdf$expanded_cgr[imp_num2] <- "amplicon_end"
        if (cnvdf$expanded_cgr[i] == "amplicon_end"){ # the case that the original amplicon_itself is already modified
          cnvdf$expanded_cgr[i] <- "no"
        } else {
          cnvdf$expanded_cgr[i] <- "amplicon_start"
        }
        rm(imp_num2)
      }
    }
  }
  colnames(cnvdf)[ncol(cnvdf)-1] <- "pre_acgr"
  colnames(cnvdf)[ncol(cnvdf)] <- "acgr"
  
  if (sum(cnvdf$cnbase >= 6) != 0){
    print(paste0("high CN baseline (>=6): ", paste(unique(cnvdf$chromosome[cnvdf$cnbase >= 6]), collapse = ";")))
  }
  write.table(cnvdf, paste0(outputdir, sumdf$icgc_donor_id[w], ".pcawg.cnv.somatic.acgr.tsv"), row.names=F, col.names=T, quote=F, sep="\t")
}
