library(dplyr)
library(grid)
library(fields)
library(stringr)
library(plotrix)
library(mapplots)
library(viridis)
library(reshape2)
options(scipen=999)

# Collecting amplified segment information after the definition of amplified segments in each sample.
# The products of the codes below is "List.amplicons_780breast.final.txt".

sourcedir <- "" # Location of the source directory, where the output files from "HMF_definition_amp_segment.R" are located.

## Collecting amplified segment information starting from master file
infodf <- data.frame(matrix(ncol=20, nrow=0))
colnames(infodf) <- c("study_id", "chromosome", "start", "end", "mcj_s", "mcj_e", "class_s", "svclass_s", "svname_s", "svf_s", "svlen_s", "leftcn_s", "rightcn_s", "class_e", "svclass_e", "svname_e", "svf_e", "svlen_e", "leftcn_e", "rightcn_e")
for (w in which(masterdf$availability=="yes")){
  acdf <- read.csv(paste0(sourcedir, masterdf$study_id[w], ".purple.cnv.somatic.acgr.tsv"), header=T, as.is=T, sep="\t")
  df <- acdf[acdf$acgr %in% c("amplicon_start", "amplicon_end", "amplicon_itself"),]
  if (nrow(df) != 0){
    svdf <- read.csv(paste0(gsub("acgr", "rearrangement_bedpe", sourcedir), masterdf$study_id[w], ".purple.sv.vcf.gz.bedpe"), header=T, as.is=T, sep="\t")
    # Print ID
    print(paste0(masterdf$study_id[w], " ; ", w, " ; ", masterdf$file_nov26[w]))
    for (i in 1:nrow(df)){
      if (df$acgr[i] == "amplicon_itself"){
        infodf[nrow(infodf)+1,] <- NA
        infodf$study_id[nrow(infodf)] <- masterdf$study_id[w]
        infodf$chromosome[nrow(infodf)] <- df$chromosome[i]
        infodf$start[nrow(infodf)] <- df$start[i]
        infodf$end[nrow(infodf)] <- df$end[i]
        if (df$segmentStartSupport[i] != "TELOMERE"){
          if (acdf$cnstate[which(acdf$chromosome==df$chromosome[i] & acdf$start==df$start[i])-1] %in% c("baseline", "baseline+1", "below_baseline")){
            infodf$mcj_s[nrow(infodf)] <- acdf$cnstate[which(acdf$chromosome==df$chromosome[i] & acdf$start==df$start[i])-1]
          } else {
            infodf$mcj_s[nrow(infodf)] <- "no"
          }
          infodf$class_s[nrow(infodf)] <- df$segmentStartSupport[i]
          if (infodf$class_s[nrow(infodf)] %in% c("BND", "DEL", "DUP", "INV", "MULTIPLE") & nrow(svdf[svdf$chrom1 == df$chromosome[i] & svdf$start1 <= df$start[i] & svdf$end1 >= df$start[i] | svdf$chrom2 == df$chromosome[i] & svdf$start2 <= df$start[i] & svdf$end2 >= df$start[i],]) != 0){
            tempdf <- svdf[svdf$chrom1 == df$chromosome[i] & svdf$start1 <= df$start[i] & svdf$end1 >= df$start[i] | svdf$chrom2 == df$chromosome[i] & svdf$start2 <= df$start[i] & svdf$end2 >= df$start[i],]
            infodf$svclass_s[nrow(infodf)] <- tempdf$svclass[tempdf$score == max(tempdf$score, na.rm=T)]
            infodf$svname_s[nrow(infodf)] <- tempdf$name[tempdf$score == max(tempdf$score, na.rm=T)]
            infodf$svf_s[nrow(infodf)] <- tempdf$vf[tempdf$score == max(tempdf$score, na.rm=T)]
            infodf$svlen_s[nrow(infodf)] <- tempdf$svlen[tempdf$score == max(tempdf$score, na.rm=T)]
          } else {
            infodf$svclass_s[nrow(infodf)] <- "no"
            infodf$svname_s[nrow(infodf)] <- "no"
            infodf$svf_s[nrow(infodf)] <- NA
            infodf$svlen_s[nrow(infodf)] <- NA
          }
          infodf$leftcn_s[nrow(infodf)] <- acdf$copyNumber[which(acdf$chromosome==df$chromosome[i] & acdf$start==df$start[i])-1]
          infodf$rightcn_s[nrow(infodf)] <- df$copyNumber[i] 
        } else {
          infodf$mcj_s[nrow(infodf)] <- "no"
          infodf$class_s[nrow(infodf)] <- df$segmentStartSupport[i]
          infodf$svclass_s[nrow(infodf)] <- "no"
          infodf$svname_s[nrow(infodf)] <- "no"
          infodf$svf_s[nrow(infodf)] <- NA
          infodf$svlen_s[nrow(infodf)] <- NA
          infodf$leftcn_s[nrow(infodf)] <- NA
          infodf$rightcn_s[nrow(infodf)] <- df$copyNumber[i] 
        }
        if (df$segmentEndSupport[i] != "TELOMERE"){
          if (acdf$cnstate[which(acdf$chromosome==df$chromosome[i] & acdf$start==df$start[i])+1] %in% c("baseline", "baseline+1", "below_baseline")){
            infodf$mcj_e[nrow(infodf)] <- acdf$cnstate[which(acdf$chromosome==df$chromosome[i] & acdf$start==df$start[i])+1]
          } else {
            infodf$mcj_e[nrow(infodf)] <- "no"
          }
          infodf$class_e[nrow(infodf)] <- df$segmentEndSupport[i]
          if (infodf$class_e[nrow(infodf)] %in% c("BND", "DEL", "DUP", "INV", "MULTIPLE") & nrow(svdf[svdf$chrom1 == df$chromosome[i] & svdf$start1 <= df$end[i] & svdf$end1 >= df$end[i] | svdf$chrom2 == df$chromosome[i] & svdf$start2 <= df$end[i] & svdf$end2 >= df$end[i],]) != 0){
            tempdf <- svdf[svdf$chrom1 == df$chromosome[i] & svdf$start1 <= df$end[i] & svdf$end1 >= df$end[i] | svdf$chrom2 == df$chromosome[i] & svdf$start2 <= df$end[i] & svdf$end2 >= df$end[i],]
            infodf$svclass_e[nrow(infodf)] <- tempdf$svclass[tempdf$score == max(tempdf$score, na.rm=T)]
            infodf$svname_e[nrow(infodf)] <- tempdf$name[tempdf$score == max(tempdf$score, na.rm=T)]
            infodf$svf_e[nrow(infodf)] <- tempdf$vf[tempdf$score == max(tempdf$score, na.rm=T)]
            infodf$svlen_e[nrow(infodf)] <- tempdf$svlen[tempdf$score == max(tempdf$score, na.rm=T)]
          } else {
            infodf$svclass_e[nrow(infodf)] <- "no"
            infodf$svname_e[nrow(infodf)] <- "no"
            infodf$svf_e[nrow(infodf)] <- NA
            infodf$svlen_e[nrow(infodf)] <- NA
          }
          infodf$leftcn_e[nrow(infodf)] <- df$copyNumber[i]
          infodf$rightcn_e[nrow(infodf)] <- acdf$copyNumber[which(acdf$chromosome==df$chromosome[i] & acdf$start==df$start[i])+1]
        } else {
          infodf$mcj_e[nrow(infodf)] <- "no"
          infodf$class_e[nrow(infodf)] <- df$segmentEndSupport[i]
          infodf$svclass_e[nrow(infodf)] <- "no"
          infodf$svname_e[nrow(infodf)] <- "no"
          infodf$svf_e[nrow(infodf)] <- NA
          infodf$svlen_e[nrow(infodf)] <- NA
          infodf$leftcn_e[nrow(infodf)] <- df$copyNumber[i]
          infodf$rightcn_e[nrow(infodf)] <- NA
        }
      } else if (df$acgr[i] == "amplicon_start"){
        if (df$acgr[i+1] != "amplicon_end"){
          print(paste0(masterdf$study_id[w], " unexpected error in acgr"))
          break
        }
        infodf[nrow(infodf)+1,] <- NA
        infodf$study_id[nrow(infodf)] <- masterdf$study_id[w]
        infodf$chromosome[nrow(infodf)] <- df$chromosome[i]
        infodf$start[nrow(infodf)] <- df$start[i]
        infodf$end[nrow(infodf)] <- df$end[i+1]
        if (df$segmentStartSupport[i] != "TELOMERE"){
          if (acdf$cnstate[which(acdf$chromosome==df$chromosome[i] & acdf$start==df$start[i])-1] %in% c("baseline", "baseline+1", "below_baseline")){
            infodf$mcj_s[nrow(infodf)] <- acdf$cnstate[which(acdf$chromosome==df$chromosome[i] & acdf$start==df$start[i])-1]
          } else {
            infodf$mcj_s[nrow(infodf)] <- "no"
          }
          infodf$class_s[nrow(infodf)] <- df$segmentStartSupport[i]
          if (infodf$class_s[nrow(infodf)] %in% c("BND", "DEL", "DUP", "INV", "MULTIPLE") & nrow(svdf[svdf$chrom1 == df$chromosome[i] & svdf$start1 <= df$start[i] & svdf$end1 >= df$start[i] | svdf$chrom2 == df$chromosome[i] & svdf$start2 <= df$start[i] & svdf$end2 >= df$start[i],]) != 0){
            tempdf <- svdf[svdf$chrom1 == df$chromosome[i] & svdf$start1 <= df$start[i] & svdf$end1 >= df$start[i] | svdf$chrom2 == df$chromosome[i] & svdf$start2 <= df$start[i] & svdf$end2 >= df$start[i],]
            infodf$svclass_s[nrow(infodf)] <- tempdf$svclass[tempdf$score == max(tempdf$score, na.rm=T)]
            infodf$svname_s[nrow(infodf)] <- tempdf$name[tempdf$score == max(tempdf$score, na.rm=T)]
            infodf$svf_s[nrow(infodf)] <- tempdf$vf[tempdf$score == max(tempdf$score, na.rm=T)]
            infodf$svlen_s[nrow(infodf)] <- tempdf$svlen[tempdf$score == max(tempdf$score, na.rm=T)]
          } else {
            infodf$svclass_s[nrow(infodf)] <- "no"
            infodf$svname_s[nrow(infodf)] <- "no"
            infodf$svf_s[nrow(infodf)] <- NA
            infodf$svlen_s[nrow(infodf)] <- NA
          }
          infodf$leftcn_s[nrow(infodf)] <- acdf$copyNumber[which(acdf$chromosome==df$chromosome[i] & acdf$start==df$start[i])-1]
          infodf$rightcn_s[nrow(infodf)] <- df$copyNumber[i] 
        } else {
          infodf$mcj_s[nrow(infodf)] <- "no"
          infodf$class_s[nrow(infodf)] <- df$segmentStartSupport[i]
          infodf$svclass_s[nrow(infodf)] <- "no"
          infodf$svname_s[nrow(infodf)] <- "no"
          infodf$svf_s[nrow(infodf)] <- NA
          infodf$svlen_s[nrow(infodf)] <- NA
          infodf$leftcn_s[nrow(infodf)] <- NA
          infodf$rightcn_s[nrow(infodf)] <- df$copyNumber[i] 
        }
        if (df$segmentEndSupport[i+1] != "TELOMERE"){
          if (acdf$cnstate[which(acdf$chromosome==df$chromosome[i+1] & acdf$start==df$start[i+1])+1] %in% c("baseline", "baseline+1", "below_baseline")){
            infodf$mcj_e[nrow(infodf)] <- acdf$cnstate[which(acdf$chromosome==df$chromosome[i+1] & acdf$start==df$start[i+1])+1]
          } else {
            infodf$mcj_e[nrow(infodf)] <- "no"
          }
          infodf$class_e[nrow(infodf)] <- df$segmentEndSupport[i+1]
          if (infodf$class_e[nrow(infodf)] %in% c("BND", "DEL", "DUP", "INV", "MULTIPLE") & nrow(svdf[svdf$chrom1 == df$chromosome[i+1] & svdf$start1 <= df$end[i+1] & svdf$end1 >= df$end[i+1] | svdf$chrom2 == df$chromosome[i+1] & svdf$start2 <= df$end[i+1] & svdf$end2 >= df$end[i+1],]) != 0){
            tempdf <- svdf[svdf$chrom1 == df$chromosome[i+1] & svdf$start1 <= df$end[i+1] & svdf$end1 >= df$end[i+1] | svdf$chrom2 == df$chromosome[i+1] & svdf$start2 <= df$end[i+1] & svdf$end2 >= df$end[i+1],]
            infodf$svclass_e[nrow(infodf)] <- tempdf$svclass[tempdf$score == max(tempdf$score, na.rm=T)]
            infodf$svname_e[nrow(infodf)] <- tempdf$name[tempdf$score == max(tempdf$score, na.rm=T)]
            infodf$svf_e[nrow(infodf)] <- tempdf$vf[tempdf$score == max(tempdf$score, na.rm=T)]
            infodf$svlen_e[nrow(infodf)] <- tempdf$svlen[tempdf$score == max(tempdf$score, na.rm=T)]
          } else {
            infodf$svclass_e[nrow(infodf)] <- "no"
            infodf$svname_e[nrow(infodf)] <- "no"
            infodf$svf_e[nrow(infodf)] <- NA
            infodf$svlen_e[nrow(infodf)] <- NA
          }
          infodf$leftcn_e[nrow(infodf)] <- df$copyNumber[i+1]
          infodf$rightcn_e[nrow(infodf)] <- acdf$copyNumber[which(acdf$chromosome==df$chromosome[i+1] & acdf$start==df$start[i+1])+1]
        } else {
          infodf$mcj_e[nrow(infodf)] <- "no"
          infodf$class_e[nrow(infodf)] <- df$segmentEndSupport[i+1]
          infodf$svclass_e[nrow(infodf)] <- "no"
          infodf$svname_e[nrow(infodf)] <- "no"
          infodf$svf_e[nrow(infodf)] <- NA
          infodf$svlen_e[nrow(infodf)] <- NA
          infodf$leftcn_e[nrow(infodf)] <- df$copyNumber[i+1]
          infodf$rightcn_e[nrow(infodf)] <- NA
        }
      } else if (df$acgr[i] == "amplicon_end"){
        
      }
    }
  }
  
}  

## Assignment of mechanisms
infodf$mech_s <- "not_analyzed"
infodf$mech_e <- "not_analyzed"
infodf$mech_s[infodf$mcj_s != "no"] <- "analyzed"
infodf$mech_e[infodf$mcj_e != "no"] <- "analyzed"
infodf$mech_s[infodf$mcj_s != "no" & infodf$mcj_e != "no" & infodf$svclass_s == "DUP" & infodf$svclass_e == "DUP" & infodf$svname_s == infodf$svname_e] <- "tandem_duplication"
infodf$mech_e[infodf$mcj_s != "no" & infodf$mcj_e != "no" & infodf$svclass_s == "DUP" & infodf$svclass_e == "DUP" & infodf$svname_s == infodf$svname_e] <- "tandem_duplication"
infodf$mech_s[infodf$mcj_s != "no" & infodf$mcj_e != "no" & infodf$svclass_s == "DUP" & infodf$svclass_e == "DUP" & infodf$svname_s == infodf$svname_e & infodf$rightcn_s > 4*infodf$leftcn_s & infodf$leftcn_e > 4*infodf$rightcn_e] <- "td_extrachromosomal"
infodf$mech_e[infodf$mcj_s != "no" & infodf$mcj_e != "no" & infodf$svclass_s == "DUP" & infodf$svclass_e == "DUP" & infodf$svname_s == infodf$svname_e & infodf$rightcn_s > 4*infodf$leftcn_s & infodf$leftcn_e > 4*infodf$rightcn_e] <- "td_extrachromosomal"
infodf$mech_s[infodf$mcj_s != "no" & infodf$svclass_s == "t2tINV" & infodf$svlen_s <=5000 & infodf$svlen_s < (infodf$end-infodf$start) & !(infodf$svname_s %in% infodf$svname_s[infodf$mcj_s != "no" & infodf$svclass_s == "t2tINV" & infodf$svlen_s <=5000][duplicated(infodf$svname_s[infodf$mcj_s != "no" & infodf$svclass_s == "t2tINV" & infodf$svlen_s <=5000])])] <- "fold_back_inversion"
infodf$mech_e[infodf$mcj_e != "no" & infodf$svclass_e == "h2hINV" & infodf$svlen_e <=5000 & infodf$svlen_e < (infodf$end-infodf$start) & !(infodf$svname_e %in% infodf$svname_e[infodf$mcj_e != "no" & infodf$svclass_e == "h2hINV" & infodf$svlen_e <=5000][duplicated(infodf$svname_e[infodf$mcj_e != "no" & infodf$svclass_e == "h2hINV" & infodf$svlen_e <=5000])])] <- "fold_back_inversion"
infodf$mech_s[infodf$mcj_s != "no" & infodf$svclass_s == "TRA"] <- "translocation"
infodf$mech_e[infodf$mcj_e != "no" & infodf$svclass_e == "TRA"] <- "translocation"
infodf$mech_s[infodf$mech_s == "analyzed" & infodf$svclass_s != "no"] <- "intrachromosomal_complex"
infodf$mech_s[infodf$mech_s == "analyzed" & infodf$svclass_s == "no"] <- "no_supporting_sv"
infodf$mech_e[infodf$mech_e == "analyzed" & infodf$svclass_e != "no"] <- "intrachromosomal_complex"
infodf$mech_e[infodf$mech_e == "analyzed" & infodf$svclass_e == "no"] <- "no_supporting_sv"

infodf$mech_c <- NA
infodf$mech_c[infodf$mech_s == "td_extrachromosomal" & infodf$mech_e == "td_extrachromosomal"] <- "double_minute"
infodf$mech_c[infodf$mech_s == "tandem_duplication" & infodf$mech_e == "tandem_duplication"] <- "tandem_duplication"
infodf$mech_c[infodf$mech_s == "fold_back_inversion" | infodf$mech_e == "fold_back_inversion"] <- "fold_back_inversion"
infodf$mech_c[infodf$mech_s == "translocation" | infodf$mech_e == "translocation"] <- "translocation"
infodf$mech_c[(infodf$mech_s == "translocation" & infodf$mech_e == "fold_back_inversion") | (infodf$mech_s == "fold_back_inversion" & infodf$mech_e == "translocation")] <- "tra+fbi"
infodf$mech_c[is.na(infodf$mech_c)] <- "intrachromosomal_complex"
infodf$mech_c[infodf$mech_s == "no_supporting_sv" & infodf$mech_e == "no_supporting_sv"] <- "no_supporting_sv"
infodf$mech_c[infodf$mech_s == "not_analyzed" & infodf$mech_e == "not_analyzed"] <- "excluded"

write.table(infodf, "List.amplicons_780breast.final.txt", row.names=F, col.names=T, quote=F, sep="\t")
amphmf <- infodf ## Renaming for convenience.


## Collecting SVs in the major boundaries
tempdf <- amphmf[amphmf$mcj_s != "no" | amphmf$mcj_e != "no",]
tumorlist <- unique(tempdf$study_id)
for (i in tumorlist){
  print(i)
  isv <- read.csv(paste0("../27_EBI/final_data_022222/rearrangement_bedpe_linx/", i, ".purple.sv.linx.bedpe.tsv"), header=T, as.is=T, sep="\t")
  if (which(tumorlist == i) == 1){
    df <- isv[isv$name %in% unique(c(tempdf$svname_s[tempdf$study_id==i & tempdf$mcj_s != "no"], tempdf$svname_e[tempdf$study_id==i & tempdf$mcj_e != "no"])),]
    df$study_id <- i
  } else {
    temporarydf <-isv[isv$name %in% unique(c(tempdf$svname_s[tempdf$study_id==i & tempdf$mcj_s != "no"], tempdf$svname_e[tempdf$study_id==i & tempdf$mcj_e != "no"])),]
    if (nrow(temporarydf) != 0){
      temporarydf$study_id <- i
      df <- rbind(df, temporarydf)
    }
  }
}
df$arm1 <- "q"
df$arm2 <- "q"
for (i in 1:nrow(df)){
  if ((df$start1[i]+df$end1[i])/2 < hg19_coord$centro_mid[hg19_coord$chr == df$chrom1[i]]){
    df$arm1[i] <- "p"
  } 
  if ((df$start2[i]+df$end2[i])/2 < hg19_coord$centro_mid[hg19_coord$chr == df$chrom2[i]]){
    df$arm2[i] <- "p"
  } 
}
svkey <- df

## Descriptive statistics
library(plotrix)
tempdf <- svkey[!is.na(svkey$svlen),]
tempdf$svclass <- factor(tempdf$svclass, levels = c("DEL", "DUP", "h2hINV", "t2tINV"))
tempdf$loglen <- log10(tempdf$svlen)
histStack(loglen ~ svclass, data=tempdf, col=c(rgb(0/255,0/255,255/255), rgb(0/255,128/255,128/255), rgb(220/255,20/255,60/255), rgb(128/255,128/255,0/255)), breaks=seq(0,9,0.1), las=1, xlab="Size of the intrachromosomal SVs (log10)", main=paste0("Boundary SVs flanking unamplified region (n=", nrow(tempdf), ")"))

## Collecting all SVs
df <- masterdf[masterdf$availability == "yes",c("study_id", "cohort")]
for (i in 1:nrow(df)){
  print(paste0(df$study_id[i], "; ", i))
  if (i == 1){
    allsv <- read.csv(paste0("../27_EBI/final_data_022222/rearrangement_bedpe_linx/", df$study_id[i], ".purple.sv.linx.bedpe.tsv"), header=T, as.is=T, sep="\t")
    allsv$study_id <- df$study_id[i]
    allsv$boundary <- "no"
    if (nrow(svkey[svkey$study_id == df$study_id[i],]) != 0){
      allsv$boundary[allsv$name %in% svkey$name[svkey$study_id == df$study_id[i]]] <- "yes"
    }
  } else {
    tempdf <- read.csv(paste0("../27_EBI/final_data_022222/rearrangement_bedpe_linx/", df$study_id[i], ".purple.sv.linx.bedpe.tsv"), header=T, as.is=T, sep="\t")
    if (nrow(tempdf) != 0){
      tempdf$study_id <- df$study_id[i]
      tempdf$boundary <- "no"
      if (nrow(svkey[svkey$study_id == df$study_id[i],]) != 0){
        tempdf$boundary[tempdf$name %in% svkey$name[svkey$study_id == df$study_id[i]]] <- "yes"
      }
      allsv <- rbind(allsv, tempdf)
    }
  }
}
allsv$plotbound <- "no"
allsv$plotbound[allsv$boundary=="yes" & allsv$svclass!="TRA"] <- "yes_intra"
allsv$plotbound[allsv$boundary=="yes" & allsv$svclass=="TRA"] <- "yes_trans"
allsv$arm1 <- paste0(allsv$chrom1, "q")
allsv$arm2 <- paste0(allsv$chrom2, "q")
for (i in 1:nrow(hg19_coord)){
  allsv$arm1[allsv$chrom1 == hg19_coord$chr[i] & allsv$start1 < hg19_coord$centro_mid[i]] <- paste0(hg19_coord$chr[i], "p")
  allsv$arm2[allsv$chrom2 == hg19_coord$chr[i] & allsv$start2 < hg19_coord$centro_mid[i]] <- paste0(hg19_coord$chr[i], "p")
}

# Comparing read supports between non-boundary vs boundary
boxplot(allsv$vf[!is.na(allsv$vf)] ~ allsv$plotbound[!is.na(allsv$vf)], log='y', las=1, frame=F, outline=F, xlab="", ylab="Read support", ylim=c(1,4000), col=c(rgb(204/255, 204/255, 204/255, .7), rgb(67/255, 110/255, 238/255, .7), rgb(128/255, 0/255, 128/255, .7)))
stripchart(allsv$vf[!is.na(allsv$vf)] ~ allsv$plotbound[!is.na(allsv$vf)], log='y', pch=20, col=rgb(0,0,0,.01), vertical=T, add=T, method="jitter", xlab="", ylab="", frame=F, las=1, xaxt='n', yaxt='n', cex=.3)
vioplot(allsv$vf[!is.na(allsv$vf)] ~ allsv$plotbound[!is.na(allsv$vf)], ylog=T, las=1, frame=F, outline=F, xlab="", ylab="Read support", ylim=c(1,4000), col=c(rgb(204/255, 204/255, 204/255, .7), rgb(67/255, 110/255, 238/255, .7), rgb(128/255, 0/255, 128/255, .7)))
boxplot(allsv$vf[!is.na(allsv$vf)] ~ allsv$plotbound[!is.na(allsv$vf)], log='y', las=1, frame=F, outline=F, xlab="", ylab="Read support", ylim=c(1,4000), col=c(rgb(204/255, 204/255, 204/255, .7), rgb(67/255, 110/255, 238/255, .7), rgb(128/255, 0/255, 128/255, .7)))

t.test(allsv$vf[!is.na(allsv$vf) & allsv$plotbound=="no"], allsv$vf[!is.na(allsv$vf) & allsv$plotbound %in% c("yes_intra", "yes_trans")])
t.test(allsv$vf[!is.na(allsv$vf) & allsv$plotbound=="no"], allsv$vf[!is.na(allsv$vf) & allsv$plotbound=="yes_intra"])
t.test(allsv$vf[!is.na(allsv$vf) & allsv$plotbound=="no"], allsv$vf[!is.na(allsv$vf) & allsv$plotbound=="yes_trans"])
length(allsv$vf[!is.na(allsv$vf) & allsv$plotbound=="no"])
length(allsv$vf[!is.na(allsv$vf) & allsv$plotbound=="yes_intra"])
length(allsv$vf[!is.na(allsv$vf) & allsv$plotbound=="yes_trans"])

t.test(log10(allsv$vf[!is.na(allsv$vf) & allsv$plotbound=="no"]), allsv$vf[!is.na(allsv$vf) & allsv$plotbound=="yes_intra"])
t.test(log10(allsv$vf[!is.na(allsv$vf) & allsv$plotbound=="no"]), allsv$vf[!is.na(allsv$vf) & allsv$plotbound=="yes_trans"])

## Comparing length of microhomology
vioplot(allsv$homlen ~ allsv$plotbound, las=1, frame=F, outline=F, col=c(rgb(204/255, 204/255, 204/255, .7), rgb(67/255, 110/255, 238/255, .7), rgb(128/255, 0/255, 128/255, .7)))
boxplot(allsv$homlen ~ allsv$boundary, las=1, frame=F, outline=F, main=paste0("p-value = ", round(t.test(allsv$homlen ~ allsv$boundary)$p.value, 4)), col=rev(twocolors))
t.test(allsv$homlen ~ allsv$boundary)

wilcox.test(allsv$homlen[allsv$plotbound != "yes_intra"] ~ allsv$boundary[allsv$plotbound != "yes_intra"])
wilcox.test(allsv$homlen[allsv$plotbound != "yes_trans"] ~ allsv$boundary[allsv$plotbound != "yes_trans"])

## Replication timing comparison
fisher.test(c(allsv$fragileSiteStart, allsv$fragileSiteEnd), c(allsv$boundary, allsv$boundary))
boxplot(allsv$replicationTimingStart ~ allsv$boundary, las=1, frame=F, outline=F, main=paste0("p-value = ", round(t.test(allsv$replicationTimingStart ~ allsv$boundary)$p.value, 4)), col=rev(twocolors))
boxplot(allsv$replicationTimingEnd ~ allsv$boundary, las=1, frame=F, outline=F, main=paste0("p-value = ", round(t.test(allsv$replicationTimingStart ~ allsv$boundary)$p.value, 4)), col=rev(twocolors))
boxplot(allsv$replicationTimingStart ~ allsv$svclass, las=1, frame=F, outline=F)


