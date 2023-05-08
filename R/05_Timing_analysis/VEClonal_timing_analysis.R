library(dplyr)
library(grid)
library(fields)
library(stringr)
library(mapplots)
library(reshape2)
library(data.table)
library(beeswarm)

# Timing of arm-level gains based on VEClonal mutations
## Integrating VEClonal mutation to the arm-level data
tempdf <- data.frame(matrix(vector(), 0, 5, dimnames=list(c(), c("study_id", "is.wgd", "time", "d_veclonal","arm"))), stringsAsFactors=F)

lengthminimum <- 5000000
for (i in c("1q", "8q", "16p", "20q")){
  print(i)
  df <- tmdf[,c(4,1,2,3,c(5:ncol(tmdf)))]
  df$p_nseg <- NA
  df$q_nseg <- NA
  df$p_majcn <- NA
  df$q_majcn <- NA
  df$p_mincn <- NA
  df$q_mincn <- NA
  df$p_time <- NA
  df$q_time <- NA
  df$d_veclonal <- NA
  for (w in 1:nrow(df)){
    cnvdf <- read.csv(paste0("../27_EBI/final_data_022222/acgr/", df$study_id[w], ".purple.cnv.somatic.acgr.tsv"), header=T, as.is=T, sep="\t")
    timer <- read.csv(paste0("../27_EBI/final_data_022222/timeR/cna_timing/", df$study_id[w], ".cnv.timeR.txt"), header=T, as.is=T, sep="\t")
    vecdf <- read.csv(paste0("../27_EBI/final_data_022222/timing_veclonal/", df$study_id[w], ".sage30.snv.timeR.VEClonal.txt"), header=T, as.is=T, sep="\t")
    vecdf <- vecdf[vecdf$kt == "not_analyzed",]
    timer$arm <- cnvdf$arm
    timer$loh <- cnvdf$loh
    timer$seglen <- timer$end-timer$start+1
    df$p_nseg[w] <- nrow(timer[timer$arm == paste0(str_sub(i, 1, nchar(i)-1), "p") & timer$seglen > lengthminimum,])
    df$q_nseg[w] <- nrow(timer[timer$arm == paste0(str_sub(i, 1, nchar(i)-1), "q") & timer$seglen > lengthminimum,])
    df$p_majcn[w] <- weighted.mean(timer$major_cn[timer$arm == paste0(str_sub(i, 1, nchar(i)-1), "p") & timer$seglen > lengthminimum], timer$seglen[timer$arm == paste0(str_sub(i, 1, nchar(i)-1), "p") & timer$seglen > lengthminimum], na.rm=T)
    df$q_majcn[w] <- weighted.mean(timer$major_cn[timer$arm == paste0(str_sub(i, 1, nchar(i)-1), "q") & timer$seglen > lengthminimum], timer$seglen[timer$arm == paste0(str_sub(i, 1, nchar(i)-1), "q") & timer$seglen > lengthminimum], na.rm=T)
    df$p_mincn[w] <- weighted.mean(timer$minor_cn[timer$arm == paste0(str_sub(i, 1, nchar(i)-1), "p") & timer$seglen > lengthminimum], timer$seglen[timer$arm == paste0(str_sub(i, 1, nchar(i)-1), "p") & timer$seglen > lengthminimum], na.rm=T)
    df$q_mincn[w] <- weighted.mean(timer$minor_cn[timer$arm == paste0(str_sub(i, 1, nchar(i)-1), "q") & timer$seglen > lengthminimum], timer$seglen[timer$arm == paste0(str_sub(i, 1, nchar(i)-1), "q") & timer$seglen > lengthminimum], na.rm=T)
    df$p_time[w] <- weighted.mean(timer$time1[timer$arm == paste0(str_sub(i, 1, nchar(i)-1), "p") & timer$seglen > lengthminimum], timer$seglen[timer$arm == paste0(str_sub(i, 1, nchar(i)-1), "p") & timer$seglen > lengthminimum], na.rm=T)
    df$q_time[w] <- weighted.mean(timer$time1[timer$arm == paste0(str_sub(i, 1, nchar(i)-1), "q") & timer$seglen > lengthminimum], timer$seglen[timer$arm == paste0(str_sub(i, 1, nchar(i)-1), "q") & timer$seglen > lengthminimum], na.rm=T)
    ndf <- timer[timer$arm == i & timer$seglen > lengthminimum,]
    if (nrow(ndf) != 0){
      ndf$n_veclonal <- 0
      ndf$n_snv_recount <- 0
      for (j in 1:nrow(ndf)){
        ndf$n_veclonal[j] <- nrow(vecdf[vecdf$chrom == ndf$chrom[j] & vecdf$pos >= ndf$start[j] & vecdf$pos <= ndf$end[j] & vecdf$pVEClonal > vecdf$pOtherClonal & vecdf$pVEClonal > vecdf$pSClonal,])
        ndf$n_snv_recount[j] <- nrow(vecdf[vecdf$chrom == ndf$chrom[j] & vecdf$pos >= ndf$start[j] & vecdf$pos <= ndf$end[j],])
      }
      ndf$d_veclonal <- ndf$n_veclonal*1000000/ndf$seglen
      df$d_veclonal[w] <- weighted.mean(ndf$d_veclonal, ndf$seglen, na.rm = T)
      rm(ndf)
    } 
  }
  if (str_sub(i, -1) == "p"){
    xdf <- df[(df$q_majcn+1 <= df$p_majcn) & df$q_mincn < 0.5 & df$q_nseg < 10 & !is.na(df$p_majcn) & !is.na(df$q_majcn) & !is.na(df$p_mincn) & !is.na(df$q_mincn) & df$p_time != "NaN", c("study_id", "is.wgd", "p_time", "d_veclonal")]
    colnames(xdf) <- c("study_id", "is.wgd", "time", "d_veclonal")
    xdf$arm <- i
    tempdf <- rbind(tempdf, xdf)
    rm(xdf)
  } else {
    xdf <- df[(df$q_majcn >= df$p_majcn+1) & df$q_nseg < 10 & !is.na(df$p_majcn) & !is.na(df$q_majcn) & !is.na(df$p_mincn) & !is.na(df$q_mincn) & df$q_time != "NaN", c("study_id", "is.wgd", "q_time", "d_veclonal")]
    colnames(xdf) <- c("study_id", "is.wgd", "time", "d_veclonal")
    xdf$arm <- i
    tempdf <- rbind(tempdf, xdf)
    rm(xdf)
  }
}
tempdf$arm <- factor(tempdf$arm, levels=c("1q", "16p", "8q", "20q"))
table(tempdf$arm, tempdf$is.wgd) ## Important information, consistent with previous papers indicating that 16q loss is associated with low-grade tumors.

df <- tmdf[tmdf$is.wgd=="WGD" & tmdf$timingclass=="sync", c("study_id", "is.wgd", "time.wgd")]
colnames(df) <- c("study_id", "is.wgd", "time")
df$d_veclonal <- NA
df$arm <- "WGD"
df$time_reassessed <- NA
for (w in 1:nrow(df)){
  cnvdf <- read.csv(paste0("../27_EBI/final_data_022222/acgr/", df$study_id[w], ".purple.cnv.somatic.acgr.tsv"), header=T, as.is=T, sep="\t")
  timer <- read.csv(paste0("../27_EBI/final_data_022222/timeR/cna_timing/", df$study_id[w], ".cnv.timeR.txt"), header=T, as.is=T, sep="\t")
  vecdf <- read.csv(paste0("../27_EBI/final_data_022222/timing_veclonal/", df$study_id[w], ".sage30.snv.timeR.VEClonal.txt"), header=T, as.is=T, sep="\t")
  vecdf <- vecdf[vecdf$kt == "not_analyzed",]
  timer$arm <- cnvdf$arm
  timer$loh <- cnvdf$loh
  timer$seglen <- timer$end-timer$start+1
  ndf <- timer[timer$seglen > lengthminimum,]
  if (nrow(ndf) != 0){
    ndf$n_veclonal <- 0
    ndf$n_snv_recount <- 0
    for (j in 1:nrow(ndf)){
      ndf$n_veclonal[j] <- nrow(vecdf[vecdf$chrom == ndf$chrom[j] & vecdf$pos >= ndf$start[j] & vecdf$pos <= ndf$end[j] & vecdf$pVEClonal > vecdf$pOtherClonal & vecdf$pVEClonal > vecdf$pSClonal,])
      ndf$n_snv_recount[j] <- nrow(vecdf[vecdf$chrom == ndf$chrom[j] & vecdf$pos >= ndf$start[j] & vecdf$pos <= ndf$end[j],])
    }
    ndf$d_veclonal <- ndf$n_veclonal*1000000/ndf$seglen
    df$time_reassessed[w] <- weighted.mean(ndf$time1[ndf$type == "Bi-allelic Gain (WGD)" & is.na(ndf$time2)], ndf$seglen[ndf$type == "Bi-allelic Gain (WGD)" & is.na(ndf$time2)], na.rm = T)
    df$d_veclonal[w] <- weighted.mean(ndf$d_veclonal[ndf$type == "Bi-allelic Gain (WGD)" & is.na(ndf$time2)], ndf$seglen[ndf$type == "Bi-allelic Gain (WGD)" & is.na(ndf$time2)], na.rm = T)
    rm(ndf)
  }
}

df <- df[,colnames(df) != "time_reassessed"]
df$d_veclonal_diploid <- df$d_veclonal
tempdf$d_veclonal_diploid <- 2*tempdf$d_veclonal
tempdf <- rbind(df, tempdf)
tempdf$arm <- factor(tempdf$arm, levels=c("WGD", "1q", "16p", "8q", "20q"))

aneutm <- tempdf
aneutm$er_plot <- "unknown"
aneutm$er_plot <- factor(aneutm$er_plot, levels = c("positive", "negative", "unknown"))
for (i in 1:nrow(aneutm)){
  aneutm$er_plot[i] <- demdf$er_plot[demdf$study_id == aneutm$study_id[i] & !is.na(demdf$er_plot)]
}
aneutm$purity <- 0
for (i in 1:nrow(aneutm)){
  aneutm$purity[i] <- purdf$purity[purdf$study_id == aneutm$study_id[i]]
}

aneutm$time.wgd <- NA
aneutm$lo.wgd <- NA
aneutm$hi.wgd <- NA
for (i in which(aneutm$is.wgd == "WGD" & aneutm$arm != "WGD")){
  if (tmdf$timingclass[tmdf$study_id == aneutm$study_id[i]] == "sync"){
    aneutm$time.wgd[i] <- tmdf$time.wgd[tmdf$study_id == aneutm$study_id[i]]
    aneutm$lo.wgd[i] <- tmdf$time.wgd[tmdf$study_id == aneutm$study_id[i]] - 2*tmdf$sd.wgd[tmdf$study_id == aneutm$study_id[i]]
    aneutm$hi.wgd[i] <- tmdf$time.wgd[tmdf$study_id == aneutm$study_id[i]] + 2*tmdf$sd.wgd[tmdf$study_id == aneutm$study_id[i]]
  }
}
aneutm$overlap.wgd <- "not_available"
for (i in which(aneutm$is.wgd == "WGD" & !is.na(aneutm$time.wgd))){
  if (!is.na(aneutm$time[i])){
    aneutm$overlap.wgd[i] <- ifelse(aneutm$time[i] >= aneutm$lo.wgd[i] & aneutm$time[i] <= aneutm$hi.wgd[i], "yes", "no")
  }
}
aneutm$arm <- factor(aneutm$arm, levels=c("WGD", "1q", "16p", "8q", "20q"))


## Hemichromosomal timing
nrow(comparms)
length(unique(bdclust$study_id[bdclust$ct_trakey != 0 & bdclust$allsvCount >= 10 & bdclust$f_trakeyarm >= 0.7]))
nrow(bdclust[bdclust$ct_trakey != 0 & bdclust$allsvCount >= 10 & bdclust$f_trakeyarm >= 0.7,])
nrow(bdclust[bdclust$ct_trakey != 0 & bdclust$allsvCount >= 10,])
length(unique(bdclust$study_id[bdclust$ct_trakey != 0]))
nrow(bdclust[bdclust$ct_trakey != 0,])
nrow(bdclust[bdclust$ct_trakey != 0 & bdclust$allsvCount >= 10,])

df <- tbclust[tbclust$study_id %in% purdf$study_id[purdf$availability == "yes" & purdf$timer == "yes" & purdf$veclonal == "yes"],]
hemichr <- data.frame(matrix(vector(), 0, 13, dimnames = list(c(), c("study_id", "pair", "ct_tra_pair", "hemi1", "len_hemi1", "cn_hemi1", "dvec_hemi1", "time_hemi1", "hemi2", "len_hemi2", "cn_hemi2", "dvec_hemi2", "time_hemi2"))), stringsAsFactors = F)
for (w in unique(df$study_id)){
  print(paste0(w))
  svdf <- read.csv(paste0("../27_EBI/final_data_022222/rearrangement_bedpe_linx/", w, ".purple.sv.linx.bedpe.tsv"), header=T, as.is=T, sep="\t")
  svdf$arm1 <- "no"
  svdf$arm2 <- "no"
  for (i in 1:nrow(svdf)){
    svdf$arm1[i] <- paste0(svdf$chrom1[i], ifelse((svdf$start1[i]+svdf$end1[i])/2 < hg19_coord$centro_mid[hg19_coord$chr == svdf$chrom1[i]], "p", "q"))
    svdf$arm2[i] <- paste0(svdf$chrom2[i], ifelse((svdf$start2[i]+svdf$end2[i])/2 < hg19_coord$centro_mid[hg19_coord$chr == svdf$chrom2[i]], "p", "q"))
  }
  cnvdf <- read.csv(paste0("../27_EBI/final_data_022222/acgr/", w, ".purple.cnv.somatic.acgr.tsv"), header=T, as.is=T, sep="\t")
  timer <- read.csv(paste0("../27_EBI/final_data_022222/timeR/cna_timing/", w, ".cnv.timeR.txt"), header=T, as.is=T, sep="\t")
  vecdf <- read.csv(paste0("../27_EBI/final_data_022222/timing_veclonal/", w, ".sage30.snv.timeR.VEClonal.txt"), header=T, as.is=T, sep="\t")
  vecdf <- vecdf[vecdf$kt == "not_analyzed",]
  timer$arm <- cnvdf$arm
  timer$loh <- cnvdf$loh
  timer$seglen <- timer$end-timer$start+1
  timer$cn <- cnvdf$copyNumber
  timer$blacklist <- "no"
  for (i in 1:nrow(timer)){
    if (timer$start[i] == 1 & timer$end[i] >10000000 & timer$end[i] <20000000 & timer$chrom[i] %in% c("13", "14", "15", "21", "22")){
      timer$blacklist[i] <- "acrocentric"
    }
    if (timer$chrom[i] == "1" & timer$start[i] == 123035434){
      timer$blacklist[i] <- "centro_gap"
    } else if (timer$chrom[i] == "9" & (timer$end[i] == 48867678 | timer$start[i] == 48867679)){
      timer$blacklist[i] <- "centro_gap"
    } else if (timer$chrom[i] == "16" & timer$start[i] == 36835801){
      timer$blacklist[i] <- "centro_gap"
    }
  }
  timer <- timer[timer$blacklist == "no",]
  threshold <- ifelse(purdf$is.wgd[purdf$study_id == w] == "WGD", 3, 2)
  for (j in unique(sub("\\[.*", "", strsplit(paste(df$pair_trakey[df$study_id == w], collapse = ";"), "];", fixed=T)[[1]]))){
    hemichr[nrow(hemichr)+1,] <- NA
    hemichr$study_id[nrow(hemichr)] <- w
    hemichr$pair[nrow(hemichr)] <- j
    l <- strsplit(j, ";", fixed=T)[[1]]
    hemichr$ct_tra_pair[nrow(hemichr)] <- nrow(svdf[svdf$arm1 == l[1] & svdf$arm2 == l[2],]) + nrow(svdf[svdf$arm1 == l[2] & svdf$arm2 == l[1],])
    m <- gsub("a", "q", gsub("q", "p", gsub("p", "a", l)))
    hemichr$hemi1[nrow(hemichr)] <- m[1]
    hemichr$hemi2[nrow(hemichr)] <- m[2]
    #hemichr$floh_hemi1[nrow(hemichr)] <- sum(timer$seglen[timer$arm == m[1] & timer$loh == "loh"])
    #hemichr$floh_hemi2[nrow(hemichr)] <- sum(timer$seglen[timer$arm == m[2] & timer$loh == "loh"])
    if (nrow(timer[timer$arm %in% m & timer$seglen > lengthminimum,]) != 0){
      ndf <- timer[timer$arm %in% m & timer$seglen > lengthminimum,]
      ndf$n_veclonal <- 0
      for (z in 1:nrow(ndf)){
        ndf$n_veclonal[z] <- nrow(vecdf[vecdf$chrom == ndf$chrom[z] & vecdf$pos >= ndf$start[z] & vecdf$pos <= ndf$end[z] & vecdf$pVEClonal > vecdf$pOtherClonal & vecdf$pVEClonal > vecdf$pSClonal,])
      }
      ndf$d_veclonal <- ndf$n_veclonal*1000000/ndf$seglen
      hemichr$len_hemi1[nrow(hemichr)] <- sum(ndf$seglen[ndf$arm == m[1] & ndf$major_cn >= threshold])
      hemichr$cn_hemi1[nrow(hemichr)] <- weighted.mean(ndf$cn[ndf$arm == m[1] & ndf$major_cn >= threshold], ndf$seglen[ndf$arm == m[1] & ndf$major_cn >= threshold], na.rm = T)
      hemichr$dvec_hemi1[nrow(hemichr)] <- weighted.mean(ndf$d_veclonal[ndf$arm == m[1] & ndf$major_cn >= threshold], ndf$seglen[ndf$arm == m[1] & ndf$major_cn >= threshold], na.rm = T)
      hemichr$time_hemi1[nrow(hemichr)] <- weighted.mean(ndf$time1[ndf$arm == m[1] & ndf$major_cn >= threshold], ndf$seglen[ndf$arm == m[1] & ndf$major_cn >= threshold], na.rm = T)
      hemichr$len_hemi2[nrow(hemichr)] <- sum(ndf$seglen[ndf$arm == m[2] & ndf$major_cn >= threshold])
      hemichr$cn_hemi2[nrow(hemichr)] <- weighted.mean(ndf$cn[ndf$arm == m[2] & ndf$major_cn >= threshold], ndf$seglen[ndf$arm == m[2] & ndf$major_cn >= threshold], na.rm = T)
      hemichr$dvec_hemi2[nrow(hemichr)] <- weighted.mean(ndf$d_veclonal[ndf$arm == m[2] & ndf$major_cn >= threshold], ndf$seglen[ndf$arm == m[2] & ndf$major_cn >= threshold], na.rm = T)
      hemichr$time_hemi2[nrow(hemichr)] <- weighted.mean(ndf$time1[ndf$arm == m[2] & ndf$major_cn >= threshold], ndf$seglen[ndf$arm == m[2] & ndf$major_cn >= threshold], na.rm = T)
      rm(ndf)
    } else {
      hemichr$len_hemi1[nrow(hemichr)] <- NA
      hemichr$cn_hemi1[nrow(hemichr)] <- NA
      hemichr$dvec_hemi1[nrow(hemichr)] <- NA
      hemichr$time_hemi1[nrow(hemichr)] <- NA
      hemichr$len_hemi2[nrow(hemichr)] <- NA
      hemichr$cn_hemi2[nrow(hemichr)] <- NA
      hemichr$dvec_hemi2[nrow(hemichr)] <- NA
      hemichr$time_hemi2[nrow(hemichr)] <- NA
    }
  }
}
hemichr$is.wgd <- ""
for (i in 1:nrow(hemichr)){
  hemichr$is.wgd[i] <- purdf$is.wgd[purdf$study_id == hemichr$study_id[i]]
}
hemichr$time.wgd <- NA
hemichr$lo.wgd <- NA
hemichr$hi.wgd <- NA
for (i in which(hemichr$is.wgd == "WGD")){
  if (tmdf$timingclass[tmdf$study_id == hemichr$study_id[i]] == "sync"){
    hemichr$time.wgd[i] <- tmdf$time.wgd[tmdf$study_id == hemichr$study_id[i]]
    hemichr$lo.wgd[i] <- tmdf$time.wgd[tmdf$study_id == hemichr$study_id[i]] - 2*tmdf$sd.wgd[tmdf$study_id == hemichr$study_id[i]]
    hemichr$hi.wgd[i] <- tmdf$time.wgd[tmdf$study_id == hemichr$study_id[i]] + 2*tmdf$sd.wgd[tmdf$study_id == hemichr$study_id[i]]
  }
}
hemichr$syn.hemi1 <- "not_available"
hemichr$syn.hemi2 <- "not_available"
for (i in which(hemichr$is.wgd == "WGD" & !is.na(hemichr$time.wgd))){
  if (!is.na(hemichr$time_hemi1[i])){
    hemichr$syn.hemi1[i] <- ifelse(hemichr$time_hemi1[i] >= hemichr$lo.wgd[i] & hemichr$time_hemi1[i] <= hemichr$hi.wgd[i], "yes", "no")
  }
  if (!is.na(hemichr$time_hemi2[i])){
    hemichr$syn.hemi2[i] <- ifelse(hemichr$time_hemi2[i] >= hemichr$lo.wgd[i] & hemichr$time_hemi2[i] <= hemichr$hi.wgd[i], "yes", "no")
  }
}

df <- hemichr[hemichr$ct_tra_pair > 5 & (!is.na(hemichr$dvec_hemi1) | !is.na(hemichr$dvec_hemi2)),]
tbatime <- data.frame(matrix(vector(), 0, 6, dimnames=list(c(), c("study_id", "is.wgd", "time", "d_veclonal", "arm", "pair"))), stringsAsFactors=F)
for (i in 1:nrow(df)){
  if (is.na(df$time.wgd[i])){
    if (!is.na(df$dvec_hemi1[i]) & is.na(df$dvec_hemi2[i])){
      tbatime[nrow(tbatime)+1,] <- df[i,c(1,14,8,7,4,2)]
    } else if (is.na(df$dvec_hemi1[i]) & !is.na(df$dvec_hemi2[i])){
      tbatime[nrow(tbatime)+1,] <- df[i,c(1,14,13,12,9,2)]
    } else if (!is.na(df$dvec_hemi1[i]) & !is.na(df$dvec_hemi2[i]) & df$dvec_hemi1[i] < df$dvec_hemi2[i]){
      tbatime[nrow(tbatime)+1,] <- df[i,c(1,14,8,7,4,2)]
    } else if (!is.na(df$dvec_hemi1[i]) & !is.na(df$dvec_hemi2[i]) & df$dvec_hemi1[i] > df$dvec_hemi2[i]){
      tbatime[nrow(tbatime)+1,] <- df[i,c(1,14,13,12,9,2)]
    }
  } else {
    if (df$syn.hemi1[i] == "yes" & df$syn.hemi2[i] == "yes"){
      
    } else if (df$syn.hemi1[i] != "yes" & df$syn.hemi2[i] == "yes"){
      if (!is.na(df$dvec_hemi1[i])){
        tbatime[nrow(tbatime)+1,] <- df[i,c(1,14,8,7,4,2)]
      }
    } else if (df$syn.hemi1[i] == "yes" & df$syn.hemi2[i] != "yes"){
      if (!is.na(df$dvec_hemi2[i])){
        tbatime[nrow(tbatime)+1,] <- df[i,c(1,14,13,12,9,2)]
      }
    } else if (df$syn.hemi1[i] != "yes" & df$syn.hemi2[i] != "yes"){
      if (!is.na(df$dvec_hemi1[i]) & is.na(df$dvec_hemi2[i])){
        tbatime[nrow(tbatime)+1,] <- df[i,c(1,14,8,7,4,2)]
      } else if (is.na(df$dvec_hemi1[i]) & !is.na(df$dvec_hemi2[i])){
        tbatime[nrow(tbatime)+1,] <- df[i,c(1,14,13,12,9,2)]
      } else if (!is.na(df$dvec_hemi1[i]) & !is.na(df$dvec_hemi2[i]) & df$dvec_hemi1[i] < df$dvec_hemi2[i]){
        tbatime[nrow(tbatime)+1,] <- df[i,c(1,14,8,7,4,2)]
      } else if (!is.na(df$dvec_hemi1[i]) & !is.na(df$dvec_hemi2[i]) & df$dvec_hemi1[i] > df$dvec_hemi2[i]){
        tbatime[nrow(tbatime)+1,] <- df[i,c(1,14,13,12,9,2)]
      }
    }
  }
}
tbatime$d_veclonal_diploid <- tbatime$d_veclonal*2
tbatime$idseg <- paste0(tbatime$study_id, ":", tbatime$arm)
tbatime$duplicate <- 1
for (i in 1:nrow(tbatime)){
  tbatime$duplicate[i] <- nrow(tbatime[tbatime$idseg == tbatime$idseg[i],])
}
for (i in 1:nrow(tbatime)){
  if (tbatime$duplicate[i] == 2){
    tbatime$duplicate[tbatime$idseg == tbatime$idseg[i]][2] <- 1
  }
}
tbatime <- tbatime[tbatime$duplicate == 1,]

aneutm$tbamp <- "no"
for (i in 1:nrow(aneutm)){
  if (nrow(tbatime[tbatime$study_id == aneutm$study_id[i] & tbatime$arm == aneutm$arm[i],]) != 0){
    aneutm$tbamp[i] <- "yes"
  }
}

## Beeswarm plot for the timing of common arm-level aneuploidies
pdf(paste0("FIP_revision/beeswarm.absolute.timing.chroms.5.all.cases.diploid.kt.excluded.expanded.", lengthminimum/1000000, "k.no.overlap.wgd.new.pdf"), height=4.5, width=10) # Current version
beeswarm(aneutm$d_veclonal_diploid[aneutm$tbamp == "no" & aneutm$overlap.wgd != "yes"] ~ aneutm$arm[aneutm$tbamp == "no" & aneutm$overlap.wgd != "yes"], method='hex', pch=20, las=1, col=c("#FF9912", "#1874CD", "#00868B", "#27408B", "#607B8B"), ylab="Mutational burden at the time of CNA (per Mbp)", xlab="Chromosome arm", cex=0.7, ylim=c(0,2))
dev.off()

pdf(paste0("FIP_revision/beeswarm.tbatime.kt.excluded.expanded.", lengthminimum/1000000, "k.pdf"), height=4.5, width=4)
beeswarm(tbatime$d_veclonal_diploid, method='hex', pch=20, las=1, col=c("#800080"), ylab="Mutational burden at the time of CNA (per Mbp)", xlab="TB Amplifications", cex=0.7, ylim=c(0,2))
dev.off()

## Mutational signature of VEClonal mutations
msdf <- aneutm[aneutm$tbamp == "no" & aneutm$overlap.wgd != "yes" & aneutm$arm %in% c("1q", "16p") & aneutm$is.wgd == "ND",]
msdf$ct_vecmut_on_arm <- 0
msdf$ACAA <- 0
msdf$ACAC <- 0
msdf$ACAG <- 0
msdf$ACAT <- 0
msdf$CCAA <- 0
msdf$CCAC <- 0
msdf$CCAG <- 0
msdf$CCAT <- 0
msdf$GCAA <- 0
msdf$GCAC <- 0
msdf$GCAG <- 0
msdf$GCAT <- 0
msdf$TCAA <- 0
msdf$TCAC <- 0
msdf$TCAG <- 0
msdf$TCAT <- 0
msdf$ACGA <- 0
msdf$ACGC <- 0
msdf$ACGG <- 0
msdf$ACGT <- 0
msdf$CCGA <- 0
msdf$CCGC <- 0
msdf$CCGG <- 0
msdf$CCGT <- 0
msdf$GCGA <- 0
msdf$GCGC <- 0
msdf$GCGG <- 0
msdf$GCGT <- 0
msdf$TCGA <- 0
msdf$TCGC <- 0
msdf$TCGG <- 0
msdf$TCGT <- 0
msdf$ACTA <- 0
msdf$ACTC <- 0
msdf$ACTG <- 0
msdf$ACTT <- 0
msdf$CCTA <- 0
msdf$CCTC <- 0
msdf$CCTG <- 0
msdf$CCTT <- 0
msdf$GCTA <- 0
msdf$GCTC <- 0
msdf$GCTG <- 0
msdf$GCTT <- 0
msdf$TCTA <- 0
msdf$TCTC <- 0
msdf$TCTG <- 0
msdf$TCTT <- 0
msdf$ATAA <- 0
msdf$ATAC <- 0
msdf$ATAG <- 0
msdf$ATAT <- 0
msdf$CTAA <- 0
msdf$CTAC <- 0
msdf$CTAG <- 0
msdf$CTAT <- 0
msdf$GTAA <- 0
msdf$GTAC <- 0
msdf$GTAG <- 0
msdf$GTAT <- 0
msdf$TTAA <- 0
msdf$TTAC <- 0
msdf$TTAG <- 0
msdf$TTAT <- 0
msdf$ATCA <- 0
msdf$ATCC <- 0
msdf$ATCG <- 0
msdf$ATCT <- 0
msdf$CTCA <- 0
msdf$CTCC <- 0
msdf$CTCG <- 0
msdf$CTCT <- 0
msdf$GTCA <- 0
msdf$GTCC <- 0
msdf$GTCG <- 0
msdf$GTCT <- 0
msdf$TTCA <- 0
msdf$TTCC <- 0
msdf$TTCG <- 0
msdf$TTCT <- 0
msdf$ATGA <- 0
msdf$ATGC <- 0
msdf$ATGG <- 0
msdf$ATGT <- 0
msdf$CTGA <- 0
msdf$CTGC <- 0
msdf$CTGG <- 0
msdf$CTGT <- 0
msdf$GTGA <- 0
msdf$GTGC <- 0
msdf$GTGG <- 0
msdf$GTGT <- 0
msdf$TTGA <- 0
msdf$TTGC <- 0
msdf$TTGG <- 0
msdf$TTGT <- 0
for (w in 1:nrow(msdf)){
  print(paste0(msdf$study_id[w], " ", msdf$arm[w]))
  cnvdf <- read.csv(paste0("../27_EBI/final_data_022222/acgr/", msdf$study_id[w], ".purple.cnv.somatic.acgr.tsv"), header=T, as.is=T, sep="\t")
  timer <- read.csv(paste0("../27_EBI/final_data_022222/timeR/cna_timing/", msdf$study_id[w], ".cnv.timeR.txt"), header=T, as.is=T, sep="\t")
  vecdf <- read.csv(paste0("../27_EBI/final_data_022222/timing_veclonal/", msdf$study_id[w], ".sage30.snv.timeR.VEClonal.txt"), header=T, as.is=T, sep="\t")
  vecdf <- vecdf[vecdf$kt == "not_analyzed",]
  timer$arm <- cnvdf$arm
  timer$loh <- cnvdf$loh
  timer$seglen <- timer$end-timer$start+1
  ndf <- timer[timer$arm == as.character(msdf$arm[w]) & timer$seglen > lengthminimum,]
  ndf <- ndf[ndf$major_cn >=2 & ndf$minor_cn < ndf$major_cn,]
  if (nrow(ndf) != 0){
    vecdf <- vecdf[vecdf$chrom == str_sub(msdf$arm[w], 1, nchar(as.character(msdf$arm[w]))-1),]
    vecdf$inclusion <- "no"
    for (i in 1:nrow(vecdf)){
      if (nrow(ndf[ndf$chrom == vecdf$chrom[i] & ndf$start <= vecdf$pos[i] & ndf$end >= vecdf$pos[i],]) != 0){
        vecdf$inclusion[i] <- "yes"
      }
    }
    vecdf <- vecdf[vecdf$inclusion == "yes" & vecdf$pVEClonal < vecdf$pSClonal & vecdf$pOtherClonal < vecdf$pSClonal,] # subclonal
    msdf$ct_vecmut_on_arm[w] <- nrow(vecdf)
    msdf$ACAA[w] <- nrow(vecdf[vecdf$mutspectra == "A[C>A]A",]) 
    msdf$ACAC[w] <- nrow(vecdf[vecdf$mutspectra == "A[C>A]C",]) 
    msdf$ACAG[w] <- nrow(vecdf[vecdf$mutspectra == "A[C>A]G",]) 
    msdf$ACAT[w] <- nrow(vecdf[vecdf$mutspectra == "A[C>A]T",]) 
    msdf$CCAA[w] <- nrow(vecdf[vecdf$mutspectra == "C[C>A]A",]) 
    msdf$CCAC[w] <- nrow(vecdf[vecdf$mutspectra == "C[C>A]C",]) 
    msdf$CCAG[w] <- nrow(vecdf[vecdf$mutspectra == "C[C>A]G",]) 
    msdf$CCAT[w] <- nrow(vecdf[vecdf$mutspectra == "C[C>A]T",]) 
    msdf$GCAA[w] <- nrow(vecdf[vecdf$mutspectra == "G[C>A]A",]) 
    msdf$GCAC[w] <- nrow(vecdf[vecdf$mutspectra == "G[C>A]C",]) 
    msdf$GCAG[w] <- nrow(vecdf[vecdf$mutspectra == "G[C>A]G",]) 
    msdf$GCAT[w] <- nrow(vecdf[vecdf$mutspectra == "G[C>A]T",]) 
    msdf$TCAA[w] <- nrow(vecdf[vecdf$mutspectra == "T[C>A]A",]) 
    msdf$TCAC[w] <- nrow(vecdf[vecdf$mutspectra == "T[C>A]C",]) 
    msdf$TCAG[w] <- nrow(vecdf[vecdf$mutspectra == "T[C>A]G",]) 
    msdf$TCAT[w] <- nrow(vecdf[vecdf$mutspectra == "T[C>A]T",]) 
    msdf$ACGA[w] <- nrow(vecdf[vecdf$mutspectra == "A[C>G]A",]) 
    msdf$ACGC[w] <- nrow(vecdf[vecdf$mutspectra == "A[C>G]C",]) 
    msdf$ACGG[w] <- nrow(vecdf[vecdf$mutspectra == "A[C>G]G",]) 
    msdf$ACGT[w] <- nrow(vecdf[vecdf$mutspectra == "A[C>G]T",]) 
    msdf$CCGA[w] <- nrow(vecdf[vecdf$mutspectra == "C[C>G]A",]) 
    msdf$CCGC[w] <- nrow(vecdf[vecdf$mutspectra == "C[C>G]C",]) 
    msdf$CCGG[w] <- nrow(vecdf[vecdf$mutspectra == "C[C>G]G",]) 
    msdf$CCGT[w] <- nrow(vecdf[vecdf$mutspectra == "C[C>G]T",]) 
    msdf$GCGA[w] <- nrow(vecdf[vecdf$mutspectra == "G[C>G]A",]) 
    msdf$GCGC[w] <- nrow(vecdf[vecdf$mutspectra == "G[C>G]C",]) 
    msdf$GCGG[w] <- nrow(vecdf[vecdf$mutspectra == "G[C>G]G",]) 
    msdf$GCGT[w] <- nrow(vecdf[vecdf$mutspectra == "G[C>G]T",]) 
    msdf$TCGA[w] <- nrow(vecdf[vecdf$mutspectra == "T[C>G]A",]) 
    msdf$TCGC[w] <- nrow(vecdf[vecdf$mutspectra == "T[C>G]C",]) 
    msdf$TCGG[w] <- nrow(vecdf[vecdf$mutspectra == "T[C>G]G",]) 
    msdf$TCGT[w] <- nrow(vecdf[vecdf$mutspectra == "T[C>G]T",]) 
    msdf$ACTA[w] <- nrow(vecdf[vecdf$mutspectra == "A[C>T]A",]) 
    msdf$ACTC[w] <- nrow(vecdf[vecdf$mutspectra == "A[C>T]C",]) 
    msdf$ACTG[w] <- nrow(vecdf[vecdf$mutspectra == "A[C>T]G",]) 
    msdf$ACTT[w] <- nrow(vecdf[vecdf$mutspectra == "A[C>T]T",]) 
    msdf$CCTA[w] <- nrow(vecdf[vecdf$mutspectra == "C[C>T]A",]) 
    msdf$CCTC[w] <- nrow(vecdf[vecdf$mutspectra == "C[C>T]C",]) 
    msdf$CCTG[w] <- nrow(vecdf[vecdf$mutspectra == "C[C>T]G",]) 
    msdf$CCTT[w] <- nrow(vecdf[vecdf$mutspectra == "C[C>T]T",]) 
    msdf$GCTA[w] <- nrow(vecdf[vecdf$mutspectra == "G[C>T]A",]) 
    msdf$GCTC[w] <- nrow(vecdf[vecdf$mutspectra == "G[C>T]C",]) 
    msdf$GCTG[w] <- nrow(vecdf[vecdf$mutspectra == "G[C>T]G",]) 
    msdf$GCTT[w] <- nrow(vecdf[vecdf$mutspectra == "G[C>T]T",]) 
    msdf$TCTA[w] <- nrow(vecdf[vecdf$mutspectra == "T[C>T]A",]) 
    msdf$TCTC[w] <- nrow(vecdf[vecdf$mutspectra == "T[C>T]C",]) 
    msdf$TCTG[w] <- nrow(vecdf[vecdf$mutspectra == "T[C>T]G",]) 
    msdf$TCTT[w] <- nrow(vecdf[vecdf$mutspectra == "T[C>T]T",]) 
    msdf$ATAA[w] <- nrow(vecdf[vecdf$mutspectra == "A[T>A]A",]) 
    msdf$ATAC[w] <- nrow(vecdf[vecdf$mutspectra == "A[T>A]C",]) 
    msdf$ATAG[w] <- nrow(vecdf[vecdf$mutspectra == "A[T>A]G",]) 
    msdf$ATAT[w] <- nrow(vecdf[vecdf$mutspectra == "A[T>A]T",]) 
    msdf$CTAA[w] <- nrow(vecdf[vecdf$mutspectra == "C[T>A]A",]) 
    msdf$CTAC[w] <- nrow(vecdf[vecdf$mutspectra == "C[T>A]C",]) 
    msdf$CTAG[w] <- nrow(vecdf[vecdf$mutspectra == "C[T>A]G",]) 
    msdf$CTAT[w] <- nrow(vecdf[vecdf$mutspectra == "C[T>A]T",]) 
    msdf$GTAA[w] <- nrow(vecdf[vecdf$mutspectra == "G[T>A]A",]) 
    msdf$GTAC[w] <- nrow(vecdf[vecdf$mutspectra == "G[T>A]C",]) 
    msdf$GTAG[w] <- nrow(vecdf[vecdf$mutspectra == "G[T>A]G",]) 
    msdf$GTAT[w] <- nrow(vecdf[vecdf$mutspectra == "G[T>A]T",]) 
    msdf$TTAA[w] <- nrow(vecdf[vecdf$mutspectra == "T[T>A]A",]) 
    msdf$TTAC[w] <- nrow(vecdf[vecdf$mutspectra == "T[T>A]C",]) 
    msdf$TTAG[w] <- nrow(vecdf[vecdf$mutspectra == "T[T>A]G",]) 
    msdf$TTAT[w] <- nrow(vecdf[vecdf$mutspectra == "T[T>A]T",]) 
    msdf$ATCA[w] <- nrow(vecdf[vecdf$mutspectra == "A[T>C]A",]) 
    msdf$ATCC[w] <- nrow(vecdf[vecdf$mutspectra == "A[T>C]C",]) 
    msdf$ATCG[w] <- nrow(vecdf[vecdf$mutspectra == "A[T>C]G",]) 
    msdf$ATCT[w] <- nrow(vecdf[vecdf$mutspectra == "A[T>C]T",]) 
    msdf$CTCA[w] <- nrow(vecdf[vecdf$mutspectra == "C[T>C]A",]) 
    msdf$CTCC[w] <- nrow(vecdf[vecdf$mutspectra == "C[T>C]C",]) 
    msdf$CTCG[w] <- nrow(vecdf[vecdf$mutspectra == "C[T>C]G",]) 
    msdf$CTCT[w] <- nrow(vecdf[vecdf$mutspectra == "C[T>C]T",]) 
    msdf$GTCA[w] <- nrow(vecdf[vecdf$mutspectra == "G[T>C]A",]) 
    msdf$GTCC[w] <- nrow(vecdf[vecdf$mutspectra == "G[T>C]C",]) 
    msdf$GTCG[w] <- nrow(vecdf[vecdf$mutspectra == "G[T>C]G",]) 
    msdf$GTCT[w] <- nrow(vecdf[vecdf$mutspectra == "G[T>C]T",]) 
    msdf$TTCA[w] <- nrow(vecdf[vecdf$mutspectra == "T[T>C]A",]) 
    msdf$TTCC[w] <- nrow(vecdf[vecdf$mutspectra == "T[T>C]C",]) 
    msdf$TTCG[w] <- nrow(vecdf[vecdf$mutspectra == "T[T>C]G",]) 
    msdf$TTCT[w] <- nrow(vecdf[vecdf$mutspectra == "T[T>C]T",]) 
    msdf$ATGA[w] <- nrow(vecdf[vecdf$mutspectra == "A[T>G]A",]) 
    msdf$ATGC[w] <- nrow(vecdf[vecdf$mutspectra == "A[T>G]C",]) 
    msdf$ATGG[w] <- nrow(vecdf[vecdf$mutspectra == "A[T>G]G",]) 
    msdf$ATGT[w] <- nrow(vecdf[vecdf$mutspectra == "A[T>G]T",]) 
    msdf$CTGA[w] <- nrow(vecdf[vecdf$mutspectra == "C[T>G]A",]) 
    msdf$CTGC[w] <- nrow(vecdf[vecdf$mutspectra == "C[T>G]C",]) 
    msdf$CTGG[w] <- nrow(vecdf[vecdf$mutspectra == "C[T>G]G",]) 
    msdf$CTGT[w] <- nrow(vecdf[vecdf$mutspectra == "C[T>G]T",]) 
    msdf$GTGA[w] <- nrow(vecdf[vecdf$mutspectra == "G[T>G]A",]) 
    msdf$GTGC[w] <- nrow(vecdf[vecdf$mutspectra == "G[T>G]C",]) 
    msdf$GTGG[w] <- nrow(vecdf[vecdf$mutspectra == "G[T>G]G",]) 
    msdf$GTGT[w] <- nrow(vecdf[vecdf$mutspectra == "G[T>G]T",]) 
    msdf$TTGA[w] <- nrow(vecdf[vecdf$mutspectra == "T[T>G]A",]) 
    msdf$TTGC[w] <- nrow(vecdf[vecdf$mutspectra == "T[T>G]C",]) 
    msdf$TTGG[w] <- nrow(vecdf[vecdf$mutspectra == "T[T>G]G",]) 
    msdf$TTGT[w] <- nrow(vecdf[vecdf$mutspectra == "T[T>G]T",])
  }
}
sumspec <- c()
for (i in c(15:110)){
  sumspec[i-14] <- sum(msdf[,i])
}

pdf("FIP_revision/Early_mutspectra/mutspectra.1q.16p.total.veclonal.pdf", height=4, width=8)
barplot(sum_veclonal_spec, border=NA, las=1, col=sigcolor, space=0.8, main=paste0("number of mutations = ", sum(sumspec)), ylab="Number of mutations", xlab="Trinucleotide context")
dev.off()

pdf("FIP_revision/Early_mutspectra/mutspectra.1q.16p.total.otherclonal.pdf", height=4, width=8)
barplot(sum_otherclonal_spec, border=NA, las=1, col=sigcolor, space=0.8, main=paste0("number of mutations = ", sum(sumspec)), ylab="Number of mutations", xlab="Trinucleotide context")
dev.off()

pdf("FIP_revision/Early_mutspectra/mutspectra.1q.16p.total.subclonal.pdf", height=4, width=8)
barplot(sum_subclonal_spec, border=NA, las=1, col=sigcolor, space=0.8, main=paste0("number of mutations = ", sum(sumspec)), ylab="Number of mutations", xlab="Trinucleotide context")
dev.off()

## Correlation between mutation burden and age
tmdf$ct_clonal_early <- NA
tmdf$ct_clonal_late <- NA
tmdf$ct_clonal_unspec <- NA
tmdf$ct_subclonal <- NA
for (i in 1:nrow(tmdf)){
  if (file.exists(paste0("../27_EBI/final_data_022222/timing_veclonal/", tmdf$study_id[i], ".sage30.snv.timeR.VEClonal.txt"))){
    vecdf <- read.csv(paste0("../27_EBI/final_data_022222/timing_veclonal/", tmdf$study_id[i], ".sage30.snv.timeR.VEClonal.txt"), header=T, as.is=T, sep="\t")
    tmdf$ct_clonal_early[i] <- nrow(vecdf[vecdf$CLS == "clonal [early]",])
    tmdf$ct_clonal_late[i] <- nrow(vecdf[vecdf$CLS == "clonal [late]",])
    tmdf$ct_clonal_unspec[i] <- nrow(vecdf[vecdf$CLS == "clonal [NA]",])
    tmdf$ct_subclonal[i] <- nrow(vecdf[vecdf$CLS == "subclonal",])
  }
}
median(tmdf$ct_clonal_early+tmdf$ct_clonal_late+tmdf$ct_clonal_unspec)
median(tmdf$age[!is.na(tmdf$age)])
tmdf$ct_clonal_all <- tmdf$ct_clonal_early+tmdf$ct_clonal_late+tmdf$ct_clonal_unspec

tmdf$f_apobec <- 0
tmdf$f_clock <- 0
for (i in 1:nrow(tmdf)){
  tmdf$f_apobec[i] <- (sigdf$SBS2[sigdf$study_id == tmdf$study_id[i]] + sigdf$SBS13[sigdf$study_id == tmdf$study_id[i]])/sigdf$ct_snv[sigdf$study_id == tmdf$study_id[i]]
  tmdf$f_clock[i] <- (sigdf$SBS1[sigdf$study_id == tmdf$study_id[i]] + sigdf$SBS5[sigdf$study_id == tmdf$study_id[i]])/sigdf$ct_snv[sigdf$study_id == tmdf$study_id[i]]
}

## Linear model based on gradual mutagenesis in the early evolution
### Using clonal mutation count
df <- tmdf[tmdf$hr_status == "HR_proficient" & tmdf$is.wgd == "ND" & tmdf$qc == "PASS" & tmdf$f_apobec < 0.5 & tmdf$purity >= 0.6 & tmdf$ct_snv < 10000 & tmdf$msStatus == "MSS",]
plot(df$ct_clonal_all ~ df$age, pch=20, col=onecolor, las=1, xlab="Age at diagnosis", ylab="Number of clonal mutations", xlim=c(0,100), ylim=c(0,10000), frame=F)
abline(lm(df$ct_clonal_all ~ df$age))
lm(df$ct_clonal_all ~ df$age)
median(tbatime$d_veclonal_diploid) * 2859/as.vector(lm(df$ct_clonal_all ~ df$age)$coefficients[2]) # TB amplification happens before 51 years of age. -- Used in the text
lm(df$ct_clonal_all ~ df$age)$coefficients[2]
summary(lm(df$ct_clonal_all ~ df$age))

### Using total SNV count
median(tbatime$d_veclonal_diploid) * 2859/as.vector(lm(df$ct_snv ~ df$age)$coefficients[2]) # TB amplification happens before 51 years of age. -- Used in the text
lm(df$ct_snv ~ df$age)$coefficients[2]
summary(lm(df$ct_snv ~ df$age))

pdf("FIP_revision/scatter.age.clonal.mutations.limit.hypermutation.purity.0.6.pdf", height=5, width=4.5)
plot(df$ct_clonal_all ~ df$age, pch=20, col=c(rgb(219/255, 112/255, 147/255, 0.7), rgb(238/255, 213/255, 210/255, 0.8), rgb(0,0,0,.3))[df$er_plot], las=1, xlab="Age at diagnosis", ylab="Number of clonal mutations", xlim=c(0,100), ylim=c(0,10000), frame=F)
abline(lm(df$ct_clonal_all ~ df$age))
dev.off()

pdf("FIP_revision/scatter.age.all.snvs.limit.hypermutation.purity.0.6.pdf", height=5, width=4.5)
plot(df$ct_snv ~ df$age, pch=20, col=c(rgb(219/255, 112/255, 147/255, 0.7), rgb(238/255, 213/255, 210/255, 0.8), rgb(0,0,0,.3))[df$er_plot], las=1, xlab="Age at diagnosis", ylab="Number of SNVs", xlim=c(0,100), ylim=c(0,10000), frame=F)
abline(lm(df$ct_snv ~ df$age))
dev.off()



