library(dplyr)
library(grid)
library(fields)
library(stringr)
library(plotrix)
library(mapplots)
library(viridis)
library(reshape2)

# Key metadata files are available in the 'https://github.com/parklab/focal-amplification/tree/main/Data'
# Loading three key metadata files

masterdf <- read.csv("List.patients_summary.final.txt", header=T, as.is=T, sep="\t")
purdf <- read.csv("List.purple.final.txt", header=T, as.is=T, sep="\t")
hrdf <- read.csv("List.hrd_status.final.txt", header=T, as.is=T, sep="\t")

# purdf and hrdf were created based on the masterdf dataframe, by collecting information from individual variant calls of each sample, using the codes below.
# The products of the codes below are "List.purple.final.txt" and "List.hrd_status.final.txt".

## Collection of summary files
purdf <- masterdf[,c("icgc_donor_id", "age", "submitted_sample_id", "study_id", "class", "cohort", "proc_nov26", "folder_nov26", "file_nov26")]
tempdf <- read.csv(paste0("../27_EBI/purple_purity_summary/pcawg/", purdf$file_nov26[1], ".purple.purity.tsv"), header=T, as.is=T, sep="\t")
purdf[1,colnames(tempdf)] <- tempdf[1,]
for (i in 2:nrow(purdf)){
  if (purdf$file_nov26[i] != "no"){
    if (purdf$class[i] == "PCAWG"){
      tempdf <- read.csv(paste0("../27_EBI/purple_purity_summary/pcawg/", purdf$file_nov26[i], ".purple.purity.tsv"), header=T, as.is=T, sep="\t")
    } else if (purdf$class[i] == "French72"){
      tempdf <- read.csv(paste0("../27_EBI/purple_purity_summary/ferrari/", purdf$file_nov26[i], ".purple.purity.tsv"), header=T, as.is=T, sep="\t")
    } else if (purdf$class[i] == "BritCol87"){
      tempdf <- read.csv(paste0("../27_EBI/purple_purity_summary/zhao/", purdf$file_nov26[i], ".purple.purity.tsv"), header=T, as.is=T, sep="\t")
    } else if (purdf$class[i] == "Sanger560"){
      tempdf <- read.csv(paste0("../27_EBI/purple_purity_summary/basis/", purdf$file_nov26[i], ".purple.purity.tsv"), header=T, as.is=T, sep="\t")
    } else if (purdf$class[i] == "Yale20"){
      tempdf <- read.csv(paste0("../27_EBI/purple_purity_summary/yale/", purdf$file_nov26[i], ".purple.purity.tsv"), header=T, as.is=T, sep="\t")
    }
    purdf[i,colnames(tempdf)] <- tempdf[1,]
  }
}
purdf$tmbGenome <- purdf$tmbPerMb*2859

## Collection of qc reports
purdf$qc <- NA
for (i in 1:nrow(purdf)){
  if (purdf$file_nov26[i] != "no"){
    purdf$qc[i] <- strsplit(readLines(paste0("../27_EBI/purple_qc/merged/", purdf$file_nov26[i], ".purple.qc"), n=1), "\t", fixed=T)[[1]][2]
  }
}

## Updating purdf information from the final datasets
purdf$cohort <- "no"
purdf$cohort[purdf$class == "PCAWG"] <- "pcawg"
purdf$cohort[purdf$class == "French72"] <- "ferrari"
purdf$cohort[purdf$class == "BritCol87"] <- "zhao"
purdf$cohort[purdf$class == "Yale20"] <- "yale"
purdf$cohort[purdf$class == "Sanger560"] <- "basis"

## Below code runs with a premise of the same order and length or the rows in both masterdf and purdf
for (i in 1:nrow(purdf)){
  purdf$folder_nov26[i] <- masterdf$folder_nov26[masterdf$study_id==purdf$study_id[i]]
  purdf$file_nov26[i] <- masterdf$file_nov26[masterdf$study_id==purdf$study_id[i]]
  purdf$availability[i] <- masterdf$availability[masterdf$study_id==purdf$study_id[i]]
}

for (i in which(masterdf$availability == "yes")){
  #print(i)
  tempdf <- read.csv(paste0("../27_EBI/final_data_022222/purity/", purdf$study_id[i], ".purple.purity.tsv"), header=T, as.is=T, sep="\t")
  purdf[i, colnames(tempdf)] <- tempdf[1,]
  purdf$qc[i] <- strsplit(readLines(paste0("../27_EBI/final_data_022222/qc/", purdf$study_id[i], ".purple.qc"), n=1), "\t", fixed=T)[[1]][2]
  hcall <- read.csv(paste0("../27_EBI/final_data_022222/rearrangement_bedpe/", purdf$study_id[i], ".purple.sv.vcf.gz.bedpe"), header=T, as.is=T, sep="\t")
  if (nrow(hcall) != 0){
    hcall <- read.csv(paste0("../27_EBI/final_data_022222/rearrangement_bedpe/", purdf$study_id[i], ".purple.sv.vcf.gz.bedpe"), header=T, as.is=T, sep="\t")
    hcall$chr_num1 <- hcall$chrom1
    hcall$chr_num1[hcall$chr_num1 == "X"] <- "23"
    hcall$chr_num1[hcall$chr_num1 == "Y"] <- "24"
    hcall$chr_num1[hcall$chr_num1 == "MT"] <- "25"
    hcall$chr_num1 <- as.numeric(as.character(hcall$chr_num1))
    hcall$chr_num2 <- hcall$chrom2
    hcall$chr_num2[hcall$chr_num2 == "X"] <- "23"
    hcall$chr_num2[hcall$chr_num2 == "Y"] <- "24"
    hcall$chr_num2[hcall$chr_num2 == "MT"] <- "25"
    hcall$chr_num2 <- as.numeric(as.character(hcall$chr_num2))
    hcall$chord <- "normal"
    for (j in 1:nrow(hcall)){
      if (hcall$svclass[j] == "TRA"){
        if (hcall$chr_num2[j] < hcall$chr_num1[j]){
          hcall$chord[j] <- "abnormal"
        }
      } else {
        if (hcall$start2[j] < hcall$start1[j]){
          hcall$chord[j] <- "abnormal"
        }
      }
    }
    hcall$c1 <- 0
    hcall$p1 <- 0
    hcall$c2 <- 0
    hcall$p2 <- 0
    for (j in 1:nrow(hcall)){
      if (hcall$chord[j] == "normal"){
        hcall$c1[j] <- hcall$chr_num1[j]
        hcall$p1[j] <- round(mean(c(hcall$start1[j], hcall$end1[j])))
        hcall$c2[j] <- hcall$chr_num2[j]
        hcall$p2[j] <- round(mean(c(hcall$start2[j], hcall$end2[j])))
      } else {
        hcall$c1[j] <- hcall$chr_num2[j]
        hcall$p1[j] <- round(mean(c(hcall$start2[j], hcall$end2[j])))
        hcall$c2[j] <- hcall$chr_num1[j]
        hcall$p2[j] <- round(mean(c(hcall$start1[j], hcall$end1[j])))
      }
    }
    purdf$ct_sv_hmf[i] <- nrow(hcall)
    purdf$ct_DEL[i] <- nrow(hcall[hcall$svclass == "DEL",])
    purdf$ct_DUP[i] <- nrow(hcall[hcall$svclass == "DUP",])
    purdf$ct_h2hINV[i] <- nrow(hcall[hcall$svclass == "h2hINV",])
    purdf$ct_t2tINV[i] <- nrow(hcall[hcall$svclass == "t2tINV",])
    purdf$ct_TRA[i] <- nrow(hcall[hcall$svclass == "TRA",])
    purdf$ori_abn[i] <- nrow(hcall[hcall$chord == "abnormal",])
    purdf$inv_rate[i] <- (purdf$ct_h2hINV[i]+purdf$ct_t2tINV[i])/(purdf$ct_sv_hmf[i]-purdf$ct_TRA[i])
  } else {
    print(paste0("no SV file for ", purdf$study_id[i]))
    next
  }
}
purdf$tmbGenome <- purdf$tmbPerMb*2859

## Number of linx-SVs in each sample
purdf$ct_linxsv <- NA
for (i in which(masterdf$availability=="yes")){
  tempdf <- read.csv(paste0("../27_EBI/final_data_022222/linx_sv_annotations/", masterdf$study_id[i], ".linx.svs.tsv"), header=T, as.is=T, sep="\t")
  purdf$ct_linxsv[i] <- nrow(tempdf)
}
hist((purdf$ct_sv_hmf/purdf$ct_linxsv), las=1) # linx analysis needs to be restricted to those with QC == PASS

## Number of SNVs and Indels in each sample
purdf$ct_snv <- NA
purdf$ct_indel <- NA
for (i in which(purdf$availability=="yes")){
  print(i)
  snvdf <- read.csv(paste0("../27_EBI/final_data_022222/snv_dataframe/", purdf$study_id[i], ".sage30.snv.txt"), header=T, as.is=T, sep="\t")
  indeldf <- read.csv(paste0("../27_EBI/final_data_022222/indel_dataframe/", purdf$study_id[i], ".sage30.indels.txt"), header=T, as.is=T, sep="\t")
  purdf$ct_snv[i] <- nrow(snvdf)
  purdf$ct_indel[i] <- nrow(indeldf)
  rm(snvdf)
  rm(indeldf)
}

purdf$haploidProportion <- NA
for (i in which(purdf$availability=="yes")){
  print(i)
  cnvdf <- read.csv(paste0("../27_EBI/final_data_022222/acgr/", purdf$study_id[i], ".purple.cnv.somatic.acgr.tsv"), header=T, as.is=T, sep="\t")
  purdf$haploidProportion[i] <- sum(cnvdf$seglen[cnvdf$loh=="loh"])/sum(cnvdf$seglen)
  rm(cnvdf)
}
purdf$is.wgd <- ifelse(2.85 - 1.6*purdf$haploidProportion <= purdf$ploidy, "WGD", "ND")

purdf$hr_status <- "not_tested"
for (i in which(purdf$availability == "yes")){
  purdf$hr_status[i] <- hrdf$hr_status[hrdf$study_id == purdf$study_id[i]]
}

write.table(purdf, "List.purple.final.txt", row.names=F, col.names=T, quote=F, sep="\t")

## Integrating CHORD result to the samples
hrdf <- purdf
hrdf$p_BRCA1 <- NA
hrdf$p_BRCA2 <- NA
hrdf$p_hrd <- NA
hrdf$hr_status <- "not_analyzed"
hrdf$hrd_type <- "not_analyzed"
hrdf <- hrdf[hrdf$availability == "yes",]

for (k in unique(hrdf$cohort)){
  tempdf <- read.csv(paste0("../27_EBI/chord/", k, ".chord_prediction.txt"), header=T, as.is=T, sep="\t")
  for (i in which(hrdf$cohort == k)){
    hrdf$p_BRCA1[i] <- tempdf$p_BRCA1[tempdf$sample == hrdf$study_id[i]]
    hrdf$p_BRCA2[i] <- tempdf$p_BRCA2[tempdf$sample == hrdf$study_id[i]]
    hrdf$p_hrd[i] <- tempdf$p_hrd[tempdf$sample == hrdf$study_id[i]]
    hrdf$hr_status[i] <- tempdf$hr_status[tempdf$sample == hrdf$study_id[i]]
    hrdf$hrd_type[i] <- tempdf$hrd_type[tempdf$sample == hrdf$study_id[i]]
  }
}

hrdf$er_plot <- "unknown"
hrdf$pr_plot <- "unknown"
hrdf$her2_plot <- "unknown"
hrdf$pam50_plot <- "unknown"
hrdf$er_plot <- factor(hrdf$er_plot, levels = c("positive", "negative", "unknown"))
hrdf$pr_plot <- factor(hrdf$pr_plot, levels = c("positive", "negative", "unknown"))
hrdf$her2_plot <- factor(hrdf$her2_plot, levels = c("positive", "equivocal", "negative", "unknown"))
hrdf$pam50_plot <- factor(hrdf$pam50_plot, levels = c("LumA", "LumB", "Her2", "Basal", "Normal", "unknown"))
hrdf$age_group <- NA
hrdf$histo_plot <- "no"

for (i in 1:nrow(hrdf)){
  hrdf$er_plot[i] <- demdf$er_plot[demdf$study_id == hrdf$study_id[i]]
  hrdf$pr_plot[i] <- demdf$pr_plot[demdf$study_id == hrdf$study_id[i]]
  hrdf$her2_plot[i] <- demdf$her2_plot[demdf$study_id == hrdf$study_id[i]]
  hrdf$pam50_plot[i] <- demdf$pam50_plot[demdf$study_id == hrdf$study_id[i]]
  hrdf$age_group[i] <- demdf$age_group[demdf$study_id == hrdf$study_id[i]]
  hrdf$histo_plot[i] <- demdf$histo_plot[demdf$study_id == hrdf$study_id[i]]
}

hrdf$germline <- "no"
hrdf$mut_germline <- 0
for (i in 1:nrow(hrdf)){
  tempdf <- read.csv(paste0("../27_EBI/final_data_022222/driver_germline/", hrdf$study_id[i], ".driver.catalog.germline.tsv"), header=T, as.is=T, sep="\t")
  if (nrow(tempdf) != 0){
    hrdf$germline[i] <- paste(tempdf$gene, collapse=";")
    if ("BRCA1" %in% tempdf$gene){
      hrdf$mut_germline[i] <- 1
    }
    if ("BRCA2" %in% tempdf$gene){
      hrdf$mut_germline[i] <- 2
    }
    if ("PALB2" %in% tempdf$gene){
      hrdf$mut_germline[i] <- 3
    }
  }
}

hrdf$mut_somatic <- 0
hrdf$mut_tp53 <- 0
for (i in 1:nrow(hrdf)){
  tempdf <- read.csv(paste0("../27_EBI/final_data_022222/driver_somatic/", hrdf$study_id[i], ".driver.catalog.somatic.tsv"), header=T, as.is=T, sep="\t")
  if (nrow(tempdf) != 0){
    if ("BRCA1" %in% tempdf$gene){
      if (tempdf$biallelic[tempdf$gene == "BRCA1"] == "true"){
        hrdf$mut_somatic[i] <- 1
      }
    }
    if ("BRCA2" %in% tempdf$gene){
      if (tempdf$biallelic[tempdf$gene == "BRCA2"] == "true"){
        hrdf$mut_somatic[i] <- 2
      }
    }
    if ("TP53" %in% tempdf$gene){
      hrdf$mut_tp53[i] <- 1
    }
  }
}

write.table(hrdf, "List.hrd_status.final.txt", row.names = F, col.names = T, quote = F, sep = "\t")



