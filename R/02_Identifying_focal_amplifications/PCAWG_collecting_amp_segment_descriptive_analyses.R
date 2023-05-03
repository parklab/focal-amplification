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
# The products of the codes below is "List.amplicons_PCAWG.final.txt".

sourcedir <- "" # Location of the source directory, where the output files from "PCAWG_definition_amp_segment.R" are located.
svdir <- "" # Location of the SV directory. SV files are available "SV_BEDPE_annotated.zip" in the Data folder.

infodf <- data.frame(matrix(ncol=21, nrow=0))
colnames(infodf) <- c("icgc_donor_id", "project_code", "chromosome", "start", "end", "mcj_s", "mcj_e", "class_s", "svclass_s", "svname_s", "svf_s", "svlen_s", "leftcn_s", "rightcn_s", "class_e", "svclass_e", "svname_e", "svf_e", "svlen_e", "leftcn_e", "rightcn_e")

## Starting from sumdf dataframe that was created in "PCAWG_definition_amp_segment.R"
for (w in which(!is.na(sumdf$SV.events))){
  acdf <- read.csv(paste0(sourcedir, sumdf$icgc_donor_id[w], ".pcawg.cnv.somatic.acgr.tsv"), header=T, as.is=T, sep="\t")
  acdf$segmentStartSupport[acdf$segmentStartSupport=="t2tINV" | acdf$segmentStartSupport=="h2hINV"] <- "INV"
  acdf$segmentEndSupport[acdf$segmentEndSupport=="t2tINV" | acdf$segmentEndSupport=="h2hINV"] <- "INV"
  df <- acdf[acdf$acgr %in% c("amplicon_start", "amplicon_end", "amplicon_itself"),]
  if (nrow(df) != 0){
    svdf <- read.csv(paste0(svdir, sumdf$tumor[w], ".pcawg_consensus_1.6.161022.somatic.sv.bedpe"), header=T, as.is=T, sep="\t")
    svdf$svlen[svdf$svclass != "TRA"] <- round(abs((svdf$start2[svdf$svclass != "TRA"] + svdf$end2[svdf$svclass != "TRA"])/2 - (svdf$start1[svdf$svclass != "TRA"] + svdf$end1[svdf$svclass != "TRA"])/2))
    # Print ID
    print(paste0(sumdf$icgc_donor_id[w], " ; ", w, " ; ", sumdf$tumor[w]))
    for (i in 1:nrow(df)){
      if (df$acgr[i] == "amplicon_itself"){
        infodf[nrow(infodf)+1,] <- NA
        infodf$icgc_donor_id[nrow(infodf)] <- sumdf$icgc_donor_id[w]
        infodf$project_code[nrow(infodf)] <- sumdf$project_code[w]
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
            infodf$svclass_s[nrow(infodf)] <- tempdf$svclass[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
            infodf$svname_s[nrow(infodf)] <- tempdf$sv_id[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
            infodf$svf_s[nrow(infodf)] <- tempdf$pe_support[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
            infodf$svlen_s[nrow(infodf)] <- tempdf$svlen[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
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
          if (infodf$class_e[nrow(infodf)] %in% c("BND", "DEL", "DUP", "INV", "MULTIPLE") & nrow(svdf[svdf$chrom1 == df$chromosome[i] & svdf$start1-1 <= df$end[i] & svdf$end1-1 >= df$end[i] | svdf$chrom2 == df$chromosome[i] & svdf$start2-1 <= df$end[i] & svdf$end2-1 >= df$end[i],]) != 0){
            tempdf <- svdf[svdf$chrom1 == df$chromosome[i] & svdf$start1-1 <= df$end[i] & svdf$end1-1 >= df$end[i] | svdf$chrom2 == df$chromosome[i] & svdf$start2-1 <= df$end[i] & svdf$end2-1 >= df$end[i],]
            infodf$svclass_e[nrow(infodf)] <- tempdf$svclass[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
            infodf$svname_e[nrow(infodf)] <- tempdf$sv_id[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
            infodf$svf_e[nrow(infodf)] <- tempdf$pe_support[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
            infodf$svlen_e[nrow(infodf)] <- tempdf$svlen[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
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
        infodf$icgc_donor_id[nrow(infodf)] <- sumdf$icgc_donor_id[w]
        infodf$project_code[nrow(infodf)] <- sumdf$project_code[w]
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
            infodf$svclass_s[nrow(infodf)] <- tempdf$svclass[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
            infodf$svname_s[nrow(infodf)] <- tempdf$sv_id[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
            infodf$svf_s[nrow(infodf)] <- tempdf$pe_support[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
            infodf$svlen_s[nrow(infodf)] <- tempdf$svlen[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
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
          if (infodf$class_e[nrow(infodf)] %in% c("BND", "DEL", "DUP", "INV", "MULTIPLE") & nrow(svdf[svdf$chrom1 == df$chromosome[i+1] & svdf$start1-1 <= df$end[i+1] & svdf$end1-1 >= df$end[i+1] | svdf$chrom2 == df$chromosome[i+1] & svdf$start2-1 <= df$end[i+1] & svdf$end2-1 >= df$end[i+1],]) != 0){
            tempdf <- svdf[svdf$chrom1 == df$chromosome[i+1] & svdf$start1-1 <= df$end[i+1] & svdf$end1-1 >= df$end[i+1] | svdf$chrom2 == df$chromosome[i+1] & svdf$start2-1 <= df$end[i+1] & svdf$end2-1 >= df$end[i+1],]
            infodf$svclass_e[nrow(infodf)] <- tempdf$svclass[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
            infodf$svname_e[nrow(infodf)] <- tempdf$sv_id[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
            infodf$svf_e[nrow(infodf)] <- tempdf$pe_support[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
            infodf$svlen_e[nrow(infodf)] <- tempdf$svlen[tempdf$pe_support == max(tempdf$pe_support, na.rm=T)][1]
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

## Assignment of mechanisms -- order of the commands matters
infodf$mech_s <- "not_analyzed"
infodf$mech_e <- "not_analyzed"
infodf$mech_s[infodf$mcj_s != "no"] <- "analyzed"
infodf$mech_e[infodf$mcj_e != "no"] <- "analyzed"
infodf$mech_s[infodf$mcj_s != "no" & infodf$mcj_e != "no" & infodf$svclass_s == "DUP" & infodf$svclass_e == "DUP" & infodf$svname_s == infodf$svname_e] <- "tandem_duplication"
infodf$mech_e[infodf$mcj_s != "no" & infodf$mcj_e != "no" & infodf$svclass_s == "DUP" & infodf$svclass_e == "DUP" & infodf$svname_s == infodf$svname_e] <- "tandem_duplication"
infodf$mech_s[infodf$mcj_s != "no" & infodf$mcj_e != "no" & infodf$svclass_s == "DUP" & infodf$svclass_e == "DUP" & infodf$svname_s == infodf$svname_e & (infodf$rightcn_s > 3*infodf$leftcn_s | infodf$leftcn_e > 3*infodf$rightcn_e)] <- "td_extrachromosomal" # too strict criteria >4x to >3x
infodf$mech_e[infodf$mcj_s != "no" & infodf$mcj_e != "no" & infodf$svclass_s == "DUP" & infodf$svclass_e == "DUP" & infodf$svname_s == infodf$svname_e & (infodf$rightcn_s > 3*infodf$leftcn_s | infodf$leftcn_e > 3*infodf$rightcn_e)] <- "td_extrachromosomal"
infodf$mech_s[infodf$mcj_s != "no" & infodf$svclass_s == "t2tINV" & infodf$svlen_s <=5000] <- "fold_back_inversion"
infodf$mech_e[infodf$mcj_e != "no" & infodf$svclass_e == "h2hINV" & infodf$svlen_e <=5000] <- "fold_back_inversion"
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
infodf$mech_c[is.na(infodf$mech_c)] <- "intrachromosomal_complex"
infodf$mech_c[infodf$mech_s == "no_supporting_sv" & infodf$mech_e == "no_supporting_sv"] <- "no_supporting_sv"
infodf$mech_c[infodf$mech_s == "not_analyzed" & infodf$mech_e == "not_analyzed"] <- "excluded"

write.table(infodf, "List.amplicons_PCAWG.final.txt", row.names=F, col.names=T, quote=F, sep="\t")
ampcawg <- infodf



# Codes for the downstream analyses in the pan-cancer cohort 
## Clustering by mechanism
ampcawg$histology_abbreviation <- "no"
for (i in unique(ampcawg$icgc_donor_id)){
  ampcawg$histology_abbreviation[ampcawg$icgc_donor_id == i] <- sumdf$histology_abbreviation[sumdf$icgc_donor_id == i]
}
ampcawg$histology_abbreviation[ampcawg$histology_abbreviation=="Breast-LobularCA"] <- "Breast-AdenoCA"
ampcawg$histology_abbreviation[ampcawg$histology_abbreviation=="Breast-DCIS"] <- "Breast-AdenoCA"

typdf <- as.data.frame(matrix(NA, ncol=8, nrow=24))
colnames(typdf) <- c("histology_abbreviation", "n_cgr", "n_dn", "n_td", "n_fbi", "n_tra", "n_intcomp", "n_nosv")
typdf$histology_abbreviation <- names(sort(table(ampcawg$histology_abbreviation), decreasing = T))[1:24]
for (i in 1:nrow(typdf)){
  typdf$n_cgr[i] <- nrow(ampcawg[ampcawg$histology_abbreviation == typdf$histology_abbreviation[i] & ampcawg$mech_c != "excluded",])
  typdf$n_dn[i] <- nrow(ampcawg[ampcawg$histology_abbreviation == typdf$histology_abbreviation[i] & ampcawg$mech_c != "excluded" & ampcawg$mech_c == "double_minute",])
  typdf$n_td[i] <- nrow(ampcawg[ampcawg$histology_abbreviation == typdf$histology_abbreviation[i] & ampcawg$mech_c != "excluded" & ampcawg$mech_c == "tandem_duplication",])
  typdf$n_fbi[i] <- nrow(ampcawg[ampcawg$histology_abbreviation == typdf$histology_abbreviation[i] & ampcawg$mech_c != "excluded" & ampcawg$mech_c == "fold_back_inversion",])
  typdf$n_tra[i] <- nrow(ampcawg[ampcawg$histology_abbreviation == typdf$histology_abbreviation[i] & ampcawg$mech_c != "excluded" & ampcawg$mech_c == "translocation",])
  typdf$n_intcomp[i] <- nrow(ampcawg[ampcawg$histology_abbreviation == typdf$histology_abbreviation[i] & ampcawg$mech_c != "excluded" & ampcawg$mech_c == "intrachromosomal_complex",])
  typdf$n_nosv[i] <- nrow(ampcawg[ampcawg$histology_abbreviation == typdf$histology_abbreviation[i] & ampcawg$mech_c != "excluded" & ampcawg$mech_c == "no_supporting_sv",])
}
typdf$f_dn <- typdf$n_dn/(typdf$n_cgr-typdf$n_nosv)
typdf$f_td <- typdf$n_td/(typdf$n_cgr-typdf$n_nosv)
typdf$f_fbi <- typdf$n_fbi/(typdf$n_cgr-typdf$n_nosv)
typdf$f_tra <- typdf$n_tra/(typdf$n_cgr-typdf$n_nosv)
typdf$f_intcomp <- typdf$n_intcomp/(typdf$n_cgr-typdf$n_nosv)

## Dendrogram
selected_num <- c(1:22) 
data <- typdf[selected_num,c(9,11,12)] # only using td, fbi, and tra
data <- typdf[selected_num,c(9,10,11,12,13)]
row.names(data) <- typdf$histology_abbreviation[selected_num]
d <- dist(data, method="euclidean")
hc <- hclust(d, method = "complete")

pdf("FIP_revision/dendrogram.hclust.pcawg.052622.pdf", height=5, width=5)
plot(hc, hang = -1, xlab = "Histologic types")
rect.hclust(hc, k = 4)
dev.off()

# Sex annotation
ampcawg$gender <- "unknown"
tumorlist <- unique(ampcawg$icgc_donor_id)
for (i in tumorlist){
  ampcawg$gender[ampcawg$icgc_donor_id == i] <- sumdf$gender[sumdf$icgc_donor_id == i]
}

## Pie graph for FBI-TRA classes (Fig. 5b)
library(reshape2)
library(mapplots)
classcolors <- c("#ADFF2F", "#3D9140", "#DC143C", "#551A8B", "#C1CDCD", "#5A5A5A")
histolist <- c("ColoRect-AdenoCA", "Cervix-SCC", "Head-SCC", "Lung-SCC", "Panc-AdenoCA", "Liver-HCC", "Eso-AdenoCA", "Stomach-AdenoCA", "Bladder-TCC", "Ovary-AdenoCA", "Uterus-AdenoCA", "Prost-AdenoCA", "Biliary-AdenoCA", "Lymph-BNHL", "CNS-GBM", "CNS-Medullo", "Kidney-RCC", "SoftTissue-Liposarc", "Bone-Osteosarc", "Skin-Melanoma", "Lung-AdenoCA", "Breast-AdenoCA")
numberlist <- NA
for (i in 1:length(histolist)){
  numberlist[i] <- length(unique(ampcawg$icgc_donor_id[ampcawg$histology_abbreviation == histolist[i]]))
}
cgrdf <- as.data.frame(matrix(NA, ncol=12, nrow=24*nrow(data))) # 21 --> 22
colnames(cgrdf) <- c("histology_abbreviation", "xvalue", "yvalue", "chrom", "n_cgr", "n_dn", "n_td", "n_fbi", "n_tra", "n_intcomp", "n_nosv")
for (i in 1:length(histolist)){
  cgrdf$histology_abbreviation[(i*24-23):(i*24)] <- histolist[i]
  cgrdf$xvalue[(i*24-23):(i*24)] <- c(1:24)
  cgrdf$yvalue[(i*24-23):(i*24)] <- nrow(data)+2-i
  cgrdf$chrom[(i*24-23):(i*24)] <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
  cgrdf$n_cgr[(i*24-23):(i*24)] <- rep(0,24)
  for (j in (i*24-23):(i*24)){
    cgrdf$n_cgr[j] <- nrow(ampcawg[ampcawg$chromosome == cgrdf$chrom[j] & ampcawg$histology_abbreviation == histolist[i],])
    cgrdf$n_dn[j] <- nrow(ampcawg[ampcawg$chromosome == cgrdf$chrom[j] & ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "double_minute",])
    cgrdf$n_td[j] <- nrow(ampcawg[ampcawg$chromosome == cgrdf$chrom[j] & ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "tandem_duplication",])
    cgrdf$n_fbi[j] <- nrow(ampcawg[ampcawg$chromosome == cgrdf$chrom[j] & ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "fold_back_inversion",])
    cgrdf$n_tra[j] <- nrow(ampcawg[ampcawg$chromosome == cgrdf$chrom[j] & ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "translocation",])
    cgrdf$n_intcomp[j] <- nrow(ampcawg[ampcawg$chromosome == cgrdf$chrom[j] & ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "intrachromosomal_complex",])
    cgrdf$n_nosv[j] <- nrow(ampcawg[ampcawg$chromosome == cgrdf$chrom[j] & ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "no_supporting_sv",])
  }
}
cgr_short <- cgrdf[,c("xvalue", "yvalue", "n_dn", "n_td", "n_fbi", "n_tra", "n_intcomp", "n_nosv")]
cgr_melt <- melt(cgr_short, id.vars = c("xvalue", "yvalue"))

xyz <- make.xyz(cgr_melt$xvalue, cgr_melt$yvalue, cgr_melt$value, cgr_melt$variable)

pdf("FIP_revision/ampseg_boundaries_bychrom_tumortype_class.052622.pdf", height=7, width=8)
par(mar=c(4.1, 12.1, 4.1, 2.1))
plot(NA, xlim = c(0,24), ylim = c(-1,23), frame = F, xlab = "Chromosome", yaxt='n', xaxt='n', ylab='', main = paste0("Number of amplified segments (total n = ", nrow(ampcawg[ampcawg$mech_c != "excluded",]), ")"))
axis(2, at=c(1:(nrow(data)+2)), las=1, labels=c("", paste0(rev(histolist), " (", rev(numberlist),")"), ""))
axis(1, at=c(0:24), las=1, labels=c("", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))
for (i in 1:23){
  abline(h = i, col=rgb(232/255, 232/255, 232/255))
}
for (i in 1:24){
  abline(v = i, col=rgb(232/255, 232/255, 232/255))
}
draw.pie(xyz$x, xyz$y, xyz$z, radius = 1.4, col=classcolors, border = NA)
legend.pie(16, -0.4, labels=c("Double Minutes", "Tandem duplication", "Fold-back inversion", "Interchromosomal", "Intrachromosomal complex", "No SV support"), radius = 0.6, bty="n", col=classcolors, border=NA, label.dist=1.3)
legend.z <- round(max(rowSums(xyz$z,na.rm=TRUE)))
legend.bubble(3, -0.4, z=legend.z, round=0, maxradius=1.4, bty="n", txt.cex=0.6)
dev.off()

## Male and female separately
histolist <- c("ColoRect-AdenoCA", "Cervix-SCC", "Head-SCC", "Lung-SCC", "Panc-AdenoCA", "Liver-HCC", "Eso-AdenoCA", "Stomach-AdenoCA", "Bladder-TCC", "Ovary-AdenoCA", "Uterus-AdenoCA", "Prost-AdenoCA", "Biliary-AdenoCA", "Lymph-BNHL", "CNS-GBM", "CNS-Medullo", "Kidney-RCC", "SoftTissue-Liposarc", "Bone-Osteosarc", "Skin-Melanoma", "Lung-AdenoCA", "Breast-AdenoCA")
for (w in c("male", "female")){
  numberlist <- NA
  for (i in 1:length(histolist)){
    numberlist[i] <- length(unique(ampcawg$icgc_donor_id[ampcawg$histology_abbreviation == histolist[i] & ampcawg$gender == w]))
  }
  cgrdf <- as.data.frame(matrix(NA, ncol=12, nrow=24*nrow(data))) # 21 --> 22
  colnames(cgrdf) <- c("histology_abbreviation", "xvalue", "yvalue", "chrom", "n_cgr", "n_dn", "n_td", "n_fbi", "n_tra", "n_intcomp", "n_nosv")
  for (i in 1:length(histolist)){
    cgrdf$histology_abbreviation[(i*24-23):(i*24)] <- histolist[i]
    cgrdf$xvalue[(i*24-23):(i*24)] <- c(1:24)
    cgrdf$yvalue[(i*24-23):(i*24)] <- nrow(data)+2-i
    cgrdf$chrom[(i*24-23):(i*24)] <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
    cgrdf$n_cgr[(i*24-23):(i*24)] <- rep(0,24)
    for (j in (i*24-23):(i*24)){
      cgrdf$n_cgr[j] <- nrow(ampcawg[ampcawg$chromosome == cgrdf$chrom[j] & ampcawg$histology_abbreviation == histolist[i] & ampcawg$gender == w,])
      cgrdf$n_dn[j] <- nrow(ampcawg[ampcawg$chromosome == cgrdf$chrom[j] & ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "double_minute" & ampcawg$gender == w,])
      cgrdf$n_td[j] <- nrow(ampcawg[ampcawg$chromosome == cgrdf$chrom[j] & ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "tandem_duplication" & ampcawg$gender == w,])
      cgrdf$n_fbi[j] <- nrow(ampcawg[ampcawg$chromosome == cgrdf$chrom[j] & ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "fold_back_inversion" & ampcawg$gender == w,])
      cgrdf$n_tra[j] <- nrow(ampcawg[ampcawg$chromosome == cgrdf$chrom[j] & ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "translocation" & ampcawg$gender == w,])
      cgrdf$n_intcomp[j] <- nrow(ampcawg[ampcawg$chromosome == cgrdf$chrom[j] & ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "intrachromosomal_complex" & ampcawg$gender == w,])
      cgrdf$n_nosv[j] <- nrow(ampcawg[ampcawg$chromosome == cgrdf$chrom[j] & ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "no_supporting_sv" & ampcawg$gender == w,])
    }
  }
  cgr_short <- cgrdf[,c("xvalue", "yvalue", "n_dn", "n_td", "n_fbi", "n_tra", "n_intcomp", "n_nosv")]
  cgr_melt <- melt(cgr_short, id.vars = c("xvalue", "yvalue"))
  cgr_melt$value[cgr_melt$xvalue==24 & cgr_melt$yvalue==2 & cgr_melt$variable == "n_nosv"] <- 600
  #cgr_melt[nrow(cgr_melt)+1,] <- c(1,1,"n_nosv", 600)
  xyz <- make.xyz(cgr_melt$xvalue, cgr_melt$yvalue, cgr_melt$value, cgr_melt$variable)
  pdf(paste0("FIP_revision/ampseg_boundaries_bychrom_tumortype_class.update.", w, ".052722.pdf"), height=7, width=8)
  par(mar=c(4.1, 12.1, 4.1, 2.1))
  plot(NA, xlim = c(0,24), ylim = c(-1,23), frame = F, xlab = "Chromosome", yaxt='n', xaxt='n', ylab='', main = paste0("Number of amplified segments (total n = ", nrow(ampcawg[ampcawg$mech_c != "excluded" & ampcawg$gender == w,]), ")"))
  axis(2, at=c(1:(nrow(data)+2)), las=1, labels=c("", paste0(rev(histolist), " (", rev(numberlist),")"), ""))
  axis(1, at=c(0:24), las=1, labels=c("", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))
  for (i in 1:23){
    abline(h = i, col=rgb(232/255, 232/255, 232/255))
  }
  for (i in 1:24){
    abline(v = i, col=rgb(232/255, 232/255, 232/255))
  }
  draw.pie(xyz$x, xyz$y, xyz$z, radius = 1.4, col=classcolors, border = NA)
  legend.pie(16, -0.4, labels=c("Double Minutes", "Tandem duplication", "Fold-back inversion", "Interchromosomal", "Intrachromosomal complex", "No SV support"), radius = 0.6, bty="n", col=classcolors, border=NA, label.dist=1.3)
  legend.z <- round(max(rowSums(xyz$z,na.rm=TRUE)))
  legend.bubble(3, -0.4, z=legend.z, round=0, maxradius=1.4, bty="n", txt.cex=0.6)
  dev.off()
}

gdrdf <- as.data.frame(matrix(NA, ncol=10, nrow=0))
colnames(gdrdf) <- c("histology_abbreviation", "gender", "n_patient", "n_cgr", "n_dn", "n_td", "n_fbi", "n_tra", "n_intcomp", "n_nosv")
for (i in 1:length(histolist)){
  gdrdf[(nrow(gdrdf)+1),] <- NA
  gdrdf$histology_abbreviation[nrow(gdrdf)] <- histolist[i]
  gdrdf$gender[nrow(gdrdf)] <- "male"
  gdrdf$n_patient[nrow(gdrdf)] <- length(unique(ampcawg$icgc_donor_id[ampcawg$histology_abbreviation == histolist[i] & ampcawg$gender == "male"]))
  gdrdf$n_cgr[nrow(gdrdf)] <- nrow(ampcawg[ampcawg$histology_abbreviation == histolist[i] & ampcawg$gender == "male",])
  gdrdf$n_dn[nrow(gdrdf)] <- nrow(ampcawg[ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "double_minute" & ampcawg$gender == "male",])
  gdrdf$n_td[nrow(gdrdf)] <- nrow(ampcawg[ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "tandem_duplication" & ampcawg$gender == "male",])
  gdrdf$n_fbi[nrow(gdrdf)] <- nrow(ampcawg[ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "fold_back_inversion" & ampcawg$gender == "male",])
  gdrdf$n_tra[nrow(gdrdf)] <- nrow(ampcawg[ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "translocation" & ampcawg$gender == "male",])
  gdrdf$n_intcomp[nrow(gdrdf)] <- nrow(ampcawg[ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "intrachromosomal_complex" & ampcawg$gender == "male",])
  gdrdf$n_nosv[nrow(gdrdf)] <- nrow(ampcawg[ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "no_supporting_sv" & ampcawg$gender == "male",])
  gdrdf[(nrow(gdrdf)+1),] <- NA
  gdrdf$histology_abbreviation[nrow(gdrdf)] <- histolist[i]
  gdrdf$gender[nrow(gdrdf)] <- "female"
  gdrdf$n_patient[nrow(gdrdf)] <- length(unique(ampcawg$icgc_donor_id[ampcawg$histology_abbreviation == histolist[i] & ampcawg$gender == "female"]))
  gdrdf$n_cgr[nrow(gdrdf)] <- nrow(ampcawg[ampcawg$histology_abbreviation == histolist[i] & ampcawg$gender == "female",])
  gdrdf$n_dn[nrow(gdrdf)] <- nrow(ampcawg[ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "double_minute" & ampcawg$gender == "female",])
  gdrdf$n_td[nrow(gdrdf)] <- nrow(ampcawg[ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "tandem_duplication" & ampcawg$gender == "female",])
  gdrdf$n_fbi[nrow(gdrdf)] <- nrow(ampcawg[ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "fold_back_inversion" & ampcawg$gender == "female",])
  gdrdf$n_tra[nrow(gdrdf)] <- nrow(ampcawg[ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "translocation" & ampcawg$gender == "female",])
  gdrdf$n_intcomp[nrow(gdrdf)] <- nrow(ampcawg[ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "intrachromosomal_complex" & ampcawg$gender == "female",])
  gdrdf$n_nosv[nrow(gdrdf)] <- nrow(ampcawg[ampcawg$histology_abbreviation == histolist[i] & ampcawg$mech_c == "no_supporting_sv" & ampcawg$gender == "female",])
}
gdrdf$f_dn <- gdrdf$n_dn/(gdrdf$n_dn + gdrdf$n_td + gdrdf$n_fbi + gdrdf$n_tra + gdrdf$n_intcomp + gdrdf$n_nosv)
gdrdf$f_td <- gdrdf$n_td/(gdrdf$n_dn + gdrdf$n_td + gdrdf$n_fbi + gdrdf$n_tra + gdrdf$n_intcomp + gdrdf$n_nosv)
gdrdf$f_fbi <- gdrdf$n_fbi/(gdrdf$n_dn + gdrdf$n_td + gdrdf$n_fbi + gdrdf$n_tra + gdrdf$n_intcomp + gdrdf$n_nosv)
gdrdf$f_tra <- gdrdf$n_tra/(gdrdf$n_dn + gdrdf$n_td + gdrdf$n_fbi + gdrdf$n_tra + gdrdf$n_intcomp + gdrdf$n_nosv)
gdrdf$f_intcomp <- gdrdf$n_intcomp/(gdrdf$n_dn + gdrdf$n_td + gdrdf$n_fbi + gdrdf$n_tra + gdrdf$n_intcomp + gdrdf$n_nosv)
gdrdf$f_nosv <- gdrdf$n_nosv/(gdrdf$n_dn + gdrdf$n_td + gdrdf$n_fbi + gdrdf$n_tra + gdrdf$n_intcomp + gdrdf$n_nosv)

tempdf <- as.data.frame(matrix(NA, ncol=9, nrow=0))
colnames(tempdf) <- c("histology_abbreviation", "n_male", "n_female", "cgr_male", "cgr_female", "fbi_male", "fbi_female", "tra_male", "tra_female")
for (i in 1:length(histolist)){
  tempdf[i,] <- NA
  tempdf$histology_abbreviation[i] <- histolist[i]
  tempdf$n_male[i] <- gdrdf$n_patient[gdrdf$histology_abbreviation == histolist[i] & gdrdf$gender == "male"]
  tempdf$n_female[i] <- gdrdf$n_patient[gdrdf$histology_abbreviation == histolist[i] & gdrdf$gender == "female"]
  tempdf$cgr_male[i] <- gdrdf$n_cgr[gdrdf$histology_abbreviation == histolist[i] & gdrdf$gender == "male"]
  tempdf$cgr_female[i] <- gdrdf$n_cgr[gdrdf$histology_abbreviation == histolist[i] & gdrdf$gender == "female"]
  tempdf$fbi_male[i] <- gdrdf$n_fbi[gdrdf$histology_abbreviation == histolist[i] & gdrdf$gender == "male"]
  tempdf$fbi_female[i] <- gdrdf$n_fbi[gdrdf$histology_abbreviation == histolist[i] & gdrdf$gender == "female"]
  tempdf$tra_male[i] <- gdrdf$n_tra[gdrdf$histology_abbreviation == histolist[i] & gdrdf$gender == "male"]
  tempdf$tra_female[i] <- gdrdf$n_tra[gdrdf$histology_abbreviation == histolist[i] & gdrdf$gender == "female"]
}
tempdf$p.cgr <- NA
tempdf$p.fbi <- NA
tempdf$o.fbi <- NA
tempdf$p.tra <- NA
tempdf$o.tra <- NA
tempdf$lci.tra <- NA
tempdf$uci.tra <- NA

for (i in which(tempdf$n_male != 0 & tempdf$n_female != 0)){
  tempdf$p.cgr[i] <- fisher.test(matrix(c(tempdf$cgr_female[i], tempdf$cgr_male[i], tempdf$n_female[i], tempdf$n_male[i]), nrow = 2, ncol = 2))$p.value
  tempdf$p.fbi[i] <- fisher.test(matrix(c(tempdf$fbi_female[i], tempdf$fbi_male[i], tempdf$cgr_female[i], tempdf$cgr_male[i]), nrow = 2, ncol = 2))$p.value
  tempdf$o.fbi[i] <- fisher.test(matrix(c(tempdf$fbi_female[i], tempdf$fbi_male[i], tempdf$cgr_female[i], tempdf$cgr_male[i]), nrow = 2, ncol = 2))$estimate
  tempdf$p.tra[i] <- fisher.test(matrix(c(tempdf$tra_female[i], tempdf$tra_male[i], tempdf$cgr_female[i], tempdf$cgr_male[i]), nrow = 2, ncol = 2))$p.value
  tempdf$o.tra[i] <- fisher.test(matrix(c(tempdf$tra_female[i], tempdf$tra_male[i], tempdf$cgr_female[i], tempdf$cgr_male[i]), nrow = 2, ncol = 2))$estimate
  tempdf$lci.tra[i] <- fisher.test(matrix(c(tempdf$tra_female[i], tempdf$tra_male[i], tempdf$cgr_female[i], tempdf$cgr_male[i]), nrow = 2, ncol = 2))$conf.int[1]
  tempdf$uci.tra[i] <- fisher.test(matrix(c(tempdf$tra_female[i], tempdf$tra_male[i], tempdf$cgr_female[i], tempdf$cgr_male[i]), nrow = 2, ncol = 2))$conf.int[2]
}

fisher.test(matrix(c(sum(tempdf$cgr_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$n_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$n_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))])), nrow = 2, ncol = 2))$p.value
fisher.test(matrix(c(sum(tempdf$fbi_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$fbi_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))])), nrow = 2, ncol = 2))$p.value
fisher.test(matrix(c(sum(tempdf$tra_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$tra_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))])), nrow = 2, ncol = 2))$p.value
fisher.test(matrix(c(sum(tempdf$tra_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$tra_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))])), nrow = 2, ncol = 2))$estimate
fisher.test(matrix(c(sum(tempdf$tra_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$tra_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))])), nrow = 2, ncol = 2))$conf.int[1]
fisher.test(matrix(c(sum(tempdf$tra_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$tra_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))])), nrow = 2, ncol = 2))$conf.int[2]

tempdf <- tempdf[tempdf$tra_male > 2 & tempdf$tra_female > 2,]
tempdf <- tempdf[order(tempdf$o.tra, decreasing = T),]

## Forest plot (Extended Data Fig. 10d)
pdf("FIP_revision/forest.sex.boundary.translocations.052722.pdf", height=5, width=5)
par(mar=c(5.1, 9.1, 4.1, 2.1))
plot(NULL, log='x', xlim=c(0.05, 20), ylim=c(0,nrow(tempdf)+2), frame=F, las=1, yaxt='n', ylab="", xlab="Odds ratio")
segments(1,0,1,12, col="gray")
for (i in nrow(tempdf):1){
  segments(tempdf$lci.tra[i], 12-i, tempdf$uci.tra[i], 12-i, col=colordf$color[colordf$histology_abbreviation==tempdf$histology_abbreviation[i]], lend="butt")
  points(tempdf$o.tra[i], 12-i, pch=20, col=colordf$color[colordf$histology_abbreviation==tempdf$histology_abbreviation[i]])
}
segments(fisher.test(matrix(c(sum(tempdf$tra_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$tra_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))])), nrow = 2, ncol = 2))$conf.int[1], 1, fisher.test(matrix(c(sum(tempdf$tra_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$tra_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))])), nrow = 2, ncol = 2))$conf.int[2], 1, col="dark red", lwd=2, lend="butt")
points(fisher.test(matrix(c(sum(tempdf$tra_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$tra_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_female[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))]), sum(tempdf$cgr_male[!(tempdf$histology_abbreviation %in% c("Breast-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Cervix-SCC", "Prost-AdenoCA"))])), nrow = 2, ncol = 2))$estimate, 1, pch=20, col="dark red")
axis(2, at=seq(1,11,1), labels=rev(c(tempdf$histology_abbreviation, "All")), las=1)
dev.off()

## Identification of TB amplification cases in the PCAWG cohort
df <- sumdf[,c("icgc_donor_id", "project_code", "gender", "age", "SV.events", "ploidy", "histology_abbreviation")]
df$histology_abbreviation[df$histology_abbreviation=="Breast-LobularCA"] <- "Breast-AdenoCA"
df$histology_abbreviation[df$histology_abbreviation=="Breast-DCIS"] <- "Breast-AdenoCA"
df <- df[!is.na(df$SV.events) & df$SV.events != 0,] # n = 2439

cutval <- 3
df$num.clustra <- 0
df$list.clustra <- "no"
for (i in 1:nrow(df)){
  print(paste0("Processing ", i))
  svdf <- read.csv(paste0("~/Documents/ICGC_PCAWG_calls/SV_calls_1.6/BEDPE/", sumdf$tumor[sumdf$icgc_donor_id == df$icgc_donor_id[i]], ".pcawg_consensus_1.6.161022.somatic.sv.bedpe"), header=T, as.is=T, sep="\t")
  svdf <- svdf[svdf$svclass == "TRA",]
  if (nrow(svdf) != 0){
    svdf$arm1 <- "not_available"
    svdf$arm2 <- "not_available"
    for (j in 1:nrow(svdf)){
      if ((svdf$start1[j]+svdf$end1[j])/2 <= hg19_coord$centro_mid[hg19_coord$chr == svdf$chrom1[j]]){
        svdf$arm1[j] <- paste0(svdf$chrom1[j], "p")
      } else {
        svdf$arm1[j] <- paste0(svdf$chrom1[j], "q")
      }
      if ((svdf$start2[j]+svdf$end2[j])/2 <= hg19_coord$centro_mid[hg19_coord$chr == svdf$chrom2[j]]){
        svdf$arm2[j] <- paste0(svdf$chrom2[j], "p")
      } else {
        svdf$arm2[j] <- paste0(svdf$chrom2[j], "q")
      }
    }
    svdf$boundary <- "no"
    infolist <- unique(c(ampcawg$svname_s[ampcawg$icgc_donor_id == df$icgc_donor_id[i] & ampcawg$mcj_s != "no" & ampcawg$svclass_s == "TRA"], ampcawg$svname_e[ampcawg$icgc_donor_id == df$icgc_donor_id[i] & ampcawg$mcj_e != "no" & ampcawg$svclass_e == "TRA"]))
    svdf$boundary[svdf$sv_id %in% infolist] <- "yes"
    svdf$merged <- paste(svdf$arm1, svdf$arm2, sep = ";")
    for (j in unique(svdf$merged[svdf$boundary=="yes"])){
      if (nrow(svdf[svdf$merged == j,]) >= cutval){
        df$num.clustra[i] <- df$num.clustra[i] + 1
        if (df$list.clustra[i] == "no"){
          df$list.clustra[i] <- paste0(j, "[", nrow(svdf[svdf$merged == j,]), "]")
        } else {
          df$list.clustra[i] <- paste0(df$list.clustra[i], ",", j, "[", nrow(svdf[svdf$merged == j,]), "]")
        }
      }
    }
  }
}

df$color = "no"
for (i in 1:nrow(df)){
  df$color[i] <- colordf$color[colordf$histology_abbreviation == df$histology_abbreviation[i]]
}

df <- df[order(df$icgc_donor_id, decreasing = F),]
pctba <- df[order(df$histology_abbreviation, decreasing = F),]

## ESR1/ESR2 mRNA expression analysis using PCAWG RNA-seq data
df <- gexp.tumor[gexp.tumor$geneName == "ESR1",]
df <- rbind(df, gexp.tumor[gexp.tumor$geneName == "ESR2",])
df <- t(df[c(1,2),c(2:ncol(df))])
colnames(df) <- c("ESR1", "ESR2")
df <- as.data.frame(df)
df$icgc_donor_id <- rownames(df)
df$histology_abbreviation <- "not_available"
for (i in 1:nrow(df)){
  if (nrow(sumdf[sumdf$icgc_donor_id == df$icgc_donor_id[i],]) != 0){
    df$histology_abbreviation[i] <- sumdf$histology_abbreviation[sumdf$icgc_donor_id == df$icgc_donor_id[i]]
  }
}
df <- df[df$histology_abbreviation != "not_available",]

pctba$ESR1 <- NA
pctba$ESR2 <- NA
for (i in 1:nrow(pctba)){
  if (nrow(df[df$icgc_donor_id == pctba$icgc_donor_id[i],]) != 0){
    pctba$ESR1[i] <- df$ESR1[df$icgc_donor_id == pctba$icgc_donor_id[i]]
    pctba$ESR2[i] <- df$ESR2[df$icgc_donor_id == pctba$icgc_donor_id[i]]
  }
}

survdf <- read.csv("../../ICGC_PCAWG_calls/pcawg_donor_clinical_August2016_v9.tsv", header=T, as.is=T, sep="\t") # Used PCAWG clinical data available at dcc.icgc.org
pctba$donor_survival_time <- NA
pctba$donor_event <- NA
pctba$surv.event <- 0
for (i in 1:nrow(pctba)){
  if (!is.na(survdf$donor_survival_time[survdf$icgc_donor_id == pctba$icgc_donor_id[i]])){
    pctba$donor_survival_time[i] <- survdf$donor_survival_time[survdf$icgc_donor_id == pctba$icgc_donor_id[i]]
    pctba$donor_event[i] <- survdf$donor_vital_status[survdf$icgc_donor_id == pctba$icgc_donor_id[i]]
  }
}
pctba$surv.event[pctba$donor_event == "deceased"] <- 1
pctba$tba <- 1
pctba$tba[pctba$num.clustra == 0] <- 0

## Correlation with age and gene expression
library(survival)
tempdf <- as.data.frame(matrix(NA, ncol=14, nrow=0))
colnames(tempdf) <- c("histology_abbreviation", "no_tba", "yes_tba", "age_no_tba", "age_yes_tba", "age.p.val", "esr1_no_tba", "esr1_yes_tba", "esr1.p.val", "num_survinfo", "num_surv_yes_tba", "msurv_no_tba", "msurv_yes_tba", "msurv.log.rank.p.val")
tempdf[1:length(c("Cervix-SCC", "Head-SCC", "Lung-SCC", "Panc-AdenoCA", "Eso-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Bladder-TCC", "Prost-AdenoCA", "ColoRect-AdenoCA", "Liver-HCC", "Biliary-AdenoCA", "Lymph-BNHL", "SoftTissue-Liposarc", "Bone-Osteosarc", "Lung-AdenoCA", "Stomach-AdenoCA", "Skin-Melanoma", "Breast-AdenoCA", "CNS-GBM", "CNS-Medullo")),] <- NA
tempdf$histology_abbreviation <- c("Cervix-SCC", "Head-SCC", "Lung-SCC", "Panc-AdenoCA", "Eso-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Bladder-TCC", "Prost-AdenoCA", "ColoRect-AdenoCA", "Liver-HCC", "Biliary-AdenoCA", "Lymph-BNHL", "SoftTissue-Liposarc", "Bone-Osteosarc", "Lung-AdenoCA", "Stomach-AdenoCA", "Skin-Melanoma", "Breast-AdenoCA", "CNS-GBM", "CNS-Medullo")
for (i in 1:nrow(tempdf)){
  tempdf$no_tba[i] <- nrow(pctba[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & pctba$tba == 0,])
  tempdf$yes_tba[i] <- nrow(pctba[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & pctba$tba != 0,])
  if (nrow(pctba[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & pctba$tba != 0,]) > 1 & nrow(pctba[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & pctba$tba == 0,]) > 1){
    tempdf$age_no_tba[i] <- mean(pctba$age[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & pctba$tba == 0], na.rm=T)
    tempdf$age_yes_tba[i] <- mean(pctba$age[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & pctba$tba != 0], na.rm=T)
    tempdf$age.p.val[i] <- t.test(pctba$age[pctba$histology_abbreviation == tempdf$histology_abbreviation[i]] ~ ifelse(pctba$tba[pctba$histology_abbreviation == tempdf$histology_abbreviation[i]] == 0, "no_tba", "yes_tba"))$p.value
  }
  if (nrow(pctba[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$ESR1) & pctba$tba != 0,]) > 1 & nrow(pctba[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$ESR1) & pctba$tba == 0,]) > 1){
    tempdf$esr1_no_tba[i] <- mean(pctba$ESR1[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$ESR1) & pctba$tba == 0], na.rm=T)
    tempdf$esr1_yes_tba[i] <- mean(pctba$ESR1[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$ESR1) & pctba$tba != 0], na.rm=T)
    tempdf$esr1.p.val[i] <- t.test(pctba$ESR1[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$ESR1)] ~ ifelse(pctba$tba[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$ESR1)] == 0, "no_tba", "yes_tba"))$p.value
  }
  if (nrow(pctba[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$donor_survival_time) & pctba$tba != 0,]) & nrow(pctba[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$donor_survival_time) & pctba$tba == 0,]) > 1){
    tempdf$num_survinfo[i] <- sum(!is.na(pctba$donor_survival_time[pctba$histology_abbreviation == tempdf$histology_abbreviation[i]]))
    tempdf$num_surv_yes_tba[i] <- sum(!is.na(pctba$donor_survival_time[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & pctba$tba != 0]))
    tempdf$msurv_no_tba[i] <- median(pctba$donor_survival_time[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$donor_survival_time) & pctba$tba == 0], na.rm=T)
    tempdf$msurv_yes_tba[i] <- median(pctba$donor_survival_time[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$donor_survival_time) & pctba$tba != 0], na.rm=T)
    print(tempdf$histology_abbreviation[i])
    print(survdiff(Surv(pctba$donor_survival_time[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$donor_survival_time)], pctba$surv.event[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$donor_survival_time)]) ~ ifelse(pctba$tba[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$donor_survival_time)] == 0, "NO_TBA", "YES_TBA")))
    print(survfit(Surv(pctba$donor_survival_time[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$donor_survival_time)], pctba$surv.event[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$donor_survival_time)]) ~ ifelse(pctba$tba[pctba$histology_abbreviation == tempdf$histology_abbreviation[i] & !is.na(pctba$donor_survival_time)] == 0, "NO_TBA", "YES_TBA")))
  }
}
tempdf$esr1.fdr <- p.adjust(tempdf$esr1.p.val, method="fdr", length(tempdf$esr1.p.val))
tempdf$age.fdr <- p.adjust(tempdf$age.p.val, method="fdr", length(tempdf$age.p.val))

## Age comparison
tempdf[tempdf$yes_tba >=5 & tempdf$no_tba >=5,]
cutoff_value <- 5
infolist <- tempdf$histology_abbreviation[tempdf$yes_tba >=cutoff_value & tempdf$no_tba >=cutoff_value]
infolist <- infolist[infolist != "Breast-AdenoCA"]
df <- pctba[pctba$histology_abbreviation %in% infolist,]
df$merged <- paste0(df$histology_abbreviation, "_", as.character(df$tba))
df <- df[order(df$merged, decreasing = F),]

pdf("FIP_revision/boxplot.age.tbamp.vs.no.112522.pdf", height=5, width=4)
par(mar=c(5.1, 9.1, 4.1, 2.1))
boxplot(df$age ~ df$merged, horizontal=T, las=1, ylim=c(0,100), outline=F, frame=F, ylab="", xlab="Age", col=rep(unique(df$color), each=2))
stripchart(df$age ~ df$merged, pch=20, col=rgb(0,0,0,.1), vertical=F, method="jitter", add=T)
dev.off()

## ESR1 comparison
infolist <- infolist[infolist != "Eso-AdenoCA"]
infolist <- infolist[infolist != "Bone-Osteosarc"]
df <- pctba[pctba$histology_abbreviation %in% infolist,]
df$merged <- paste0(df$histology_abbreviation, "_", as.character(df$tba))
df <- df[order(df$merged, decreasing = F),]

pdf("FIP_revision/boxplot.esr1.tbamp.vs.no.112522.pdf", height=5, width=4)
par(mar=c(5.1, 9.1, 4.1, 2.1))
boxplot(log10(df$ESR1+1) ~ df$merged, horizontal=T, las=1, outline=F, frame=F, ylab="", xlab="ESR1 mRNA expression", col=rep(unique(df$color), each=2))
stripchart(log10(df$ESR1+1) ~ df$merged, pch=20, col=rgb(0,0,0,.1), vertical=F, method="jitter", add=T)
dev.off()

## Survival difference depending on the TB amplification status in each cancer type
for (i in tempdf$histology_abbreviation[!is.na(tempdf$msurv_no_tba) & tempdf$num_survinfo >50]){
  pdf(paste0("FIP_revision/Survival/survival.", i, ".tbamp.vs.no.112522.pdf"), height=5, width=5)
  plot(survfit(Surv(pctba$donor_survival_time[pctba$histology_abbreviation == i & !is.na(pctba$donor_survival_time)], pctba$surv.event[pctba$histology_abbreviation == i & !is.na(pctba$donor_survival_time)]) ~ ifelse(pctba$tba[pctba$histology_abbreviation == i & !is.na(pctba$donor_survival_time)] == 0, "NO_TBA", "YES_TBA")), frame=F, col=c("dark blue", "purple"), lwd=2, xaxt='n', yaxt='n', xlab="Time (years)", ylab="Overall survival (%)")
  axis(1, at=seq(0,10950,by=365), labels = seq(0,30,1))
  axis(2, at=seq(0,1,by=0.1), labels=seq(0,100,10), las=1)
  legend(2000, 1, c("without TB amp", "with TB amp"), lwd=2, col=c("dark blue", "purple"))
  dev.off()
}

