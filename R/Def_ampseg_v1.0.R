library(dplyr)
library(grid)
library(fields)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
tumortype <- args
sumdf <- read.csv("Summaryinfo.table.1.19.txt", sep="\t", as.is=T, header=T)

mainDir <- getwd()
subDir <- "Result_1mb"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
ifelse(!dir.exists(file.path(mainDir, subDir, tumortype)), dir.create(file.path(mainDir, subDir, tumortype)), FALSE)
ifelse(!dir.exists(file.path(mainDir, subDir, tumortype, "Bychr_copynumber")), dir.create(file.path(mainDir, subDir, tumortype, "Bychr_copynumber")), FALSE)
ifelse(!dir.exists(file.path(mainDir, subDir, tumortype, "CNV_amplified_regions")), dir.create(file.path(mainDir, subDir, tumortype, "CNV_amplified_regions")), FALSE)
ifelse(!dir.exists(file.path(mainDir, subDir, tumortype, "SV_amp_boundaries")), dir.create(file.path(mainDir, subDir, tumortype, "SV_amp_boundaries")), FALSE)
ifelse(!dir.exists(file.path(mainDir, subDir, tumortype, "SV_breakpoints_long")), dir.create(file.path(mainDir, subDir, tumortype, "SV_breakpoints_long")), FALSE)

if (tumortype == "Breast"){
  df <- sumdf[sumdf$histology_abbreviation %in% c("Breast-AdenoCA", "Breast-LobularCA", "Breast-DCIS"),]
} else {
  df <- sumdf[sumdf$histology_abbreviation == tumortype,]
}

df <- df[!is.na(df$SV.events),]
df <- df[df$SV.events != 0,]
df$ampchrs <- 0
df$chr_list <- "no"

cadf <- data.frame(matrix(NA, ncol = 10, nrow=0))
colnames(cadf) <- c("icgc_donor_id", "amp_id", "chr", "s", "e", "adj_s", "adj_e", "seglen", "sv_s", "sv_e")
cadf$icgc_donor_id <- as.character(cadf$icgc_donor_id)
cadf$amp_id <- as.character(cadf$amp_id)
cadf$chr <- as.character(cadf$chr)
cadf$s <- as.numeric(as.character(cadf$s))
cadf$e <- as.numeric(as.character(cadf$e))
cadf$adj_s <- as.numeric(as.character(cadf$adj_s))
cadf$adj_e <- as.numeric(as.character(cadf$adj_e))
cadf$seglen <- as.numeric(as.character(cadf$seglen))
cadf$sv_s <- as.character(cadf$sv_s)
cadf$sv_e <- as.character(cadf$sv_e)

hg19_coord <- read.csv("hg19_coord.txt", sep="\t", header=T, as.is=T)

for (z in 1:nrow(df)){
  tumorname <- df$icgc_donor_id[z]
  print(paste0("Processing ", tumorname, " ", z))
  
  cnvdf <- read.csv(paste0("/CNA_absCN_GC/", df$icgc_donor_id[z], ".somatic.cna.abscn.txt"), header=T, as.is=T, sep="\t")
  
  # Choose CNV File to Analyze (preferably JaBbA for its tight correlation with SV bkpts)
  cnvdf$chromosome <- as.character(cnvdf$chromosome)
  colnames(cnvdf) <- c("chr", "s", "e", "total_cn", "major_cn", "minor_cn", "depth")
  if (df$gender[z] == "female"){
    cnvdf <- cnvdf[cnvdf$chr != "Y",]
  }
  
  # Additional clearing of the dataset for chr arm annotation (not applicable to PCAWG consensus)
  #for (i in 1:nrow(cnvdf)){
  #  if (cnvdf$e[i] > hg19_coord$centro_end[hg19_coord$chr == cnvdf$chr[i]] & cnvdf$s[i] == 1){
  #    cnvdf[nrow(cnvdf)+1,] <- NA
  #    cnvdf$chr[nrow(cnvdf)] <- cnvdf$chr[i]
  #    cnvdf$s[nrow(cnvdf)] <- 1
  #    cnvdf$e[nrow(cnvdf)] <- hg19_coord$centro_end[hg19_coord$chr == cnvdf$chr[i]] - 1
  #    cnvdf$total_cn[nrow(cnvdf)] <- cnvdf$total_cn[i]
  #    cnvdf$major_cn[nrow(cnvdf)] <- cnvdf$major_cn[i]
  #    cnvdf$minor_cn[nrow(cnvdf)] <- cnvdf$minor_cn[i]
  #    cnvdf$depth[nrow(cnvdf)] <- cnvdf$depth[i]
  #    cnvdf$s[i] <- hg19_coord$centro_end[hg19_coord$chr == cnvdf$chr[i]]
  #  }
  #}
  #cnvdf$chr_num <- cnvdf$chr
  #cnvdf$chr_num[cnvdf$chr_num == "X"] <- "23"
  #cnvdf$chr_num[cnvdf$chr_num == "Y"] <- "24"
  #cnvdf$chr_num <- as.numeric(cnvdf$chr_num)
  #cnvdf <- cnvdf[order(cnvdf$s, decreasing = F),]
  #cnvdf <- cnvdf[order(cnvdf$chr_num, decreasing = F),]
  #cnvdf <- cnvdf[,c(1:(ncol(cnvdf)-1))]
  #
  #for (i in 1:nrow(cnvdf)){
  #  if (cnvdf$total_cn[i] == 1 & is.na(cnvdf$major_cn[i])){
  #    cnvdf$major_cn[i] <- 1
  #    cnvdf$minor_cn[i] <- 0
  #  } 
  #}
  #
  #for (i in 1:nrow(cnvdf)){
  #  if (cnvdf$total_cn[i] > 1 & cnvdf$major_cn[i] == 0 & cnvdf$minor_cn[i] == 0 & !is.na(cnvdf$major_cn[i])){
  #    cnvdf$major_cn[i] <- NA
  #    cnvdf$minor_cn[i] <- NA
  #  }
  #}
  
  # Annotation of chromosome arms
  cnvdf$arm = "no"
  for (i in 1:nrow(cnvdf)){
    if (cnvdf$e[i] < hg19_coord$centro_end[hg19_coord$chr == cnvdf$chr[i]]){
      cnvdf$arm[i] <- paste0(cnvdf$chr[i], "p")
    } else {
      cnvdf$arm[i] <- paste0(cnvdf$chr[i], "q")
    }
  }
  cnvdf$seglen <- cnvdf$e - cnvdf$s
  #cnvdf <- cnvdf[cnvdf$chr != "Y",] ## Pan-cancer cohort does have male patients
  
  for (i in unique(cnvdf$arm)){
    infodf1 = NULL
    if (length(sort(unique(cnvdf$total_cn[cnvdf$arm == i]))) != 0){
      for (k in sort(unique(cnvdf$total_cn[cnvdf$arm == i]))){
        lenvalue <- sum(cnvdf$seglen[cnvdf$arm == i & cnvdf$total_cn == k], na.rm = T)
        infoline1 <- cbind(i, k, lenvalue)
        infodf1 <- rbind(infodf1, infoline1)
      }
      infodf1 <- as.data.frame(infodf1)
      colnames(infodf1) <- c("arm", "total_cn", "lenvalue")
      infodf1$arm <- as.character(infodf1$arm)
      infodf1$total_cn <- as.numeric(as.character(infodf1$total_cn))
      infodf1$lenvalue <- as.numeric(as.character(infodf1$lenvalue))
      if (length(infodf1$total_cn[infodf1$lenvalue == max(infodf1$lenvalue)]) == 2){
        #cnvdf$chrbase[cnvdf$arm == i] <- infodf1$total_cn[order(infodf1$lenvalue, decreasing = T)][2]  # second largest!
        cnvdf$chrbase[cnvdf$arm == i] <- infodf1$total_cn[order(infodf1$lenvalue, decreasing = T)][2]
      } else {
        if (infodf1$total_cn[infodf1$lenvalue == max(infodf1$lenvalue)] == 0){
          #cnvdf$chrbase[cnvdf$arm == i] <- infodf1$total_cn[order(infodf1$lenvalue, decreasing = T)][2]  # second largest!
          cnvdf$chrbase[cnvdf$arm == i] <- infodf1$total_cn[order(infodf1$lenvalue, decreasing = T)][2]
        } else {
          #cnvdf$chrbase[cnvdf$arm == i] <- infodf1$total_cn[infodf1$lenvalue == max(infodf1$lenvalue)]  # reconsider!
          cnvdf$chrbase[cnvdf$arm == i] <- infodf1$total_cn[infodf1$lenvalue == max(infodf1$lenvalue)]
        }
      }
      
      infodf2 = NULL
      for (k in sort(unique(cnvdf$major_cn[cnvdf$arm == i]))){
        lenvalue <- sum(cnvdf$seglen[cnvdf$arm == i & cnvdf$major_cn == k], na.rm = T)
        infoline2 <- cbind(i, k, lenvalue)
        infodf2 <- rbind(infodf2, infoline2)
      }
      infodf2 <- as.data.frame(infodf2)
      if (nrow(infodf2) != 0){
        colnames(infodf2) <- c("arm", "major_cn", "lenvalue")
        infodf2$arm <- as.character(infodf2$arm)
        infodf2$major_cn <- as.numeric(as.character(infodf2$major_cn))
        infodf2$lenvalue <- as.numeric(as.character(infodf2$lenvalue))
        infodf3 = NULL
        for (k in sort(unique(cnvdf$minor_cn[cnvdf$arm == i]))){
          lenvalue <- sum(cnvdf$seglen[cnvdf$arm == i & cnvdf$minor_cn == k], na.rm = T)
          infoline3 <- cbind(i, k, lenvalue)
          infodf3 <- rbind(infodf3, infoline3)
        }
        infodf3 <- as.data.frame(infodf3)
        colnames(infodf3) <- c("arm", "minor_cn", "lenvalue")
        infodf3$arm <- as.character(infodf3$arm)
        infodf3$minor_cn <- as.numeric(as.character(infodf3$minor_cn))
        infodf3$lenvalue <- as.numeric(as.character(infodf3$lenvalue))
        infodf4 <- as.data.frame(sort(unique(c(unique(cnvdf$major_cn[cnvdf$arm == i]), unique(cnvdf$minor_cn[cnvdf$arm == i])))))
        colnames(infodf4)[1] = "copy_integer"
        infodf4$copy_integer <- as.numeric(as.character(infodf4$copy_integer))
        infodf4$major_length = NA
        infodf4$minor_length = NA
        for (k in 1:nrow(infodf4)){
          if (length(infodf2$lenvalue[infodf2$major_cn == infodf4$copy_integer[k]]) == 0){
            infodf4$major_length[k] <- 0
          } else {
            infodf4$major_length[k] <- infodf2$lenvalue[infodf2$major_cn == infodf4$copy_integer[k]]
          }
          if (length(infodf3$lenvalue[infodf3$minor_cn == infodf4$copy_integer[k]]) == 0){
            infodf4$minor_length[k] <- 0
          } else {
            infodf4$minor_length[k] <- infodf3$lenvalue[infodf3$minor_cn == infodf4$copy_integer[k]]
          }
        }
        infodf4$added_length <- infodf4$major_length + infodf4$minor_length
        if (length(infodf4$copy_integer[infodf4$added_length == max(infodf4$added_length)]) == 2){
          cnvdf$allelebase[cnvdf$arm == i] <- infodf4$copy_integer[order(infodf4$added_length, decreasing = T)][2]  # second largest!
        } else {
          if (infodf4$copy_integer[infodf4$added_length == max(infodf4$added_length)] == 0){
            cnvdf$allelebase[cnvdf$arm == i] <- infodf4$copy_integer[order(infodf4$added_length, decreasing = T)][2]  # second largest!
          } else {
            cnvdf$allelebase[cnvdf$arm == i] <- infodf4$copy_integer[infodf4$added_length == max(infodf4$added_length)] # reconsider!
          }
        }
        for (k in which(cnvdf$arm == i)){
          if (infodf4$minor_length[1] > 0.99*sum(infodf4$minor_length)){ # the first item usually has cn==0, meaning the LOH. But in case for complete ROH, this will assign those major and minor alleles appropriately.
            cnvdf$alt_cn[k] <- cnvdf$major_cn[k]
            cnvdf$unalt_cn[k] <- cnvdf$minor_cn[k]
          } else {
            if (!is.na(cnvdf$total_cn[k]) & !is.na(cnvdf$major_cn[k])){
              majorcopy <- as.numeric(cnvdf$major_cn[k])
              minorcopy <- as.numeric(cnvdf$minor_cn[k])
              allelic_baseline <- as.numeric(cnvdf$allelebase[k])
              if (abs(majorcopy-allelic_baseline) > abs(minorcopy-allelic_baseline)){
                cnvdf$alt_cn[k] <- majorcopy
                cnvdf$unalt_cn[k] <- minorcopy
              } else if (abs(majorcopy-allelic_baseline) < abs(minorcopy-allelic_baseline)){
                cnvdf$alt_cn[k] <- minorcopy
                cnvdf$unalt_cn[k] <- majorcopy
              } else if (abs(majorcopy-allelic_baseline) == abs(minorcopy-allelic_baseline) & majorcopy == minorcopy){
                cnvdf$alt_cn[k] <- majorcopy
                cnvdf$unalt_cn[k] <- minorcopy
              } else if (abs(majorcopy-allelic_baseline) == abs(minorcopy-allelic_baseline) & majorcopy != minorcopy){
                cnvdf$alt_cn[k] <- min(majorcopy, minorcopy) # reconsider!
                cnvdf$unalt_cn[k] <- max(majorcopy, minorcopy) # reconsider!
              } else {
                cnvdf$alt_cn[k] <- NA
                cnvdf$unalt_cn[k] <- NA
              }
            } else {
              cnvdf$alt_cn[k] <- NA # added to deal with the cases with first row NA values
              cnvdf$unalt_cn[k] <- NA
            }
          }
        }
        # Plotting Each Chromosomal CNs
        pdf(paste0(subDir, "/", tumortype, "/Bychr_copynumber/", tumorname, ".arm", i, ".cnhisto.pdf"), height=6, width=10)
        par(mfrow=c(1,2))
        plot(infodf4$copy_integer, log(infodf4$added_length), ylim = c(0, max(log(infodf4$added_length))), type='h', col="blue", lwd=2, xlab = paste0("Copy number in chr", i), ylab = "Log(genomic length)", frame = FALSE, las = 1, main = paste0(tumorname, ", chr", i))
        plot(0, type = 'n', xlim = c(min(cnvdf$s[cnvdf$arm == i]), max(cnvdf$e[cnvdf$arm == i])), ylim = c(0, 10), xlab = paste0("Positions in chr", i), ylab = "Total copy number", frame = FALSE, las = 1, main = paste0(tumorname, ", chr", i))
        for (w in 1:nrow(cnvdf[cnvdf$arm == i,])){
          segments(cnvdf$s[cnvdf$arm == i][w], cnvdf$major_cn[cnvdf$arm == i][w]+0.05, cnvdf$e[cnvdf$arm == i][w], cnvdf$major_cn[cnvdf$arm == i][w]+0.05, lty = 1, lwd = 1, col = "red")
          segments(cnvdf$s[cnvdf$arm == i][w], cnvdf$minor_cn[cnvdf$arm == i][w], cnvdf$e[cnvdf$arm == i][w], cnvdf$minor_cn[cnvdf$arm == i][w], lty = 1, lwd = 1, col = "blue")
        }
        dev.off()
      }
    } 
  }
  
  # Adjustment of chrbase values around centromere
  for (i in unique(cnvdf$chr)){
    elp <- rev(which(cnvdf$arm == paste0(i, "p")))
    elq <- which(cnvdf$arm == paste0(i, "q"))
    pbase <- cnvdf$chrbase[elp[1]]
    qbase <- cnvdf$chrbase[elq[1]]
    if (!is.na(pbase) & !is.na(qbase)){
      if (pbase == qbase){
        next
      } else {
        if (!is.na(cnvdf$total_cn[elp[1]]) & !is.na(cnvdf$total_cn[elq[1]])){
          if (abs(cnvdf$total_cn[elp[1]] - pbase) > abs(cnvdf$total_cn[elp[1]] - qbase) & abs(cnvdf$total_cn[elq[1]] - pbase) > abs(cnvdf$total_cn[elq[1]] - qbase)){ ## when p arm is misassigned
            for (j in 1:length(elp)){
              if (!is.na(cnvdf$total_cn[elp[j]]) & abs(cnvdf$total_cn[elp[j]] - pbase) > abs(cnvdf$total_cn[elp[j]] - qbase)){
                cnvdf$chrbase[elp[j]] <- qbase
                cnvdf$allelebase[elp[j]] <- cnvdf$allelebase[elq[1]]
              } else {
                break
              }
            }
          } else if (abs(cnvdf$total_cn[elp[1]] - pbase) < abs(cnvdf$total_cn[elp[1]] - qbase) & abs(cnvdf$total_cn[elq[1]] - pbase) < abs(cnvdf$total_cn[elq[1]] - qbase)){ ## when q arm is misassigned
            for (j in 1:length(elq)){
              if (!is.na(cnvdf$total_cn[elq[j]]) & abs(cnvdf$total_cn[elq[j]] - qbase) > abs(cnvdf$total_cn[elq[j]] - pbase)){
                cnvdf$chrbase[elq[j]] <- pbase
                cnvdf$allelebase[elq[j]] <- cnvdf$allelebase[elp[1]]
              } else {
                break
              }
            }
          }
        }
      }
    }
  }
  
  # Initial annotation of amplified regions
  cnvdf$ampf = "no"
  if (exists("modeval")){
    rm(modeval)
  }
  for (i in 1:nrow(cnvdf)){
    if (!is.na(cnvdf$chrbase[i])){
      if (cnvdf$depth[i] > 3*cnvdf$chrbase[i] & cnvdf$depth[i] >= 6 & !is.na(cnvdf$depth[i])){
        modeval <- "amp"
        cnvdf$ampf[i] <- modeval
      } else if (cnvdf$depth[i] <= 3*cnvdf$chrbase[i] & !is.na(cnvdf$depth[i])){
        if (cnvdf$depth[i] >= cnvdf$chrbase[i]+6 & !is.na(cnvdf$depth[i])){
          modeval <- "amp"
          cnvdf$ampf[i] <- modeval
        } else {
          modeval <- "no"
          cnvdf$ampf[i] <- modeval
        }
      } else {
        if (exists("modeval")){
          cnvdf$ampf[i] <- modeval
        } else {
          modeval <- "no"
          cnvdf$ampf[i] <- modeval
        }
      }
    }
  }
  cnvdf$ampr = "no" # reverse annotation
  if (exists("modeval")){
    rm(modeval)
  }
  for (i in nrow(cnvdf):1){
    if (!is.na(cnvdf$chrbase[i])){
      if (cnvdf$depth[i] > 3*cnvdf$chrbase[i] & cnvdf$depth[i] >= 6 & !is.na(cnvdf$depth[i])){
        modeval <- "amp"
        cnvdf$ampr[i] <- modeval
      } else if (cnvdf$depth[i] <= 3*cnvdf$chrbase[i] & !is.na(cnvdf$depth[i])){
        if (cnvdf$depth[i] >= cnvdf$chrbase[i]+6 & !is.na(cnvdf$depth[i])){
          modeval <- "amp"
          cnvdf$ampr[i] <- modeval
        } else {
          modeval <- "no"
          cnvdf$ampr[i] <- modeval
        }
      } else {
        if (exists("modeval")){
          cnvdf$ampr[i] <- modeval
        } else {
          modeval <- "no"
          cnvdf$ampr[i] <- modeval
        }
      }
    }
  }
  cnvdf$amp = "no"
  for (i in 1:nrow(cnvdf)){
    if (cnvdf$ampf[i] == "amp" | cnvdf$ampr[i] == "amp"){
      cnvdf$amp[i] <- "amp"
    }
  }
  cnvdf <- cnvdf[,c(1:(ncol(cnvdf)-3), ncol(cnvdf))]
  

  # Annotate Boundaries of Amplified Regions (top of the mountains)
  cnvdf$cgr = "no"
  prev_chr <- cnvdf$chr[1]
  prev_s <- cnvdf$s[1]
  prev_e <- cnvdf$e[1]
  prev_amp <- cnvdf$amp[1]
  if (prev_amp == "amp" & cnvdf$amp[2] == "amp"){
    cnvdf$cgr[1] <- "amplicon_start"
  } else if (prev_amp == "amp" & cnvdf$amp[2] == "no"){
    cnvdf$cgr[1] <- "amplicon_itself"
  }
  #if (cnvdf$amp[nrow(cnvdf)] == "amp" & cnvdf$amp[nrow(cnvdf)-1] == "amp"){
  #  cnvdf$cgr[nrow(cnvdf)] <- "amplicon_end"
  #} else if (cnvdf$amp[nrow(cnvdf)] == "amp" & cnvdf$amp[nrow(cnvdf)-1] == "no"){
  #  cnvdf$cgr[nrow(cnvdf)] <- "amplicon_itself"
  #}
  for (i in 2:nrow(cnvdf)){
    curr_chr <- cnvdf$chr[i]
    curr_s <- cnvdf$s[i]
    curr_e <- cnvdf$e[i]
    curr_amp <- cnvdf$amp[i]
    if (curr_amp == "amp" & prev_amp == "no" & i != nrow(cnvdf)){
      cnvdf$cgr[i] <- paste0("amplicon_start")
      prev_chr <- curr_chr
      prev_s <- curr_s
      prev_e <- curr_e
      prev_amp <- curr_amp
    } else if (curr_amp == "amp" & prev_amp == "no" & i == nrow(cnvdf)){
      cnvdf$cgr[i] <- paste0("amplicon_itself")
      prev_chr <- curr_chr
      prev_s <- curr_s
      prev_e <- curr_e
      prev_amp <- curr_amp
    } else if (curr_amp == "amp" & prev_amp == "amp" & i == nrow(cnvdf)){
      cnvdf$cgr[i] <- paste0("amplicon_end")
      prev_chr <- curr_chr
      prev_s <- curr_s
      prev_e <- curr_e
      prev_amp <- curr_amp
    } else if (curr_amp == "amp" & prev_amp == "amp" & prev_chr != curr_chr & cnvdf$cgr[i-1] == "amplicon_start"){
      cnvdf$cgr[i] <- paste0("amplicon_start")
      cnvdf$cgr[i-1] <- paste0("amplicon_itself")
      prev_chr <- curr_chr
      prev_s <- curr_s
      prev_e <- curr_e
      prev_amp <- curr_amp
    } else if (curr_amp == "amp" & prev_amp == "amp" & prev_chr != curr_chr & cnvdf$cgr[i-1] != "amplicon_start"){
      cnvdf$cgr[i] <- paste0("amplicon_start")
      cnvdf$cgr[i-1] <- paste0("amplicon_end")
      prev_chr <- curr_chr
      prev_s <- curr_s
      prev_e <- curr_e
      prev_amp <- curr_amp
    } else if (curr_amp == "no" & prev_amp == "amp" & cnvdf$cgr[i-1] != "amplicon_start" & cnvdf$cgr[i-1] != "amplicon_itself"){ ## corrected
      cnvdf$cgr[i-1] <- paste0("amplicon_end")
      prev_chr <- curr_chr
      prev_s <- curr_s
      prev_e <- curr_e
      prev_amp <- curr_amp
    } else if (curr_amp == "no" & prev_amp == "amp" & cnvdf$cgr[i-1] == "amplicon_start"){
      cnvdf$cgr[i-1] <- paste0("amplicon_itself")
      prev_chr <- curr_chr
      prev_s <- curr_s
      prev_e <- curr_e
      prev_amp <- curr_amp
    } else {
      prev_chr <- curr_chr
      prev_s <- curr_s
      prev_e <- curr_e
      prev_amp <- curr_amp
    }
  }
  
  
  # Merge Contiguous Amplified Regions (connect tops of mountains if they are sharing the bases of less than 3Mbp)
  cnvdf$merged_cgr <- cnvdf$cgr  # we can copy the cgr values for this
  amp_mode <- "nothing"
  for (i in 2:(nrow(cnvdf)-1)){
    if ((cnvdf$cgr[i] == "amplicon_itself" | cnvdf$cgr[i] == "amplicon_end") & cnvdf$cgr[i+1] == "no" & amp_mode == "nothing"){
      interval_start <- i+1
      amp_mode <- "interval"
    } else if (cnvdf$cgr[i] == "no" & (cnvdf$cgr[i+1] == "amplicon_itself" | cnvdf$cgr[i+1] == "amplicon_start") & amp_mode == "interval"){
      interval_end <- i
      if (cnvdf$chr[interval_start] == cnvdf$chr[interval_end] & sum(cnvdf$seglen[interval_start:interval_end], na.rm = T) < 3000000 & all(cnvdf$depth[interval_start:interval_end][!is.na(cnvdf$depth[interval_start:interval_end])] >= 2*cnvdf$chrbase[interval_end]) & all(cnvdf$depth[interval_start:interval_end][!is.na(cnvdf$depth[interval_start:interval_end])] >= 4)){
#      if (cnvdf$chr[interval_start] == cnvdf$chr[interval_end] & sum(cnvdf$seglen[interval_start:interval_end], na.rm = T) < 3000000 & weighted.mean(cnvdf$depth[interval_start:interval_end], cnvdf$seglen[interval_start:interval_end], na.rm = T) > 2*cnvdf$chrbase[interval_end] & all(cnvdf$depth[interval_start:interval_end][!is.na(cnvdf$depth[interval_start:interval_end])] > cnvdf$chrbase[interval_end])){
        #cnvdf$merged_cgr[interval_start:interval_end] <- "merged"
        if (cnvdf$merged_cgr[interval_start-1] == "amplicon_end"){
          cnvdf$merged_cgr[interval_start-1] <- "no"
        } else if (cnvdf$merged_cgr[interval_start-1] == "amplicon_itself"){
          cnvdf$merged_cgr[interval_start-1] <- "amplicon_start"
        }
        if (cnvdf$merged_cgr[interval_end+1] == "amplicon_start"){
          cnvdf$merged_cgr[interval_end+1] <- "no"
        } else if (cnvdf$merged_cgr[interval_end+1] == "amplicon_itself"){
          cnvdf$merged_cgr[interval_end+1] <- "amplicon_end"
        }
        #sum(cnvdf$seglen[interval_start:interval_end], na.rm = T)
        #weighted.mean(cnvdf$alt_cn[interval_start:interval_end], cnvdf$seglen[interval_start:interval_end], na.rm = T)
      }
      amp_mode <- "nothing"
    }
  }
  
  # Merge Adjacent Amplicons Bordered by Balanced Bkpts (<1000 bp)
  gap_mode <- "nothing"
  for (i in 2:(nrow(cnvdf)-1)){
    if ((cnvdf$merged_cgr[i] == "amplicon_end" | cnvdf$merged_cgr[i] == "amplicon_itself") & gap_mode == "nothing"){
      interval_start <- i
      gap_mode <- "gap"
    } else if ((cnvdf$merged_cgr[i] == "amplicon_start" | cnvdf$merged_cgr[i] == "amplicon_itself") & gap_mode == "gap"){
      interval_end <- i
      if ((cnvdf$s[interval_end] - cnvdf$e[interval_start]) < 1000 & cnvdf$chr[interval_end] == cnvdf$chr[interval_start]){
        if (cnvdf$merged_cgr[interval_start] == "amplicon_end"){
          cnvdf$merged_cgr[interval_start] <- "no"
        } else if (cnvdf$merged_cgr[interval_start] == "amplicon_itself"){
          cnvdf$merged_cgr[interval_start] <- "amplicon_start"
        }
        if (cnvdf$merged_cgr[interval_end] == "amplicon_start"){
          cnvdf$merged_cgr[interval_end] <- "no"
        } else if (cnvdf$merged_cgr[interval_end] == "amplicon_itself"){
          cnvdf$merged_cgr[interval_end] <- "amplicon_end"
        }
      }
      gap_mode <- "nothing"
    }
  }
  
  
  # Filter out small peaks
  cnvdf$filter <- cnvdf$merged_cgr
  for (i in 2:(nrow(cnvdf)-1)){
    if (cnvdf$merged_cgr[i] == "amplicon_itself" & cnvdf$chr[i] == cnvdf$chr[i-1] & cnvdf$chr[i] == cnvdf$chr[i+1] & !is.na(cnvdf$depth[i-1]) & !is.na(cnvdf$depth[i+1])){
      if ((cnvdf$depth[i] - cnvdf$depth[i-1]) < 3 & (cnvdf$depth[i] - cnvdf$depth[i+1]) < 3 & (cnvdf$seglen[i-1] > 3000000 | cnvdf$seglen[i+1] > 3000000)){
        cnvdf$filter[i] <- "no"
      }
    }
  }
  
  # Filter out global elevations
  modeval <- "no"
  for (i in 1:nrow(cnvdf)){
    if (cnvdf$merged_cgr[i] == "amplicon_start" & modeval == "no"){
      modeval <- "amp"
      if (i == which(cnvdf$chr == cnvdf$chr[i])[1]){
        startval <- "editable"
        start_num <- i
      } else {
        if (cnvdf$depth[i] - cnvdf$depth[i-1] < 3 & !is.na(cnvdf$depth[i]) & !is.na(cnvdf$depth[i-1]) & cnvdf$seglen[i-1] > 3000000){
          startval <- "editable"
          start_num <- i
        } else if (is.na(cnvdf$depth[i]) | is.na(cnvdf$depth[i-1])){
          startval <- "editable"
          start_num <- i
        } else {
          startval <- "okay"
          start_num <- i
        }
      }
    } else if (cnvdf$merged_cgr[i] == "amplicon_end" & modeval == "amp"){
      if (i == rev(which(cnvdf$chr == cnvdf$chr[i]))[1]){
        endval <- "editable"
        end_num <- i
      } else {
        if (cnvdf$depth[i] - cnvdf$depth[i+1] < 3 & !is.na(cnvdf$depth[i]) & !is.na(cnvdf$depth[i+1]) & cnvdf$seglen[i+1] > 3000000){
          endval <- "editable"
          end_num <- i
        } else if (is.na(cnvdf$depth[i]) | is.na(cnvdf$depth[i+1])){
          endval <- "editable"
          end_num <- i
        } else {
          endval <- "okay"
          end_num <- i
        }
      }
      if (startval == "editable" & endval == "editable"){
        cnvdf$filter[start_num] <- "no"
        cnvdf$filter[end_num] <- "no"
      }
      modeval <- "no"
    }
  }
  
  
  # Expand the Merged Amplified Regions to the CN Boundaries (Borders between LOH/ROH or CN down to the allelic baseline)
  exprange <- 1000000
  cnvdf$expanded_cgr <- cnvdf$filter
  for (i in 2:(nrow(cnvdf)-1)){
    if (cnvdf$filter[i] == "amplicon_start"){
      curr_chr <- cnvdf$chr[i]
      curr_s <- cnvdf$s[i]
      if (curr_s > exprange){
        if (curr_s - exprange < cnvdf$s[which(cnvdf$chr == curr_chr)[1]]){ # exceptional case of no report of telomeric segment
          anchor_val <- which(cnvdf$chr == curr_chr)[1]
        } else {
          anchor_val <- which(cnvdf$chr == curr_chr & cnvdf$e > (curr_s-exprange) & cnvdf$s <= (curr_s-exprange))
        }
      } else {
        #anchor_val <- which(cnvdf$chr == curr_chr & cnvdf$e > 1 & cnvdf$s <= 1) # corrected d/t Error in (i - 1):anchor_val : argument of length 0
        anchor_val <- which(cnvdf$chr == curr_chr)[1]
      }
      for (j in (i-1):anchor_val){
        if (cnvdf$filter[j] != "no"){
          imp_num1 <- i
          break
        } else if (cnvdf$depth[j] <= cnvdf$chrbase[j] & !is.na(cnvdf$depth[j]) & !is.na(cnvdf$chrbase[j])){
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
    } else if (cnvdf$filter[i] == "amplicon_end"){
      curr_chr <- cnvdf$chr[i]
      curr_e <- cnvdf$e[i]
      if ((max(cnvdf$e[cnvdf$chr == curr_chr])-curr_e) > exprange){
        if (curr_e + exprange > cnvdf$e[which(cnvdf$chr == curr_chr)[length(which(cnvdf$chr == curr_chr))]]){ # exceptional case of no report of telomeric segment
          anchor_val <- which(cnvdf$chr == curr_chr)[length(which(cnvdf$chr == curr_chr))]
        } else {
          anchor_val <- which(cnvdf$chr == curr_chr & cnvdf$e >= (curr_e+exprange) & cnvdf$s <= (curr_e+exprange))
        }
      } else {
        #anchor_val <- which(cnvdf$chr == curr_chr & cnvdf$e == max(cnvdf$e[cnvdf$chr == curr_chr]))
        anchor_val <- which(cnvdf$chr == curr_chr)[length(which(cnvdf$chr == curr_chr))]
      }
      for (j in (i+1):anchor_val){
        if (cnvdf$filter[j] != "no"){
          imp_num2 <- i
          break
        } else if (cnvdf$depth[j] <= cnvdf$chrbase[j] & !is.na(cnvdf$depth[j]) & !is.na(cnvdf$chrbase[j])){
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
    } else if (cnvdf$filter[i] == "amplicon_itself"){
      curr_chr <- cnvdf$chr[i]
      curr_s <- cnvdf$s[i]
      if (curr_s > exprange){
        if (curr_s - exprange < cnvdf$s[which(cnvdf$chr == curr_chr)[1]]){ # exceptional case of no report of telomeric segment
          anchor_val <- which(cnvdf$chr == curr_chr)[1]
        } else {
          anchor_val <- which(cnvdf$chr == curr_chr & cnvdf$e >= (curr_s-exprange) & cnvdf$s <= (curr_s-exprange))
        }
      } else {
        #anchor_val <- which(cnvdf$chr == curr_chr & cnvdf$e > 1 & cnvdf$s <= 1) # corrected d/t Error in (i - 1):anchor_val : argument of length 0
        anchor_val <- which(cnvdf$chr == curr_chr)[1]
      }
      for (j in (i-1):anchor_val){
        if (cnvdf$filter[j] != "no"){
          imp_num1 <- i
          break
        } else if (cnvdf$depth[j] <= cnvdf$chrbase[j] & !is.na(cnvdf$depth[j]) & !is.na(cnvdf$chrbase[j])){
          imp_num1 <- j+1
          break
        } else {
          imp_num1 <- i
        }
      }
      curr_e <- cnvdf$e[i]
      if ((max(cnvdf$e[cnvdf$chr == curr_chr])-curr_e) > exprange){
        if (curr_e + exprange > cnvdf$e[which(cnvdf$chr == curr_chr)[length(which(cnvdf$chr == curr_chr))]]){ # exceptional case of no report of telomeric segment
          anchor_val <- which(cnvdf$chr == curr_chr)[length(which(cnvdf$chr == curr_chr))]
        } else {
          anchor_val <- which(cnvdf$chr == curr_chr & cnvdf$e >= (curr_e+exprange) & cnvdf$s <= (curr_e+exprange))
        }
      } else {
        #anchor_val <- which(cnvdf$chr == curr_chr & cnvdf$e == max(cnvdf$e[cnvdf$chr == curr_chr]))
        anchor_val <- which(cnvdf$chr == curr_chr)[length(which(cnvdf$chr == curr_chr))]
      }
      for (j in (i+1):anchor_val){
        if (cnvdf$filter[j] != "no"){
          imp_num2 <- i
          break
        } else if (cnvdf$depth[j] <= cnvdf$chrbase[j] & !is.na(cnvdf$depth[j]) & !is.na(cnvdf$chrbase[j])){
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
  
  #write.csv(cnvdf, paste0(df$icgc_donor_id[z], ".temporary.csv"), row.names=F)
  write.table(cnvdf, paste0(subDir, "/", tumortype, "/CNV_amplified_regions/", tumorname, ".amplified_regions.txt"), row.names=F, quote=F, sep="\t", col.names=T)
  
  
  # Summary of Merged/Expanded Amplified Regions
  ampdf <- data.frame(matrix(NA, ncol = 7, nrow=0))
  colnames(ampdf) <- c("icgc_donor_id", "amp_id", "chr", "s", "e", "adj_s", "adj_e")
  ampdf$icgc_donor_id <- as.character(ampdf$icgc_donor_id)
  ampdf$amp_id <- as.character(ampdf$amp_id)
  ampdf$chr <- as.character(ampdf$chr)
  ampdf$s <- as.numeric(as.character(ampdf$s))
  ampdf$e <- as.numeric(as.character(ampdf$e))
  ampdf$adj_s <- as.numeric(as.character(ampdf$adj_s))
  ampdf$adj_e <- as.numeric(as.character(ampdf$adj_e))
  
  imp_num <- 1
  for (i in 1:nrow(cnvdf)){
    if (cnvdf$expanded_cgr[i] == "amplicon_start"){
      ampdf[nrow(ampdf)+1,] <- NA
      ampdf$icgc_donor_id[imp_num] <- tumorname
      ampdf$amp_id[imp_num] <- paste0("amplicon_", imp_num)
      ampdf$chr[imp_num] <- cnvdf$chr[i]
      ampdf$s[imp_num] <- cnvdf$s[i]
    } else if (cnvdf$expanded_cgr[i] == "amplicon_end"){
      ampdf$e[imp_num] <- cnvdf$e[i]
      imp_num <- imp_num+1
    } else if (cnvdf$expanded_cgr[i] == "amplicon_itself"){
      ampdf[nrow(ampdf)+1,] <- NA
      ampdf$icgc_donor_id[imp_num] <- tumorname
      ampdf$amp_id[imp_num] <- paste0("amplicon_", imp_num)
      ampdf$chr[imp_num] <- cnvdf$chr[i]
      ampdf$s[imp_num] <- cnvdf$s[i]
      ampdf$e[imp_num] <- cnvdf$e[i]
      imp_num <- imp_num+1
    }
  }
  ampdf$seglen <- ampdf$e-ampdf$s
  if (nrow(ampdf) == 0){
    df$ampchrs[z] <- 0
    df$chr_list[z] <- "no"
  } else {
    df$ampchrs[z] <- length(unique(ampdf$chr))
    df$chr_list[z] <- paste(unique(ampdf$chr), collapse = ";")
    
    # Loading SV Datasets and Creating a Long Format
    svdf <- read.csv(paste0("/SV_BEDPE_annotated/", sumdf$tumor[sumdf$icgc_donor_id == df$icgc_donor_id[z]], ".pcawg_consensus_1.6.161116.somatic.sv.bedpe.anv.txt"), header=T, as.is=T, sep="\t")
    svdf$chrom1 <- as.character(svdf$chrom1)
    svdf$chrom2 <- as.character(svdf$chrom2)
    svdf$seglen = NA
    svdf$fbi = "no"
    svdf$cn = "no"
    for (i in 1:nrow(svdf)){
      if ((svdf$svclass[i] == "h2hINV" | svdf$svclass[i] == "t2tINV") & svdf$start2[i]-svdf$start1[i] < 5000){
        #if ((svdf$svclass[i] == "h2hINV" | svdf$svclass[i] == "t2tINV") & svdf$start2[i]-svdf$start1[i] < 5000 & length(unique(c(cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i]-5000 & cnvdf$e >= svdf$start1[i]-5000], cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i] & cnvdf$e >= svdf$start1[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i] & cnvdf$e >= svdf$start2[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i]+5000 & cnvdf$e >= svdf$start2[i]+5000])))!=1){
        if (nrow(cnvdf[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i]-5000,])==0 | nrow(cnvdf[cnvdf$chr==svdf$chrom2[i] & cnvdf$e >= svdf$start2[i]+5000,])==0){
          svdf$seglen[i] <- svdf$start2[i]-svdf$start1[i]
          svdf$cn[i] <- paste(cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i]-5000 & cnvdf$e >= svdf$start1[i]-5000], cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i] & cnvdf$e >= svdf$start1[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i] & cnvdf$e >= svdf$start2[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i]+5000 & cnvdf$e >= svdf$start2[i]+5000], sep = ";")
          if (length(unique(c(cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i]-5000 & cnvdf$e >= svdf$start1[i]-5000], cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i] & cnvdf$e >= svdf$start1[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i] & cnvdf$e >= svdf$start2[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i]+5000 & cnvdf$e >= svdf$start2[i]+5000]))) > 1){
            svdf$fbi[i] = "possible"
          }
        } else if (is.na(cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i]-5000 & cnvdf$e >= svdf$start1[i]-5000]) & is.na(cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i] & cnvdf$e >= svdf$start1[i]]) & is.na(cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i] & cnvdf$e >= svdf$start2[i]]) & is.na(cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i]+5000 & cnvdf$e >= svdf$start2[i]+5000])){
          svdf$seglen[i] <- svdf$start2[i]-svdf$start1[i]
          svdf$cn[i] <- paste(cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i]-5000 & cnvdf$e >= svdf$start1[i]-5000], cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i] & cnvdf$e >= svdf$start1[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i] & cnvdf$e >= svdf$start2[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i]+5000 & cnvdf$e >= svdf$start2[i]+5000], sep = ";")
          if (length(unique(c(cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i]-5000 & cnvdf$e >= svdf$start1[i]-5000], cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i] & cnvdf$e >= svdf$start1[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i] & cnvdf$e >= svdf$start2[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i]+5000 & cnvdf$e >= svdf$start2[i]+5000]))) > 1){
            svdf$fbi[i] = "possible"
          }
        } else if (length(unique(c(cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i]-5000 & cnvdf$e >= svdf$start1[i]-5000], cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i] & cnvdf$e >= svdf$start1[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i] & cnvdf$e >= svdf$start2[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i]+5000 & cnvdf$e >= svdf$start2[i]+5000])))!=1){
          svdf$seglen[i] <- svdf$start2[i]-svdf$start1[i]
          svdf$cn[i] <- paste(cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i]-5000 & cnvdf$e >= svdf$start1[i]-5000], cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i] & cnvdf$e >= svdf$start1[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i] & cnvdf$e >= svdf$start2[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i]+5000 & cnvdf$e >= svdf$start2[i]+5000], sep = ";")
          if (length(unique(c(cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i]-5000 & cnvdf$e >= svdf$start1[i]-5000], cnvdf$depth[cnvdf$chr==svdf$chrom1[i] & cnvdf$s <= svdf$start1[i] & cnvdf$e >= svdf$start1[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i] & cnvdf$e >= svdf$start2[i]], cnvdf$depth[cnvdf$chr==svdf$chrom2[i] & cnvdf$s <= svdf$start2[i]+5000 & cnvdf$e >= svdf$start2[i]+5000]))) > 1){
            svdf$fbi[i] = "possible"
          }
        }
      }
    }
    
    tempdf1 <- svdf[,c("chrom1", "start1", "strand1", "sv_id", "pe_support", "svclass", "fbi")]
    tempdf1$sv_id2 <- paste0(tempdf1$sv_id, "-1")
    colnames(tempdf1) <- c("chr", "s", "ori", "sv_id", "pe_support", "svclass", "fbi", "sv_id2")
    tempdf2 <- svdf[,c("chrom2", "start2", "strand2", "sv_id", "pe_support", "svclass", "fbi")]
    tempdf2$sv_id2 <- paste0(tempdf2$sv_id, "-2")
    colnames(tempdf2) <- c("chr", "s", "ori", "sv_id", "pe_support", "svclass", "fbi", "sv_id2")
    bpdf <- rbind(tempdf1, tempdf2)
    rm(tempdf1, tempdf2)
    bpdf$chr_num <- bpdf$chr
    bpdf$chr_num[bpdf$chr_num == "X"] <- "23"
    bpdf$chr_num[bpdf$chr_num == "Y"] <- "24"
    bpdf$chr_num[bpdf$chr_num == "MT"] <- "25"
    bpdf$chr_num <- as.numeric(bpdf$chr_num)
    bpdf <- bpdf[order(bpdf$s, decreasing = F),]
    bpdf <- bpdf[order(bpdf$chr_num, decreasing = F),]
    bpdf <- bpdf[,c(1:(ncol(bpdf)-1))]
    
    # Correlating Boundaries of Amplified Regions with SVs
    svmargin = 10000
    ampdf$sv_s <- "no"
    ampdf$sv_e <- "no"
    svdf$amp_boundary <- "no"
    bpdf$amp_boundary <- "no"
    for (i in 1:nrow(ampdf)){
      tempdf <- bpdf[bpdf$chr == ampdf$chr[i],]
      if (nrow(tempdf) != 0){
        # Left boundary of the amplicon
        if (length(which(abs(tempdf$s - ampdf$s[i]) == min(abs(tempdf$s - ampdf$s[i])))) == 1){
          if (tempdf$ori[which(abs(tempdf$s - ampdf$s[i]) == min(abs(tempdf$s - ampdf$s[i])))] == "-" & abs(tempdf$s[which(abs(tempdf$s - ampdf$s[i]) == min(abs(tempdf$s - ampdf$s[i])))] - ampdf$s[i]) < svmargin){
            j <- which(abs(tempdf$s - ampdf$s[i]) == min(abs(tempdf$s - ampdf$s[i])))
            ampdf$sv_s[i] <- tempdf$sv_id[j]
            ampdf$adj_s[i] <- tempdf$s[j]
            svdf$amp_boundary[svdf$sv_id==tempdf$sv_id[j]] <- "boundary"
            bpdf$amp_boundary[bpdf$sv_id==tempdf$sv_id[j]] <- "boundary"
          }
        } else if (length(which(abs(tempdf$s - ampdf$s[i]) == min(abs(tempdf$s - ampdf$s[i])) & abs(tempdf$s - ampdf$s[i]) < svmargin)) > 1 & ("-" %in% tempdf$ori[which(abs(tempdf$s - ampdf$s[i]) == min(abs(tempdf$s - ampdf$s[i])))])){
          first_value <- 0
          for (j in which(abs(tempdf$s - ampdf$s[i]) == min(abs(tempdf$s - ampdf$s[i])) & tempdf$ori == "-" & abs(tempdf$s - ampdf$s[i]) < svmargin)){
            next_value <- tempdf$pe_support[j]
            if (first_value < next_value){
              first_value <- next_value
              imp_num <- j
            }
          }
          ampdf$sv_s[i] <- tempdf$sv_id[j]
          ampdf$adj_s[i] <- tempdf$s[j]
          svdf$amp_boundary[svdf$sv_id==tempdf$sv_id[j]] <- "boundary"
          bpdf$amp_boundary[bpdf$sv_id==tempdf$sv_id[j]] <- "boundary"
        }
        # Right boundary of the amplicon
        if (length(which(abs(tempdf$s - ampdf$e[i]) == min(abs(tempdf$s - ampdf$e[i])))) == 1){
          if (tempdf$ori[which(abs(tempdf$s - ampdf$e[i]) == min(abs(tempdf$s - ampdf$e[i])))] == "+" & abs(tempdf$s[which(abs(tempdf$s - ampdf$e[i]) == min(abs(tempdf$s - ampdf$e[i])))] - ampdf$e[i]) < svmargin){
            j <- which(abs(tempdf$s - ampdf$e[i]) == min(abs(tempdf$s - ampdf$e[i])))
            ampdf$sv_e[i] <- tempdf$sv_id[j]
            ampdf$adj_e[i] <- tempdf$s[j]
            svdf$amp_boundary[svdf$sv_id==tempdf$sv_id[j]] <- "boundary"
            bpdf$amp_boundary[bpdf$sv_id==tempdf$sv_id[j]] <- "boundary"
          }
        } else if (length(which(abs(tempdf$s - ampdf$e[i]) == min(abs(tempdf$s - ampdf$e[i])) & abs(tempdf$s - ampdf$e[i]) < svmargin)) > 1 & ("+" %in% tempdf$ori[which(abs(tempdf$s - ampdf$e[i]) == min(abs(tempdf$s - ampdf$e[i])))])){
          first_value <- 0
          for (j in which(abs(tempdf$s - ampdf$e[i]) == min(abs(tempdf$s - ampdf$e[i])) & tempdf$ori == "+" & abs(tempdf$s - ampdf$e[i]) < svmargin)){
            next_value <- tempdf$pe_support[j]
            if (first_value < next_value){
              first_value <- next_value
              imp_num <- j
            }
          }
          ampdf$sv_e[i] <- tempdf$sv_id[j]
          ampdf$adj_e[i] <- tempdf$s[j]
          svdf$amp_boundary[svdf$sv_id==tempdf$sv_id[j]] <- "boundary"
          bpdf$amp_boundary[bpdf$sv_id==tempdf$sv_id[j]] <- "boundary"
        }
      } 
    }
    
    colnames(bpdf)[2] <- "pos"
    bpdf <- bpdf[,c("chr", "pos", "ori", "svclass", "pe_support", "sv_id", "fbi", "amp_boundary", "sv_id2")]
    
    ampdf$num_bkpt = NA
    ampdf$num_fbi = NA
    ampdf$num_tra = NA
    ampdf$cn_peak = NA
    ampdf$cn_avg = NA
    
    # Calculate density
    bpdf$amp_id <- "no"
    for (j in 1:nrow(ampdf)){
      if (nrow(bpdf[bpdf$chr == ampdf$chr[j] & bpdf$pos >= ampdf$s[j] & bpdf$pos <= (ampdf$e[j]+1),]) > 0){
        ampdf$num_bkpt[j] <- nrow(bpdf[bpdf$chr == ampdf$chr[j] & bpdf$pos >= ampdf$s[j] & bpdf$pos <= (ampdf$e[j]+1),])
        ampdf$num_fbi[j] <- nrow(bpdf[bpdf$chr == ampdf$chr[j] & bpdf$pos >= ampdf$s[j] & bpdf$pos <= (ampdf$e[j]+1) & bpdf$fbi == "possible",])
        ampdf$num_tra[j] <- nrow(bpdf[bpdf$chr == ampdf$chr[j] & bpdf$pos >= ampdf$s[j] & bpdf$pos <= (ampdf$e[j]+1) & bpdf$svclass == "TRA",])
        bpdf[bpdf$chr == ampdf$chr[j] & bpdf$pos >= ampdf$s[j] & bpdf$pos <= (ampdf$e[j]+1),]$amp_id <- ampdf$amp_id[j]
      } else {
        ampdf$num_bkpt[j] <- nrow(bpdf[bpdf$chr == ampdf$chr[j] & bpdf$pos >= ampdf$s[j] & bpdf$pos <= (ampdf$e[j]+1),])
        ampdf$num_fbi[j] <- nrow(bpdf[bpdf$chr == ampdf$chr[j] & bpdf$pos >= ampdf$s[j] & bpdf$pos <= (ampdf$e[j]+1) & bpdf$fbi == "possible",])
        ampdf$num_tra[j] <- nrow(bpdf[bpdf$chr == ampdf$chr[j] & bpdf$pos >= ampdf$s[j] & bpdf$pos <= (ampdf$e[j]+1) & bpdf$svclass == "TRA",])
      }
      if (length(cnvdf$depth[cnvdf$chr == ampdf$chr[j] & cnvdf$s >= ampdf$s[j] & cnvdf$e <= ampdf$e[j] & !is.na(cnvdf$depth) & cnvdf$depth != 500]) == 0){
        ampdf$cn_peak[j] <- NA
        ampdf$cn_avg[j] <- NA
      } else {
        ampdf$cn_peak[j] <- max(cnvdf$depth[cnvdf$chr == ampdf$chr[j] & cnvdf$s >= ampdf$s[j] & cnvdf$e <= ampdf$e[j] & !is.na(cnvdf$depth) & cnvdf$depth != 500]) # newly added
        ampdf$cn_avg[j] <- weighted.mean(cnvdf$depth[cnvdf$chr == ampdf$chr[j] & cnvdf$s >= ampdf$s[j] & cnvdf$e <= ampdf$e[j] & !is.na(cnvdf$depth) & cnvdf$depth != 500], cnvdf$seglen[cnvdf$chr == ampdf$chr[j] & cnvdf$s >= ampdf$s[j] & cnvdf$e <= ampdf$e[j] & !is.na(cnvdf$depth) & cnvdf$depth != 500])  # newly added
      }
    }
    
    # Identify connections
    ampdf$connection = NA
    bpdf$connection <- "no"
    for (j in 1:nrow(bpdf)){
      bpdf$connection[j] <- bpdf$amp_id[bpdf$sv_id == bpdf$sv_id[j] & bpdf$sv_id2 != bpdf$sv_id2[j]]
    }
    
    for (j in 1:nrow(ampdf)){
      if (length(setdiff(unique(bpdf$connection[bpdf$amp_id == ampdf$amp_id[j]]), "no")) == 0){
        ampdf$connection[j] <- "no"
      } else {
        ampdf$connection[j] <- paste(as.character(sort(as.numeric(gsub("amplicon_", "", setdiff(unique(bpdf$connection[bpdf$amp_id == ampdf$amp_id[j]]), "no"))))), collapse = "-")
      }
    }
    
    
    # Per-boundary annotation
    ampdf$bd_s = "no"
    ampdf$binfo_s = "no"
    ampdf$bd_e = "no"
    ampdf$binfo_e = "no"
    
    svlist <- c(ampdf$sv_s, ampdf$sv_e)
    for (j in 1:nrow(ampdf)){
      if (ampdf$sv_s[j] != "no"){
        if (bpdf$amp_id[bpdf$sv_id==ampdf$sv_s[j] & bpdf$pos!=ampdf$adj_s[j]] == "no"){
          ampdf$bd_s[j] <- "e2o"
        } else {
          if (sum(svlist == ampdf$sv_s[j]) == 2){
            ampdf$bd_s[j] <- "e2e"
          } else {
            ampdf$bd_s[j] <- "e2i"
          }
        }
      }
      if (ampdf$sv_e[j] != "no"){
        if (bpdf$amp_id[bpdf$sv_id==ampdf$sv_e[j] & bpdf$pos!=ampdf$adj_e[j]] == "no"){
          ampdf$bd_e[j] <- "e2o"
        } else {
          if (sum(svlist == ampdf$sv_e[j]) == 2){
            ampdf$bd_e[j] <- "e2e"
          } else {
            ampdf$bd_e[j] <- "e2i"
          }
        }
      }
    }
    for (j in 1:nrow(ampdf)){
      if (ampdf$bd_s[j] == "e2e"){
        if (ampdf$sv_s[j] == ampdf$sv_e[j]){
          ampdf$binfo_s[j] <- "TD"
        } else if (unique(bpdf$fbi[bpdf$sv_id == ampdf$sv_s[j]]) == "possible"){
          ampdf$binfo_s[j] <- "FBI"
        } else if (unique(bpdf$svclass[bpdf$sv_id == ampdf$sv_s[j]]) == "TRA"){
          ampdf$binfo_s[j] <- "OA-TRA"
        } else {
          ampdf$binfo_s[j] <- "OA-INT"
        }
      } else if (ampdf$bd_s[j] == "e2i"){
        if (bpdf$amp_id[bpdf$sv_id==ampdf$sv_s[j] & bpdf$pos==ampdf$adj_s[j]] == bpdf$connection[bpdf$sv_id==ampdf$sv_s[j] & bpdf$pos==ampdf$adj_s[j]]){
          ampdf$binfo_s[j] <- "Self"
        } else if (bpdf$amp_id[bpdf$sv_id==ampdf$sv_s[j] & bpdf$pos==ampdf$adj_s[j]] != bpdf$connection[bpdf$sv_id==ampdf$sv_s[j] & bpdf$pos==ampdf$adj_s[j]] & bpdf$svclass[bpdf$sv_id==ampdf$sv_s[j] & bpdf$pos==ampdf$adj_s[j]] == "TRA"){
          ampdf$binfo_s[j] <- "OA-TRA"
        } else {
          ampdf$binfo_s[j] <- "OA-INT"
        }
      } else if (ampdf$bd_s[j] == "e2o"){
        if (bpdf$svclass[bpdf$sv_id==ampdf$sv_s[j] & bpdf$pos==ampdf$adj_s[j]] == "TRA"){
          if (bpdf$ori[bpdf$sv_id==ampdf$sv_s[j] & bpdf$pos!=ampdf$adj_s[j]] == "+"){
            imp_num <- which(bpdf$sv_id==ampdf$sv_s[j] & bpdf$pos!=ampdf$adj_s[j])
            if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-1) & cnvdf$e>=(bpdf$pos[imp_num]-1)]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-1) & cnvdf$e>=(bpdf$pos[imp_num]-1)]) && length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)])){
              if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)]) != 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)] == 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-1) & cnvdf$e>=(bpdf$pos[imp_num]-1)] != 0){
                ampdf$binfo_s[j] <- "LOH-TRA"
              } else {
                ampdf$binfo_s[j] <- "Unspec1-TRA" # Translocation to the amplified chromosome, but the CN pattern is not consistent with LOH
              }
            } else {
              ampdf$binfo_s[j] <- "Unspec2-TRA" # Translocation to the unamplified chromosome OR amplified chromosome but no CN info due to NA values OR telomeric rearrangement
            }
          } else { ## ori == "-"
            imp_num <- which(bpdf$sv_id==ampdf$sv_s[j] & bpdf$pos!=ampdf$adj_s[j])
            if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=bpdf$pos[imp_num] & cnvdf$e>=bpdf$pos[imp_num]]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=bpdf$pos[imp_num] & cnvdf$e>=bpdf$pos[imp_num]]) && length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)])){
              if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)]) != 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)] == 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=bpdf$pos[imp_num] & cnvdf$e>=bpdf$pos[imp_num]] != 0){
                ampdf$binfo_s[j] <- "LOH-TRA"
              } else {
                ampdf$binfo_s[j] <- "Unspec1-TRA" # Translocation to the amplified chromosome, but the CN pattern is not consistent with LOH
              }
            } else {
              ampdf$binfo_s[j] <- "Unspec2-TRA" # Translocation to the unamplified chromosome OR amplified chromosome but no CN info due to NA values OR telomeric rearrangement
            }
          }
        } else { ## Non-TRA
          if (bpdf$ori[bpdf$sv_id==ampdf$sv_s[j] & bpdf$pos!=ampdf$adj_s[j]] == "+"){
            imp_num <- which(bpdf$sv_id==ampdf$sv_s[j] & bpdf$pos!=ampdf$adj_s[j])
            if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-1) & cnvdf$e>=(bpdf$pos[imp_num]-1)]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-1) & cnvdf$e>=(bpdf$pos[imp_num]-1)]) && length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)])){
              if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)]) != 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)] == 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-1) & cnvdf$e>=(bpdf$pos[imp_num]-1)] != 0){
                ampdf$binfo_s[j] <- "LOH-INT"
              } else {
                ampdf$binfo_s[j] <- "Unspec1-INT" # Rearrangement to the amplified chromosome, but the CN pattern is not consistent with LOH
              }
            } else {
              ampdf$binfo_s[j] <- "Unspec2-INT" # Rearrangement to the unamplified chromosome OR amplified chromosome but no CN info due to NA values OR telomeric rearrangement
            }
          } else { ## ori == "-"
            imp_num <- which(bpdf$sv_id==ampdf$sv_s[j] & bpdf$pos!=ampdf$adj_s[j])
            if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=bpdf$pos[imp_num] & cnvdf$e>=bpdf$pos[imp_num]]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=bpdf$pos[imp_num] & cnvdf$e>=bpdf$pos[imp_num]]) && length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)])){
              if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)]) != 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)] == 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=bpdf$pos[imp_num] & cnvdf$e>=bpdf$pos[imp_num]] != 0){
                ampdf$binfo_s[j] <- "LOH-INT"
              } else {
                ampdf$binfo_s[j] <- "Unspec1-INT" # Rearrangement to the amplified chromosome, but the CN pattern is not consistent with LOH
              }
            } else {
              ampdf$binfo_s[j] <- "Unspec2-INT" # Rearrangement to the unamplified chromosome OR amplified chromosome but no CN info due to NA values OR telomeric rearrangement
            }
          }
        }
      }
      
      if (ampdf$bd_e[j] == "e2e"){
        if (ampdf$sv_s[j] == ampdf$sv_e[j]){
          ampdf$binfo_e[j] <- "TD"
        } else if (unique(bpdf$fbi[bpdf$sv_id == ampdf$sv_e[j]]) == "possible"){
          ampdf$binfo_e[j] <- "FBI"
        } else if (unique(bpdf$svclass[bpdf$sv_id == ampdf$sv_e[j]]) == "TRA"){
          ampdf$binfo_e[j] <- "OA-TRA"
        } else {
          ampdf$binfo_e[j] <- "OA-INT"
        }
      } else if (ampdf$bd_e[j] == "e2i"){
        if (bpdf$amp_id[bpdf$sv_id==ampdf$sv_e[j] & bpdf$pos==ampdf$adj_e[j]] == bpdf$connection[bpdf$sv_id==ampdf$sv_e[j] & bpdf$pos==ampdf$adj_e[j]]){
          ampdf$binfo_e[j] <- "Self"
        } else if (bpdf$amp_id[bpdf$sv_id==ampdf$sv_e[j] & bpdf$pos==ampdf$adj_e[j]] != bpdf$connection[bpdf$sv_id==ampdf$sv_e[j] & bpdf$pos==ampdf$adj_e[j]] & bpdf$svclass[bpdf$sv_id==ampdf$sv_e[j] & bpdf$pos==ampdf$adj_e[j]] == "TRA"){
          ampdf$binfo_e[j] <- "OA-TRA"
        } else {
          ampdf$binfo_e[j] <- "OA-INT"
        }
      } else if (ampdf$bd_e[j] == "e2o"){
        if (bpdf$svclass[bpdf$sv_id==ampdf$sv_e[j] & bpdf$pos==ampdf$adj_e[j]] == "TRA"){
          if (bpdf$ori[bpdf$sv_id==ampdf$sv_e[j] & bpdf$pos!=ampdf$adj_e[j]] == "+"){
            imp_num <- which(bpdf$sv_id==ampdf$sv_e[j] & bpdf$pos!=ampdf$adj_e[j])
            if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-1) & cnvdf$e>=(bpdf$pos[imp_num]-1)]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-1) & cnvdf$e>=(bpdf$pos[imp_num]-1)]) && length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)])){
              if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)]) != 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)] == 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-1) & cnvdf$e>=(bpdf$pos[imp_num]-1)] != 0){
                ampdf$binfo_e[j] <- "LOH-TRA"
              } else {
                ampdf$binfo_e[j] <- "Unspec1-TRA" # Translocation to the amplified chromosome, but the CN pattern is not consistent with LOH
              }
            } else {
              ampdf$binfo_e[j] <- "Unspec2-TRA" # Translocation to the unamplified chromosome OR amplified chromosome but no CN info due to NA values OR telomeric rearrangement
            }
          } else { ## ori == "-"
            imp_num <- which(bpdf$sv_id==ampdf$sv_e[j] & bpdf$pos!=ampdf$adj_e[j])
            if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=bpdf$pos[imp_num] & cnvdf$e>=bpdf$pos[imp_num]]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=bpdf$pos[imp_num] & cnvdf$e>=bpdf$pos[imp_num]]) && length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)])){
              if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)]) != 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)] == 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=bpdf$pos[imp_num] & cnvdf$e>=bpdf$pos[imp_num]] != 0){
                ampdf$binfo_e[j] <- "LOH-TRA"
              } else {
                ampdf$binfo_e[j] <- "Unspec1-TRA" # Translocation to the amplified chromosome, but the CN pattern is not consistent with LOH
              }
            } else {
              ampdf$binfo_e[j] <- "Unspec2-TRA" # Translocation to the unamplified chromosome OR amplified chromosome but no CN info due to NA values OR telomeric rearrangement
            }
          }
        } else { ## Non-TRA
          if (bpdf$ori[bpdf$sv_id==ampdf$sv_e[j] & bpdf$pos!=ampdf$adj_e[j]] == "+"){
            imp_num <- which(bpdf$sv_id==ampdf$sv_e[j] & bpdf$pos!=ampdf$adj_e[j])
            if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-1) & cnvdf$e>=(bpdf$pos[imp_num]-1)]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-1) & cnvdf$e>=(bpdf$pos[imp_num]-1)]) && length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)])){
              if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)]) != 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]+1) & cnvdf$e>=(bpdf$pos[imp_num]+1)] == 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-1) & cnvdf$e>=(bpdf$pos[imp_num]-1)] != 0){
                ampdf$binfo_e[j] <- "LOH-INT"
              } else {
                ampdf$binfo_e[j] <- "Unspec1-INT" # Rearrangement to the amplified chromosome, but the CN pattern is not consistent with LOH
              }
            } else {
              ampdf$binfo_e[j] <- "Unspec2-INT" # Rearrangement to the unamplified chromosome OR amplified chromosome but no CN info due to NA values OR telomeric rearrangement
            }
          } else { ## ori == "-"
            imp_num <- which(bpdf$sv_id==ampdf$sv_e[j] & bpdf$pos!=ampdf$adj_e[j])
            if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=bpdf$pos[imp_num] & cnvdf$e>=bpdf$pos[imp_num]]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=bpdf$pos[imp_num] & cnvdf$e>=bpdf$pos[imp_num]]) && length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)]) != 0 && !is.na(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)])){
              if (length(cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)]) != 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=(bpdf$pos[imp_num]-2) & cnvdf$e>=(bpdf$pos[imp_num]-2)] == 0 & cnvdf$alt_cn[cnvdf$chr==bpdf$chr[imp_num] & cnvdf$s<=bpdf$pos[imp_num] & cnvdf$e>=bpdf$pos[imp_num]] != 0){
                ampdf$binfo_e[j] <- "LOH-INT"
              } else {
                ampdf$binfo_e[j] <- "Unspec1-INT" # Rearrangement to the amplified chromosome, but the CN pattern is not consistent with LOH
              }
            } else {
              ampdf$binfo_e[j] <- "Unspec2-INT" # Rearrangement to the unamplified chromosome OR amplified chromosome but no CN info due to NA values OR telomeric rearrangement
            }
          }
        }
      }
      
    }
    for (j in 1:nrow(ampdf)){
      if (ampdf$sv_s[j] != "no"){
        if (unique(bpdf$fbi[bpdf$sv_id == ampdf$sv_s[j]]) == "possible"){
          ampdf$bd_s[j] <- "e2e"
          ampdf$binfo_s[j] <- "FBI"
        }
      }
      if (ampdf$sv_e[j] != "no"){
        if (unique(bpdf$fbi[bpdf$sv_id == ampdf$sv_e[j]]) == "possible"){
          ampdf$bd_e[j] <- "e2e"
          ampdf$binfo_e[j] <- "FBI"
        }
      }
    }
    ampdf$den_bkpt <- ampdf$num_bkpt*1000000/ampdf$seglen
    
    bpdf$balanced = "no"
    svdf$balanced = "no"
    for (j in 1:(nrow(bpdf)-1)){
      if (bpdf$chr[j] == bpdf$chr[j+1] & (bpdf$pos[j+1] - bpdf$pos[j]) < 1000 & bpdf$ori[j] == "+" & bpdf$ori[j+1] == "-" & bpdf$balanced[j] == "no"){
        bpdf$balanced[j] <- "yes"
        bpdf$balanced[j+1] <- "yes"
        svdf$balanced[svdf$sv_id == bpdf$sv_id[j]] <- "yes"
        svdf$balanced[svdf$sv_id == bpdf$sv_id[j+1]] <- "yes"
      }
    }
    
    bpdf$interd <- NA
    for (j in 1:(nrow(bpdf)-1)){
      if (bpdf$chr[j] == bpdf$chr[j+1]){
        bpdf$interd[j] <- bpdf$pos[j+1] - bpdf$pos[j]
      } 
    }
    
    write.table(svdf, paste0(subDir, "/", tumortype, "/SV_amp_boundaries/", tumorname, ".amp_boundaries.txt"), row.names=F, quote=F, sep="\t", col.names=T)
    write.table(bpdf, paste0(subDir, "/", tumortype, "/SV_breakpoints_long/", tumorname, ".breakpoints_long.txt"), row.names=F, quote=F, sep="\t", col.names=T)
    cadf <- rbind(cadf, ampdf)
  }
}

# Size filter (removing amplified segments larger than 30 Mb, 5 out of 1188, 0.4% was removed)
cadf <- cadf[cadf$seglen < 30000000,]

write.table(cadf, paste0(subDir, "/", tumortype, "/ampdf.", tumortype, ".txt"), row.names=F, quote=F, sep="\t", col.names=T)
write.table(df, paste0(subDir, "/", tumortype, "/", tumortype, ".df.txt"), row.names=F, quote=F, sep="\t", col.names=T)


