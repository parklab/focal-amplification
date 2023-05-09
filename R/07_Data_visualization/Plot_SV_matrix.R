library(dplyr)
library(grid)
library(fields)
library(stringr)
library(plotrix)
library(mapplots)
library(viridis)
library(reshape2)
library(viridis)
library(RColorBrewer)

# Plotting SV matrix (Figure 1)
## Creating a dataframe with all genomic bins
bin_size = 5000000 # 5MB
bindf = NULL
options(scipen=999)
chromosome <- hg19_coord$chr[1]
num_bins <- ceiling(hg19_coord$length[1]/bin_size)
boundaries_start <- NA
boundaries_end <- NA
for (j in 1:num_bins){
  boundaries_start[j] <- 1 + bin_size*(j-1)
  boundaries_end[j] <- bin_size*j
}
df <- data.frame(rep(chromosome, num_bins), boundaries_start, boundaries_end, stringsAsFactors = F)
colnames(df) <- c("chr", "start", "end")
bindf <- df
for (i in 2:nrow(hg19_coord)){
  chromosome <- hg19_coord$chr[i]
  num_bins <- ceiling(hg19_coord$length[i]/bin_size)
  boundaries_start <- NA
  boundaries_end <- NA
  for (j in 1:num_bins){
    boundaries_start[j] <- 1 + bin_size*(j-1)
    boundaries_end[j] <- bin_size*j
  }
  df <- data.frame(rep(chromosome, num_bins), boundaries_start, boundaries_end, stringsAsFactors = F)
  colnames(df) <- c("chr", "start", "end")
  bindf <- rbind(bindf, df)
}
write.csv(bindf, paste0("genomicbins_GRCh37_decoy_5mb.csv"), row.names=F) # now 10kb (way too small), 100kb, 1mb, 5mb, and 10mb (used in the initial submission) available

## SV matrix based on breast cancer SVs (Fig1A)
mb1 <- read.csv("genomicbins_GRCh37_decoy_5mb.csv", header=T, as.is=T)
mb1$id <- c(1:nrow(mb1))
cmtx <- matrix(0, nrow(mb1), nrow(mb1))

df <- purdf[!is.na(purdf$ct_sv_hmf) & purdf$availability=="yes",]
df <- df[df$ct_sv_hmf != 0,]
for (i in 1:nrow(df)){
  print(paste0("Processing ", df$study_id[i]))
  svdf <- read.csv(paste0(df$study_id[i], ".purple.sv.vcf.gz.bedpe"), header=T, as.is=T, sep="\t")
  svdf$ch1_num <- svdf$chrom1
  svdf$ch1_num[svdf$ch1_num == "X"] <- "23"
  svdf$ch1_num[svdf$ch1_num == "Y"] <- "24"
  svdf$ch1_num[svdf$ch1_num == "MT"] <- "25"
  svdf$ch1_num <- as.numeric(svdf$ch1_num)
  svdf$ch2_num <- svdf$chrom2
  svdf$ch2_num[svdf$ch2_num == "X"] <- "23"
  svdf$ch2_num[svdf$ch2_num == "Y"] <- "24"
  svdf$ch2_num[svdf$ch2_num == "MT"] <- "25"
  svdf$ch2_num <- as.numeric(svdf$ch2_num)
  svdf$correct <- ""
  svdf$loc <- ""
  for (j in 1:nrow(svdf)){
    if (svdf$ch1_num[j] == svdf$ch2_num[j] && svdf$start1[j] <= svdf$start2[j]){
      svdf$correct[j] <- "yes"
      svdf$loc[j] <- paste0(mb1$id[mb1$chr == svdf$chrom1[j] & mb1$start <= (svdf$start1[j]+svdf$end1[j])/2 & mb1$end >= (svdf$start1[j]+svdf$end1[j])/2], ";", mb1$id[mb1$chr == svdf$chrom2[j] & mb1$start <= (svdf$start2[j]+svdf$end2[j])/2 & mb1$end >= (svdf$start2[j]+svdf$end2[j])/2])
    } else if (svdf$ch1_num[j] < svdf$ch2_num[j]){
      svdf$correct[j] <- "yes"
      svdf$loc[j] <- paste0(mb1$id[mb1$chr == svdf$chrom1[j] & mb1$start <= (svdf$start1[j]+svdf$end1[j])/2 & mb1$end >= (svdf$start1[j]+svdf$end1[j])/2], ";", mb1$id[mb1$chr == svdf$chrom2[j] & mb1$start <= (svdf$start2[j]+svdf$end2[j])/2 & mb1$end >= (svdf$start2[j]+svdf$end2[j])/2])
    } else {
      svdf$correct[j] <- "no"
      svdf$loc[j] <- paste0(mb1$id[mb1$chr == svdf$chrom2[j] & mb1$start <= (svdf$start2[j]+svdf$end2[j])/2 & mb1$end >= (svdf$start2[j]+svdf$end2[j])/2], ";", mb1$id[mb1$chr == svdf$chrom1[j] & mb1$start <= (svdf$start1[j]+svdf$end1[j])/2 & mb1$end >= (svdf$start1[j]+svdf$end1[j])/2])
    }
  }
  loclist <- unique(svdf$loc)
  for (j in 1:length(loclist)){
    num1 <- as.numeric(strsplit(loclist[j], ";", fixed=T)[[1]][1])
    num2 <- as.numeric(strsplit(loclist[j], ";", fixed=T)[[1]][2])
    if (num1 >= num2){
      cmtx[num2, num1] <- cmtx[num2, num1] + 1
    } else {
      cmtx[num1, num2] <- cmtx[num1, num2] + 1
    }
  }
}

range(cmtx)
table(cmtx)
cutlevel <- 7 ## Truncation of the values for visualization.

temporary_matrix <- cmtx
temporary_matrix[temporary_matrix > cutlevel] <- cutlevel
brewer.pal(n=9, "Blues")

pdf(paste0("heatmap.5mb.allbreast.", nrow(df), ".cut.", cutlevel, ".final.pdf"), height=9, width=9)
heatmap(temporary_matrix, Rowv=NA, Colv=NA, col = brewer.pal(9, "Blues"), scale = "none", xlab = "", ylab = "", labRow = FALSE, labCol = FALSE)
dev.off()

pdf(paste0("heatmap.5mb.allbreast.", nrow(df), ".cut", cutlevel, ".final.legend.pdf"), height=4, width=4)
image.plot(temporary_matrix, col=brewer.pal(9, "Blues"), xaxt = "n", yaxt = "n", legend.only = T)
dev.off()


