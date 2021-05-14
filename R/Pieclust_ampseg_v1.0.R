library(dplyr)
library(grid)
library(fields)
library(stringr)
library(plotrix)
library(mapplots)
library(viridis)
library(reshape2)
options(scipen=999)

mainDir <- getwd()
subDir <- "Figures"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
ifelse(!dir.exists(file.path(mainDir, subDir, tumortype)), dir.create(file.path(mainDir, subDir, tumortype)), FALSE)
ifelse(!dir.exists(file.path(mainDir, subDir, tumortype, "Bychr_copynumber")), dir.create(file.path(mainDir, subDir, tumortype, "Bychr_copynumber")), FALSE)

# Merging ampdf files
bdcolors <- c(rgb(255/255, 153/255, 18/255), rgb(39/255, 64/255, 139/255), rgb(205/255, 50/255, 120/255), rgb(200/255, 200/255, 200/255))
classcolors <- c(rgb(61/255, 145/255, 64/255), rgb(220/255, 20/255, 60/255), rgb(139/255, 28/255, 98/255), rgb(85/255, 26/255, 139/255), rgb(193/255, 205/255, 205/255), rgb(90/255, 90/255, 90/255))

sumdf <- read.csv("Summaryinfo.table.1.19.txt", sep="\t", as.is=T, header=T)
df <- read.csv("BreastCancer278.Ampdf.v2.full.txt", header=T, as.is=T, sep="\t")
ampdf <- as.data.frame(matrix(NA, ncol=22, nrow=0))
colnames(ampdf) <- c(colnames(df), "histology_abbreviation")

histolist <- list.dirs(path = "Result_1mb", recursive = FALSE)
histolist <- gsub("Result_1mb/", "", histolist)
histolist <- setdiff(histolist, c("Bone-Benign", "Myeloid-MDS"))
for (i in 1:length(histolist)){
  if (histolist[i] != "Breast-AdenoCA"){
    df <- read.csv(paste0("Result_1mb/", histolist[i], "/ampdf.", histolist[i], ".txt"), header=T, as.is=T, sep="\t")
  } else {
    df <- read.csv("BreastCancer278.Ampdf.v2.full.txt", header=T, as.is=T, sep="\t")
  }
  if (nrow(df) > 0){
    df$histology_abbreviation <- histolist[i]
    ampdf <- rbind(ampdf, df)
  }
}

ampdf <- ampdf[,c(1,ncol(ampdf),c(2:(ncol(ampdf)-1)))]
for (j in 1:nrow(ampdf)){
  ampdf$number[j] <- strsplit(ampdf$amp_id[j], "_", fixed=T)[[1]][2]
}
for (j in 1:nrow(ampdf)){
  if (ampdf$connection[j] == ampdf$number[j]){
    ampdf$n_conn[j] <- 0
  } else if (ampdf$connection[j] == "no"){
    ampdf$n_conn[j] <- 0
  } else {
    ampdf$n_conn[j] <- length(strsplit(ampdf$connection[j], "-", fixed=T)[[1]])
  }
}

ampdf$class[ampdf$binfo_s == "TD" & ampdf$binfo_e == "TD" & ampdf$sv_s == ampdf$sv_e] <- "TD"
ampdf$class[ampdf$binfo_s == "FBI" | ampdf$binfo_e == "FBI"] <- "FBI"
ampdf$class[ampdf$binfo_s %in% c("LOH-TRA", "OA-TRA", "Unspec1-TRA", "Unspec2-TRA") | ampdf$binfo_e %in% c("LOH-TRA", "OA-TRA", "Unspec1-TRA", "Unspec2-TRA")] <- "TRA"
ampdf$class[ampdf$binfo_s == "FBI" & ampdf$binfo_e %in% c("LOH-TRA", "OA-TRA", "Unspec1-TRA", "Unspec2-TRA")] <- "FBI+TRA"
ampdf$class[ampdf$binfo_s %in% c("LOH-TRA", "OA-TRA", "Unspec1-TRA", "Unspec2-TRA") & ampdf$binfo_e == "FBI"] <- "FBI+TRA"
ampdf$class[ampdf$binfo_s == "no" & ampdf$binfo_e == "no"] <- "NO"
ampdf$class[is.na(ampdf$class)] <- "INT"
ampdf$class <- factor(ampdf$class, levels = c("TD", "FBI", "FBI+TRA", "TRA", "INT", "NO"))

ampdf <- ampdf[!is.na(ampdf$cn_avg),]
write.table(ampdf, "ampdf.focal.amplifications.entire.breast.and.pcawg.txt", row.names=F, quote=F, col.names=T, sep="\t")

# Simple pie charts for amplicon classes and boundary types
pdf("Figures/pie_boundary_types.pdf", height=5, width=10)
par(mfrow=c(1,2))
pie(table(ampdf$bd_s), col=bdcolors, border=NA, main = "Left border")
pie(table(ampdf$bd_e), col=bdcolors, border=NA, main = "Right border")
dev.off()

pdf("Figures/pie_ampseg_classes.pdf", height=5, width=5)
pie(table(ampdf$class), col=classcolors, border=NA, main = "Classes")
dev.off()

hfreq <- as.data.frame(matrix(NA, ncol=2, nrow=length(histolist)))
colnames(hfreq) <- c("histology_abbreviation", "n_tumor")
hfreq$histology_abbreviation <- histolist
for (i in 1:nrow(hfreq)){
  hfreq$n_tumor[i] <- length(sumdf$icgc_donor_id[sumdf$SV.events != 0 & !is.na(sumdf$SV.events) & sumdf$histology_abbreviation == hfreq$histology_abbreviation[i]])
}
hfreq$n_tumor[hfreq$histology_abbreviation=="Breast-AdenoCA"] <- 278

hfreq$num <- 0
for (i in 1:nrow(hfreq)){
  hfreq$num[i] <- nrow(ampdf[ampdf$histology_abbreviation == hfreq$histology_abbreviation[i],])
}

# Ampseg information by tumor type
histolist <- hfreq$histology_abbreviation[hfreq$num > 30]
typdf <- as.data.frame(matrix(NA, ncol=8, nrow=length(histolist)))
colnames(typdf) <- c("histology_abbreviation", "n_cgr", "n_td", "n_fbi", "n_fbi_tra", "n_tra", "n_int", "n_nosv")
typdf$histology_abbreviation <- histolist
for (i in 1:nrow(typdf)){
  typdf$n_cgr[i] <- nrow(ampdf[ampdf$histology_abbreviation == typdf$histology_abbreviation[i],])
  typdf$n_td[i] <- nrow(ampdf[ampdf$histology_abbreviation == typdf$histology_abbreviation[i] & ampdf$class == "TD",])
  typdf$n_fbi[i] <- nrow(ampdf[ampdf$histology_abbreviation == typdf$histology_abbreviation[i] & ampdf$class == "FBI",])
  typdf$n_fbi_tra[i] <- nrow(ampdf[ampdf$histology_abbreviation == typdf$histology_abbreviation[i] & ampdf$class == "FBI+TRA",])
  typdf$n_tra[i] <- nrow(ampdf[ampdf$histology_abbreviation == typdf$histology_abbreviation[i] & ampdf$class == "TRA",])
  typdf$n_int[i] <- nrow(ampdf[ampdf$histology_abbreviation == typdf$histology_abbreviation[i] & ampdf$class == "INT",])
  typdf$n_nosv[i] <- nrow(ampdf[ampdf$histology_abbreviation == typdf$histology_abbreviation[i] & ampdf$class == "NO",])
}

typdf$f_td <- typdf$n_td/typdf$n_cgr
typdf$f_fbi <- typdf$n_fbi/typdf$n_cgr
typdf$f_fbi_tra <- typdf$n_fbi_tra/typdf$n_cgr
typdf$f_tra <- typdf$n_tra/typdf$n_cgr
typdf$f_int <- typdf$n_int/typdf$n_cgr
typdf$f_nosv <- typdf$n_nosv/typdf$n_cgr

#data <- typdf[,c(9:14)]
data <- typdf[,c(9:13)]
row.names(data) <- typdf$histology_abbreviation
d <- dist(data, method="euclidean")
hc <- hclust(d, method = "complete")

pdf("Figures/hclust_dendrogram.pdf", height=5, width=5)
plot(hc, hang = -1, xlab = "Histologic types")
rect.hclust(hc, k = 3)
dev.off()


# Pie graph for FBI-TRA classes
histolist <- hfreq$histology_abbreviation[hfreq$num > 30]
histolist <- c("Cervix-SCC", "Head-SCC", "Lung-SCC", "Panc-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA", "Uterus-AdenoCA", "Liver-HCC", "Biliary-AdenoCA", "ColoRect-AdenoCA", "Bladder-TCC", "Eso-AdenoCA", "Lymph-BNHL", "CNS-GBM", "CNS-Medullo", "SoftTissue-Liposarc", "Bone-Osteosarc", "Lung-AdenoCA", "Stomach-AdenoCA", "Skin-Melanoma", "Breast-AdenoCA") # List ordered by dendrogram
numberlist <- NA
for (i in 1:length(histolist)){
  numberlist[i] <- hfreq$n_tumor[hfreq$histology_abbreviation == histolist[i]]
}
cgrdf <- as.data.frame(matrix(NA, ncol=11, nrow=24*21)) # 21 --> 22
colnames(cgrdf) <- c("histology_abbreviation", "xvalue", "yvalue", "chrom", "n_cgr", "n_td", "n_fbi", "n_fbi_tra", "n_tra", "n_int", "n_nosv")
for (i in 1:length(histolist)){
  cgrdf$histology_abbreviation[(i*24-23):(i*24)] <- histolist[i]
  cgrdf$xvalue[(i*24-23):(i*24)] <- c(1:24)
  cgrdf$yvalue[(i*24-23):(i*24)] <- 23-i
  cgrdf$chrom[(i*24-23):(i*24)] <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
  cgrdf$n_cgr[(i*24-23):(i*24)] <- rep(0,24)
  for (j in (i*24-23):(i*24)){
    cgrdf$n_cgr[j] <- length(ampdf$amp_id[ampdf$chr == cgrdf$chrom[j] & ampdf$histology_abbreviation == histolist[i]])
    cgrdf$n_td[j] <- length(ampdf$amp_id[ampdf$chr == cgrdf$chrom[j] & ampdf$histology_abbreviation == histolist[i] & ampdf$class == "TD"])
    cgrdf$n_fbi[j] <- length(ampdf$amp_id[ampdf$chr == cgrdf$chrom[j] & ampdf$histology_abbreviation == histolist[i] & ampdf$class == "FBI"])
    cgrdf$n_fbi_tra[j] <- length(ampdf$amp_id[ampdf$chr == cgrdf$chrom[j] & ampdf$histology_abbreviation == histolist[i] & ampdf$class == "FBI+TRA"])
    cgrdf$n_tra[j] <- length(ampdf$amp_id[ampdf$chr == cgrdf$chrom[j] & ampdf$histology_abbreviation == histolist[i] & ampdf$class == "TRA"])
    cgrdf$n_int[j] <- length(ampdf$amp_id[ampdf$chr == cgrdf$chrom[j] & ampdf$histology_abbreviation == histolist[i] & ampdf$class == "INT"])
    cgrdf$n_nosv[j] <- length(ampdf$amp_id[ampdf$chr == cgrdf$chrom[j] & ampdf$histology_abbreviation == histolist[i] & ampdf$class == "NO"])
  }
}

cgr_short <- cgrdf[,c("xvalue", "yvalue", "n_td", "n_fbi", "n_fbi_tra", "n_tra", "n_int", "n_nosv")]
cgr_melt <- melt(cgr_short, id.vars = c("xvalue", "yvalue"))
xyz <- make.xyz(cgr_melt$xvalue, cgr_melt$yvalue, cgr_melt$value, cgr_melt$variable)

pdf("Figures/ampseg_boundaries_bychrom_tumortype_class.pdf", height=8, width=12)
par(mar=c(4.1, 12.1, 4.1, 2.1))
plot(NA, xlim = c(0,24), ylim = c(-1,22), frame = F, xlab = "Chromosome", yaxt='n', xaxt='n', ylab='', main = paste0("Number of amplified segments (total n = ", nrow(ampdf), ")"))
axis(2, at=c(1:23), las=1, labels=c("", paste0(rev(histolist), " (", rev(numberlist),")"), ""))
axis(1, at=c(0:24), las=1, labels=c("", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))
for (i in 1:22){
  abline(h = i, col=rgb(232/255, 232/255, 232/255))
}
for (i in 1:24){
  abline(v = i, col=rgb(232/255, 232/255, 232/255))
}
draw.pie(xyz$x, xyz$y, xyz$z, radius = 1.1, col=classcolors, border = NA)
legend.pie(19, -0.6, labels=c("Tandem duplication", "Fold-back inversion", "FBI+TRA", "Interchromosomal", "Other Intrachromosomal", "No SV support"), radius = 0.6, bty="n", col=classcolors, border=NA, label.dist=1.3)
legend.z <- round(max(rowSums(xyz$z,na.rm=TRUE)))
legend.bubble(12, -0.6, z=legend.z, round=0, maxradius=1.1, bty="n", txt.cex=0.6)
dev.off()


