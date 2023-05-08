library(dplyr)
library(grid)
library(fields)
library(stringr)
library(mapplots)
library(reshape2)

# GSEA figures
df <- read.csv("GSEA.report_for_na_pos_1650575918912.tsv", header=T, as.is=T, sep="\t") # GSEA output file. Available in the Data folder.
df <- df[c(1:8),]

pdf("FIP_revision/dotplot.gsea.q.values.updated.pdf", height=4.5, width=4.5)
par(mar=c(5.1, 12.1, 4.1, 2.1))
dfx = data.frame(ev1=-log10(df$FDR.q.val), ev2=c(8:1), ev3=pi*(df$ES)^2)
with(dfx, symbols(x=ev1, y=ev2, circles=ev3, inches=1/7, ann=F, bg="#27408B", fg=NULL, xaxt = 'n', yaxt = 'n', frame=F))
axis(1)
axis(2, at=c(1:8), labels=rev(gsub("_", " ", tolower(gsub("HALLMARK_", "", df$NAME)))), las=1)
dev.off()

# Number of junctions in each library
libdf <- as.data.frame(matrix(NA, ncol = 5, nrow = 24))
colnames(libdf) <- c("cell", "cut", "e2", "library", "ct_junctions")
libdf$cut <- c(rep("SHANK2", 12), rep("RARA", 12))
libdf$cell <- c(rep("MCF7", 6), rep("T47D", 6), rep("MCF7", 6), rep("T47D", 6))
libdf$e2 <- c(rep("e2_treated", 3), rep("control", 3), rep("e2_treated", 3), rep("control", 3), rep("e2_treated", 3), rep("control", 3), rep("e2_treated", 3), rep("control", 3))
libdf$library <- c(rep(c(1, 2, 3), 8))
tempdf <- read.csv("../28_HTGTS/library_junctions", header=T, as.is=T, sep="\t")
libdf$name <- tempdf$sample
libdf$ct_junctions <- tempdf$total
libdf$ct_unique <- tempdf$unique
libdf$order <- c(2,2,2,1,1,1,4,4,4,3,3,3,6,6,6,5,5,5,8,8,8,7,7,7)

pdf("FIP_revision/boxplot.htgts.library.junction.count.pdf", height=4, width=6)
boxplot(libdf$ct_unique ~ libdf$order, frame=F, las=1, xaxt='n', ylim=c(0,25000), xlab="", ylab="# unique translocation junctions", col=rep(c("#FFF0F5", "#CD3278"), 4))
stripchart(libdf$ct_unique ~ libdf$order, pch=19, col=rgb(0,0,0,.2), vertical=T, add=T, method="jitter")
axis(1, at = c(1:8), labels = c(rep("", 8)))
dev.off()

t.test(libdf$ct_unique[libdf$cell == "MCF7" & libdf$cut == "SHANK2"] ~ libdf$e2[libdf$cell == "MCF7" & libdf$cut == "SHANK2"])
t.test(libdf$ct_unique[libdf$cell == "MCF7" & libdf$cut == "RARA"] ~ libdf$e2[libdf$cell == "MCF7" & libdf$cut == "RARA"])
t.test(libdf$ct_unique[libdf$cell == "T47D" & libdf$cut == "SHANK2"] ~ libdf$e2[libdf$cell == "T47D" & libdf$cut == "SHANK2"])
t.test(libdf$ct_unique[libdf$cell == "T47D" & libdf$cut == "RARA"] ~ libdf$e2[libdf$cell == "T47D" & libdf$cut == "RARA"])


# Annotation of HTGTS junctions
## Collecting all the junctions into a single file
library(stringr)
filelist <- list.files("bedfiles/", pattern = "_withoutE2.hg19.bed", full.names = F)
for (w in filelist){
  df <- read.csv(paste0("bedfiles/", w), header=F, as.is=T, sep="\t")
  df$chr <- gsub("chr", "", df$V1)
  df$start <- df$V2
  df <- df[df$chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),]
  df$chr_num <-df$chr
  df$chr_num[df$chr_num == "X"] <- "23"
  df$chr_num[df$chr_num == "Y"] <- "24"
  df$chr_num <- as.numeric(as.character(df$chr_num))
  df <- df[order(df$start, decreasing = F),]
  df <- df[order(df$chr_num, decreasing = F),]
  df <- df[,colnames(df) != "chr_num"]
  df$ori <- df$V6
  df$inclusion <- "yes"
  if (strsplit(w, "_", fixed = T)[[1]][2] == "hRARA"){
    df$inclusion[df$chr=="17" & df$start <= 38491451 + 1000000 & df$start >= 38491451 - 1000000] <- "no"
  } else if (strsplit(w, "_", fixed = T)[[1]][2] == "hSHANK2"){
    df$inclusion[df$chr=="11" & df$start <= 70658290 + 1000000 & df$start >= 70658290 - 1000000] <- "no"
  }
  df <- df[df$inclusion == "yes",]
  df$end <- df$start
  df$ref <- "-"
  df$alt <- "A"
  df <- df[,c("chr", "start", "end", "ref", "alt")]
  if (w == filelist[1]){
    tempdf <- df
  } else {
    tempdf <- rbind(tempdf, df)
  }
}
df <- tempdf
df$duplicate <- "no"
for (i in 2:nrow(df)){
  if (df$start[i] == df$start[i-1]){
    df$duplicate[i] <- "yes"
  }
}
df <- df[df$duplicate == "no"]
write.table(df[,c("chr", "start", "end", "ref", "alt")], "htgts.junctions.withoutE2.avinput", row.names = F, col.names = F, quote = F, sep = "\t")

## Annotated files for general statistics
ctdf <- read.csv("../28_HTGTS/htgts.junctions.withoutE2.avinput.bkpt.anv.hg19_multianno.txt", header=T, as.is=T, sep="\t")
e2df <- read.csv("../28_HTGTS/htgts.junctions.withE2.avinput.bkpt.anv.hg19_multianno.txt", header=T, as.is=T, sep="\t")

table(ctdf$Func.refGene)*100/nrow(ctdf)
table(e2df$Func.refGene)*100/nrow(e2df)
sum(genedf$end-genedf$start)
sum(hg19_coord$length)

nrow(ctdf[ctdf$Func.refGene %in% c("intergenic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing", "ncRNA_UTR5"),]) # intergenic region, control
nrow(ctdf) - nrow(ctdf[ctdf$Func.refGene %in% c("intergenic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing", "ncRNA_UTR5"),]) # genic region, control
nrow(e2df[e2df$Func.refGene %in% c("intergenic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing", "ncRNA_UTR5"),]) # intergenic region, e2
nrow(e2df) - nrow(e2df[e2df$Func.refGene %in% c("intergenic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing", "ncRNA_UTR5"),]) # genic region, e2

## Density ratio
((nrow(ctdf) - nrow(ctdf[ctdf$Func.refGene %in% c("intergenic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing", "ncRNA_UTR5"),]))/sum(genedf$end-genedf$start))/(nrow(ctdf[ctdf$Func.refGene %in% c("intergenic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing", "ncRNA_UTR5"),])/(sum(hg19_coord$length)-sum(genedf$end-genedf$start))) # control
nrow(ctdf[ctdf$Func.refGene %in% c("intergenic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing", "ncRNA_UTR5"),])/(sum(hg19_coord$length)-sum(genedf$end-genedf$start))
((nrow(e2df) - nrow(e2df[e2df$Func.refGene %in% c("intergenic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing", "ncRNA_UTR5"),]))/sum(genedf$end-genedf$start))/(nrow(e2df[e2df$Func.refGene %in% c("intergenic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing", "ncRNA_UTR5"),])/(sum(hg19_coord$length)-sum(genedf$end-genedf$start)))

## Fisher's test for control
fisher.test(matrix(c(nrow(ctdf) - nrow(ctdf[ctdf$Func.refGene %in% c("intergenic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing", "ncRNA_UTR5"),]), nrow(ctdf[ctdf$Func.refGene %in% c("intergenic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing", "ncRNA_UTR5"),]), sum(genedf$end-genedf$start), (sum(hg19_coord$length)-sum(genedf$end-genedf$start))), nrow = 2, ncol = 2))
## Fisher's test for e2
fisher.test(matrix(c(nrow(e2df) - nrow(e2df[e2df$Func.refGene %in% c("intergenic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing", "ncRNA_UTR5"),]), nrow(e2df[e2df$Func.refGene %in% c("intergenic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing", "ncRNA_UTR5"),]), sum(genedf$end-genedf$start), (sum(hg19_coord$length)-sum(genedf$end-genedf$start))), nrow = 2, ncol = 2))

## Genome-wide average
df <- as.data.frame(table(ctdf$Func.refGene)/nrow(ctdf))
colnames(df) <- c("class", "fraction")
df$fraction[df$class == "exonic"] <- df$fraction[df$class == "exonic"] + df$fraction[df$class == "exonic;splicing"]
df$fraction[df$class == "ncRNA_exonic"] <- df$fraction[df$class == "ncRNA_exonic"] + df$fraction[df$class == "ncRNA_exonic;splicing"]
df$fraction[df$class == "UTR5"] <- df$fraction[df$class == "UTR5"] + df$fraction[df$class == "UTR5;UTR3"]
df$fraction[df$class == "upstream"] <- df$fraction[df$class == "upstream"] + df$fraction[df$class == "upstream;downstream"]
df <- df[!(df$class %in% c("exonic;splicing", "ncRNA_exonic;splicing", "UTR5;UTR3", "upstream;downstream")),]

df$class <- factor(df$class, levels = c("upstream", "UTR5", "exonic", "splicing", "intronic", "UTR3", "downstream", "ncRNA_UTR5", "ncRNA_exonic", "ncRNA_splicing", "ncRNA_intronic", "intergenic"), )
pie(df$fraction)
df <- df[order(df$class, decreasing = F),]

pdf("FIP_revision/pie.htgts.junction.distribution.control.pdf", height=4.5, width=5)
pie(df$fraction, init.angle=90, labels=paste0(df$class, ",\n", round(df$fraction*100), "%"), col=c("#FFFACD", "#FFEC8B", "#FFC125", "#EE9A00", "#FF4500", "#FF6347", "#F08080", "#EEB4B4", "#FFC1C1", "#D15FEE", "#D8BFD8", "#C6E2FF"), main=paste0("Control, # of unique junctions = ", nrow(ctdf)))
par(new=TRUE)
symbols(x=0,y=0,circles=0.4, inches=FALSE, add=TRUE, bg="white")
dev.off()

df <- as.data.frame(table(e2df$Func.refGene)/nrow(e2df))
colnames(df) <- c("class", "fraction")
df$fraction[df$class == "exonic"] <- df$fraction[df$class == "exonic"] + df$fraction[df$class == "exonic;splicing"]
df$fraction[df$class == "ncRNA_exonic"] <- df$fraction[df$class == "ncRNA_exonic"] + df$fraction[df$class == "ncRNA_exonic;splicing"]
df$fraction[df$class == "UTR5"] <- df$fraction[df$class == "UTR5"] + df$fraction[df$class == "UTR5;UTR3"]
df$fraction[df$class == "upstream"] <- df$fraction[df$class == "upstream"] + df$fraction[df$class == "upstream;downstream"]
df <- df[!(df$class %in% c("exonic;splicing", "ncRNA_exonic;splicing", "UTR5;UTR3", "upstream;downstream")),]

df$class <- factor(df$class, levels = c("upstream", "UTR5", "exonic", "splicing", "intronic", "UTR3", "downstream", "ncRNA_UTR5", "ncRNA_exonic", "ncRNA_splicing", "ncRNA_intronic", "intergenic"), )
pie(df$fraction)
df <- df[order(df$class, decreasing = F),]

pdf("FIP_revision/pie.htgts.junction.distribution.e2treated.pdf", height=4.5, width=5)
pie(df$fraction, init.angle=90, labels=paste0(df$class, ",\n", round(df$fraction*100), "%"), col=c("#FFFACD", "#FFEC8B", "#FFC125", "#EE9A00", "#FF4500", "#FF6347", "#F08080", "#EEB4B4", "#FFC1C1", "#D15FEE", "#D8BFD8", "#C6E2FF"), main=paste0("E2 treated, # of unique junctions = ", nrow(e2df)))
par(new=TRUE)
symbols(x=0,y=0,circles=0.4, inches=FALSE, add=TRUE, bg="white")
dev.off()


# Circos links
df <- read.csv("HTGTS.count.lograt.average.allcells.bait.marked.txt", header=T, as.is=T, sep="\t")

## MCF7 cells
tempdf <- df[df$msr >=2 & !is.na(df$msr),]
tempdf$V1 <- paste0("hs", tempdf$chr)
tempdf$V2 <- round((tempdf$start+tempdf$end)/2)
tempdf$V3 <- tempdf$V2
tempdf$V4 <- "hs11"
tempdf$V5 <- 70658290
tempdf$V6 <- 70658290

df <- df[df$mrr >=2 & !is.na(df$mrr),]
df$V1 <- paste0("hs", df$chr)
df$V2 <- round((df$start+df$end)/2)
df$V3 <- df$V2
df$V4 <- "hs17"
df$V5 <- 38491451
df$V6 <- 38491451

df <- rbind(tempdf, df)
df$V7[df$V4 == "hs11"] <- "color=vdgreen_a3"
df$V7[df$V4 == "hs17"] <- "color=vdpurple_a3"
df$chr_num <-df$chr
df$chr_num[df$chr_num == "X"] <- "23"
df$chr_num[df$chr_num == "Y"] <- "24"
df$chr_num <- as.numeric(df$chr_num)
df <- df[order(df$start, decreasing = F),]
df <- df[order(df$chr_num, decreasing = F),]
df <- df[,colnames(df) != "chr_num"]
write.table(df[,c("V1", "V2", "V3", "V4", "V5", "V6", "V7")], ".MCF7.svlinks.circos", row.names = F, col.names = F, quote = F, sep = "\t")

## T47D cells
df <- read.csv("HTGTS.count.lograt.average.allcells.bait.marked.txt", header=T, as.is=T, sep="\t")
tempdf <- df[df$tsr >=2 & !is.na(df$tsr),]
tempdf$V1 <- paste0("hs", tempdf$chr)
tempdf$V2 <- round((tempdf$start+tempdf$end)/2)
tempdf$V3 <- tempdf$V2
tempdf$V4 <- "hs11"
tempdf$V5 <- 70658290
tempdf$V6 <- 70658290

df <- df[df$trr >=2 & !is.na(df$trr),]
df$V1 <- paste0("hs", df$chr)
df$V2 <- round((df$start+df$end)/2)
df$V3 <- df$V2
df$V4 <- "hs17"
df$V5 <- 38491451
df$V6 <- 38491451

df <- rbind(tempdf, df)
df$V7[df$V4 == "hs11"] <- "color=vdgreen_a3"
df$V7[df$V4 == "hs17"] <- "color=vdpurple_a3"
df$chr_num <-df$chr
df$chr_num[df$chr_num == "X"] <- "23"
df$chr_num[df$chr_num == "Y"] <- "24"
df$chr_num <- as.numeric(df$chr_num)
df <- df[order(df$start, decreasing = F),]
df <- df[order(df$chr_num, decreasing = F),]
df <- df[,colnames(df) != "chr_num"]
write.table(df[,c("V1", "V2", "V3", "V4", "V5", "V6", "V7")], ".T47D.svlinks.circos", row.names = F, col.names = F, quote = F, sep = "\t")

