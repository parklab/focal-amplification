library(dplyr)
library(grid)
library(fields)
library(stringr)
library(plotrix)
library(mapplots)
library(viridis)
library(reshape2)
library(devtools)
library(ComplexHeatmap)

# Driver OncoPrint
drvp <- read.csv("List.drivers.final.txt", header=T, as.is=T, sep="\t")
oncdf <- purdf[purdf$availability=="yes",]

## Number of genes to be plotted from each class
sort(table(drvp$gene[drvp$driver=="MUTATION"]), decreasing = T) # top 20
sort(table(drvp$gene[drvp$driver=="GERMLINE"]), decreasing = T) # top 4
sort(table(drvp$gene[drvp$driver=="AMP"]), decreasing = T) 
sort(table(drvp$gene[drvp$driver=="DEL"]), decreasing = T)

sort(table(drvp$protein[drvp$driver=="MUTATION" & drvp$gene=="PIK3CA"]), decreasing = T)
sort(table(drvp$protein[drvp$driver=="MUTATION" & drvp$gene=="TP53"]), decreasing = T)

genelist <- unique(c(names(sort(table(drvp$gene[drvp$driver=="AMP"]), decreasing = T))[1:10], names(sort(table(drvp$gene[drvp$driver=="MUTATION"]), decreasing = T))[1:20], names(sort(table(drvp$gene[drvp$driver=="GERMLINE"]), decreasing = T))[1:4], names(sort(table(drvp$gene[drvp$driver=="DEL"]), decreasing = T))[1:10]))
genelist <- c("ERBB2", "CCND1", "FGF3", "RSF1", "PAK1", "ZNF703", "FGFR1", "MYC", "TRPS1", "PPM1D", "USP32", "ZNF217", "MDM2", "MCL1", "PTEN", "CDKN2A", "RB1", "TP53", "BRCA1", "BRCA2", "PALB2", "MAP3K1", "MAP2K4", "PIK3CA", "PIK3R1", "AKT1", "GATA3", "RUNX1", "CBFB", "SF3B1", "CDH1", "KMT2C", "KMT2D", "ARID1A", "SETD2") # Reordered gene list
length(genelist)

for (i in 1:length(genelist)){
  print(paste0("processing ", genelist[i]))
  oncdf[,ncol(oncdf)+1] <- ""
  colnames(oncdf)[ncol(oncdf)] <- genelist[i]
  for (j in 1:nrow(oncdf)){
    if (nrow(drvp[drvp$study_id == oncdf$study_id[j] & drvp$gene == genelist[i] & drvp$driver != "PARTIAL_AMP",]) != 0){
      tempdf <- drvp[drvp$study_id == oncdf$study_id[j] & drvp$gene == genelist[i] & drvp$driver != "PARTIAL_AMP",]
      for (k in 1:nrow(tempdf)){
        infoline <- c()
        if (tempdf$missense[k] == 1){
          infoline <- c(infoline, "8missense")
        }
        if (tempdf$nonsense[k] == 1){
          infoline <- c(infoline, "6nonsense")
        }
        if (tempdf$splice[k] == 1){
          infoline <- c(infoline, "5splice")
        }
        if (tempdf$inframe[k] == 1){
          infoline <- c(infoline, "7inframe")
        }
        if (tempdf$frameshift[k] == 1){
          infoline <- c(infoline, "4frameshift")
        }
        if (tempdf$biallelic[k] == "true"){
          infoline <- c(infoline, "3biallelic")
        }
        if (tempdf$driver[k] == "AMP"){
          infoline <- c(infoline, "9amp")
        }
        if (tempdf$driver[k] == "DEL"){
          infoline <- c(infoline, "2del")
        }
        if (tempdf$driver[k] == "GERMLINE"){
          infoline <- c(infoline, "1germline")
        }
      }
      oncdf[j,ncol(oncdf)] <- paste(unique(infoline), collapse="::")
      rm(infoline)
    }
  }
}

oncdf$er_plot <- "no"
oncdf$pr_plot <- "no"
oncdf$her2_plot <- "no"
oncdf$pam50_plot <- "no"
for (i in 1:nrow(oncdf)){
  oncdf$er_plot[i] <- demdf$er_plot[demdf$study_id==oncdf$study_id[i]]
  oncdf$pr_plot[i] <- demdf$pr_plot[demdf$study_id==oncdf$study_id[i]]
  oncdf$her2_plot[i] <- demdf$her2_plot[demdf$study_id==oncdf$study_id[i]]
  oncdf$pam50_plot[i] <- demdf$pam50_plot[demdf$study_id==oncdf$study_id[i]]
}

oncdf$p_BRCA1 <- NA
oncdf$p_BRCA2 <- NA
oncdf$p_hrd <- NA
oncdf$hr_status <- "not_analyzed"
for (i in 1:nrow(oncdf)){
  oncdf$p_BRCA1[i] <- hrdf$p_BRCA1[hrdf$study_id == oncdf$study_id[i]]
  oncdf$p_BRCA2[i] <- hrdf$p_BRCA2[hrdf$study_id == oncdf$study_id[i]]
  oncdf$p_hrd[i] <- hrdf$p_hrd[hrdf$study_id == oncdf$study_id[i]]
  oncdf$hr_status[i] <- hrdf$hr_status[hrdf$study_id == oncdf$study_id[i]]
}

oncdf$tba_num <- 0
for (i in 1:nrow(oncdf)){
  oncdf$tba_num[i] <- nrow(hemichr[hemichr$study_id == oncdf$study_id[i],])
}
oncdf$tba_num[oncdf$tba_num > 8] <- 8
oncdf$hrdicho <- "no"
oncdf$hrdicho[oncdf$hr_status == "HR_deficient"] <- "2hrd"
oncdf$hrdicho[oncdf$hr_status == "cannot_be_determined"] <- "1hrp"
oncdf$hrdicho[oncdf$hr_status == "HR_proficient"] <- "1hrp"
oncdf <- oncdf[order(oncdf$ct_sv_hmf, decreasing = T),]
for (i in rev(genelist)){
  oncdf <- oncdf[order(oncdf[,i], decreasing = T),]
}
oncdf <- oncdf[order(oncdf$tba_num, decreasing = F),]
oncdf <- oncdf[order(oncdf$er_plot, decreasing = F),]
oncdf <- oncdf[order(oncdf$hrdicho, decreasing = F),]
oncdf$ord_num_plot <- 1:nrow(oncdf)

mutix <- oncdf[,c(which(colnames(oncdf)=="ERBB2"):which(colnames(oncdf)=="SETD2"))]
mutix <- as.matrix(mutix)
mutix <- t(mutix)
mutix <- gsub("9amp", "amp", mutix)
mutix <- gsub("8missense", "missense", mutix)
mutix <- gsub("7inframe", "inframe", mutix)
mutix <- gsub("6nonsense", "nonsense", mutix)
mutix <- gsub("5splice", "splice", mutix)
mutix <- gsub("4frameshift", "frameshift", mutix)
mutix <- gsub("3biallelic", "biallelic", mutix)
mutix <- gsub("2del", "del", mutix)
mutix <- gsub("1germline", "germline", mutix)

col = c(missense = "#006400", nonsense = "#000000", splice = "#FF34B3", inframe = "#CD0000", frameshift = "#0000CD", biallelic = "#BABABA", amp = "#EE2C2C", del = "#00688B", germline = "#FF9912")

test.oncoprint = oncoPrint(mutix, get_type = function(x) strsplit(x, "::")[[1]],
                           alter_fun = list(
                             background = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#E0E0E0", col = NA)),
                             biallelic = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["biallelic"], col = NA)),
                             amp = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["amp"], col = NA)),
                             del = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["del"], col = NA)),
                             missense = function(x, y, w, h) grid.rect(x, y, w*0.7, h*0.35, gp = gpar(fill = col["missense"], col = NA)),
                             nonsense = function(x, y, w, h) grid.rect(x, y, w*0.7, h*0.35, gp = gpar(fill = col["nonsense"], col = NA)),
                             splice = function(x, y, w, h) grid.rect(x, y, w*0.7, h*0.35, gp = gpar(fill = col["splice"], col = NA)),
                             inframe = function(x, y, w, h) grid.rect(x, y, w*0.7, h*0.35, gp = gpar(fill = col["inframe"], col = NA)),
                             frameshift = function(x, y, w, h) grid.rect(x, y, w*0.7, h*0.35, gp = gpar(fill = col["frameshift"], col = NA)),
                             germline = function(x, y, w, h) grid.rect(x, y, w*0.7, h*0.6, gp = gpar(fill = col["germline"], col = NA))
                           ), col = col, row_order = 1:nrow(mutix), column_order = 1:ncol(mutix), remove_empty_columns = FALSE, column_title = paste0("Breast cancer (n = ", nrow(oncdf), ")")
                           , heatmap_legend_param = list(title = "Alterations", nrow = 2, title_position = "leftcenter"))

## Plotting the oncoprint plot
pdf("oncoprint.pdf", height=8, width=9)
draw(test.oncoprint, heatmap_legend_side = "bottom")
dev.off()

## Ribbon plot for each clinical feature
pdf("image.er.status.oncoprint.pdf", height=2, width=7)
image(matrix(as.numeric(oncdf$er_plot)), col=c(rgb(219/255, 112/255, 147/255, 1), rgb(255/255, 240/255, 245/255, 1), rgb(200/255, 200/255, 200/255, 1)), xaxt = "n", yaxt = "n")
dev.off()

pdf("image.pr.status.oncoprint.pdf", height=2, width=7)
image(matrix(as.numeric(oncdf$pr_plot)), col=c(rgb(154/255, 50/255, 205/255, 1), rgb(255/255, 240/255, 245/255, 1), rgb(200/255, 200/255, 200/255, 1)), xaxt = "n", yaxt = "n")
dev.off()

pdf("image.her2.status.oncoprint.pdf", height=2, width=7)
image(matrix(as.numeric(oncdf$her2_plot)), col=c(rgb(220/255, 20/255, 60/255, 1), rgb(255/255, 106/255, 106/255, 1), rgb(255/255, 240/255, 245/255, 1), rgb(200/255, 200/255, 200/255, 1)), xaxt = "n", yaxt = "n")
dev.off()

pdf("image.pam50.status.oncoprint.pdf", height=2, width=7)
image(matrix(as.numeric(oncdf$pam50_plot)), col=c("#FFB5C5", "#CD6889", "#DC143C", "#FF9912", "#63B8FF", rgb(200/255, 200/255, 200/255, 1)), xaxt = "n", yaxt = "n")
dev.off()

pdf("image.wgd.status.oncoprint.pdf", height=2, width=7)
image(matrix(as.numeric(factor(oncdf$is.wgd, levels = c("ND", "WGD")))), col=rev(twocolors), xaxt = "n", yaxt = "n")
dev.off()

pdf("image.sex.status.oncoprint.pdf", height=2, width=7)
image(matrix(as.numeric(factor(oncdf$gender, levels = c("FEMALE", "MALE")))), col=c("#FFF0F5", "#8B8386"), xaxt = "n", yaxt = "n")
dev.off()

pdf("image.age.level.oncoprint.pdf", height=2, width=7)
image(matrix(as.numeric(oncdf$age)), col=gray.colors(20, start = 0.1, end = 0.9), xaxt = "n", yaxt = "n")
dev.off()

pdf("image.age.level.oncoprint.legend.pdf", height=4, width=4)
image.plot(matrix(as.numeric(oncdf$age)), col=gray.colors(20, start = 0.1, end = 0.9), xaxt = "n", yaxt = "n", legend.only = T, horizontal = T)
dev.off()

pdf("image.ploidy.level.oncoprint.pdf", height=2, width=7)
image(matrix(as.numeric(oncdf$ploidy)), col=RColorBrewer::brewer.pal(9, "YlOrRd"), xaxt = "n", yaxt = "n")
dev.off()

pdf("image.ploidy.level.oncoprint.legend.pdf", height=4, width=4)
image.plot(matrix(as.numeric(oncdf$ploidy)), col=RColorBrewer::brewer.pal(9, "YlOrRd"), xaxt = "n", yaxt = "n", legend.only = T, horizontal = T)
dev.off()

pdf("image.tbamp.level.oncoprint.pdf", height=2, width=7)
image(matrix(as.numeric(oncdf$tba_num)), col=RColorBrewer::brewer.pal(9, "YlGnBu"), xaxt = "n", yaxt = "n")
dev.off()

pdf("image.tbamp.level.oncoprint.legend.pdf", height=4, width=4)
image.plot(matrix(as.numeric(oncdf$tba_num)), col=RColorBrewer::brewer.pal(9, "YlGnBu"), xaxt = "n", yaxt = "n", legend.only = T, horizontal = T)
dev.off()

pdf("barplot.svtype.oncoprint.pdf", height=4, width=7)
barplot(t(as.matrix(oncdf[,c("ct_DEL", "ct_DUP", "ct_h2hINV", "ct_t2tINV", "ct_TRA")])), col=c(rgb(0/255,0/255,255/255,.7), rgb(0/255,128/255,128/255,.7), rgb(220/255,20/255,60/255,.7), rgb(128/255,128/255,0/255,.7), rgb(85/255,26/255,139/255,.7)), border=NA, space=0, las=1)
dev.off()

pdf("image.hrd.status.oncoprint.pdf", height=2, width=7)
image(matrix(as.numeric(factor(oncdf$hrdicho, levels = c("1hrp", "2hrd")))), col=c("#FFDAB9", "#FF9912"), xaxt = "n", yaxt = "n")
dev.off()

pdf("image.ct_snv.level.oncoprint.pdf", height=2, width=7)
image(matrix(as.numeric(log10(oncdf$ct_snv))), col=RColorBrewer::brewer.pal(9, "RdPu"), xaxt = "n", yaxt = "n")
dev.off()

pdf("image.ct_snv.level.oncoprint.legend.pdf", height=4, width=4)
image.plot(matrix(as.numeric(log10(oncdf$ct_snv))), col=RColorBrewer::brewer.pal(9, "RdPu"), xaxt = "n", yaxt = "n", legend.only = T, horizontal = T)
dev.off()

pdf("image.ct_indel.level.oncoprint.pdf", height=2, width=7)
image(matrix(as.numeric(log10(oncdf$ct_indel))), col=RColorBrewer::brewer.pal(9, "PuBuGn"), xaxt = "n", yaxt = "n")
dev.off()

pdf("image.ct_indel.level.oncoprint.legend.pdf", height=4, width=4)
image.plot(matrix(as.numeric(log10(oncdf$ct_indel))), col=RColorBrewer::brewer.pal(9, "PuBuGn"), xaxt = "n", yaxt = "n", legend.only = T, horizontal = T)
dev.off()

sigdf <- read.csv("List.SBS.ID.signatures.final.txt", header=T, as.is=T, sep="\t")
df <- sigdf[sigdf$availability == "yes",]
df$ord_num_plot <- NA
for (i in 1:nrow(df)){
  df$ord_num_plot[i] <- oncdf$ord_num_plot[oncdf$study_id == df$study_id[i]]
}
df <- df[order(df$ord_num_plot, decreasing = F),]
df$tba_num <- NA
for (i in 1:nrow(df)){
  df$tba_num[i] <- oncdf$tba_num[oncdf$study_id == df$study_id[i]]
}

df$SBS_msi <- df$SBS15 + df$SBS20	+ df$SBS21 + df$SBS26 + df$SBS95 + df$SBS96 + df$SBS97 + df$SBS98
df$SBS_other <- df$SBS12 + df$SBS36 + df$SBS41 + df$SBS84 + df$SBS85 + df$SBS93 + df$SBS31 + df$SBS35 + df$SBS86
tempdf <- df[,c("SBS1", "SBS2", "SBS3", "SBS5", "SBS8", "SBS13", "SBS17a", "SBS17b", "SBS18", "SBS30", "SBS_msi", "SBS_other", "ct_musical_snv", "study_id", "tba_num")]
for (i in 1:nrow(tempdf)){
  tempdf[i,c(1:12)] <- tempdf[i,c(1:12)]/tempdf$ct_musical_snv[i]
}
pdf("barplot.SBS.spectra.pdf", height=2.8, width=7)
barplot(as.matrix(t(tempdf[,c(1:12)])), border=NA, las=1, col=snvspeccolor, space = 0)
dev.off()

## Comparisons of contribution of signatures
t.test(tempdf$SBS18 ~ ifelse(tempdf$tba_num == 0, "no_tba", "tba"))
t.test(tempdf$SBS18[tempdf$SBS18 != 0] ~ ifelse(tempdf$tba_num[tempdf$SBS18 != 0] == 0, "no_tba", "tba"))
t.test(tempdf$SBS18 ~ ifelse(tempdf$tba_num == 0, "no_tba", "tba"))
t.test(tempdf$SBS2 ~ ifelse(tempdf$tba_num == 0, "no_tba", "tba"))
t.test(tempdf$SBS13 ~ ifelse(tempdf$tba_num == 0, "no_tba", "tba"))
t.test(tempdf$SBS2[tempdf$SBS2 != 0] ~ ifelse(tempdf$tba_num[tempdf$SBS2 != 0] == 0, "no_tba", "tba"))
t.test(tempdf$SBS13[tempdf$SBS13 != 0] ~ ifelse(tempdf$tba_num[tempdf$SBS13 != 0] == 0, "no_tba", "tba"))
