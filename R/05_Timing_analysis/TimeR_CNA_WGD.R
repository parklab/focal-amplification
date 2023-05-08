library(dplyr)
library(grid)
library(fields)
library(stringr)
library(mapplots)
library(reshape2)
library(data.table)
library(MutationTimeR)

# Running Mutation TimeR
args <- commandArgs(trailingOnly = TRUE)

purity <- as.numeric(as.character(args[2]))
vcf <- readVcf(paste("/n/data1/hms/dbmi/park/jake/TBAmp/TimeR/rerun/Renamed_edited_VCF/", args[1], ".purple.somatic.filtered.updated.vcf.gz", sep = ""))
cn_read <- read.table(paste("/n/data1/hms/dbmi/park/jake/TBAmp/TimeR/rerun/Renamed_edited_CNA/", args[1], ".purple.cnv.somatic.edited.tsv", sep = ""), header=T)
cn_read$majorAlleleCopyNumber_round[cn_read$majorAlleleCopyNumber_round <0] <- 0
cn_read$clonal_frequency <- purity
bb <- GRanges(cn_read, major_cn=cn_read$majorAlleleCopyNumber_round , minor_cn=cn_read$minorAlleleCopyNumber_round, clonal_frequency=purity)
mt <- mutationTime(vcf, bb, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)
writeVcf(vcf, paste("/n/data1/hms/dbmi/park/jake/TBAmp/TimeR/rerun/Timer_annotated_VCF/", args[1], ".vcf", sep=""))
pdf(paste("/n/data1/hms/dbmi/park/jake/TBAmp/TimeR/rerun/Timer_PDF/", args[1], ".pdf", sep=""))
plotSample(vcf,bb)
dev.off()

df <- data.frame(chrom=seqnames(bb), start=start(bb), end=end(bb))
df$major_cn <- bb$major_cn
df$minor_cn <- bb$minor_cn
df$type <- bb$type
df$time1 <- bb$time
df$time1_lo <- bb$time.lo
df$time1_up <- bb$time.up
df$time2 <- bb$time.2nd
df$time2_lo <- bb$time.2nd.lo
df$time2_up <- bb$time.2nd.up
df$time_star <- bb$time.star
df$n_snv <- bb$n.snv_mnv

write.table(df, paste0("/n/data1/hms/dbmi/park/jake/TBAmp/TimeR/rerun/Timer_annotated_CNA/", args[1], ".cnv.timeR.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

min.dist=0.05
w <- which(bb$n.snv_mnv > 20 & !is.na(bb$time))
s <- seq(0,1,0.01)
l2 <- pmin(bb$time.lo, bb$time - min.dist)[w]
u2 <- pmax(bb$time.up, bb$time + min.dist)[w]
l1 <- (l2 +  bb$time[w])/2
u1 <- (u2+  bb$time[w])/2
wd <- as.numeric(width(bb)[w])
o <- sapply(s, function(i) sum(wd * ( (l2 <= i & u2 >=i) + (l1 <= i & u1 >= i))))
m <- s[which.max(o)]
l <- pmin(bb$time.lo, bb$time - min.dist)
u <- pmax(bb$time.up, bb$time + min.dist)
w <- which(l <= m & u >= m)
avgCi <- weighted.mean(bb$time.up- bb$time.lo, width(bb), na.rm=TRUE)
sd.wgd <- sqrt(weighted.mean((bb$time[w] - m)^2, width(bb)[w], na.rm=TRUE))
sd.all <- sqrt(weighted.mean((bb$time - m)^2, width(bb), na.rm=TRUE))
write.table(c(nt.wgd=sum(as.numeric(width(bb))[w]), nt.total=sum(as.numeric(width(bb))[!is.na(bb$time)]), time.wgd=m, n.wgd=length(w), n.all = sum(!is.na(bb$time)), chr.wgd = length(unique(seqnames(bb)[w])), chr.all = length(unique(seqnames(bb)[!is.na(bb$time)])), sd.wgd=sd.wgd, avg.ci=avgCi, sd.all=sd.all), paste0("/n/data1/hms/dbmi/park/jake/TBAmp/TimeR/rerun/Timer_WGD_assessment/", args[1], ".wgd.synchronicity.txt"), sep="\t", quote=F)
