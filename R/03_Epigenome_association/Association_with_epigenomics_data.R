
# codes blow in R ver 4.1.2

association.with.chromatin.features<-function(file="Data/epigenomics.data.boundary.number.250kb.RData"){
  library(lasso2)
  load(file) #load the file contains values of number of amplication boundaries and number of binding for each factor in 250kb bins. Values in centromere and peri-centric regions were filtered out.
  
  #all
  #this is for all amplication boudaries  
  d2=d2.all
  id=which(rowSums(d2[,2:12])>0) #bins where there was no binding from any of the epigenetic features were further filtered out.
  d3=d2[id,]
   
  lm.lasso=l1ce(b.n ~ ER.E2 + CTCF + TOP2B + DHS + R.loop + H3K9me3 + H3K27ac + H3K4me3 + Pol2 + H3K36me3, data = d3, standardize=T, absolute.t=T, bound=0.39, trace=T) 
  #lasso regression. panalty parameters were chosen so that the least two important features were not selected.   
  
  pval=list();
  type="all";pval[[type]]=summary(lm.lasso)$coefficients[-1,4]
  
  #ER positive
  #this is for amplicon boundaries from ER positive samples  
  d2=d2.er.pos
     
  id=which(rowSums(d2[,2:12])>0)
  length(id)
  d3=d2[id,]
  lm.lasso<-l1ce(b.n ~ ER.E2 + CTCF + TOP2B + DHS + R.loop + H3K9me3 + H3K27ac + H3K4me3 + Pol2 + H3K36me3, data = d3, standardize=T, absolute.t=T, trace=T, bound=0.29)
  type="er.pos";pval[[type]]=summary(lm.lasso)$coefficients[-1,4]
  pval
  
#ER negative
#this is for amplicon boundaries from ER negative samples  
  d2=d2.er.neg  
  id=which(rowSums(d2[,2:12])>0)
  length(id)
  d3=d2[id,]
  lm.lasso<-l1ce(b.n ~ ER.E2 + CTCF + TOP2B+ DHS + R.loop + H3K9me3 + H3K27ac + H3K4me3 + Pol2 + H3K36me3, data = d3, standardize=T, absolute.t=T, trace=T, bound=0.08)
  type="er.neg";pval[[type]]=summary(lm.lasso)$coefficients[-1,4]
  pval 
  return(pval)
}

comparison.er.e2.control<-function(file="Data/epigenomics.data.breast.RData"){
load(file)
t=er.intensity.1mb
d1=density(t[,4],na.rm=T)
plot(d1$x,d1$y,type="l",xlim=c(0,14),lwd=2,col=rgb(0, 0, 1,0.5),xlab="ER binding intensity",ylab="Density",main="Control vs E2",axes=F)
axis(1)
axis(2,las=2)
d2=density(t[,5],na.rm=T)
lines(d2$x,d2$y,type="l",col=rgb(1, 0, 0,0.5),lwd=2)
polygon(c(d1$x, rev(d1$x)), c(d1$y ,rev(rep(0,times=length(d1$y)))), col = rgb(0, 0, 1,0.3) )
polygon(c(d2$x, rev(d2$x)), c(d2$y ,rev(rep(0,times=length(d2$y)))), col = rgb(1, 0, 0,0.3) )
id1=which(t$chr=="chr11" & t$end==71e6)
id2=which(t$chr=="chr11" & t$end==79e6)
id3=which(t$chr=="chr8" & t$end==38e6)
id4=which(t$chr=="chr17" & t$end==39e6)
ids=c(id1,id2,id3,id4);
points(t[ids,5],rep(0,times=4),col=rgb(1, 0, 0,0.5),cex=1,pch=16)
points(t[ids,4],rep(0,times=4),col=rgb(0, 0, 1,0.5),cex=1,pch=16)
text(c(t[ids,5]),c(0.15,0.1,0.1,0.1),labels=c("SHANK2","TENM4","ZNF703","RARA"))
}  

association.recurrence.e2.er.intensity<-function(rec.file="Data/all.boundary.patient.recurrence.breast.100kb.txt",er.file="Data/epigenomics.data.breast.RData"){
load(er.file);
ta=read.table(file=rec.file,sep="\t",stringsAsFactor=F);
t=e2.er.acc.peak.int.100kb;
cols=c("#eff3ff","#bdd7e7","#6baed6","#3182bd","#08519c");
boxplot(t$V4[ta$V4==0],t$V4[ta$V4>=1 & ta$V4<=2],t$V4[ta$V4>=3 & ta$V4<=4],t$V4[ta$V4>=5 & ta$V4<=6],t$V4[ta$V4>=7],outline=F,axes=F,ylab="Accumulated E2-responsive ERa binding intensity in 100kb bin",col=cols)
axis(1,at=1:5,label=c("Recurrence=0","1-2","3-4","5-6",">=7"))
mtext(1,at=1:5,text=c(paste0("(n=",sum(ta$V4==0),")"),paste0("(",sum(ta$V4>=1 & ta$V4<=2),")"),paste0("(",sum(ta$V4>=3 & ta$V4<=4),")"),paste0("(",sum(ta$V4>=5 & ta$V4<=6),")"),paste0("(",sum(ta$V4>=7),")")),line=2.2)
axis(2)
}  

association.recurrence.3d.contact.t47d<-function(sv.file="Data/breast.cancer.278.boundary.sv.bedpe",3d.file.folder="Data/contact"){
 chrs=c(1:22,"X");
 library(gtools)
 pairs=combinations(23,2,chrs,repeats.allowed=F)
 files=paste0("Data/contact/",pairs[,1],"_",pairs[,2],"_oe_2.5mb.txt");
 all=unlist(lapply(files,function(file) {
   t=read.table(file,sep="",stringsAsFactor=F)
   t$V3}))

t=read.table(file=sv.file,sep="",stringsAsFactor=F)
bin=2.5e6;
d=unlist(lapply(i:dim(t)[1],function(i){
t.chr1=gsub("chr","",t$V1[i])
t.p1=t$V2[i];
t.chr2=gsub("chr","",t$V4[i])
t.p2=t$V5[i];
if(t.chr1!="X" & t.chr2!="X"){
if(as.numeric(t.chr1)<as.numeric(t.chr2)){chr1=t.chr1;p1=t.p1;chr2=t.chr2;p2=t.p2}
if(as.numeric(t.chr1)>as.numeric(t.chr2)){chr1=t.chr2;p1=t.p2;chr2=t.chr1;p2=t.p1}}
if(t.chr2=="X"){chr1=t.chr1;p1=t.p1;chr2=t.chr2;p2=t.p2}
if(t.chr1=="X"){chr1=t.chr2;p1=t.p2;chr2=t.chr1;p2=t.p1}
chr1;p1;chr2;p2
options(scipen=999)
p1=floor(p1/bin)*bin;p1
p2=floor(p2/bin)*bin;p2
paste0(chr1,",",p1,",",chr2,",",p2)}))
s=table(d)
all.amp=unique(d)
rec.amp=names(s[s>=2])

  chr1=strsplit(a,",")[[1]][1];
all.amp.contact=unlist(lapply(all.amp,function(a){
  chr1=strsplit(a,",")[[1]][1];
  p1=as.numeric(strsplit(a,",")[[1]][2]);
  chr2=strsplit(a,",")[[1]][3];
  p2=as.numeric(strsplit(a,",")[[1]][4]);
  t=read.table(file=paste0("~/groups/chromoplexy/data/hic/contact/",chr1,"_",chr2,"_oe_",bin/  p2=as.numeric(strsplit(a,",")[[1]][4]);
  t=read.table(file=paste0("~/groups/chromoplexy/data/hic/contact/",chr1,"_",chr2,"_oe_",bin/1e6,"mb_hg38.txt"),sep="",stringsAsFactor=F)
  id=which(t$V1==p1 & t$V2==p2);
  t$V3[id]
}))

names(rec.amp)=rec.amp
rec.amp.contact=unlist(lapply(rec.amp,function(a){
  chr1=strsplit(a,",")[[1]][1];
  p1=as.numeric(strsplit(a,",")[[1]][2]);
  chr2=strsplit(a,",")[[1]][3];
  p2=as.numeric(strsplit(a,",")[[1]][4]);
  t=read.table(file=paste0("~/groups/chromoplexy/data/hic/contact/",chr1,"_",chr2,"_oe_",bin/1e6,"mb_hg38.txt"),sep="",stringsAsFactor=F)
  id=which(t$V1==p1 & t$V2==p2);
  t$V3[id]
}))
wilcox.test(all,all.amp.contact)
wilcox.test(all.amp.contact,rec.amp.contact) 
  
par(mar=c(6,5,4,5))
ddf=data.frame(Contact=c(all,all.amp.contact,rec.amp.contact),Freq=c(rep("all",times=length(all)),rep("amp",times=length(all.amp.contact)),rep("rec.amp",times=length(rec.amp.contact))))
cols1=c(rgb(0,0,0,0.3),rgb(0.0,0,0.6,0.3),rgb(0.4,0,0.6,0.3));
cols2=c(rgb(0,0,0,0.0),rgb(0.0,0,0.6,0.1),rgb(0.4,0,0.6,0.5));
boxplot(Contact ~ Freq, data = ddf, lwd = 1, ylab = 'Contact frequencies (Observed/Expected)',ylim=c(0,3),outline=F,col=cols1,axes=F)
stripchart(Contact ~ Freq, vertical = TRUE, data = ddf, 
    method = "jitter", add = TRUE, pch = 19, col = cols2,cex=1.5)
axis(1,at=1:3,label=c("All pairs","Amp boundaries","Rec"),line=-0.5)
axis(1,at=1:3,label=c(paste0("(n=",length(all),")"),paste0("(",length(all.amp.contact),")"),paste0("(",length(rec.amp.contact),")")),line=0.5,tick=F)
mtext(3,at=c(1.5,2.5),text=c("***","*"))
axis(2,las=2)  
}
  
translocation.network<-function(file="Data/trans.frequencies.byarm.breast.cancer.278.txt"){
t=read.table(file=file,sep="",header=T,stringsAsFactor=F)
thr=3;
id=which(t$amptracount>=thr)
links=data.frame(from=t$arm1[id],to=t$arm2[id],weight=t$amptracount[id],type="chr.arm")
library(igraph);
chrs=as.character(nodes$id);
names(chrs)=chrs;
count=unlist(lapply(chrs,function(chr){
  co=sum(t$amptracount[t$arm1==chr])+sum(t$amptracount[t$arm2==chr])
  return(co)
}))
nodes2=data.frame(id=chrs[rev(order(count))],type="chr.arm")
count2=rev(sort(count))
library(RColorBrewer)
thr=6
cols1=unlist(lapply(count2,function(cnt){
  n=round(cnt/10);
  if(n>=thr){n=thr}
  return(brewer.pal(n=thr,name="Blues")[n])
  }))
  net2 <- graph.data.frame(links, nodes2, directed=F)
  x11();
col2=rgb(0,0,0,0.3) #edge color
V(net2)$color=cols1;
E(net2)$color=col2;
plot(net2,edge.width=E(net2)$weight, vertex.label.color="black",vertex.label.family="Helvetica",vertex.size=20,vertex.frame.color="black");
}  
