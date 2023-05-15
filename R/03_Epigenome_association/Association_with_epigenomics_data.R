
# R ver 4.1.2

association.with.chromatin.features.htgts<-function(file="Data/epigenomics.data.htgts.data.250kb.RData"){
  library(lasso2)
  load(file)
  #tht: HTGTS data
  #t: recurrence of amp boundaries in 250kb
  #t2: recurrence of amp boundaries in 250kb for non-pericentric/centromere regions
  ids=match(paste0(t2$V1,":",t2$V2),paste0(t$V1,":",t$V2)) # remaining bins index after filtering

  #ratio in mcf7 
  sum.e2=tht$mcf7.shank2.e2+tht$mcf7.rara.e2 #breakpoints in E2 treated MCF7
  sum.ctrl=tht$mcf7.shank2.ctrl+tht$mcf7.rara.ctrl #breakpoints in contral MCF7
  logr=log2(sum.e2/(sum.ctrl+1)+1) #log ratio

  mcf=unlist(lapply(1:dim(t)[1],function(i){
  id=which(t$V1[i]==paste0("chr",tht$chr) & t$V2[i]==tht$start);
  if(length(id)>0){return(logr[id])}
  if(length(id)==0){return(0)}
  }))
  # ratios in filtered regions

  #d: number of binding of epigenomic marks in 250 kb
  d=d[ids,] #filter out pericentric/cetromic regions
  d$ratio=mcf[ids] #b.n: ratio of E2/control

  id=which(rowSums(d[,2:12])>0)
  d2=d[id,] #filter out bins with no epigenetic marks

  lm.lasso<-l1ce(ratio ~ ER.E2 + CTCF + TOP2B+ DHS + R.loop + H3K9me3 + H3K27ac + H3K4me3 + Pol2 + H3K36me3, data = d2, standardize=T, absolute.t=T, trace=T, bound=0.215)
  #lasso regression. a penalized parameter was chosen to remove two least important factors

  pval=summary(lm.lasso)$coefficients[-1,4]
  pval
  return(pval)
}

association.with.chromatin.features<-function(file="Data/epigenomics.data.boundary.number.250kb.RData"){
  library(lasso2)
  load(file) #load the file contains values of number of amplifon boundaries and number of binding for each factor in 250kb bins. Values in centromere and peri-centric regions were filtered out.
  
  #all
  #this is for all amplication boudaries  
  d2=d2.all
  id=which(rowSums(d2[,2:12])>0) #bins where there was no binding from any of the epigenetic features were further filtered out.
  d3=d2[id,]
   
  lm.lasso=l1ce(b.n ~ ER.E2 + CTCF + TOP2B + DHS + R.loop + H3K9me3 + H3K27ac + H3K4me3 + Pol2 + H3K36me3, data = d3, standardize=T, absolute.t=T, bound=0.39, trace=T) 
  #lasso regression. panalty parameters were chosen so that two least important features were not selected.   
  
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
  text(c(t[ids,5]),c(0.15,0.1,0.1,0.1),labels=c("SHANK2","TENM4","KCNU1","RARA"))
}  

association.recurrence.e2.er.non.amp<-function(file="Data/sample.number.er.e2.binding.non.amp.100kb.RData"){
  load(file)
  #load SV breakpoints (number of samples) and number of ERa binding in E2 treated cells and info of ampin 100kb bins

  id=which(t1$"amp"==0) #unamplified regions
  t2=t1[id,] 

  n=t2$num_ERa_E2 #number of ERa binding

  a=list();b=list();
  for(i in 0:10){
  a[[i+1]]=sum(n==0 & t2$num_sample==i) #bins without ERa binding
  b[[i+1]]=sum(n>=1 & t2$num_sample==i) #bins with ERa binding
  }
  d=data.frame(recurrence=0:10,fraction=unlist(b)/(unlist(a)+unlist(b)))
  d #recurrence vs percentage of bins with ERa binding
  cor(1:10,(unlist(b)/(unlist(a)+unlist(b)))[2:11])
  cor.test(1:10,(unlist(b)/(unlist(a)+unlist(b)))[2:11])
  return(d)
}
 
