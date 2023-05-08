library(dplyr)
library(grid)
library(fields)
library(stringr)
library(mapplots)
library(reshape2)
library(data.table)

# Probability of mutations preceding the regional copy number gain
## Very early clonal: mutations amplified up to the maximal copy of the region (preceding multi-copy gains)
## Early clonal: amplified mutations but not up to the maximal copy of the region
args <- commandArgs(trailingOnly = TRUE)

threshold <- 0.75
df <- read.csv(paste0(args[1], ".sage30.snv.timeR.txt"), header=T, as.is=T, sep="\t")
purity_val <- as.numeric(as.character(args[2]))

df$pVEClonal <- 0
df$pOtherClonal <- 0
df$pSClonal <- 0

for (i in 1:nrow(df)){
  dp <- df$t_alt_count[i]+df$t_ref_count[i]
  var <- df$t_alt_count[i]
  if (is.na(df$MajCN[i])){
    df$pVEClonal[i] <- NA
    df$pOtherClonal[i] <- NA
    df$pSClonal[i] <- NA
  } else {
    tCN <- df$MajCN[i]+df$MinCN[i]
    mjCN <- df$MajCN[i]
    
    tfraction <- tCN*purity_val/(tCN*purity_val + 2*(1-purity_val))
    cutoff <- ceiling(mjCN*threshold)
    
    nEarly <- 0
    nLate <- 0
    nSubcl <- 0
    
    for (j in c(cutoff:mjCN)){
      if (dp*tfraction*j/tCN < var){
        nEarly <- nEarly + sum(dbinom(c(var:dp), dp, tfraction*j/tCN))
      } else if (dp*tfraction*j/tCN >= var){
        nEarly <- nEarly + sum(dbinom(c(0:var), dp, tfraction*j/tCN))
      }
    }
    for (j in c(1:(cutoff-1))){
      if (dp*tfraction*j/tCN < var){
        nLate <- nLate + sum(dbinom(c(var:dp), dp, tfraction*j/tCN))
      } else if (dp*tfraction*j/tCN >= var){
        nLate <- nLate + sum(dbinom(c(0:var), dp, tfraction*j/tCN))
      }
    }
    if (dp*tfraction*0.5/tCN < var){
      nSubcl <- nSubcl + sum(dbinom(c(var:dp), dp, tfraction*0.5/tCN))
    } else if (dp*tfraction*0.5/tCN >= var){
      nSubcl <- nSubcl + sum(dbinom(c(0:var), dp, tfraction*0.5/tCN))
    }
    df$pVEClonal[i] <- round(nEarly/(nEarly+nLate+nSubcl), 9)
    df$pOtherClonal[i] <- round(nLate/(nEarly+nLate+nSubcl), 9)
    df$pSClonal[i] <- round(nSubcl/(nEarly+nLate+nSubcl), 9)
    rm(nEarly)
    rm(nLate)
    rm(nSubcl)
  }
}
write.table(df, paste0(args[1], ".sage30.snv.timeR.VEClonal.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

