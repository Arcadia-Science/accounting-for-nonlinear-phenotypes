library("RColorBrewer")
library('gplots')
library('missMDA')
library('FactoMineR')
library('RNOmni')
library('umap')
library('energy')
library('mgcv')
library(lmtest)
library(scales)

#######################################
#####Load and clean phenotype data#####
#######################################
#Load
dat = read.csv('~/Desktop/strainmeans.csv')

#Get all strains
strains = unique(dat$strainid)

#Split on phenotype
x = split(dat, dat$varname)

#Match on strains
for(i in 1:length(x)){
  x[[i]] = x[[i]][match(strains, x[[i]]$strainid),]
}

#Combine
phenos = as.data.frame(do.call(cbind, lapply(x, function(y) y$zscore)))
rownames(phenos) = strains

#Filter phenotypes
phenos = phenos[,apply(phenos, 2, function(x) sum(!is.na(x)))>50]

#Filter strains
phenos = phenos[apply(phenos, 1, function(x) sum(!is.na(x)))>100,]


