library(RColorBrewer)
library(gplots)
library(missMDA)
library(FactoMineR)
library(RNOmni)
library(umap)
library(energy)
library(mgcv)
library(lmtest)
library(scales)

##############################################################
#####Various pre-processing (just needed to be done once)#####
##############################################################
# dat = readxl::read_xlsx('~/Desktop/dgrp_phenotypes/lafuente_2018_thermalplasticity-bodysize/pgen.1007686.s012.xlsx', skip = 1)
# dat = as.data.frame(dat)
# dat$id = paste(dat$`Body part`, dat$Temperature, sep = '_')
# dat = split(dat, dat$id)
# dat = lapply(dat, function(x) split(x, x$`DGRP line`))
# for(i in 1:length(dat)){
#   for(j in 1:length(dat[[i]])){
#     x = dat[[i]][[j]][1,]
#     y = mean(as.numeric(dat[[i]][[j]]$`Length (mm)`), na.rm = TRUE)
#     x$`Length (mm)` = y
#     dat[[i]][[j]] = x
#   }
# }
# 
# for(i in 1:length(dat)){
#   x = do.call(rbind, dat[[i]])
#   n = gsub('_', '\\.', x$id)
#   n = paste(n, '.temperature.length', sep = '')
#   x = cbind(x$`DGRP line`, x$`Length (mm)`)
#   colnames(x) = c('line', n[1])
#   x[,1] = gsub('line_', '', x[,1])
#   dat[[n[1]]] = x
# }
# 
# for(i in 5:8){
#   write.csv(dat[[i]], paste('~/Desktop/dgrp_phenotypes/lafuente_2018_thermalplasticity-bodysize/', names(dat)[i], '.csv', sep = ''))
# }
# 
# 
# 
# dat = readxl::read_xlsx('~/Desktop/dgrp_phenotypes/gaertner_2015_courtship-behavior/014811_files1.xlsx')
# dat = split(dat, dat$Line)
# 
# for(i in 1:length(dat)){
#   dat[[i]] = colMeans(dat[[i]][,7:ncol(dat[[i]])])
# }
# 
# d = cbind(names(dat), do.call(rbind, dat))
# colnames(d)[1] = 'line'
# write.csv(d, '~/Desktop/dgrp_phenotypes/gaertner_2015_courtship-behavior/courtship.behavior.csv', row.names = FALSE)
# 
# 
# setwd('~/Desktop/dgrp_phenotypes/arya_2015_olfactory-behavior/')
# files = list.files()
# for(i in 1:length(files)){
#   x = read.csv(files[i], header = FALSE)
#   colnames(x) = c('line', gsub('.csv', '', files[i]))
#   x$line = gsub('line_', '', x$line)
#   write.csv(x, files[i])
# }
# 
# setwd('~/Desktop/nonlinear_phenotype_review/dgrp_phenotypes/arya_2015_olfactory-behavior/')
# files = list.files()
# for(i in 1:length(files)){
#   x = read.csv(files[i])
#   write.csv(x[,2:3], files[i], row.names = FALSE)
# }

# 
# setwd('~/Desktop/dgrp_phenotypes/morozova_2015_alchohol-sensitivity/')
# files = list.files()
# for(i in 1:length(files)){
#   x = read.csv(files[i], header = FALSE)
#   colnames(x) = c('line', gsub('.csv', '', files[i]))
#   x$line = gsub('line_', '', x$line)
#   write.csv(x, files[i], row.names = FALSE)
# }
# 
# setwd('~/Desktop/dgrp_phenotypes/mackay_2012_starvation-startle/')
# files = list.files()
# for(i in 1:length(files)){
#   x = read.csv(files[i], header = FALSE)
#   colnames(x) = c('line', gsub('.csv', '', files[i]))
#   x$line = gsub('line_', '', x$line)
#   write.csv(x, files[i], row.names = FALSE)
# }
# 
# setwd('~/Desktop/dgrp_phenotypes/jordan_2012_oxidative-stress/')
# files = list.files()
# for(i in 1:length(files)){
#   x = read.csv(files[i], header = FALSE)
#   colnames(x) = c('line', gsub('.csv', '', files[i]))
#   x$line = gsub('line_', '', x$line)
#   write.csv(x, files[i], row.names = FALSE)
# }

###################################################
#####Pre-process gene expression (compute pcs)#####
###################################################
#Set working directory to repository
setwd('~/Documents/Research/github/accounting-for-nonlinear-phenotypes/')

#Load male
male = read.table('00_data/dgrp_phenotypes/huang_2015_gene-expression/dgrp.array.exp.male.txt', header = TRUE, row.names = 1)

#Load female
female = read.table('00_data/dgrp_phenotypes/huang_2015_gene-expression/dgrp.array.exp.female.txt', header = TRUE, row.names = 1)

#Calculate mean expression per line
n = unique(unlist(lapply(strsplit(colnames(female), "\\."), function(v){v[1]})))
f2 = as.data.frame(matrix(nrow = nrow(female), ncol = 0))
for(i in 1:length(n)){
  x = rowMeans(as.matrix(female[,grep(n[i], colnames(female))]))
  f2[,n[i]] = x
}

n = unique(unlist(lapply(strsplit(colnames(male), "\\."), function(v){v[1]})))
m2 = as.data.frame(matrix(nrow = nrow(male), ncol = 0))
for(i in 1:length(n)){
  x = rowMeans(as.matrix(male[,grep(n[i], colnames(male))]))
  m2[,n[i]] = x
}

#Calculate PCs
pca_f = prcomp(t(f2))
pca_m = prcomp(t(m2))

xf = cbind(rownames(pca_f$x), pca_f$x[,1:30])
xm = cbind(rownames(pca_m$x), pca_m$x[,1:30])

colnames(xf)[1] = 'line'
colnames(xm)[1] = 'line'

colnames(xf)[2:ncol(xf)] = paste('dgrp.array.exp.female.', colnames(xf)[2:ncol(xf)], sep = '')
colnames(xm)[2:ncol(xm)] = paste('dgrp.array.exp.male.', colnames(xm)[2:ncol(xm)], sep = '')

#Save first n pcs
write.csv(xf, '00_data/dgrp_phenotypes/huang_2015_gene-expression/dgrp.array.exp.female.30pcs.csv', row.names = FALSE)
write.csv(xm, '00_data/dgrp_phenotypes/huang_2015_gene-expression/dgrp.array.exp.male.30pcs.csv', row.names = FALSE)

#########################################
#####Load and compile all phenotypes#####
#########################################
#Set working directory
setwd('00_data/dgrp_phenotypes/')

#List directories
files = list.files()

#Load
dgrp = list()
for(i in 1:length(files)){
  setwd(files[i])
  x = list.files()
  x = x[grep('.csv', x)]
  if(length(x)>0){
    for(j in 1:length(x)){
      d = as.data.frame(read.csv(x[j], row.names = NULL))
      if(ncol(d)>2){
        lines = d[,1]
        n = colnames(d[,2:ncol(d)])
        d = split(t(d[,2:ncol(d)]), 2:ncol(d))
        for(k in 1:length(d)){
          dgrp[[paste(n[k], files[i], sep = '-')]] = as.data.frame(cbind(lines, d[[k]]))
        }
      }else{
        dgrp[[paste(colnames(d)[2], files[i], sep = '-')]] = d
      }
    }
  }
  setwd('../')
}

#Get all unique lines
lines = unique(unlist(lapply(dgrp, function(x) x$line)))
lines = lines[!is.na(lines)]

#Combine
tmp = do.call(cbind, lapply(dgrp, function(x) x[match(lines, x$line),2:ncol(x)]))
rownames(tmp) = lines

#Save
write.csv(tmp, '02_output/all_dgrp_phenotypes_091222.csv')

##################
#####PCA etc.#####
##################
#Read in
n = readLines('02_output/all_dgrp_phenotypes_091222.csv', n = 1)
n = lapply(strsplit(n, ','), function(x) substring(x, 2))
n = lapply(n, function(x) substr(x, 1, nchar(x)-1))
n = unlist(n)[-1]

dgrp = as.data.frame(data.matrix(read.csv('02_output/all_dgrp_phenotypes_091222.csv',
                                          row.names = 1)))
colnames(dgrp) = n

#Filter
dgrp = dgrp[apply(dgrp, 1, function(x) sum(is.na(x)))<50,]
dgrp = dgrp[,apply(dgrp, 2, function(x) sum(is.na(x)))<50]
dgrp = dgrp[,apply(dgrp, 2, function(x) length(unique(x)))>5]

#Impute for PCA (using column mean)
for(i in 1:ncol(dgrp)){
  x = which(is.na(dgrp[,i]))
  dgrp[x,i] = mean(dgrp[,i], na.rm = TRUE)
}

#Scale?
#dgrp_s = as.data.frame(apply(dgrp, 2, function(x) scale(x)))

#Rank normalization
dgrp_s = as.data.frame(apply(dgrp, 2, function(x) RankNorm(x)))

#Correlate
corr = cor(data.matrix(dgrp_s), use = 'complete.obs')
heatmap.2(corr,
          col=rev(brewer.pal(11,"RdBu")),
          scale="row", 
          trace="none")

#Run PCA
pca = prcomp(t(dgrp_s), scale. = TRUE, center = FALSE)

#Plot by study
s = unlist(lapply(strsplit(colnames(dgrp), "-"), function(v){v[2]}))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:length(unique(s))]
cols = WGCNA::standardColors(n = length(unique(s)))
names(cols) = unique(s)
cols = cols[match(s, names(cols))]

plot(pca$x[,1:2], pch = 20, col = cols)

#UMAP
u = umap(t(dgrp_s), verbose = TRUE)

plot(u$layout, 
     pch = 20, 
     cex = 0.8,
     col = alpha(cols, 0.75),
     xlab = 'UMAP 1',
     ylab = 'UMAP 2',
     cex.axis = 1.5,
     cex.lab = 1.5)

y = seq(-2, -8, -0.25)
for(i in 1:length(unique(s))){
  text(-1, y[i], unique(s)[i], cex = 0.5, col = unique(cols)[i])
}

#############################################################
#####Testing linear vs. non-linear components via models#####
#############################################################
#Get all comparisons 
all = expand.grid(1:ncol(dgrp_s), 1:ncol(dgrp_s))

#Convert data to numeric matrix
z = data.matrix(dgrp_s)

#Progress bar
pb <- txtProgressBar(min = 0,      
                     max = nrow(all), 
                     style = 3,    
                     width = 50,
                     char = ".")

#Test linear vs. nonlinear models via AIC for all traits
res = list()
for(i in 1:nrow(all)){
  
  #Update progress bar
  setTxtProgressBar(pb, i)
  
  mod1 = lm(z[,all[i,1]]~z[,all[i,2]])
  mod2 = gam(z[,all[i,1]]~s(z[,all[i,2]], k = length(unique(z[,all[i,2]]))-1))
  
  out = lrtest(mod1, mod2)
  a = AIC(mod1, mod2)
  
  l = list(out, a)
  names(l) = c('lrtest', 'AIC')
  
  res[[paste(colnames(dgrp_s)[all[i,1]],
             colnames(dgrp_s)[all[i,2]], 
             sep = '_')]] = l
}

#Compare with AIC
x = unlist(lapply(res, function(x) x$AIC[2,2]-x$AIC[1,2]))
hist(x[x>(-20)&x<20], breaks = 100)

#Compare with likelihood
x = unlist(lapply(res, function(x) x$lrtest$LogLik[1]-x$lrtest$LogLik[2]))
hist(x[x>(-50)&x<50], breaks = 100)

d = density(x[x>(-50)&x<50])




