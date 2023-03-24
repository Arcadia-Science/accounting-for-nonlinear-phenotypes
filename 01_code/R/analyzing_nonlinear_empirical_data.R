#Set working directory
setwd('~/Documents/Research/github/accounting-for-nonlinear-phenotypes/')

#Source utility functions
source('01_code/R/nonlinear_phenotype_utils.R')

#########################
#####Clean datasets######
#########################
#####Yeast#####
#Load
yeast = read.delim('~/Desktop/nonlinear_phenotype_review/00_data/empirical/yeast/rau_et_al_2020.tsv')[,-1]

#Rank normalize
yeast = apply(yeast, 2, function(x) RankNorm(x))

#Save
saveRDS(yeast, '02_output/yeast_cleaned_phenos.RDS')

#####Cowpea#####
#Load
cowpea = data.matrix(t(readxl::read_xlsx('~/Desktop/nonlinear_phenotype_review/00_data/empirical/cowpea/huynh_et_al_2018.xlsx')[,-1]))

#Convert to numeric
cowpea = apply(cowpea, 2, function(x) as.numeric(x))

#Impute
for(i in 1:ncol(cowpea)){
  x = which(is.na(cowpea[,i]))
  cowpea[x,i] = mean(cowpea[,i], na.rm = TRUE)
}
cowpea = as.data.frame(cowpea)

#Remove duplicated rows
cowpea = cowpea[!duplicated(cowpea), ]

#Rank normalize
cowpea = apply(cowpea, 2, function(x) RankNorm(x))

#Save
saveRDS(cowpea, '02_output/cowpea_cleaned_phenos.RDS')

#####C. elegans#####
#Load
nematode = as.data.frame(readxl::read_xlsx('~/Desktop/nonlinear_phenotype_review/00_data/empirical/nematode/snoek_et_al_2019.xlsx'))

#Select phenotype columns
nematode = nematode[,6:ncol(nematode)]

#Impute
for(i in 1:ncol(nematode)){
  x = which(is.na(nematode[,i]))
  nematode[x,i] = mean(nematode[,i], na.rm = TRUE)
}

#Rank normalize
nematode = apply(nematode, 2, function(x) RankNorm(x))

#Select variable columns
nematode = nematode[,apply(nematode, 2, function(x) length(unique(x)))>5]

#Save
saveRDS(nematode, '02_output/nematode_cleaned_phenos.RDS')

#####AIL (mouse)#####
ail = read.delim('~/Desktop/nonlinear_phenotype_review/00_data/empirical/mouse/ail_phenotypes/ail.phenos.final.txt')[,-1]
ail = ail[,13:ncol(ail)]
for(i in 1:ncol(ail )){
  x = which(is.na(ail [,i]))
  ail [x,i] = mean(ail [,i], na.rm = TRUE)
}
ail = apply(ail, 2, function(x) RankNorm(x))
ail = ail[,apply(ail, 2, function(x) length(unique(x)))>5]
saveRDS(ail, '02_output/ail_cleaned_phenos.RDS')

#####QTL archive (mouse and rat)#####
setwd('~/Desktop/nonlinear_phenotype_review/00_data/empirical/mouse/qtlarchive_phenotypes/')
files = list.files()
qtlarchive = list()
for(i in 1:length(files)){
  tmp = read.csv(files[i])
  tmp = apply(tmp, 2, function(x) gsub('-', 'NA', x))
  tmp = as.data.frame(tmp)
  for(j in 1:ncol(tmp)){
    tmp[,j] = as.numeric(tmp[,j])
    x = which(is.na(as.numeric(tmp [,j])))
    tmp [x,j] = mean(as.numeric(tmp [,j]), na.rm = TRUE)
  }
  tmp = apply(tmp, 2, function(x) RankNorm(x))
  tmp = tmp[,apply(tmp, 2, function(x) length(unique(x)))>5]
  qtlarchive[[files[i]]] = tmp
}

setwd('~/Documents/Research/github/accounting-for-nonlinear-phenotypes/')
saveRDS(qtlarchive, '02_output/qtlarchive_cleaned_phenos.RDS')

#########################
#####Load clean data#####
#########################
#Yeast
yeast = readRDS('02_output/empirical-phenotypes/yeast_cleaned_phenos.RDS')

#Cowpea
cowpea = readRDS('02_output/empirical-phenotypes/cowpea_cleaned_phenos.RDS')

#C. elegans
nematode = readRDS('02_output/empirical-phenotypes/nematode_cleaned_phenos.RDS')

#Arabadopsis
ara = readRDS('02_output/empirical-phenotypes/arapheno_cleaned_phenos.RDS')

#DGRP
dgrp = readRDS('02_output/empirical-phenotypes/dgrp_cleaned_phenos.RDS')

#JAX
jax = readRDS('02_output/empirical-phenotypes/jax_cleaned_phenos.RDS')

#AIL
ail = readRDS('02_output/empirical-phenotypes/ail_cleaned_phenos.RDS')

#QTL archive
qtlarch = readRDS('02_output/empirical-phenotypes/qtlarchive_cleaned_phenos.RDS')

#Compile
phenos = list(yeast,
              cowpea,
              nematode,
              ara[[3]],
              dgrp,
              jax,
              ail)
names(phenos) = c('yeast', 'cowpea', 'nematode', 'ara', 'dgrp', 'jax', 'ail')

################################
#####Calculate nonlinearity#####
################################
#Run permutation test for each dataset (skipping cowpea, not enough phenos to permute on)
toTest = c(1, 3:length(phenos))
nonlinear_permutations = list()
for(i in toTest){
  print(paste(i, 'out of', length(phenos)))
  nonlinear_permutations[[names(phenos)[i]]] = nonlinear.permutation(phenos[[i]], data_proportion = 0.25, permutation_number = 100)
}

#Calculate random distribution
dat = list()
for(i in 1:20){
  dat[[i]] = sample(1:1000, 600, replace = TRUE)
}
dat = do.call(cbind, dat)

#Calculate nonlinearity and entropy of random
random = nonlinear.permutation(dat, data_proportion = 0.25, permutation_number = 100)

#Add
nonlinear_permutations$random = random

#Plot
nonlinear_permutations = nonlinear_permutations[order(unlist(lapply(nonlinear_permutations, function(x) mean(x))), decreasing = TRUE)]
cols = c(arcadia.pal(n = 6, name = 'Accent'), 'gray50')
names(cols) = names(nonlinear_permutations)

vioplot::vioplot(nonlinear_permutations[[1]],
                 nonlinear_permutations[[2]],
                 nonlinear_permutations[[3]],
                 nonlinear_permutations[[4]],
                 nonlinear_permutations[[5]],
                 nonlinear_permutations[[6]],
                 nonlinear_permutations[[7]],
                 col = cols,
                 side = "right",
                 border = darken_color(cols),
                 ylab = '% nonlinear',
                 xlab = '', 
                 las = 2, 
                 names = names(nonlinear_permutations),
                 cex.axis = 1.5,
                 cex.lab = 1.5)

stripchart(nonlinear_permutations,
           col = cols,
           at = seq(0.8, (length(nonlinear_permutations)-1)+0.8, 1), 
           jitter = 0.1,
           method = "jitter", 
           vertical = TRUE, 
           cex = 1,
           pch = 20, 
           add = TRUE)

#Calculate overall for each
nonlinear = list()
for(i in toTest){
  print(paste(i, 'out of', length(phenos)))
  nonlinear[[names(phenos)[i]]] = compute.phenotype.stats(list(phenos[[i]]), run_nonlinear = TRUE, sample_sizes = seq(0.1, 0.9, 0.1))
}

#Calculate for random distribution
random_stats = compute.phenotype.stats(list(dat), run_nonlinear = TRUE, sample_sizes = seq(0.1, 0.9, 0.1))

#Plot entropy distributions
plot(seq(0.1, 0.9, 0.1),
     colMeans(do.call(cbind, nonlinear$yeast$subsampled_entropy_results[[1]]$entropies)), 
     type = 'l',
     ylim = c(0, 7),
     col = cols[1],
     lwd = 1.5,
     ylab = 'Entropy (bits)',
     xlab = 'Proportion of phenotypes',
     cex.axis = 1.5,
     cex.lab = 1.5)
for(i in 2:length(nonlinear)){
  lines(seq(0.1, 0.9, 0.1),
        colMeans(do.call(cbind, nonlinear[[i]]$subsampled_entropy_results[[1]]$entropies)), col = cols[i], lwd = 1.5)
}

#Plot entropy as a function of phenotype number
e = unlist(lapply(nonlinear, function(x) x$subsampled_entropy_slopes))
names(e) = names(nonlinear)
o = order(e, decreasing = TRUE)
e = e[o]

plot(e,
     pch = 20,
     cex = 3,
     col = cols[match(names(e), names(cols))],
     ylim = c(0, 0.35),
     cex.lab = 1.5,
     cex.axis = 1.5,
     xaxt = 'n',
     xlab = '',
     ylab = 'Entropy fit (slope)',
     bty = 'n')
axis(1, 
     1:length(e), names(e),
     cex.axis = 1.5,
     las = 2)
for(i in 1:length(e)){
  segments(i, 0, i, e[i], col = cols[match(names(e), names(cols))][i], lwd = 1.5)
}



#Plot entropy compared to nonlinearity
plot(unlist(lapply(nonlinear, function(x) x$nonlinearity_aic_ratios)),
     unlist(lapply(nonlinear, function(x) x$subsampled_entropy_slopes)))

###################################
#####Compare to synthetic data#####
###################################
#Load
all_stats = readRDS('~/Desktop/nonlinear_phenotype_review/01_output/all_synthetic_phenos_all_stats_no_models_030923.RDS')
p_int = paste(rep('int', length(seq(0, 1, 0.01))), seq(0, 1, 0.01), sep = '')

#Split
all_stats_int = split(all_stats, p_int)

#Smooth AIC ratios
a_smooth = lapply(all_stats_int, function(x) {
  z = unlist(lapply(x, function(y) y$nonlinearity_aic_ratios))
  smooth.spline(1:length(z), z, spar = 0.75)$y
})

#Plot
slopes = unlist(lapply(all_stats_int, function(x) lapply(x, function(y) y$subsampled_entropy_slopes)))

plot(unlist(a_smooth),
     unlist(slopes),
     xlab = '% nonlinear',
     ylab = 'Entropy fit (slope)',
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex = 0.75,
     pch = 20,
     col = 'gray80')
points(unlist(lapply(nonlinear, function(x) x$nonlinearity_aic_ratios)),
       unlist(lapply(nonlinear, function(x) x$subsampled_entropy_slopes)),
       bg = cols[order(unlist(lapply(nonlinear, function(x) x$nonlinearity_aic_ratios)), decreasing = TRUE)],
       col = darken_color(cols[order(unlist(lapply(nonlinear, function(x) x$nonlinearity_aic_ratios)), decreasing = TRUE)]),
       cex = 2,
       pch = 21)

#Number of populations/crosses
#dgrp = 1
#jax = 1?
#yeast = 15
#ail = 2
#nematode = 4
#arabidopsis = 517?

plot(log10(c(15, 4, 517, 1, 1, 2, 1)),
     unlist(lapply(nonlinear, function(x) x$nonlinearity_aic_ratios)),
     ylab = '% nonlinear',
     xlab = 'log10 (n populations)',
     cex.axis = 1.5,
     cex.lab = 1.5,
     pch = 21,
     bg = cols[order(unlist(lapply(nonlinear, function(x) x$nonlinearity_aic_ratios)), decreasing = TRUE)],
     col = darken_color(cols[order(unlist(lapply(nonlinear, function(x) x$nonlinearity_aic_ratios)), decreasing = TRUE)]),
     cex = 2)
abline(lm(unlist(lapply(nonlinear, function(x) x$nonlinearity_aic_ratios))~log10(c(15, 4, 517, 1, 1, 2, 1))),
       lty = 'dashed')











#Bacteria
#setwd('~/Desktop/nonlinear_phenotype_review/00_data/bacteria_phenotypes/fitness_browser_data/')
#files = list.files()
#phenos = list()
#for(i in 1:length(files)){
  print(paste(i, 'out of', length(files)))
  x = read.delim(files[i])
  x = x[x$used == TRUE,]
  image(cov(x[,6:ncol(x)]))
  phenos[[as.character(files[i])]] = as.data.frame(x[,7:ncol(x)])
}

#Remove huge experiment (very slow to calculate nonlinearity on)
#phenos = phenos[-5]

#bacteria = list()
#for(i in 1:length(phenos)){
  print(paste(i, 'out of', length(phenos)))
  bacteria[[names(phenos)[i]]] = compute.phenotype.stats(list(phenos[[i]][sample(1:nrow(phenos[[i]]), 1000),]), run_nonlinear = TRUE)
}

#Plot
plot(unlist(a_smooth),
     unlist(slopes),
     xlab = '% nonlinear',
     ylab = 'Entropy fit (slope)',
     cex.lab = 1.5,
     cex.axis = 1.5,
     xlim = c(0.2, 1),
     ylim = c(0, 0.4),
     cex = 0.5,
     pch = 20,
     col = alpha('gray50', 0.25))
points(unlist(lapply(ara, function(x) x$nonlinearity_aic_ratios)),
       unlist(lapply(ara, function(x) x$subsampled_entropy_slopes)),
       pch = 20, cex = 2, col = 'gold4')
points(dgrp$nonlinearity_aic_ratios,
       dgrp$subsampled_entropy_slopes,
       pch = 20, cex = 2, col = 'darkred')
points(jax$nonlinearity_aic_ratios,
       jax$subsampled_entropy_slopes,
       pch = 20, cex = 2, col = 'cyan4')
points(unlist(lapply(bacteria, function(x) x$nonlinearity_aic_ratios)),
       unlist(lapply(bacteria, function(x) x$subsampled_entropy_slopes)),
       pch = 20, cex = 2, col = 'darkgreen')

plot(unlist(a_smooth),
     unlist(mi),
     xlab = '% nonlinear',
     ylab = 'Mean mutual information',
     cex.lab = 1.5,
     cex.axis = 1.5,
     xlim = c(0.2, 1),
     cex = 0.5,
     pch = 20,
     col = alpha('gray50', 0.25))
points(unlist(lapply(ara, function(x) x$nonlinearity_aic_ratios)),
       unlist(lapply(ara, function(x) x$mutual_information)),
       pch = 20, cex = 2, col = 'gold4')
points(dgrp$nonlinearity_aic_ratios,
       dgrp$mutual_information,
       pch = 20, cex = 2, col = 'darkred')
points(jax$nonlinearity_aic_ratios,
       jax$mutual_information,
       pch = 20, cex = 2, col = 'cyan4')
points(unlist(lapply(bacteria, function(x) x$nonlinearity_aic_ratios)),
       unlist(lapply(bacteria, function(x) x$mutual_information)),
       pch = 20, cex = 2, col = 'darkgreen')

