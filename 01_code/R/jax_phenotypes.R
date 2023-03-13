setwd('~/Documents/Research/github/accounting-for-nonlinear-phenotypes/')
source('01_code/R/nonlinear_phenotype_utils.R')

#######################################
#####Load and clean phenotype data#####
#######################################
#Set working directory
setwd('~/Desktop/accounting-for-nonlinear-phenotypes/00_data/jax_phenotypes/')

#Load
dat = read.csv('strainmeans.csv')

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

#Filter strains
phenos = phenos[apply(phenos, 1, function(x) sum(!is.na(x)))>100,]

#Filter phenotypes
phenos = phenos[,apply(phenos, 2, function(x) sum(!is.na(x)))>100]

#Filter strains
phenos = phenos[apply(phenos, 1, function(x) sum(!is.na(x)))>100,]

#Save
write.csv(phenos, '~/Documents/Research/github/accounting-for-nonlinear-phenotypes/02_output/all_jax_phenotypes.csv')

##################
#####PCA etc.#####
##################
#Read in
phenos = read.csv('~/Documents/Research/github/accounting-for-nonlinear-phenotypes/02_output/all_jax_phenotypes.csv')

#Impute for PCA (using column mean)
for(i in 1:ncol(phenos)){
  x = which(is.na(phenos[,i]))
  phenos[x,i] = mean(phenos[,i], na.rm = TRUE)
}

#Rank normalization
phenos_s = as.data.frame(apply(phenos, 2, function(x) RankNorm(x)))

#Correlate
corr = cor(data.matrix(phenos_s), use = 'complete.obs')
heatmap.2(corr,
          col=rev(brewer.pal(11,"RdBu")),
          scale="row", 
          trace="none")

#Run PCA
pca = prcomp(t(phenos_s), scale. = TRUE, center = FALSE)
plot(pca$x[,1:2], pch = 20)

#UMAP
u = umap(t(phenos_s), verbose = TRUE)

plot(u$layout, 
     pch = 20, 
     cex = 0.8,
     xlab = 'UMAP 1',
     ylab = 'UMAP 2',
     cex.axis = 1.5,
     cex.lab = 1.5)

#Save 
saveRDS(phenos_s, '02_output/jax_cleaned_phenos.RDS')

#############################################################
#####Testing linear vs. non-linear components via models#####
#############################################################
#Get all comparisons 
all = expand.grid(1:ncol(phenos_s), 1:ncol(phenos_s))

#Convert data to numeric matrix
z = data.matrix(phenos_s)

#Test linear vs. nonlinear models via AIC for all traits
res = list()

#Progress bar
pb <- txtProgressBar(min = 0,      
                     max = nrow(all), 
                     style = 3,    
                     width = 50,
                     char = ".")
for(i in 1:nrow(all)){
  
  #Update progress bar
  setTxtProgressBar(pb, i)
  
  mod1 = lm(z[,all[i,1]]~z[,all[i,2]])
  mod2 = gam(z[,all[i,1]]~s(z[,all[i,2]], k = length(unique(z[,all[i,2]]))-1))
  
  out = lrtest(mod1, mod2)
  a = AIC(mod1, mod2)
  
  l = list(out, a)
  names(l) = c('lrtest', 'AIC')
  
  res[[paste(colnames(phenos_s)[all[i,1]],
             colnames(phenos_s)[all[i,2]], 
             sep = '_')]] = l
}

#Compare with AIC
x = unlist(lapply(res, function(x) x$AIC[2,2]-x$AIC[1,2]))
hist(x[x>(-20)&x<20], breaks = 100)

#Compare with likelihood
x = unlist(lapply(res, function(x) x$lrtest$LogLik[1]-x$lrtest$LogLik[2]))
hist(x[x>(-50)&x<50], breaks = 100)

d = density(x[x>(-50)&x<50])





