#p_pleio and p_int had values of: [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#The organization of the data is: data[p_pleio][p_int][phens][individuals]

#Set working directory
setwd('~/Documents/Research/github/accounting-for-nonlinear-phenotypes/')

#Source utility functions
source('01_code/R/nonlinear_phenotype_utils.R')

#######################################
#####Load synthetic phenotype data#####
#######################################
##0-1 by 0.01
#Load pickled phenotype data
phens <- pd$read_pickle("~/Desktop/nonlinear_phenotype_review/00_data/synthetic/phen_pleio_int_01_0_1.pk")

#Combine phenotypes into matrices
phenos = list()

#Generate naming vectors
p_pleio = paste(rep('pleio', length(seq(0, 1, 0.01))), seq(0, 1, 0.01), sep = '')
p_int = paste(rep('int', length(seq(0, 1, 0.01))), seq(0, 1, 0.01), sep = '')

#Combine phenotypes into matrices
for(i in 1:length(phens)){
  for(j in 1:length(phens[[i]])){
    phenos[[paste(p_pleio[i], p_int[j], sep = '_')]] = do.call(cbind, phens[[i]][[j]])
  }
}

#######################################################
#####Calculate entropy distributions for p_int = 0#####
#######################################################
#Entropy by data fraction
p = split(phenos, p_int)

#Calculate for each study
entropies = list()
for(i in 1:length(phenos)){
  print(paste(i, 'out of', length(phenos)))
  x = cov(phenos[[i]]/max(phenos[[i]]))
  e = eigen(x)
  
  l = list(entropy::entropy.empirical(e$values, unit = 'log2'),
           entropy::entropy.empirical(e$values, unit = 'log2')/log2(length(e$values)),
           det(x),
           ncol(x))
  names(l) = c('entropy', 'normalized_entropy', 'determinant', 'n_phenotypes')
  entropies[[i]] = l
}

#Split
entropies = split(entropies, p_int)

#Plot
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))(length(entropies)))
plot(unlist(lapply(entropies[[1]], function(x) x$entropy)), 
     type = 'l', 
     col = cols[1])
for(i in 2:length(entropies)){
  lines(unlist(lapply(entropies[[i]], function(x) x$entropy)),
        col = cols[i])
}
lapply(entropies, function(x) lines(unlist(lapply(x, function(y) y$entropy)), lwd = 0.5, col = alpha('gray50', 0.5)))
  
#Entropy by data fraction
p = split(phenos, p_int)

all_entropies = list()
for(i in 1:length(p$int0)){
  print(names(p$int0)[i])
  all_entropies[[names(p$int0)[i]]] = subsample.entropy(p$int0[[i]],
                                                        numBins = 100,
                                                        permutations = 1000,
                                                        normalize_entropy = FALSE,
                                                        sample_sizes = seq(0.1, 0.9, 0.05))}

#Plot
layout(matrix(1:3,ncol=3), width = c(2,1,2),height = c(1,1,1))

cols = rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))(length(all_entropies)))
plot(colMeans(do.call(cbind, all_entropies[[1]]$entropies)), 
     type = 'l',
     ylim = c(0, 4.5),
     col = cols[1],
     ylab = 'Entropy (bits)',
     xlab = 'n phenotypes',
     cex.axis = 1.5,
     cex.lab = 1.5)
for(i in 2:length(all_entropies)){
  lines(colMeans(do.call(cbind, all_entropies[[i]]$entropies)), col = cols[i])
}

legend_image <- as.raster(matrix(cols, ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pleiotropy', font.main = 1)
text(x=1.5, y = rev(seq(0,1,l=5)), labels = seq(0,1,l=5), pos = 1)
rasterImage(legend_image, 0, 0, 1,1)

#Fit slope
slopes = c()
for(i in 1:length(all_entropies)){
  y = colMeans(do.call(cbind, all_entropies[[i]]$entropies))
  slopes = c(slopes, coef(lm(y~seq(1, length(y), 1)))[2])
}

#plot(slopes, type = 'l')

#Correlate slope and pleiotropy
a = smooth.spline(seq(0, 1, 0.01),
                  slopes, 
                  spar = 0.75)$y
plot(seq(0, 1, 0.01),
     slopes,
     pch = 20,
     cex = 3,
     col = cols,
     cex.axis = 1.5,
     cex.lab = 1.5,
     xlab = 'Probability pleiotropy',
     ylab = 'Entropy (slope)')
lines(seq(0, 1, 0.01), a, lwd = 1.5)
text(0.9, 
     0.17, 
     paste('cor =', round(cor(seq(0, 1, 0.01), slopes), 2)),
     cex = 1.5)

######################################
#####Calculate stats for all data#####
######################################
#Run
all_stats = list()
for(i in 1:length(phenos)){
  print(paste(i, 'out of', length(phenos)))
  all_stats[[names(phenos)[i]]] = compute.phenotype.stats(list(phenos[[i]]), run_nonlinear = TRUE)
}

#Save
saveRDS(all_stats, '~/Desktop/all_synthetic_phenos_all_stats_030923.RDS')
saveRDS(lapply(all_stats, function(x) x[1:7]), '~/Desktop/all_synthetic_phenos_all_stats_no_models_030923.RDS')

#Split
all_stats_int = split(lapply(all_stats, function(x) x[1:7]), p_int)
all_stats_pleio = split(lapply(all_stats, function(x) x[1:7]), unlist(lapply(strsplit(names(all_stats), "_"), function(v){v[1]})))

##########################################
#####Plot non-linearity distributions#####
##########################################
#Smooth AIC ratios
a_smooth = lapply(all_stats_int, function(x) {
  z = unlist(lapply(x, function(y) y$nonlinearity_aic_ratios))
  smooth.spline(1:length(z), z, spar = 0.75)$y
})

#Plot
layout(matrix(1:3,ncol=3), width = c(2,1,2),height = c(1,1,1))

cols = rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))(length(a_smooth)))
plot(a_smooth[[1]], 
     type = 'l',
     ylim = c(0, 1),
     col = cols[1],
     ylab = '% nonlinear',
     xlab = 'Pleiotropy (probability)',
     cex.axis = 1.5,
     cex.lab = 1.5)
abline(h = 0.5, lwd = 1.5, lty = 'dashed')
for(i in 2:length(a_smooth)){
  lines(a_smooth[[i]], col = cols[i])
}
lines(a_smooth[[1]], lwd = 2, col = cols[1])

legend_image <- as.raster(matrix(cols, ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'P (gene interaction)', font.main = 1)
text(x=1.5, y = rev(seq(0,1,l=5)), labels = seq(0,1,l=5), pos = 1)
rasterImage(legend_image, 0, 0, 1,1)

#Plot entropy and non-linearity
slopes = unlist(lapply(all_stats_int, function(x) lapply(x, function(y) y$subsampled_entropy_slopes)))

cols = rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))(length(unlist(a_smooth))))
plot(unlist(a_smooth),
     unlist(slopes),
     xlab = '% nonlinear',
     ylab = 'Entropy fit (slope)',
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex = 0.75,
     pch = 20,
     col = cols)

##############################################
#####Entropy and non-linearity landscapes#####
##############################################
#Extract pleio and int as for x and y coords
x = unlist(lapply(strsplit(names(all_stats), "_"), function(v){v[1]}))
x = gsub('pleio', '', x)
x = as.numeric(x)

y = unlist(lapply(strsplit(names(all_stats), "_"), function(v){v[2]}))
y = gsub('int', '', y)
y = as.numeric(y)

#Nonlinearity
h = unlist(lapply(all_stats, function(x) x$nonlinearity_aic_ratios))
m = gam(h ~ s(x, y, k=30))
fit <- predict(m)
vis.gam(m, 
        se=F,
        xlab = 'Prob pleiotropy',
        ylab = 'Prob gene interactions',
        plot.type = "contour", 
        cex.lab = 1.5,
        cex.axis = 1.5,
        color = 'topo',
        main = '% non-linear',
        font.main = 1,
        cex.main = 1.5,
        type="response")

vis.gam(m, 
        se=F,
        xlab = 'Prob pleiotropy',
        ylab = 'Prob gene interactions',
        plot.type = "persp", 
        cex.lab = 1.5,
        cex.axis = 1.5,
        color = 'topo',
        main = '',
        zlab = '% non-linear',
        font.main = 1,
        cex.main = 1.5,
        theta = 210, 
        n.grid = 50,
        phi=40,
        axes = TRUE, 
        ticktype = "detailed") 

#Entropy
h = unlist(lapply(all_stats, function(x) x$subsampled_entropy_slopes))
m = gam(h ~ s(x, y, k=30))
fit <- predict(m)
myvis.gam(m, 
        se=F,
        xlab = 'Prob pleiotropy',
        ylab = 'Prob gene interactions',
        plot.type = "contour",
        color = 'ylgnbu', 
        cex.lab = 1.5,
        cex.axis = 1.5,
        main = 'Entropy',
        font.main = 1,
        cex.main = 1.5,
        type="response") 

myvis.gam(m, 
        se=F,
        xlab = 'Prob pleiotropy',
        ylab = 'Prob gene interactions',
        plot.type = "persp",
        color = 'ylgnbu',
        cex.lab = 1.5,
        cex.axis = 1.5,
        main = '',
        zlab = 'Entropy',
        font.main = 1,
        cex.main = 1.5,
        theta = 140, 
        n.grid = 50,
        phi=40,
        axes = TRUE, 
        ticktype = "detailed") 

#Mutual information
h = unlist(lapply(all_stats, function(x) x$mutual_information))
m = gam(h ~ s(x, y, k=30))
fit <- predict(m)
vis.gam(m, 
        se=F,
        xlab = 'Prob pleiotropy',
        ylab = 'Prob gene interactions',
        plot.type = "contour", 
        cex.lab = 1.5,
        cex.axis = 1.5,
        color = 'topo',
        main = 'Mutual information',
        font.main = 1,
        cex.main = 1.5,
        type="response")

vis.gam(m, 
        se=F,
        xlab = 'Prob pleiotropy',
        ylab = 'Prob gene interactions',
        plot.type = "persp", 
        cex.lab = 1.5,
        cex.axis = 1.5,
        color = 'topo',
        main = 'Mutual information',
        font.main = 1,
        cex.main = 1.5,
        theta = 210,
        phi=40,
        too.far=.07,
        type="response") 

#Determinant
h = unlist(lapply(all_stats, function(x) x$determinant))
m = gam(h ~ s(x, y, k=30))
fit <- predict(m)
vis.gam(m, 
        se=F,
        xlab = 'Prob pleiotropy',
        ylab = 'Prob gene interactions',
        plot.type = "contour", 
        cex.lab = 1.5,
        cex.axis = 1.5,
        color = 'topo',
        main = 'Determinant',
        font.main = 1,
        cex.main = 1.5,
        type="response")

vis.gam(m, 
        se=F,
        xlab = 'Prob pleiotropy',
        ylab = 'Prob gene interactions',
        plot.type = "persp", 
        cex.lab = 1.5,
        cex.axis = 1.5,
        color = 'topo',
        main = 'Determinant',
        font.main = 1,
        cex.main = 1.5,
        theta = 210,
        phi=40,
        nCol = 100,
        type="response") 

############################
#####Mutual information#####
############################
#Calculate mutual information
mi = lapply(phenos, function(x){
  m = infotheo::mutinformation(infotheo::discretize(x))
  mean(m[lower.tri(m)])
})

#################################
#####Calculate non-linearity#####
#################################
#Run function on all experiments
res = list()
for(i in 1:length(phenos)){
  print(paste(i, 'out of', length(phenos)))
  res[[as.character(i)]] = compare_nonlinear(phenos[[i]], verbose = TRUE, return_models = TRUE)
  
  a = unlist(lapply(res[[as.character(i)]]$models, function(x) abs(x$AIC[2,2])-abs(x$AIC[1,2])))
  print(sum(a<0)/length(a))
}

#Save
saveRDS(res, '~/Desktop/synthetic_data_nonlinear_models.RDS')

#Calculate AIC differences
a = c()
for(i in 1:length(res)){
  tmp = unlist(lapply(res[[i]]$models, function(x) abs(x$AIC[2,2])-abs(x$AIC[1,2])))
  a = c(a, sum(tmp<0)/length(tmp))
}

#Split
a = split(a, p_int)

#Smooth
a_smooth = lapply(a, function(x) smooth.spline(1:length(x), x, spar = 0.5)$y)

#Plot
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))(length(a_smooth)))
plot(a_smooth[[1]], 
     type = 'l',
     ylim = c(0, 1),
     col = cols[1],
     ylab = '% nonlinear',
     xlab = 'Pleiotropy (probability)',
     cex.axis = 1.5,
     cex.lab = 1.5)
abline(h = 0.5, lwd = 1.5, lty = 'dashed')
for(i in 2:length(a_smooth)){
  lines(a_smooth[[i]], col = cols[i])
}
lines(a_smooth[[1]], lwd = 2, col = cols[1])

###########################################
#####Compare entropy and non-linearity#####
###########################################
#Calculate subsampled entropy
full_entropies = lapply(p, function(x){
  lapply(x, function(z){
    subsample.entropy(z,
                      numBins = 100,
                      permutations = 100,
                      normalize_entropy = FALSE,
                      sample_sizes = round(seq(0.1, 0.9, 0.1)*ncol(z)))
  })
})

#Fit slope
slopes = lapply(full_entropies, function(x){
  lapply(x, function(z){
    y = colMeans(do.call(cbind, z$entropies))
    coef(lm(y~seq(1, length(y), 1)))[2]
  })
})

#Plot
m = as.matrix(do.call(rbind, lapply(slopes, function(x) unlist(x))))
image(m)

#Compare subsampled entropy to nonlinearity
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))(length(unlist(a_smooth))))
plot(unlist(a_smooth),
     unlist(slopes),
     xlab = '% nonlinear',
     ylab = 'Entropy fit (slope)',
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex = 0.5,
     pch = 20,
     col = alpha(cols, 0.25))

#Compare total entropy to nonlinearity
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))(length(unlist(a_smooth))))
plot(unlist(a_smooth),
     unlist(lapply(entropies, function(x) lapply(x, function(y) y$entropy))),
     xlab = '% nonlinear',
     ylab = 'Entropy (bits)',
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex = 0.5,
     pch = 20,
     col = alpha(cols, 0.25))


plot(a_smooth[[1]], unlist(slopes[[1]]))
