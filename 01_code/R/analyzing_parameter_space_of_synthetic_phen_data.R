#p_pleio and p_int had values of: [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#The organization of the data is: data[p_pleio][p_int][phens][individuals]

#Set working directory
setwd('~/Documents/Research/github/accounting-for-nonlinear-phenotypes/')

#Source utility functions
source('01_code/R/nonlinear_phenotype_utils.R')

#Source python script for generating synthetic genotype/phenotype data
source_python('01_code/python/tools_for_phen_gen_creation.py')

#######################################
#####Load synthetic phenotype data#####
#######################################
#Load pickled phenotype data
library(reticulate)
pd <- import("pandas")
phens <- pd$read_pickle("02_output/ppleio_pint_sweep/phen_pleio_int.pk")

#Generate naming vectors
p_pleio = paste(rep('pleio', length(seq(0, 1, 0.1))), seq(0, 1, 0.1), sep = '')
p_int = paste(rep('int', length(seq(0, 1, 0.1))), seq(0, 1, 0.1), sep = '')

#Combine phenotypes into matrices
phenos = list()
for(i in 1:length(phens)){
  for(j in 1:length(phens[[i]])){
    phenos[[paste(p_pleio[i], p_int[j], sep = '_')]] = do.call(cbind, phens[[i]][[j]])
  }
}

########################################
#####Compare phenotype correlations#####
########################################
#Calculate mean correlations
corrs = list()
for(i in 1:length(phenos)){
  x = cor(phenos[[i]])
  x = x[!x == 1]
  corrs[[i]] = mean(x)
}

#Combine into matrix
mat = matrix(unlist(corrs), nrow = 11)

#Plot as heatmap
image(mat, col = colorRamps::matlab.like2(1000))

#Plot as lines
plot(mat[1,], type = 'l')
for(i in 2:nrow(mat)){
  lines(mat[i,])
}

###################################
#####Test linear vs. nonlinear#####
###################################
#Compare linear vs. nonlinear
counter <- 0
pb <- txtProgressBar(min = 1,      
                     max = length(phenos), 
                     style = 3,    
                     width = 50,
                     char = ".")
res = lapply(phenos, function(x){
  counter <<- counter + 1
  setTxtProgressBar(pb, counter)
  compare_nonlinear(x, return_models = TRUE)})

#Save
saveRDS(res, '02_output/ppleio_pint_sweep/ppleio_pint_sweep_linear_nonlinear_ratios_20phenos.RDS')

#Extract nonlinearity rates
lin = unlist(lapply(res, function(x) x$linear))

#Combine into matrix
mat = matrix(unlist(lin), nrow = 11)

#Plot as heatmap
image(mat, col = colorRamps::matlab.like2(1000))

#Plot as lines
par(mfrow = c(1,2))

lin = unlist(lapply(res, function(x) x$linear))
mat = matrix(unlist(lin), nrow = 11)

cols = colorRampPalette(arcadia.pal(6, 'Accent'))(length(seq(0, 1, 0.1)))
plot(mat[1,], 
     type = 'l', 
     ylim = c(0,1), 
     col = cols[1], 
     lwd = 1.5,
     cex.axis = 1.5,
     cex.lab = 1.5,
     xaxt = 'n',
     xlab = 'Probability pleiotropy',
     ylab = 'Linearity')
axis(1, seq(1, length(seq(0, 1, 0.1)), 1), seq(0, 1, 0.1), cex.axis = 1.5)
for(i in 2:nrow(mat)){
  lines(mat[i,], col = cols[i], lwd = 1.5)
}

lin = unlist(lapply(res, function(x) x$nonlinear))
mat = matrix(unlist(lin), nrow = 11)

cols = colorRampPalette(arcadia.pal(6, 'Accent'))(length(seq(0, 1, 0.1)))
plot(mat[1,], 
     type = 'l', 
     ylim = c(0,1), 
     col = cols[1], 
     lwd = 1.5,
     cex.axis = 1.5,
     cex.lab = 1.5,
     xaxt = 'n',
     xlab = 'Probability pleiotropy',
     ylab = 'Nonlinearity')
axis(1, seq(1, length(seq(0, 1, 0.1)), 1), seq(0, 1, 0.1), cex.axis = 1.5)
for(i in 2:nrow(mat)){
  lines(mat[i,], col = cols[i], lwd = 1.5)
}



#Using MIC
counter <- 0
pb <- txtProgressBar(min = 1,      
                     max = length(phenos), 
                     style = 3,    
                     width = 50,
                     char = ".")
res = lapply(phenos, function(x){
  counter <<- counter + 1
  setTxtProgressBar(pb, counter)
  compare_nonlinear_MIC(x[,1:20])})

r2 = lapply(res, function(x) unlist(x)[1])
mat = matrix(unlist(r2), nrow = 11)

#Plot as lines
plot(mat[1,], type = 'l', ylim = c(-1,1))
for(i in 2:nrow(mat)){
  lines(mat[i,])
}