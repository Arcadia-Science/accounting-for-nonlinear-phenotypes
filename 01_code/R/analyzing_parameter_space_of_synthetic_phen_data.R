#p_pleio and p_int had values of: [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#The organization of the data is: data[p_pleio][p_int][phens][individuals]

#Set working directory
setwd('~/Documents/Research/github/accounting-for-nonlinear-phenotypes/')

#Source utility functions
source('01_code/R/nonlinear_phenotype_utils.R')

#######################################
#####Load synthetic phenotype data#####
#######################################
##0.01-0.1
#Load pickled phenotype data
phens <- pd$read_pickle("02_output/ppleio_pint_sweep/0_0.1_sweep/phen_pleio_int_01.pk")

#Generate naming vectors
p_pleio = paste(rep('pleio', length(seq(0, 0.1, 0.01))), seq(0, 0.1, 0.01), sep = '')
p_int = paste(rep('int', length(seq(0, 0.1, 0.01))), seq(0, 0.1, 0.01), sep = '')

#Combine phenotypes into matrices
phenos = list()

#Combine phenotypes into matrices
for(i in 1:length(phens)){
  for(j in 1:length(phens[[i]])){
    phenos[[paste(p_pleio[i], p_int[j], sep = '_')]] = do.call(cbind, phens[[i]][[j]])
  }
}

##0.1-1
#Load pickled phenotype data
phens <- pd$read_pickle("02_output/ppleio_pint_sweep/0_1_sweep/phen_pleio_int.pk")

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

##0-1 by 0.01
#Load pickled phenotype data
phens <- pd$read_pickle("~/Desktop/nonlinear_phenotype_review/00_data/phen_pleio_int_01_0_1.pk")

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

#Filter?
x = paste(rep('pleio', length(seq(0, 1, 0.02))), seq(0, 1, 0.02), sep = '')
y = paste(rep('int', length(seq(0, 1, 0.02))), seq(0, 1, 0.02), sep = '')
toMatch = expand.grid(x, y)
toMatch = paste(toMatch[,1], toMatch[,2], sep = '_')
phenos = phenos[match(toMatch, names(phenos))]

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

#Split on pleio value
names(corrs) = names(phenos)
mat = split(corrs, unlist(lapply(strsplit(names(phenos), '_'), function(x) x[1])))

#Combine into matrix
mat = matrix(unlist(corrs), nrow = 22)

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
                     #width = 100,
                     char = ".")
system.time({
res = pbmclapply(phenos, mc.cores = detectCores()-3, function(x){
  counter <<- counter + 1
  setTxtProgressBar(pb, counter)
  compare_nonlinear(x)})
})

#Save
#saveRDS(res, '02_output/ppleio_pint_sweep/ppleio_pint_sweep_linear_nonlinear_ratios_allphenos.RDS')
#saveRDS(res, '02_output/ppleio_pint_sweep/ppleio_pint_expandedsweep_linear_nonlinear_ratios_20phenos.RDS')
#saveRDS(res, '02_output/ppleio_pint_sweep/ppleio_pint_0-1-0.02_linear_nonlinear_ratios_phenos.RDS')

#Load
res = readRDS('02_output/ppleio_pint_sweep/ppleio_pint_sweep_linear_nonlinear_ratios_allphenos.RDS')

#Extract nonlinearity rates
lin = unlist(lapply(res, function(x) x$linear))

#Split on pleio
#p = unlist(lapply(strsplit(names(lin), '_'), function(x) x[1]))
#lin = split(lin, p)

#Combine into matrix
mat = matrix(unlist(lin), nrow = 11)
#mat = do.call(rbind, lin)

#Plot as heatmap
image(mat, col = colorRamps::matlab.like2(1000))

#Plot as lines
par(mfrow = c(1,2))

lin = unlist(lapply(res, function(x) x$linear))
p = unlist(lapply(strsplit(names(lin), '_'), function(x) x[1]))
lin = split(lin, p)
mat = do.call(rbind, lin)

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
n = seq(0, 1, 0.1)
y = seq(0.45, 1, 0.05)
text(1, 1, 'P int', adj = 0, cex = 1.25)
for(i in 1:length(y)){
  text(1, y[i], rev(n)[i], col = rev(cols)[i], adj = 0, cex = 1.25)
}

lin = unlist(lapply(res, function(x) x$nonlinear))
p = unlist(lapply(strsplit(names(lin), '_'), function(x) x[1]))
lin = split(lin, p)
mat = do.call(rbind, lin)

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
n = seq(0, 1, 0.1)
y = seq(0.45, 1, 0.05)
text(1, 1, 'P int', adj = 0, cex = 1.25)
for(i in 1:length(y)){
  text(1, y[i], rev(n)[i], col = rev(cols)[i], adj = 0, cex = 1.25)
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