setwd('~/Documents/Research/github/accounting-for-nonlinear-phenotypes/')
source('01_code/R/nonlinear_phenotype_utils.R')

#######################################
#####Load and clean phenotype data#####
#######################################
#Set working directory
setwd('~/Desktop/nonlinear_phenotype_review/00_data/arapheno_phenotypes/')

#List directories
files = list.files()
files = files[1:(length(files)-1)]

#Loop through and load
phenos = list()
for(i in 1:length(files)){
  setwd(files[i])
  z = list.files()[grep('values.csv', list.files())]
  x = read.csv(z)
  phenos[[as.character(files[i])]] = as.data.frame(x[,3:ncol(x)])
  setwd('../')
}

#Count n phenotypes and samples in each study
n = unlist(lapply(phenos, function(x) ncol(x)))
s = unlist(lapply(phenos, function(x) nrow(x)))
plot(n, 
     s,
     cex = 1.5,
     xlab = 'n phenotypes',
     ylab = 'n samples',
     pch = 20,
     cex.axis = 1.5,
     cex.lab = 1.5)

#Filter on phenotype #
phenos = phenos[which(unlist(lapply(phenos, function(x) ncol(x)))>=5)]

#Filter on NAs
f = list()
for(i in 1:length(phenos)){
  
  if(ncol(phenos[[i]])>1){
  
  #Filter strains
  x = phenos[[i]][apply(phenos[[i]], 1, function(x) sum(!is.na(x)))>20,]
  
  #Filter phenotypes
  x = x[,apply(x, 2, function(x) sum(!is.na(x)))>20]
  
  #Require at least 3 unique values
  x = x[,apply(x, 2, function(x) length(unique(x)))>3]
    
  #Filter strains
  f[[names(phenos)[i]]] = x
  }
}

#Filter on phenotype #
phenos = f[unlist(lapply(f, function(x) ncol(x)))>10]

#Impute and rank normalize
phenos_s = list()
for(h in 1:length(phenos)){
  
  #Impute for PCA (using column mean)
  for(i in 1:ncol(phenos[[h]])){
    x = which(is.na(phenos[[h]][,i]))
    phenos[[h]][x,i] = mean(phenos[[h]][,i], na.rm = TRUE)
  }
  
  #Rank normalization
  phenos_s[[names(phenos)[h]]] = as.data.frame(apply(phenos[[h]], 2, function(x) RankNorm(x)))
} 

#Filter on unique measurements per phenotype
phenos_s = lapply(phenos_s, function(x) x[,apply(x, 2, function(y) length(unique(y)))>=5])

#Save
saveRDS(phenos_s, '~/Documents/Research/github/accounting-for-nonlinear-phenotypes/02_output/arapheno_cleaned_phenos.RDS')

###########################################################
#####Compare entropy distributions across experiments?#####
###########################################################
#Set working directory and load
setwd('~/Documents/Research/github/accounting-for-nonlinear-phenotypes/')
phenos_s = readRDS('02_output/arapheno_cleaned_phenos.RDS')

#Sweep through and identify optimal binning size
toTest = c(0.0001, 0.001, seq(0.01, 0.1, 0.01))
res = list()
for(k in 1:length(toTest)){
  print(paste(k, 'out of', length(toTest)))
  entropies = list()
  for(a in 1:length(phenos_s)){
    ent = c()
    for(i in 1:nrow(phenos_s[[a]])){
      x = as.numeric(phenos_s[[a]][i,]+abs(min(phenos_s[[a]][i,])))
      x = x/max(x)
      x = round(x, 2)
      h = hist(x, breaks = seq(0, 1.1, toTest[k]), plot = FALSE)
      ent = c(ent, entropy::entropy(h$density))
    }
    entropies[[names(phenos_s)[a]]] = ent
  }
  res[[as.character(toTest[k])]] = unlist(lapply(entropies, function(x) mean(x)))
}

#Variance
plot(unlist(lapply(res, function(x) var(x))))

#Correlations
cors = c()
for(i in 1:(length(res)-1)){
  cors = c(cors, cor(res[[i]], res[[i+1]]))
}

#Loop through and calculate entropy using desired bin size
entropies = list()
seqs = list()
for(a in 1:length(phenos_s)){
  ent = c()
  s = list()
  for(i in 1:nrow(phenos_s[[a]])){
    x = as.numeric(phenos_s[[a]][i,]+abs(min(phenos_s[[a]][i,])))
    x = x/max(x)
    #h = hist(x, breaks = seq(0, 1, 0.0001), plot = FALSE)
    h = discretize(x, numBins = 1000)
    
    #Calculate entropy
    e = entropy::entropy.empirical(h, unit = 'log2')
    
    #Normalize entropy by dividing by maximum entropy (log2 of n unique(measurements))
    e = e/log2(sum(h>0))
    
    #Add to results
    ent = c(ent, e)
    s[[i]] = h
  }
  entropies[[names(phenos_s)[a]]] = ent
  seqs[[names(phenos_s)[a]]] = s
}

#Plot entropy distribution
vioplot::vioplot(entropies)

#Kruskal-wallist test
kruskal.test(entropies)

#Relationship with n phenos
plot(unlist(lapply(phenos_s, function(x) ncol(x))),
     lapply(entropies, function(x) mean(x)))

#Plot bin distributions
bins = lapply(seqs, function(x) colMeans(do.call(rbind, x)))

#Kullback-leibler divergence between bins
kl = c()
for(i in 1:(length(bins)-1)){
  kl = c(kl, KL.plugin(bins[[i]], bins[[i+1]]))
}

#################################################################
#####Entropy as a function of sample size (via permutations)#####
#################################################################
#Function
subsample.entropy = function(data,
                             sample_sizes = seq(10, 100, 10),
                             permutations = 10,
                             numBins = 100,
                             normalize_entropy = FALSE,
                             plot = FALSE){
  
  out = list()
  distributions = list()
  for(a in 1:length(sample_sizes)){
    
    #print(paste(a, 'out of', length(sample_sizes)))
    
    entropies = list()
    seqs = list()
    for(j in 1:perms){
      
      dat = data[,sample(ncol(data), sample_sizes[a])]
      
      #ent = c()
      #s = list()
      # for(i in 1:nrow(dat)){
      #   x = as.numeric(dat[i,]+abs(min(dat[i,])))
      #   x = x/max(x)
      #   h = discretize(x, numBins = numBins, r = c(0,1))
      #   
      #   #Calculate entropy
      #   e = entropy::entropy.empirical(h, unit = 'log2')
      #   
      #   #Normalize entropy by dividing by maximum entropy (log2 of n unique(measurements))
      #   if(normalize_entropy == TRUE){
      #     e = e/log2(sum(h>0))
      #     
      #     ent = c(ent, e)
      #   }else{
      #     ent = c(ent, e)
      #   }
      #   s[[i]] = h
      # }
      
      x = cov(dat)
      e = eigen(x)
      if(normalize_entropy == TRUE){
        ent = entropy::entropy.empirical(e$values, unit = 'log2')/log2(length(e$values))
      }else{
        ent = entropy::entropy.empirical(e$values, unit = 'log2')
      }
      
      entropies[[as.character(j)]] = ent
      seqs[[as.character(j)]] = e$values
      #seqs[[as.character(j)]] = colMeans(do.call(rbind, s))
    }
    out[[as.character(sample_sizes[a])]] = entropies
    distributions[[as.character(sample_sizes[a])]] = seqs
  }
  
  out = lapply(out, function(x) unlist(x))
  
  #Plot if desired
  if(plot == TRUE){
    vioplot::vioplot(out,
                     ylab = 'Entropy',
                     cex.axis = 1.5,
                     cex.lab = 1.5,
                     cex.sub = 1.5,
                     xlab = 'n phenotypes')
  }
  
  #Return
  l = list(out, distributions)
  names(l) = c('entropies', 'distributions')
  return(l)
}

#Run
all_entropies = list()
for(i in 1:length(phenos_s)){
  print(names(phenos_s)[i])
  all_entropies[[names(phenos_s)[i]]] = subsample.entropy(phenos_s[[i]],
                                                          numBins = 100,
                                                          permutations = 1000,
                                                          normalize_entropy = TRUE,
                                                          sample_sizes = round(seq(0.1, 0.9, 0.05)*ncol(phenos_s[[i]])),
                                                          plot = TRUE)}

#Kullback-leibler divergence
kls = list()
for(h in 1:length(all_entropies)){
  bins = lapply(all_entropies[[h]]$distributions, function(x) colMeans(do.call(rbind, x)))
  kl = c()
  for(i in 1:(length(bins)-1)){
    kl = c(kl, KL.plugin(bins[[i]], bins[[i+1]]))
  }
  kls[[h]] = kl
  plot(kl, type = 'l')
}

plot(kls[[1]], 
     ylim =c(0,0.4), 
     type = 'l',
     lwd = 1.5)
for(i in 2:length(kls)){
  lines(kls[[i]],
        lwd = 1.5)
}

#############################################
#####Calculate phenotypic non-linearity######
#############################################
#Run function on all experiments
res = list()
for(i in 1:length(phenos_s)){
  print(paste(i, 'out of', length(phenos_s)))
  res[[as.character(i)]] = compare_nonlinear(phenos_s[[i]], verbose = TRUE)
}

#Compare to entropy
non = lapply(res, function(x) x$nonlinear)

#Non-linearity vs. entropy
plot(lapply(entropies, function(x) mean(x)), 
     non,
     ylab = '% nonlinear',
     xlab = 'Mean entropy',
     cex.axis = 1.5,
     cex.lab = 1.5,
     bty = 'n',
     pch = 20,
     cex = 1.5)

#Determinant of covariation vs. entropy
plot(lapply(entropies, function(x) mean(x)),
     log(unlist(lapply(phenos_s, function(x) det(cov(x))))), 
     ylab = 'Determinant of covariation (log)',
     xlab = 'Mean entropy',
     cex.axis = 1.5,
     cex.lab = 1.5,
     bty = 'n',
     pch = 20,
     cex = 1.5)

#Sample size vs. non-linearity
plot(lapply(phenos_s, function(x) ncol(x)),
     non, 
     ylab = '% nonlinear',
     xlab = 'n phenotypes',
     cex.axis = 1.5,
     cex.lab = 1.5,
     bty = 'n',
     pch = 20,
     cex = 1.5)


library(BioQC)
spec = list()
for(i in 1:length(phenos_s)){
  x = as.numeric(phenos_s[[i]]+abs(min(phenos_s[[i]])))
  x = x/max(x)
  spec[[i]] = mean(entropySpecificity(as.matrix(x)))))
}

plot(unlist(lapply(phenos_s, function(x) mean(entropySpecificity(as.matrix(x))))),
     lapply(phenos_s, function(x) ncol(x)))


plot(unlist(lapply(phenos_s, function(x) mean(entropyDiversity(as.matrix(x))))),
     lapply(phenos_s, function(x) ncol(x)))



x = cov(phenos_s[[1]])
e = eigen(x)
entropy::entropy.empirical(e$values, unit = 'log2')/log2(length(e$values))


out = list()
for(i in 1:length(phenos_s)){
  x = cov(phenos_s[[i]])
  e = eigen(x)
  out[[i]] = entropy::entropy.empirical(e$values, unit = 'log2')/log2(length(e$values))
}

#Sample size vs. non-linearity
plot(log(unlist(lapply(phenos_s, function(x) det(cov(x))))),
     unlist(out), 
     ylab = '% nonlinear',
     xlab = 'n phenotypes',
     cex.axis = 1.5,
     cex.lab = 1.5,
     bty = 'n',
     pch = 20,
     cex = 1.5)






