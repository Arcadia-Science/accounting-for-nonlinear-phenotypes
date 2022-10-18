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
library(reticulate)
np = import('numpy')

###################
#####Functions#####
###################
#Function to compare linear vs. non-linear models for a phenotype matrix
compare_nonlinear = function(dat, 
                             verbose = FALSE, 
                             return_models = FALSE){
  
  #Get all comparisons 
  all = expand.grid(1:ncol(dat), 1:ncol(dat))
  
  #Get just unique combos
  all = all[!all[,1] == all[,2],]
  
  #Remove reciprocal pairings
  for(i in 1:nrow(all)){
    all[i,] = sort(unlist(all[i,]))
  }
  all = all[-duplicated(all),]
  
  #Convert data to numeric matrix
  z = data.matrix(dat)
  
  #Test linear vs. nonlinear models via AIC for all traits
  res = list()
  
  if(verbose == TRUE){
    
    #Set progress bar
    pb <- txtProgressBar(min = 1,      
                         max = nrow(all), 
                         style = 3,    
                         width = 50,
                         char = ".")
  }
  for(i in 1:nrow(all)){
    
    if(verbose == TRUE){
      #Update progress bar
      setTxtProgressBar(pb, i)
    }
    
    mod1 = lm(z[,all[i,1]]~z[,all[i,2]])
    if(length(unique(z[,all[i,2]]))<=50){
      mod2 = gam(z[,all[i,1]]~s(z[,all[i,2]],
                                k = length(unique(z[,all[i,2]]))-1))
    }else{
      mod2 = gam(z[,all[i,1]]~s(z[,all[i,2]],
                                k = 50))
    }
    
    out = lrtest(mod1, mod2)
    a = AIC(mod1, mod2)
    
    l = list(out, a)
    names(l) = c('lrtest', 'AIC')
    
    res[[as.character(i)]] = l
  }
  
  #Calculate proportions
  x = unlist(lapply(res, function(x) x$lrtest$LogLik[1]-x$lrtest$LogLik[2]))
  linear = sum(x>0)/length(x)
  nonlinear = sum(x<0)/length(x)
  
  #Return
  if(return_models == TRUE){
    r = list(res, linear, nonlinear)
    names(r) = c('models', 'linear', 'nonlinear')
    return(r)
  }else{
    r = list(linear, nonlinear)
    names(r) = c('linear', 'nonlinear')
    return(r)
  }
}

#Function to run linearity test on independent random data
distribution_linearity_test_independent = function(n_reps = 10,
                                                   n_samples = 1000,
                                                   n_phenotypes = 20, 
                                                   mean = 10,
                                                   sd = 1,
                                                   distribution = c('poisson', 'uniform', 'gaussian')){
  
  #Set number of replicates
  n = n_reps
  
  #Progress bar
  pb <- txtProgressBar(min = 1,      
                       max = n, 
                       style = 3,    
                       width = 50,
                       char = ".")
  
  #Create empty list for results
  linear_nonlinear = list()
  
  #Loop through and calculate
  for(h in 1:n){
    
    #print(h)
    
    #Update progress bar
    setTxtProgressBar(pb, h)
    
    #Generate synthetic data (via Poisson distribution)
    dat = list()
    for(i in 1:n_phenotypes){
      if(distribution == 'poisson'){
        dat[[i]] = rpois(n_samples,mean)
      }else if(distribution == 'uniform'){
        dat[[i]] = runif(n_samples)
      }else if(distribution == 'gaussian'){
        dat[[i]] = rnorm(n_samples, mean, sd)
      }
    }
    dat = do.call(cbind, dat)
    colnames(dat) = as.character(seq(1, ncol(dat), 1))
    
    #Get all comparisons 
    all = expand.grid(1:ncol(dat), 1:ncol(dat))
    
    #Get just unique combos
    all = all[!all[,1] == all[,2],]
    
    #Remove reciprocal pairings
    for(i in 1:nrow(all)){
      all[i,] = sort(unlist(all[i,]))
    }
    all = all[-duplicated(all),]
    
    #Convert data to numeric matrix
    z = data.matrix(dat)
    
    #Test linear vs. nonlinear models via AIC for all traits
    res = list()
    for(i in 1:nrow(all)){
      
      mod1 = lm(z[,all[i,1]]~z[,all[i,2]])
      mod2 = gam(z[,all[i,1]]~s(z[,all[i,2]], k = length(unique(z[,all[i,2]]))-1))
      
      out = lrtest(mod1, mod2)
      a = AIC(mod1, mod2)
      
      l = list(out, a)
      names(l) = c('lrtest', 'AIC')
      
      res[[paste(colnames(dat)[all[i,1]],
                 colnames(dat)[all[i,2]], 
                 sep = '_')]] = l
    }
    
    #Calculate proportions
    x = unlist(lapply(res, function(x) x$lrtest$LogLik[1]-x$lrtest$LogLik[2]))
    linear = sum(x>0)/length(x)
    nonlinear = sum(x<0)/length(x)
    #print(c(linear, nonlinear))
    
    #Return
    r = list(linear, nonlinear)
    names(r) = c('linear', 'nonlinear')
    linear_nonlinear[[as.character(h)]] = r
  }
  
  #Return
  return(linear_nonlinear)
}

#Examples
poisson = distribution_linearity_test_independent(distribution = 'poisson', n_phenotypes = 20, n_samples = 1000)
gaussian = distribution_linearity_test_independent(distribution = 'gaussian', n_phenotypes = 20, n_samples = 1000)
uniform = distribution_linearity_test_independent(distribution = 'uniform', n_phenotypes = 20, n_samples = 1000)

#Plot
linear = lapply(poisson, function(x) x$linear)
nonlinear = lapply(poisson, function(x) x$nonlinear)

barplot(c(linear, nonlinear),
        ylab = 'Proportion',
        ylim = c(0,1),
        cex.axis = 1.5, cex.lab = 1.5, cex.names = 1.5,
        names = c('linear', 'nonlinear'),
        las = 2)
