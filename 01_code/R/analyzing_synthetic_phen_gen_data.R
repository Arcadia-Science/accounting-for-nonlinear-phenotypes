library(reticulate)
library(here)
np = import('numpy')

#Set working directory
setwd('~/Documents/Research/github/accounting-for-nonlinear-phenotypes/')

#Source python script
source_python('01_code/python/tools_for_phen_gen_creation.py')

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

##############################################
#####Sweeping number of important alleles#####
##############################################
#Sweep n of important alleles
toTest = 1:100
out = list()
pb <- txtProgressBar(min = 1,      
                     max = length(toTest), 
                     style = 3,    
                     width = 50,
                     char = ".")

for(i in 1:length(toTest)){
  setTxtProgressBar(pb, i)
  out[[as.character(toTest[i])]] = make_genotype(n_loci_ip = as.integer(i),  
                                                 n_loci = as.integer(1000), 
                                                 n_phens = as.integer(10))}

#Combine phenotypes
phenos = lapply(out, function(x) do.call(cbind, x[[3]]))

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
  #cat(counter, "\n")
  compare_nonlinear(x)})

