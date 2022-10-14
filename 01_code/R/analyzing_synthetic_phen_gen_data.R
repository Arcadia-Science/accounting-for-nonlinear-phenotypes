#Set working directory
setwd('~/Documents/Research/github/accounting-for-nonlinear-phenotypes/')

#Source utility functions
source('01_code/R/nonlinear_phenotype_utils.R')

#Source python script for generating synthetic genotype/phenotype data
source_python('01_code/python/tools_for_phen_gen_creation.py')

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

